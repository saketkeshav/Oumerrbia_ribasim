# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 15:46:10 2020

@author: keshav
"""

import glob
import os
import random
import shutil
import string
import pandas as pd
import numpy as np


import subprocess as subprocess

import struct as struct
import math as math

from datetime import datetime, timedelta

import time


from ema_workbench.em_framework.optimization import (HyperVolume,
                                                     EpsilonProgress)
from ema_workbench import (RealParameter, ScalarOutcome,
                           MultiprocessingEvaluator, ema_logging,
                           SequentialEvaluator)
from ema_workbench.em_framework.model import AbstractModel, SingleReplication
from ema_workbench.util.ema_logging import method_logger, get_module_logger


_logger = get_module_logger(__name__)


class HisFile:
    def __init__(self, fname):
        self.fname = fname
        self.sysnames = list()
        self.segnames = list()
        self.segnums = list()
        self.datetime = list()
        self.nstep=0

    def read(self):
        f = open(self.fname, 'rb')
        self.moname = f.read(160)
        year = int(self.moname[124:128])
        month = int(self.moname[129:131])
        day = int(self.moname[132:134])
        hour = int(self.moname[135:137])
        minute = int(self.moname[138:140])
        second = int(self.moname[141:143])
        self.startdate = datetime(year, month, day, hour, minute, second)
        self.scu = int(self.moname[150:158])
        #print(self.scu)
        self.nsys, = struct.unpack('i', f.read(4))
        self.nseg, = struct.unpack('i', f.read(4))
        #return self.nsys, self.nseg
        #print self.nsys, self.nseg
        for isys in range(self.nsys):
            self.sysnames.append(str.rstrip(str(f.read(20))))
            #print isys,self.sysnames[isys]
        #return self.sysnames
        for iseg in range(self.nseg):
            self.segnums.append(struct.unpack('i', f.read(4))[0])
            self.segnames.append(str.rstrip(str(f.read(20))))
        #return self.segnames, self.sysnames#print iseg,self.segnums[iseg],self.segnames[iseg]
        size = os.path.getsize(self.fname)
        self.nstep += (size - 168 - 20 * self.nsys - 24 * self.nseg) / (4 * (self.nsys * self.nseg + 1))
        #print self.nstep
        
        self.data = np.zeros((self.nsys, self.nseg, int(self.nstep)), 'f')
        for lstep in range(int(self.nstep)):
            fdays = struct.unpack('i', f.read(4))[0] * (self.scu / 86400.)
            iday = math.trunc(fdays)
            ihour = math.trunc((fdays - iday) * 24.)
            iminute = math.trunc((fdays - iday - ihour / 24.) * 60.)
            isecond = fdays * 86400 - 60 * math.trunc(fdays * 86400 / 60.)
            dt = self.startdate + timedelta(iday, ihour, iminute, isecond)
            self.datetime.append(dt)
            for iseg in range(self.nseg):
                for isys in range(self.nsys):
                    self.data[isys, iseg, lstep] = struct.unpack('f', f.read(4))[0]
        f.close()
    

    def getseries(self, sysnam, segnam):
        isys = self.sysnames.index(sysnam)
        iseg = self.segnames.index(segnam)
        #print isys,iseg
        res = np.zeros(int(self.nstep))
        for lstep in range(int(self.nstep)):
            res[lstep] = self.data[isys, iseg, lstep]
        return res

    def gettimeseries(self, sysnam, segnam):
        isys = self.sysnames.index(sysnam)
        iseg = self.segnames.index(segnam)
        #print isys,iseg
        resdata = np.zeros(int(self.nstep))
        resdate = list()
        for lstep in range(int(self.nstep)):
            resdate.append(self.datetime[lstep])
            resdata[lstep] = self.data[isys, iseg, lstep]
        return resdate, resdata#, self.nstep


def generate(a, maximum, k):
    b = (maximum-a)*k+a
    return b

def ribasimmodel(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,
                 fracb1,fracb2,fracb3,fracb4,fracb5,fracb6,fracb7,fracb8,fracb9,fracb10,fracb11,fracb12,
                 fracc1,fracc2,fracc3,fracc4,fracc5,fracc6,fracc7,fracc8,fracc9,fracc10,fracc11,fracc12,
                 f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,
                 frace1,frace2,frace3,frace4,frace5,frace6,frace7,frace8,frace9,frace10,frace11,frace12,
                 fracf1,fracf2,fracf3,fracf4,fracf5,fracf6,fracf7,fracf8,fracf9,fracf10,fracf11,fracf12,
                 i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,
                 frach1,frach2,frach3,frach4,frach5,frach6,frach7,frach8,frach9,frach10,frach11,frach12,
                 fraci1,fraci2,fraci3,fraci4,fraci5,fraci6,fraci7,fraci8,fraci9,fraci10,fraci11,fraci12,
                 l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,
                 frack1,frack2,frack3,frack4,frack5,frack6,frack7,frack8,frack9,frack10,frack11,frack12,
                 fracl1,fracl2,fracl3,fracl4,fracl5,fracl6,fracl7,fracl8,fracl9,fracl10,fracl11,fracl12,
                 o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,
                 fracn1,fracn2,fracn3,fracn4,fracn5,fracn6,fracn7,fracn8,fracn9,fracn10,fracn11,fracn12,
                 fraco1,fraco2,fraco3,fraco4,fraco5,fraco6,fraco7,fraco8,fraco9,fraco10,fraco11,fraco12):

    #tempfolder=tempfile.TemporaryDirectory(dir="C:\\")
    #Model_folder=os.path.join(tempfolder.name[0:3],tempfolder.name[6:])
    #src="C:\\Oer0T_1"
    #TarFolder=Model_folder
    #distutils.dir_util.copy_tree(src,TarFolder)
    #TarFolder=Model_folder+'\\OUMRBIA9.Rbn\\CMTWORK'
    #files = glob.glob(os.path.join(TarFolder,'*.*'))
    #for file in files:
        #f=open(file,'rt')
        #textin=f.readlines()
        #f.close()
        #textout=list()
        #for line in textin:
            #line=line.replace('\RIBASIM7','\Oer0T_1')
            #line=line.replace('OER0T','OER0T{}'.format(inst))
            #textout.append(line)
            #f=open(file,'wt')
            #f.writelines(textout)
            #f.close()
    #TarFolder=Model_folder+'\\OUMRBIA9.Rbn\\3'
    #files = glob.glob(os.path.join(TarFolder,'TIMESERI.*'))
    #for file in files:
        #f=open(file,'rt')
        #textin=f.readlines()
        #f.close()
        #textout=list()
        #for line in textin:
            #line=line.replace('\RIBASIM7','\Oer0T_1')
            #line=line.replace('OER0T','OER0T{}'.format(inst))
            #textout.append(line)
        #f=open(file,'wt')
        #f.writelines(textout)
        #f.close()
    
    #TarFolder=Model_folder+'\\OUMRBIA9.Rbn\\WORK'
    #files = glob.glob(os.path.join(TarFolder,'TIMESERI.*'))
    #for file in files:
     #   f=open(file,'rt')
     #   textin=f.readlines()
     #   f.close()
     #   textout=list()
     #   for line in textin:
     #       line=line.replace('\Ribasim7',Model_folder[2:])
            #line=line.replace('OER0T','OER0T{}'.format(inst))
     #       textout.append(line)
     #   f=open(file,'wt')
     #   f.writelines(textout)
     #   f.close()    
    #TarFolder=Model_folder+'\\Programs\\RL_PSTPR.INI'
    #f=open(TarFolder,'rt')
    #textin=f.readlines()
    #f.close()
    #textout=list()
    #for line in textin:
        #line=line.replace('C:\Ribasim7','C:\Oer0T_1')
        #textout.append(line)
    #f=open(TarFolder,'wt')
    #f.writelines(textout)
    #f.close()
    #TarFolder=Model_folder+'\\Programs\\rib2his4.fnm'
    #f=open(TarFolder,'rt')
    #textin=f.readlines()
    #f.close()
    #textout=list()
    #for line in textin:
        #line=line.replace('C:\Ribasim7', 'C:\Oer0T_1')
        #textout.append(line)
    #f=open(TarFolder,'wt')
    #f.writelines(textout)
    #f.close()

    #time.sleep(3)
    hisfile1 =  os.path.join('C:\\Oer0T_1', 'OUMRBIA9.Rbn', 'WORK', 'pws.his')
    hisfile2 =  os.path.join('C:\\Oer0T_1', 'OUMRBIA9.Rbn', 'WORK', 'reservoi.his')
    hisfile3 =  os.path.join('C:\\Oer0T_1', 'OUMRBIA9.Rbn', 'WORK', 'vir.his')
    hisfile4 =  os.path.join('C:\\Oer0T_1', 'OUMRBIA9.Rbn', 'WORK', 'fir.his')
    def generate(a,maximum,k):
        b=(maximum-a)*k+a
        return b

    
    part1max = 677.00
    
    
    
    part2max = 810.00
        
    
    part3max = 966.00
       
    
    part4max = 877.50
    
    
    part5max = 285.00
    i=1
    b1=generate(c1,part1max,fracb1);b2=generate(c2,part1max,fracb2);b3=generate(c3,part1max,fracb3);b4=generate(c4,part1max,fracb4);b5=generate(c5,part1max,fracb5);b6=generate(c6,part1max,fracb6);b7=generate(c7,part1max,fracb7);b8=generate(c8,part1max,fracb8);b9=generate(c9,part1max,fracb9);b10=generate(c10,part1max,fracb10);b11=generate(c11,part1max,fracb11);b12=generate(c12,part1max,fracb12)
    a1=generate(b1,part1max,fracc1);a2=generate(b2,part1max,fracc2);a3=generate(b3,part1max,fracc3);a4=generate(b4,part1max,fracc4);a5=generate(b5,part1max,fracc5);a6=generate(b6,part1max,fracc6);a7=generate(b7,part1max,fracc7);a8=generate(b8,part1max,fracc8);a9=generate(b9,part1max,fracc9);a10=generate(b10,part1max,fracc10);a11=generate(b11,part1max,fracc11);a12=generate(b12,part1max,fracc12)
    e1=generate(f1,part2max,frace1);e2=generate(f2,part2max,frace2);e3=generate(f3,part2max,frace3);e4=generate(f4,part2max,frace4);e5=generate(f5,part2max,frace5);e6=generate(f6,part2max,frace6);e7=generate(f7,part2max,frace7);e8=generate(f8,part2max,frace8);e9=generate(f9,part2max,frace9);e10=generate(f10,part2max,frace10);e11=generate(f11,part2max,frace11);e12=generate(f12,part2max,frace12)
    d1=generate(e1,part2max,fracf1);d2=generate(e2,part2max,fracf2);d3=generate(e3,part2max,fracf3);d4=generate(e4,part2max,fracf4);d5=generate(e5,part2max,fracf5);d6=generate(e6,part2max,fracf6);d7=generate(e7,part2max,fracf7);d8=generate(e8,part2max,fracf8);d9=generate(e9,part2max,fracf9);d10=generate(e10,part2max,fracf10);d11=generate(e11,part2max,fracf11);d12=generate(e12,part2max,fracf12)
    h1=generate(i1,part3max,frach1);h2=generate(i2,part3max,frach2);h3=generate(i3,part3max,frach3);h4=generate(i4,part3max,frach4);h5=generate(i5,part3max,frach5);h6=generate(i6,part3max,frach6);h7=generate(i7,part3max,frach7);h8=generate(i8,part3max,frach8);h9=generate(i9,part3max,frach9);h10=generate(i10,part3max,frach10);h11=generate(i11,part3max,frach11);h12=generate(i12,part3max,frach12)
    g1=generate(h1,part3max,fraci1);g2=generate(h2,part3max,fraci2);g3=generate(h3,part3max,fraci3);g4=generate(h4,part3max,fraci4);g5=generate(h5,part3max,fraci5);g6=generate(h6,part3max,fraci6);g7=generate(h7,part3max,fraci7);g8=generate(h8,part3max,fraci8);g9=generate(h9,part3max,fraci9);g10=generate(h10,part3max,fraci10);g11=generate(h11,part3max,fraci11);g12=generate(h12,part3max,fraci12)
    k1=generate(l1,part4max,frack1);k2=generate(l2,part4max,frack2);k3=generate(l3,part4max,frack3);k4=generate(l4,part4max,frack4);k5=generate(l5,part4max,frack5);k6=generate(l6,part4max,frack6);k7=generate(l7,part4max,frack7);k8=generate(l8,part4max,frack8);k9=generate(l9,part4max,frack9);k10=generate(l10,part4max,frack10);k11=generate(l11,part4max,frack11);k12=generate(l12,part4max,frack12)
    j1=generate(k1,part4max,fracl1);j2=generate(k2,part4max,fracl2);j3=generate(k3,part4max,fracl3);j4=generate(k4,part4max,fracl4);j5=generate(k5,part4max,fracl5);j6=generate(k6,part4max,fracl6);j7=generate(k7,part4max,fracl7);j8=generate(k8,part4max,fracl8);j9=generate(k9,part4max,fracl9);j10=generate(k10,part4max,fracl10);j11=generate(k11,part4max,fracl11);j12=generate(k12,part4max,fracl12)
    n1=generate(o1,part5max,fracn1);n2=generate(o2,part5max,fracn2);n3=generate(o3,part5max,fracn3);n4=generate(o4,part5max,fracn4);n5=generate(o5,part5max,fracn5);n6=generate(o6,part5max,fracn6);n7=generate(o7,part5max,fracn7);n8=generate(o8,part5max,fracn8);n9=generate(o9,part5max,fracn9);n10=generate(o10,part5max,fracn10);n11=generate(o11,part5max,fracn11);n12=generate(o12,part5max,fracn12)
    m1=generate(n1,part5max,fraco1);m2=generate(n2,part5max,fraco2);m3=generate(n3,part5max,fraco3);m4=generate(n4,part5max,fraco4);m5=generate(n5,part5max,fraco5);m6=generate(n6,part5max,fraco6);m7=generate(n7,part5max,fraco7);m8=generate(n8,part5max,fraco8);m9=generate(n9,part5max,fraco9);m10=generate(n10,part5max,fraco10);m11=generate(n11,part5max,fraco11);m12=generate(n12,part5max,fraco12)
                                          
    

    textfile=os.path.join("C:\\Oer0T_1\\OUMRBIA9.Rbn\\Actions\\Measures\\R001-2015-RibasimRsvNodeOperRulesDataOER.mes")
    f = open(textfile, 'r+')
    a=f.readlines()
    for i in range(len(a)):
        if (a[i]=='Node name=RSV_AHMEDELHANSALIDAM\n'): 
            a[i+2]='Flood control storage operation rule (t) (m)='+str(a1)+','+str(a2)+','+str(a3)+','+str(a4)+','+str(a5)+','+str(a6)+','+str(a7)+','+str(a8)+','+str(a9)+','+str(a10)+','+str(a11)+','+str(a12)+'\n'
            a[i+3]='Target storage operation rule (t) (m)='+ str(b1)+','+str(b2)+','+str(b3)+','+str(b4)+','+str(b5)+','+str(b6)+','+str(b7)+','+str(b8)+','+str(b9)+','+str(b10)+','+str(b11)+','+str(b12)+'\n'
            a[i+4]='Firm storage operation rule (t) (m)='+ str(c1)+','+str(c2)+','+str(c3)+','+str(c4)+','+str(c5)+','+str(c6)+','+str(c7)+','+str(c8)+','+str(c9)+','+str(c10)+','+str(c11)+','+str(c12)+'\n'
        if (a[i]=='Node name=RSV_BINELOUIDANEDAM\n'):
            a[i+2]='Flood control storage operation rule (t) (m)='+str(d1)+','+str(d2)+','+str(d3)+','+str(d4)+','+str(d5)+','+str(d6)+','+str(d7)+','+str(d8)+','+str(d9)+','+str(d10)+','+str(d11)+','+str(d12)+'\n'
            a[i+3]='Target storage operation rule (t) (m)='+ str(e1)+','+str(e2)+','+str(e3)+','+str(e4)+','+str(e5)+','+str(e6)+','+str(e7)+','+str(e8)+','+str(e9)+','+str(e10)+','+str(e11)+','+str(e12)+'\n'
            a[i+4]='Firm storage operation rule (t) (m)='+ str(f1)+','+str(f2)+','+str(f3)+','+str(f4)+','+str(f5)+','+str(f6)+','+str(f7)+','+str(f8)+','+str(f9)+','+str(f10)+','+str(f11)+','+str(f12)+'\n'
        if (a[i]=='Node name=RSV_HASSAN1ERDAM\n'):
            a[i+2]='Flood control storage operation rule (t) (m)='+str(g1)+','+str(g2)+','+str(g3)+','+str(g4)+','+str(g5)+','+str(g6)+','+str(g7)+','+str(g8)+','+str(g9)+','+str(g10)+','+str(g11)+','+str(g12)+'\n'
            a[i+3]='Target storage operation rule (t) (m)='+ str(h1)+','+str(h2)+','+str(h3)+','+str(h4)+','+str(h5)+','+str(h6)+','+str(h7)+','+str(h8)+','+str(h9)+','+str(h10)+','+str(h11)+','+str(h12)+'\n'
            a[i+4]='Firm storage operation rule (t) (m)='+ str(i1)+','+str(i2)+','+str(i3)+','+str(i4)+','+str(i5)+','+str(i6)+','+str(i7)+','+str(i8)+','+str(i9)+','+str(i10)+','+str(i11)+','+str(i12)+'\n'
        if (a[i]=='Node name=RSV_MOULAYYOUSSEFDAM\n'):
            a[i+2]='Flood control storage operation rule (t) (m)='+str(j1)+','+str(j2)+','+str(j3)+','+str(j4)+','+str(j5)+','+str(j6)+','+str(j7)+','+str(j8)+','+str(j9)+','+str(j10)+','+str(j11)+','+str(j12)+'\n'
            a[i+3]='Target storage operation rule (t) (m)='+ str(k1)+','+str(k2)+','+str(k3)+','+str(k4)+','+str(k5)+','+str(k6)+','+str(k7)+','+str(k8)+','+str(k9)+','+str(k10)+','+str(k11)+','+str(k12)+'\n'
            a[i+4]='Firm storage operation rule (t) (m)='+ str(l1)+','+str(l2)+','+str(l3)+','+str(l4)+','+str(l5)+','+str(l6)+','+str(l7)+','+str(l8)+','+str(l9)+','+str(l10)+','+str(l11)+','+str(l12)+'\n'
        if (a[i]=='Node name=RSV_ALMASSIRADAM\n'):
            a[i+2]='Flood control storage operation rule (t) (m)='+str(m1)+','+str(m2)+','+str(m3)+','+str(m4)+','+str(m5)+','+str(m6)+','+str(m7)+','+str(m8)+','+str(m9)+','+str(m10)+','+str(m11)+','+str(m12)+'\n'
            a[i+3]='Target storage operation rule (t) (m)='+ str(n1)+','+str(n2)+','+str(n3)+','+str(n4)+','+str(n5)+','+str(n6)+','+str(n7)+','+str(n8)+','+str(n9)+','+str(n10)+','+str(n11)+','+str(n12)+'\n'
            a[i+4]='Firm storage operation rule (t) (m)='+ str(o1)+','+str(o2)+','+str(o3)+','+str(o4)+','+str(o5)+','+str(o6)+','+str(o7)+','+str(o8)+','+str(o9)+','+str(o10)+','+str(o11)+','+str(o12)+'\n'
    f = open(textfile, 'wt')
    f.seek(0)
    f.writelines(a)

    f.close()
#break point
    os.chdir("C:\\Oer0T_1\\OUMRBIA9.Rbn\\CMTWORK")
    resrb = subprocess.call("C:\\Oer0T_1\\Programs\\Ribasim\\System\\Bin2prt.exe bin2prt.fnm") #Doubt
    if resrb <0 and resrb> 0:
        print('Error in RIBASIM:')
        print(resrb)
        exit()
    
    resrb = subprocess.call("C:\\Oer0T_1\\Programs\\Ribasim\\System\\Simproc.exe simproc.fnm") #Doubt
    if resrb <0 and resrb> 0:
        print('Error in RIBASIM:')
        print(resrb)
        exit()
    resrb = subprocess.call("C:\\Oer0T_1\\Programs\\Ribasim\\System\\Rib2His.exe rib2his4.fnm")
    if resrb <0 and resrb> 0:
        print('Error in postprocessing RIBASIM')
        print(resrb)
        exit()
        
    reshis = HisFile(hisfile1)
    reshis.read()
    segnam = "b'Pws_Marrakech:P2(CR)'"
    sysnam = "b'Shortage from networ'"
    rdate,dataMarrakech = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_ElJadida:P5(DDao'"
    rdate,dataEljadida = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_Azzemour(avalUsi'"
    rdate,dataazzemour = reshis.gettimeseries(sysnam, segnam)
    segnam =  "b'Pws_RuralDoukala(DIm'"
    rdate,dataruraldoukal = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_Safi:P3(DBSImfou'"
    rdate,datasafi = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Pws_Casablanca1(DSSM'"
    rdate,datacas=reshis.gettimeseries(sysnam,segnam)
    segnam = "b'Pws_GWTaval_Kel\\xe2aDes'"
    rdate,datagwtaval = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_MarrakechTransfe'"
    rdate,datamarrakechtransfe = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_Casablanca2_Sett'"
    rdate,datacasablanca = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_Ramnha_Benguerir'"
    rdate,databenguerir = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_Bejaad-Khourib-O'"
    rdate,databejaad = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_Azillal:P7+Kalaa'"
    rdate,dataazillal = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_Khenifra:P0(AH) '"
    rdate,datakhenifra = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_GW_BeniAmir     '"
    rdate,databeniamir = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_GW_BeniMoussa   '"
    rdate,databenimoussa = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_GW_Tadla(OCP)   '"
    rdate,datagwtadla = reshis.gettimeseries(sysnam, segnam)
    segnam =  "b'Pws_GWChaouiacotiere'"
    rdate,datagwchaouiacotiere = reshis.gettimeseries(sysnam, segnam)
    segnam =  "b'Pws_BeniMellal:P8(DA'"
    rdate,databenimellal = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_GWBahira        '"
    rdate,datagwbahira = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_GWDoukkalaSahel '"
    rdate,datagwdoukkalasahel = reshis.gettimeseries(sysnam, segnam)
    segnam = "b'Pws_GW_Dir          '"
    rdate,datagwdir = reshis.gettimeseries(sysnam, segnam)
    
    pws1=dataMarrakech
    #pws1avg=sum(pws1[880:-1])/len(pws1[880:-1])
    pws2=dataEljadida
    #pws2avg=sum(pws2[880:-1])/len(pws2[880:-1])
    pws3=dataazzemour
    #pws3avg=sum(pws3[880:-1])/len(pws3[880:-1])
    pws4=dataruraldoukal
    #pws4avg=sum(pws4[880:-1])/len(pws4[880:-1])
    pws21=datacas
    #pws21avg=sum(pws21[880:-1])/len(pws21[880:-1])
    pws5=datasafi
    #pws5avg=sum(pws5[880:-1])/len(pws5[880:-1])
    pws6=datagwtaval
    #pws6avg=sum(pws6[880:-1])/len(pws6[880:-1])
    pws7=datamarrakechtransfe
    #pws7avg=sum(pws7[880:-1])/len(pws7[880:-1])
    pws8=datacasablanca
    #pws8avg=sum(pws8[880:-1])/len(pws8[880:-1])
    pws9=databenguerir
    #pws9avg=sum(pws9[880:-1])/len(pws9[880:-1])
    pws10=databejaad
    #pws10avg=sum(pws10[880:-1])/len(pws10[880:-1])
    pws11=dataazillal
    #pws11avg=sum(pws11[880:-1])/len(pws11[880:-1])
    pws12=datakhenifra
    #pws12avg=sum(pws12[880:-1])/len(pws12[880:-1])
    pws13=databeniamir
    #pws13avg=sum(pws13[880:-1])/len(pws13[880:-1])
    pws14=databenimoussa
    #pws14avg=sum(pws14[880:-1])/len(pws14[880:-1])
    pws15=datagwtadla
    #pws15avg=sum(pws15[880:-1])/len(pws15[880:-1])
    pws16=datagwchaouiacotiere
    #pws16avg=sum(pws16[880:-1])/len(pws16[880:-1])
    pws17=databenimellal
    #pws17avg=sum(pws17[880:-1])/len(pws17[880:-1])
    pws18=datagwbahira
    #pws18avg=sum(pws18[880:-1])/len(pws18[880:-1])
    pws19=datagwdoukkalasahel
    #pws19avg=sum(pws19[880:-1])/len(pws19[880:-1])
    pws20=datagwdir
    #pws20avg=sum(pws20[880:-1])/len(pws20[880:-1])
    #pwsshortageavg=(pws1avg+pws2avg+pws3avg+pws4avg+pws5avg+pws6avg+pws7avg+pws8avg+pws9avg+pws10avg+pws11avg+pws12avg+pws13avg+pws14avg+pws15avg+pws16avg+pws17avg+pws18avg+pws19avg+pws20avg+pws21avg)/21
    pwssum=pws1+pws2+pws3+pws4+pws5+pws6+pws7+pws8+pws9+pws10+pws11+pws12+pws13+pws14+pws15+pws16+pws17+pws18+pws19+pws20+pws21
    pwssorted=sorted(pwssum)
    pwsobjective=sum(pwssorted[880:-1])
    
    
    reshis1 = HisFile(hisfile2)
    reshis1.read()
    segnam = "b'RSV_TAGZIRTDAM(P)   '"
    sysnam = "b'Energy: Generated (G'"
    rdate,datatagzirt = reshis1.gettimeseries(sysnam, segnam)
    segnam = "b'RSV_HASSAN1ERDAM    '"
    rdate,datahassan = reshis1.gettimeseries(sysnam, segnam)
    segnam =  "b'RSV_TYOUGHZADAM(P)  '"
    rdate,datatyoughza = reshis1.gettimeseries(sysnam, segnam)
    segnam =  "b'RSV_MOULAYYOUSSEFDAM'"
    rdate,datamoulayyoussef = reshis1.gettimeseries(sysnam, segnam)
    segnam =   "b'RSV_ALMASSIRADAM    '"
    rdate,dataalmassira = reshis1.gettimeseries(sysnam, segnam)
    segnam = "b'RSV_BINELOUIDANEDAM '"
    rdate,databinelouidane = reshis1.gettimeseries(sysnam, segnam)
    segnam = "b'RSV_AHMEDELHANSALIDA'"
    rdate,dataahmedhansal = reshis1.gettimeseries(sysnam, segnam)
    segnam = "b'Rsv_SIDIDRISSDAM    '"
    rdate,datasididriss = reshis1.gettimeseries(sysnam, segnam)
    
    
    energy1=datatagzirt
    energy2=datahassan
    energy3=datatyoughza
    energy4=datamoulayyoussef
    energy5=dataalmassira
    energy6=databinelouidane
    energy7=dataahmedhansal
    energy8=datasididriss
    energytotal=(energy1+energy2+energy3+energy4+energy5+energy6+energy7+energy8)
    energysorted=sorted(energytotal)
    energyobjective=sum(energysorted[880:-1])
        
    reshis = HisFile(hisfile3)
    reshis.read()
    segnam="b'Vir_PMHTessaoutAval:'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataPMHT = reshis.gettimeseries(sysnam, segnam)
    segnam="b'Vir_OeRavalMassira:I'"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataOeRaval = reshis.gettimeseries(sysnam, segnam)
    segnam="b'Vir_MoyenOeR:I4     '"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataMoyen = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_Haouz:I7        '"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataHaouz= reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_GH-TessaoutAvald'"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHTessaout = reshis.gettimeseries(sysnam, segnam)
    segnam="b'Vir_PMH:AmontTadla:I'"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataPMHAmont= reshis.gettimeseries(sysnam, segnam)
    segnam="b'Vir_PMH:I11         '"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataPMHI11 = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_GH_Haut_Doukkala'"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHHaut = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_GH_BeniAmir:I1  '"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataBeni = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_GHBeniMoussa:I5A'"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataMoussa = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_GH-TessaoutAvaln'"
   # sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHTessain = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_GH-TessaoutAmont'"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHTessaamont = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_GH_BasDoukkala:I'"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHBasDoukkala = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_Shtouka         '"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,datashtouka = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_IPAmontTadlaIp24'"
   # sysnam="b'Shortage (m3/s)     '"
    rdate,dataIp24 = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_IPMoyenOeRIp4   '"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataIp4moyen= reshis.gettimeseries(sysnam, segnam)
    segnam=  'b"Vir_Droitd\'eauLakhda"'
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataeau = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_IPLakhdaravalSid'"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,dataiplakhdasid = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_IPTadlaOeRIp3   '"
   # sysnam="b'Shortage (m3/s)     '"
    rdate,dataiptadlaoer = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_PMH_TessaoutInte'"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,datapmhtessa = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_Dir:I2          '"
   # sysnam="b'Shortage (m3/s)     '"
    rdate,dataDir = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_PMH:OERTadla:I3 '"
    #sysnam="b'Shortage (m3/s)     '"
    rdate,datapmhoertadla = reshis.gettimeseries(sysnam, segnam)
    vir1=dataPMHT
    vir2=dataOeRaval
    vir3=dataMoyen
    vir4=dataHaouz
    vir5=dataGHTessaout
    vir6=dataPMHAmont
    vir7=dataPMHI11
    vir8=dataGHHaut
    vir9=dataBeni
    vir10=dataMoussa
    vir11=dataGHTessain
    vir12=dataGHTessaamont 
    vir13=dataGHBasDoukkala 
    vir14=datashtouka 
    vir15=dataIp24 
    vir16=dataIp4moyen
    vir17=dataeau
    vir18=dataiplakhdasid
    vir19=dataiptadlaoer
    vir20=datapmhtessa 
    vir21=dataDir
    vir22=datapmhoertadla
    virsum=vir1+vir2+vir3+vir4+vir5+vir6+vir7+vir8+vir9+vir10+vir11+vir12+vir13+vir14+vir15+vir16+vir17+vir18+vir19+vir20+vir21+vir22
    virobjective=sorted(virsum)
    virobjective=sum(virobjective[880:-1])
    
    reshis = HisFile(hisfile4)
    reshis.read()
    sysnam= "b'Shortage (m3/s)     '"
    segnam= "b'Fir_AmontTyoughza:I2'"
    rdate, dataamonttyou=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_AmontElHansali:I'"
    rdate, dataamontelhansal=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_AmontTagzirt:I22'"
    rdate, dataamonttagzirt=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_AmontBinElOuidan'"
    rdate, dataamontbin=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_AmontHassan1er:I'"
    rdtae, dataamonthassan=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeSahel    '"
    rdate, dataIpnappe=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_Chaouiac\\xf4ti\\xe8re  '"
    rdate, datachaouiac=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeDoukkala '"
    rdate, dataIPnappeDoukkala=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeB\\xe9niAmir '"
    rdate, dataIPnappeBamir=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeB\\xe9niMouss'"
    rdate, dataIPnappeBMouss=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPTuroniendelaTa'"
    rdate, dataIPTuron=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeprofondeT'"
    rdate, dataIPnappepro=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeBMoussapr'"
    rdate, dataIPnappeBmoussa=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeDirprofon'"
    rdate, dataIPnappeDirprofon=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeBahira   '"
    rdate, dataIPnappeBahira=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeDir      '"
    rdate, dataIPnappeDir=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_IPnappeB\\xe9niAmirp'"
    rdate, dataIPnappeAmirp=reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Fir_PMHGhazef:I25   '"
    rdate, dataPMHGhazef=reshis.gettimeseries(sysnam, segnam)
    fir1=dataamonttyou
    fir2=dataamontelhansal
    fir3=dataamonttagzirt
    fir4=dataamontbin
    fir5= dataamonthassan
    fir6= dataIpnappe
    fir7=datachaouiac
    fir8=dataIPnappeDoukkala
    fir9=dataIPnappeBamir
    fir10=dataIPnappeBMouss
    fir11=dataIPTuron
    fir12=dataIPnappepro
    fir13=dataIPnappeBmoussa
    fir14=dataIPnappeDirprofon
    fir15=dataIPnappeBahira
    fir16=dataIPnappeDir
    fir17=dataIPnappeAmirp
    fir18=dataPMHGhazef
    firsum=fir1+fir2+fir3+fir4+fir5+fir6+fir7+fir8+fir9+fir10+fir11+fir12+fir13+fir14+fir15+fir16+fir17+fir18
    firobjective=sorted(firsum)
    firobjective=sum(firobjective[880:-1])
    
  
    return pwsobjective, energyobjective, virobjective


parameters=pd.read_csv("C:\\Users\\keshav\\OneDrive - Stichting Deltares\\Desktop\\results_250k_ribasim_final_second_run.csv")
pwsobjective_final=[]
energyobjective_final=[]
virobjective_final=[]
for i in range(len(parameters)):
    pws,energy,vir=ribasimmodel(parameters.c1[i],parameters.c2[i],parameters.c3[i],parameters.c4[i],parameters.c5[i],parameters.c6[i],parameters.c7[i],parameters.c8[i],parameters.c9[i],parameters.c10[i],parameters.c11[i],parameters.c12[i],
                                parameters.fracb1[i],parameters.fracb2[i],parameters.fracb3[i],parameters.fracb4[i],parameters.fracb5[i],parameters.fracb6[i],parameters.fracb7[i],parameters.fracb8[i],parameters.fracb9[i],parameters.fracb10[i],parameters.fracb11[i],parameters.fracb12[i],
                                parameters.fracc1[i],parameters.fracc2[i],parameters.fracc3[i],parameters.fracc4[i],parameters.fracc5[i],parameters.fracc6[i],parameters.fracc7[i],parameters.fracc8[i],parameters.fracc9[i],parameters.fracc10[i],parameters.fracc11[i],parameters.fracc12[i],
                                parameters.f1[i],parameters.f2[i],parameters.f3[i],parameters.f4[i],parameters.f5[i],parameters.f6[i],parameters.f7[i],parameters.f8[i],parameters.f9[i],parameters.f10[i],parameters.f11[i],parameters.f12[i],
                                parameters.frace1[i],parameters.frace2[i],parameters.frace3[i],parameters.frace4[i],parameters.frace5[i],parameters.frace6[i],parameters.frace7[i],parameters.frace8[i],parameters.frace9[i],parameters.frace10[i],parameters.frace11[i],parameters.frace12[i],
                                parameters.fracf1[i],parameters.fracf2[i],parameters.fracf3[i],parameters.fracf4[i],parameters.fracf5[i],parameters.fracf6[i],parameters.fracf7[i],parameters.fracf8[i],parameters.fracf9[i],parameters.fracf10[i],parameters.fracf11[i],parameters.fracf12[i],
                                parameters.i1[i],parameters.i2[i],parameters.i3[i],parameters.i4[i],parameters.i5[i],parameters.i6[i],parameters.i7[i],parameters.i8[i],parameters.i9[i],parameters.i10[i],parameters.i11[i],parameters.i12[i],
                                parameters.frach1[i],parameters.frach2[i],parameters.frach3[i],parameters.frach4[i],parameters.frach5[i],parameters.frach6[i],parameters.frach7[i],parameters.frach8[i],parameters.frach9[i],parameters.frach10[i],parameters.frach11[i],parameters.frach12[i],
                                parameters.fraci1[i],parameters.fraci2[i],parameters.fraci3[i],parameters.fraci4[i],parameters.fraci5[i],parameters.fraci6[i],parameters.fraci7[i],parameters.fraci8[i],parameters.fraci9[i],parameters.fraci10[i],parameters.fraci11[i],parameters.fraci12[i],
                                parameters.l1[i],parameters.l2[i],parameters.l3[i],parameters.l4[i],parameters.l5[i],parameters.l6[i],parameters.l7[i],parameters.l8[i],parameters.l9[i],parameters.l10[i],parameters.l11[i],parameters.l12[i],
                                parameters.frack1[i],parameters.frack2[i],parameters.frack3[i],parameters.frack4[i],parameters.frack5[i],parameters.frack6[i],parameters.frack7[i],parameters.frack8[i],parameters.frack9[i],parameters.frack10[i],parameters.frack11[i],parameters.frack12[i],
                                parameters.fracl1[i],parameters.fracl2[i],parameters.fracl3[i],parameters.fracl4[i],parameters.fracl5[i],parameters.fracl6[i],parameters.fracl7[i],parameters.fracl8[i],parameters.fracl9[i],parameters.fracl10[i],parameters.fracl11[i],parameters.fracl12[i],
                                parameters.o1[i],parameters.o2[i],parameters.o3[i],parameters.o4[i],parameters.o5[i],parameters.o6[i],parameters.o7[i],parameters.o8[i],parameters.o9[i],parameters.o10[i],parameters.o11[i],parameters.o12[i],
                                parameters.fracn1[i],parameters.fracn2[i],parameters.fracn3[i],parameters.fracn4[i],parameters.fracn5[i],parameters.fracn6[i],parameters.fracn7[i],parameters.fracn8[i],parameters.fracn9[i],parameters.fracn10[i],parameters.fracn11[i],parameters.fracn12[i],
                                parameters.fraco1[i],parameters.fraco2[i],parameters.fraco3[i],parameters.fraco4[i],parameters.fraco5[i],parameters.fraco6[i],parameters.fraco7[i],parameters.fraco8[i],parameters.fraco9[i],parameters.fraco10[i],parameters.fraco11[i],parameters.fraco12[i])
    pwsobjective_final.append(pws)
    energyobjective_final.append(energy)
    virobjective_final.append(vir)
    
np.savetxt("C:\\Users\\keshav\\OneDrive - Stichting Deltares\\Desktop\\Thesis report\\pwsobjective_final_remoteLAS_to_localLAS.csv",pwsobjective_final)
np.savetxt("C:\\Users\\keshav\\OneDrive - Stichting Deltares\\Desktop\\Thesis report\\energyobjective_final_remoteLAS_to_localLAS.csv",energyobjective_final)
np.savetxt("C:\\Users\\keshav\\OneDrive - Stichting Deltares\\Desktop\\Thesis report\\virobjective_final_remoteLAS_to_localLAS.csv",virobjective_final)