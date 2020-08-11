# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:49:25 2020

@author: keshav
"""
import numpy as np
import random as rd
import os as os
import subprocess as subprocess
import struct as struct
import math as math
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import multiprocessing as mp
from ema_workbench import Constraint
from ema_workbench.em_framework.optimization import (HyperVolume,
                                                     EpsilonProgress)
from ema_workbench import (Model, RealParameter, ScalarOutcome,
                           MultiprocessingEvaluator, ema_logging,
                           Constant)
import os
import glob
import distutils.dir_util
import sys
import shutil
from backports import tempfile
import operator
import time




#part1min=617.70
#part1max=677.00
#diff1=part1max-part1min

#part2min=741.55
#part2max=810.00
#diff2=part2max-part2min

#part3min=885.50
#part3max=966.00
#diff3=part3max-part3min

#part4min=815.50
#part4max=877.50
#diff4=part4max-part4min

#part5min=240.50 #230.70
#part5max=285.00
#diff5=part5max-part5min



#initialize paramaters
#graph=0

#ngen=100
#npart=200
#ngroup=20

#procs=4

#nparam=108
#minval=np.zeros(nparam)
#maxval=np.zeros(nparam)
#Ahmedhansal


#Bineloudine
#for iparam in range(36):
    #minval[iparam]=623.#min and maximum for the three curves, i.e. Flood curve, Target level
                       #Firm curve
    #maxval[iparam]=643.
#Hassan1er
#for iparam in range(37,72,1):
    #minval[iparam]=201.
    #maxval[iparam]=220.
#MoullayYousef
#for iparam in range(73,108,1):
    #minval[iparam]=45.
    #maxval[iparam]=106.89
#Almassira
    
#constrict=0.9
#inert1=1.5
#inert0=0.1
#c1=1.
#c2=0.5
#value power M$/GWh
#valpeakpower=65.85*1e-3
#valrestpower=31.59*1e-3
#FrPeakSag=1.
#FrPeakCir=1.
#FrPeakJat=0.21
#damage flooding (M$) if discharge Jatiluhur above floodstart (m3/s)
#floodstart=320.
#flooddamage=14.
#value irrigation water (M$/Mm3) and PWS demand (m3/s)
#agrival=0.02
#pwsdem=31.
#benscale=1.
#for elitist-mutation
#fracel=0.05
#pmut=0.9
#mutmag=0.2
#elparts=int(fracel*npart)

# define class to read his files

#if model_folder != :
   # tempfolder=tempfile.TemporaryDirectory()
   # modelfolder=tempfolder.name + ""
   # src="C:\\OER0T"
   # dst=model_folder
   # shutil.copytree(src,dst)
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
   
    #path=rd.randint(0,100)
    #if os.path.exists("C:\\OerTemp"+str(path)):
        #shutil.rmtree("C:\\OerTemp"+str(path))
    tempfolder=tempfile.TemporaryDirectory(dir="C:\\")
    Model_folder=tempfolder.name[0:3]+tempfolder.name[6:]#+"R"
    src="C:\\Oer0T"
    dst=Model_folder
    distutils.dir_util.copy_tree(src,dst)
   # os.rename(Model_folder,"C:\\OerTemp"+str(path))
    #Model_folder="C:\\OerTemp"+str(path)
    #TarFolder=Model_folder
    #distutils.dir_util.copy_tree(src,TarFolder)
    TarFolder=Model_folder+'\\OUMRBIA9.Rbn\\CMTWORK'
    print(TarFolder)
    files = glob.glob(os.path.join(TarFolder,'*.*'))
    for file in files:
        f=open(file,'rt')
        textin=f.readlines()
        f.close()
        textout=list()
        for line in textin:
            line=line.replace('\Oer0T',Model_folder[2:])
                #line=line.replace('OER0T','OER0T{}'.format(inst))
            textout.append(line)
            f=open(file,'wt')
            f.writelines(textout)
            f.close()
    TarFolder=Model_folder+'\\OUMRBIA9.Rbn\\3'
    files = glob.glob(os.path.join(TarFolder,'TIMESERI.*'))
    for file in files:
        f=open(file,'rt')
        textin=f.readlines()
        f.close()
        textout=list()
        for line in textin:
            line=line.replace('\Oer0T',Model_folder[2:])
            #line=line.replace('OER0T','OER0T{}'.format(inst))
            textout.append(line)
        f=open(file,'wt')
        f.writelines(textout)
        f.close()
    
    TarFolder=Model_folder+'\\OUMRBIA9.Rbn\\WORK'
    files = glob.glob(os.path.join(TarFolder,'TIMESERI.*'))
    for file in files:
        f=open(file,'rt')
        textin=f.readlines()
        f.close()
        textout=list()
        for line in textin:
            line=line.replace('\OER0T',Model_folder[2:])
            #line=line.replace('OER0T','OER0T{}'.format(inst))
            textout.append(line)
        f=open(file,'wt')
        f.writelines(textout)
        f.close()    
    TarFolder=Model_folder+'\\Programs\\RL_PSTPR.INI'
    f=open(TarFolder,'rt')
    textin=f.readlines()
    f.close()
    textout=list()
    for line in textin:
        line=line.replace('C:\\Oer0T', Model_folder)
        textout.append(line)
    f=open(TarFolder,'wt')
    f.writelines(textout)
    f.close()
    TarFolder=Model_folder+'\\Programs\\rib2his4.fnm'
    f=open(TarFolder,'rt')
    textin=f.readlines()
    f.close()
    textout=list()
    for line in textin:
        line=line.replace('C:\\Oer0T', Model_folder)
        textout.append(line)
    f=open(TarFolder,'wt')
    f.writelines(textout)
    f.close()
    time.sleep(3)
    def generate(a,maximum,k):
        b=(maximum-a)*k+a
        return b
    #ribdir='c:\\OER0T{}'.format(inst) 
    part1min=617.70

    part2min=741.55

    part3min=885.50

    part4min=815.50

    part5min=240.50 #230.70
    part1max=677.00

    part2max=810.00
 
    part3max=966.00

    part4max=877.50
 #230.70
    part5max=285.00 
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
    #b1=rd.uniform(part1min, a1);b2=rd.uniform(part1min, a2);b3=rd.uniform(part1min, a3);b4=rd.uniform(part1min, a4);b5=rd.uniform(part1min, a5);b6=rd.uniform(part1min, a6);b7=rd.uniform(part1min, a7);b8=rd.uniform(part1min, a8);b9=rd.uniform(part1min, a9);b10=rd.uniform(part1min, a10);b11=rd.uniform(part1min, a11);b12=rd.uniform(part1min, a12)
    #c1=rd.uniform(part1min, b1);c2=rd.uniform(part1min, b2);c3=rd.uniform(part1min, b3);c4=rd.uniform(part1min, b4);c5=rd.uniform(part1min, b5);c6=rd.uniform(part1min, b6);c7=rd.uniform(part1min, b7);c8=rd.uniform(part1min, b8);c9=rd.uniform(part1min, b9);c10=rd.uniform(part1min, b10);c11=rd.uniform(part1min, b11);c12=rd.uniform(part1min, b12)
    #e1=rd.uniform(part2min, d1);e2=rd.uniform(part2min, d2);e3=rd.uniform(part2min, d3);e4=rd.uniform(part2min, d4);e5=rd.uniform(part2min, d5);e6=rd.uniform(part2min, d6);e7=rd.uniform(part2min, d7);e8=rd.uniform(part2min, d8);e9=rd.uniform(part2min, d9);e10=rd.uniform(part2min, d10);e11=rd.uniform(part2min, d11);e12=rd.uniform(part2min, d12)
    #f1=rd.uniform(part2min, e1);f2=rd.uniform(part2min, e2);f3=rd.uniform(part2min, e3);f4=rd.uniform(part2min, e4);f5=rd.uniform(part2min, e5);f6=rd.uniform(part2min, e6);f7=rd.uniform(part2min, e7);f8=rd.uniform(part2min, e8);f9=rd.uniform(part2min, e9);f10=rd.uniform(part2min, e10);f11=rd.uniform(part2min, e11);f12=rd.uniform(part2min, e12)
    #h1=rd.uniform(part3min, g1);h2=rd.uniform(part3min, g2);h3=rd.uniform(part3min, g3);h4=rd.uniform(part3min, g4);h5=rd.uniform(part3min, g5);h6=rd.uniform(part3min, g6);h7=rd.uniform(part3min, g7);h8=rd.uniform(part3min, g8);h9=rd.uniform(part3min, g9);h10=rd.uniform(part3min, g10);h11=rd.uniform(part3min, g11);h12=rd.uniform(part3min, g12)
    #i1=rd.uniform(part3min, h1);i2=rd.uniform(part3min, h2);i3=rd.uniform(part3min, h3);i4=rd.uniform(part3min, h4);i5=rd.uniform(part3min, h5);i6=rd.uniform(part3min, h6);i7=rd.uniform(part3min, h7);i8=rd.uniform(part3min, h8);i9=rd.uniform(part3min, h9);i10=rd.uniform(part3min, h10);i11=rd.uniform(part3min, h11);i12=rd.uniform(part3min, h12)
    #k1=rd.uniform(part4min, j1);k2=rd.uniform(part4min, j2);k3=rd.uniform(part4min, j3);k4=rd.uniform(part4min, j4);k5=rd.uniform(part4min, j5);k6=rd.uniform(part4min, j6);k7=rd.uniform(part4min, j7);k8=rd.uniform(part4min, j8);k9=rd.uniform(part4min, j9);k10=rd.uniform(part4min, j10);k11=rd.uniform(part4min, j11);k12=rd.uniform(part4min, j12)
    #l1=rd.uniform(part4min, k1);l2=rd.uniform(part4min, k2);l3=rd.uniform(part4min, k3);l4=rd.uniform(part4min, k4);l5=rd.uniform(part4min, k5);l6=rd.uniform(part4min, k6);l7=rd.uniform(part4min, k7);l8=rd.uniform(part4min, k8);l9=rd.uniform(part4min, k9);l10=rd.uniform(part4min, k10);l11=rd.uniform(part4min, k11);l12=rd.uniform(part4min, k12)
    #n1=rd.uniform(part5min, m1);n2=rd.uniform(part5min, m2);n3=rd.uniform(part5min, m3);n4=rd.uniform(part5min, m4);n5=rd.uniform(part5min, m5);n6=rd.uniform(part5min, m6);n7=rd.uniform(part5min, m7);n8=rd.uniform(part5min, m8);n9=rd.uniform(part5min, m9);n10=rd.uniform(part5min, m10);n11=rd.uniform(part5min, m11);n12=rd.uniform(part5min, m12)
    #o1=rd.uniform(part5min, n1);o2=rd.uniform(part5min, n2);o3=rd.uniform(part5min, n3);o4=rd.uniform(part5min, n4);o5=rd.uniform(part5min, n5);o6=rd.uniform(part5min, n6);o7=rd.uniform(part5min, n7);o8=rd.uniform(part5min, n8);o9=rd.uniform(part5min, n9);o10=rd.uniform(part5min, n10);o11=rd.uniform(part5min, n11);o12=rd.uniform(part5min, n12)
    #ahmed=np.array([[a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12],
                    #[b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12],
                    #[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12]])
    #bine=np.array([[d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12],
                   #[e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12], 
                   #[f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12]])
    #hassan=np.array([[g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12],
                     #[h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12],
                     #[i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12]])
    #moullay=np.array([[j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12],
                      #[k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12],
                      #[l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12]])
    #almassira=np.array([[m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12],
                        #[n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12],
                        #[o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12]])
                        

    part1max=677.00

    part2max=810.00
 
    part3max=966.00

    part4max=877.50
 #230.70
    part5max=285.00                    
    
    #if checkconstraint(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,
                       #d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,
                       #g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,
                       #j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,
                       #m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,
                       #part1max,part2max,part3max,part4max,part5max):
    f = open(Model_folder+"\\OUMRBIA9.Rbn\\Actions\\Measures\\R001-2015-RibasimRsvNodeOperRulesDataOER.mes", 'r+')
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
    f = open(Model_folder+"\\OUMRBIA9.Rbn\\Actions\\Measures\\R001-2015-RibasimRsvNodeOperRulesDataOER.mes", 'wt')
    f.seek(0)
    f.writelines(a)

    f.close()
#break point

    os.chdir(Model_folder+"\\OUMRBIA9.Rbn\\CMTWORK")
    resrb = subprocess.call(Model_folder+"\\Programs\\Ribasim\\System\\Bin2prt.exe BIN2PRT.fnm") #Doubt
    if resrb <0 and resrb> 0:
        print('Error in RIBASIM:')
        print(resrb)
        exit()
    resrb = subprocess.call(Model_folder+"\\Programs\\Ribasim\\System\\Simproc.exe simproc.fnm") #Doubt
    if resrb <0 and resrb> 0:
        print('Error in RIBASIM:')
        print(resrb)
        exit()
    resrb = subprocess.call(Model_folder+"\\Programs\\Runlist\\RUNLIST.EXE " + Model_folder+"\\Programs\\RL_PSTPR.INI")
    if resrb <0 and resrb> 0:
        print('Error in postprocessing RIBASIM')
        print(resrb)
        exit()
    reshis = HisFile(Model_folder+"\\OUMRBIA9.Rbn\\WORK\\pws.his")
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
    reshis1 = HisFile(Model_folder+"\\OUMRBIA9.Rbn\\WORK\\reservoi.his")
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
    pws1=dataMarrakech
    pws1avg=sum(pws1[880:-1])/len(pws1[880:-1])
    pws2=dataEljadida
    pws2avg=sum(pws2[880:-1])/len(pws2[880:-1])
    pws3=dataazzemour
    pws3avg=sum(pws3[880:-1])/len(pws3[880:-1])
    pws4=dataruraldoukal
    pws4avg=sum(pws4[880:-1])/len(pws4[880:-1])
    pws21=datacas
    pws21avg=sum(pws21[880:-1])/len(pws21[880:-1])
    pws5=datasafi
    pws5avg=sum(pws5[880:-1])/len(pws5[880:-1])
    pws6=datagwtaval
    pws6avg=sum(pws6[880:-1])/len(pws6[880:-1])
    pws7=datamarrakechtransfe
    pws7avg=sum(pws7[880:-1])/len(pws7[880:-1])
    pws8=datacasablanca
    pws8avg=sum(pws8[880:-1])/len(pws8[880:-1])
    pws9=databenguerir
    pws9avg=sum(pws9[880:-1])/len(pws9[880:-1])
    pws10=databejaad
    pws10avg=sum(pws10[880:-1])/len(pws10[880:-1])
    pws11=dataazillal
    pws11avg=sum(pws11[880:-1])/len(pws11[880:-1])
    pws12=datakhenifra
    pws12avg=sum(pws12[880:-1])/len(pws12[880:-1])
    pws13=databeniamir
    pws13avg=sum(pws13[880:-1])/len(pws13[880:-1])
    pws14=databenimoussa
    pws14avg=sum(pws14[880:-1])/len(pws14[880:-1])
    pws15=datagwtadla
    pws15avg=sum(pws15[880:-1])/len(pws15[880:-1])
    pws16=datagwchaouiacotiere
    pws16avg=sum(pws16[880:-1])/len(pws16[880:-1])
    pws17=databenimellal
    pws17avg=sum(pws17[880:-1])/len(pws17[880:-1])
    pws18=datagwbahira
    pws18avg=sum(pws18[880:-1])/len(pws18[880:-1])
    pws19=datagwdoukkalasahel
    pws19avg=sum(pws19[880:-1])/len(pws19[880:-1])
    pws20=datagwdir
    pws20avg=sum(pws20[880:-1])/len(pws20[880:-1])
    pwsshortageavg=(pws1avg+pws2avg+pws3avg+pws4avg+pws5avg+pws6avg+pws7avg+pws8avg+pws9avg+pws10avg+pws11avg+pws12avg+pws13avg+pws14avg+pws15avg+pws16avg+pws17avg+pws18avg+pws19avg+pws20avg+pws21avg)/21
    pwssum=pws1+pws2+pws3+pws4+pws5+pws6+pws7+pws8+pws9+pws10+pws11+pws12+pws13+pws14+pws15+pws16+pws17+pws18+pws19+pws20+pws21
    pwssorted=sorted(pwssum)
    pwsobjective=sum(pwssorted[880:-1])
    energy1=datatagzirt
    #energy1=sum(energy1[880:-1])/len(energy1[880:-1])
    energy2=datahassan
    #energy2=sum(energy2[880:-1])/len(energy2[880:-1])
    energy3=datatyoughza
    #energy3=sum(energy3[880:-1])/len(energy3[880:-1])
    energy4=datamoulayyoussef
    #energy4=sum(energy4[880:-1])/len(energy4[880:-1])
    energy5=dataalmassira
    #energy5=sum(energy5[880:-1])/len(energy5[880:-1])
    energy6=databinelouidane
    #energy6=sum(energy6[880:-1])/len(energy6[880:-1])
    energy7=dataahmedhansal
    #energy7=sum(energy7[880:-1])/len(energy7[880:-1])
    energy8=datasididriss
    #energy8=sum(energy8[880:-1])/len(energy8[880:-1])
    energytotal=(energy1+energy2+energy3+energy4+energy5+energy6+energy7+energy8)
    energyobjective=sum(energytotal[880:-1])
    reshis = HisFile(Model_folder+"\\OUMRBIA9.Rbn\\WORK\\vir.his")
    reshis.read()
    segnam="b'Vir_PMHTessaoutAval:'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataPMHT = reshis.gettimeseries(sysnam, segnam)
    segnam="b'Vir_OeRavalMassira:I'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataOeRaval = reshis.gettimeseries(sysnam, segnam)
    segnam="b'Vir_MoyenOeR:I4     '"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataMoyen = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_Haouz:I7        '"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataHaouz= reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_GH-TessaoutAvald'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHTessaout = reshis.gettimeseries(sysnam, segnam)
    segnam="b'Vir_PMH:AmontTadla:I'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataPMHAmont= reshis.gettimeseries(sysnam, segnam)
    segnam="b'Vir_PMH:I11         '"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataPMHI11 = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_GH_Haut_Doukkala'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHHaut = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_GH_BeniAmir:I1  '"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataBeni = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_GHBeniMoussa:I5A'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataMoussa = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_GH-TessaoutAvaln'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHTessain = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_GH-TessaoutAmont'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHTessaamont = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_GH_BasDoukkala:I'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataGHBasDoukkala = reshis.gettimeseries(sysnam, segnam)
    segnam= "b'Vir_Shtouka         '"
    sysnam="b'Shortage (m3/s)     '"
    rdate,datashtouka = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_IPAmontTadlaIp24'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataIp24 = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_IPMoyenOeRIp4   '"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataIp4moyen= reshis.gettimeseries(sysnam, segnam)
    segnam=  'b"Vir_Droitd\'eauLakhda"'
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataeau = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_IPLakhdaravalSid'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataiplakhdasid = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_IPTadlaOeRIp3   '"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataiptadlaoer = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_PMH_TessaoutInte'"
    sysnam="b'Shortage (m3/s)     '"
    rdate,datapmhtessa = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_Dir:I2          '"
    sysnam="b'Shortage (m3/s)     '"
    rdate,dataDir = reshis.gettimeseries(sysnam, segnam)
    segnam=  "b'Vir_PMH:OERTadla:I3 '"
    sysnam="b'Shortage (m3/s)     '"
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
    os.chdir("C:\\Users\\keshav\\OneDrive - Stichting Deltares\\Desktop\\vanMarnix")
    #shutil.rmtree("C:\\OerTemp"+str(path))
    tempfolder.cleanup()
    #shutil.rmtree("C:\\OerTemp"+str(path))
    #time.sleep(3)
    return pwsobjective, energyobjective,virobjective

   # lowhis = HisFile("C:\\Ribasim7\\OUMRBIA9.Rbn\\WORK\\lowflow.his")
    #lowhis.read()
    #sysnam = "b'Realised flow(m3/s) '"
    #segnam = "b'Lfl_Irrigation      '"
    #rdate,irrisup = lowhis.gettimeseries(sysnam, segnam)
   # part.HydroPSag = 0.
   # part.HydroPCir = 0.
   ##part.irrival=0.
    #part.flooddam=0.
    #Mm3fac=86400*30.5*1e-6
    #for lstep in range(lowhis.nstep):
        #imonth=rdate[lstep].month-1
        #print(lstep,rdata[lstep],imonth,peak[imonth],rdate[lstep],min(rdate[lstep],peak[imonth]))
        #part.HydroPSag += rdataSag[lstep+1]*(FrPeakSag*valpeakpower+(1.-FrPeakSag)*valrestpower)
        #part.HydroPCir += rdataCir[lstep+1]*(FrPeakCir*valpeakpower+(1.-FrPeakCir)*valrestpower)
        #part.HydroPJat += rdataJat[lstep+1]*(FrPeakJat*valpeakpower+(1.-FrPeakJat)*valrestpower)
        #part.irrival += max(0., irrisup[lstep]-pwsdem)*agrival*Mm3fac
        #if lstep < lowhis.nstep -1:
            #if (rflowJat[lstep+1] >= floodstart) and (rflowJat[lstep+2] < floodstart) : part.flooddam += flooddamage
        #else:
            #if (rflowJat[lstep+1] >= floodstart)  : part.flooddam += flooddamage

    #part.value = (part.HydroPSag+part.HydroPCir+part.HydroPJat+part.irrival-part.flooddam)*benscale
    #if part.value > part.bestvalue:
        #part.bestvalue=part.value
        #for iparam in range(nparam): part.bestparam[iparam]=part.param[iparam]
    #output.put((ipart,part.value,part.HydroPSag,part.HydroPCir,part.HydroPJat,part.irrival,part.flooddam))
from ema_workbench import (Model, RealParameter, ScalarOutcome,
                           MultiprocessingEvaluator, ema_logging,
                           Constant, SequentialEvaluator)


if __name__ == '__main__':
    ema_logging.log_to_stderr(ema_logging.INFO)

#instantiate the model
    basin_model = Model('ribasim', function=ribasimmodel)
#basin_model.time_horizon = 1 # used to specify the number of t/imesteps

#specify uncertainties
#for i in range(180):
   # basin_model.uncertainties = [RealParameter('part1', part1min, part1max),
                               #  RealParameter('part2', part2min, part2max),
                                # RealParameter('part3', part3min, part3max),
                                # RealParameter('part4', part4min, part4max),
                               # RealParameter('part5', part5min, part5max)]
    #ahmed=np.zeros((3,12))
    #bine=np.zeros((3,12))
    #hassan=np.zeros((3,12))
    #moullay=np.zeros((3,12))
    #almassira=np.zeros((3,12))
    ahmed=np.array([ ['c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12'],
                    ['fracb1','fracb2','fracb3','fracb4','fracb5','fracb6','fracb7','fracb8','fracb9','fracb10','fracb11','fracb12'],
                    ['fracc1','fracc2','fracc3','fracc4','fracc5','fracc6','fracc7','fracc8','fracc9','fracc10','fracc11','fracc12']])
    bine=np.array([['f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f11','f12'],
                   ['frace1','frace2','frace3','frace4','frace5','frace6','frace7','frace8','frace9','frace10','frace11','frace12'], 
                   ['fracf1','fracf2','fracf3','fracf4','fracf5','fracf6','fracf7','fracf8','fracf9','fracf10','fracf11','fracf12']])
    hassan=np.array([['i1','i2','i3','i4','i5','i6','i7','i8','i9','i10','i11','i12'],
                     ['frach1','frach2','frach3','frach4','frach5','frach6','frach7','frach8','frach9','frach10','frach11','frach12'],
                     ['fraci1','fraci2','fraci3','fraci4','fraci5','fraci6','fraci7','fraci8','fraci9','fraci10','fraci11','fraci12']])
    moullay=np.array([['l1','l2','l3','l4','l5','l6','l7','l8','l9','l10','l11','l12'],
                      ['frack1','frack2','frack3','frack4','frack5','frack6','frack7','frack8','frack9','frack10','frack11','frack12'],
                      ['fracl1','fracl2','fracl3','fracl4','fracl5','fracl6','fracl7','fracl8','fracl9','fracl10','fracl11','fracl12']])
    almassira=np.array([['o1','o2','o3','o4','o5','o6','o7','o8','o9','o10','o11','o12'],
                        ['fracn1','fracn2','fracn3','fracn4','fracn5','fracn6','fracn7','fracn8','fracn9','fracn10','fracn11','fracn12'],
                        ['fraco1','fraco2','fraco3','fraco4','fraco5','fraco6','fraco7','fraco8','fraco9','fraco10','fraco11','fraco12']])
    part1min=617.70
    part1max=677.00
    diff1=part1max-part1min

    part2min=741.55
    part2max=810.00
    diff2=part2max-part2min

    part3min=885.50
    part3max=966.00
    diff3=part3max-part3min

    part4min=815.50
    part4max=877.50
    diff4=part4max-part4min
    
    part5min=240.50 #230.70
    part5max=285.00
    diff5=part5max-part5min
                           

    #ahmed=np.array([['a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12'],
                    #['b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','b11','b12'],
                    #['c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12']])
    #bine=np.array([['d1','d2','d3','d4','d5','d6','d7','d8','d9','d10','d11','d12'],
                   #['e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12'], 
                   #['f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f11','f12']])
    #hassan=np.array([['g1','g2','g3','g4','g5','g6','g7','g8','g9','g10','g11','g12'],
                     #['h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12'],
                     #['i1','i2','i3','i4','i5','i6','i7','i8','i9','i10','i11','i12']])
    #moullay=np.array([['j1','j2','j3','j4','j5','j6','j7','j8','j9','j10','j11','j12'],
                      #['k1','k2','k3','k4','k5','k6','k7','k8','k9','k10','k11','k12'],
                      #['l1','l2','l3','l4','l5','l6','l7','l8','l9','l10','l11','l12']])
    #almassira=np.array([['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12'],
                        #['n1','n2','n3','n4','n5','n6','n7','n8','n9','n10','n11','n12'],
                        #['o1','o2','o3','o4','o5','o6','o7','o8','o9','o10','o11','o12']])
    dams=[ahmed,bine,hassan,moullay,almassira]                    
    maximum=[part1max,part2max,part3max,part4max,part5max]
    minimum=[part1min,part2min,part3min,part4min,part5min]
    diff=[diff1,diff2,diff3,diff4,diff5]
    #part=[part1,part2,part3,part4,part5,part6,part7,part8,part9,part10,part11,part12,part13,part14,part15]
    #minimum=[part1min,part1min,part1min,part2min,part2min,part2min,part3min,part3min,part3min,part4min,part4min,part4min,part5min,part5min,part5min]
    #maximum=[part1max,part1max,part1max,part2max,part2max,part2max,part3max,part3max,part3max,part4max,part4max,part4max,part5max,part5max,part5max]
    l=[]

#set levers, one for each time step
    #n=1
#for i in range(basin_model.time_horizon):
    for k in range(len(dams)):
            for j in range(12):
                a=RealParameter(dams[k][0][j],minimum[k],maximum[k])
                l.append(a)
                a=RealParameter(dams[k][1][j],0,1)
                l.append(a)
                a=RealParameter(dams[k][2][j],0,1)
                l.append(a)
            #basin_model.levers.append(a)
            #l.append(a)
    print(l)
    basin_model.levers=l
        
            #basin_model.levers = [RealParameter(part1[j], part1min, part1max),
                                  #RealParameter(part2[j], part1min, part1max),
                                  #RealParameter(part3[j], part1min, part1max),
                                  #RealParameter(pat4[j], part2min, part2max),
                                  #RealParameter(part5[j], part2min, part2max),
                                  #RealParameter(part6[j], part2min, part2max),
                                  #RealParameter(part7[j], part3min, part3max),
                                  #RealParameter(part8[j], part3min, part3max),
                                  #RealParameter(part9[j], part3min, part3max),
                                  #RealParameter(part10[j], part4min, part4max),
                                  #RealParameter(part11[j], part4min, part4max),
                                  #RealParameter(part12[j], part4min, part4max),
                                  #RealParameter(part13[j], part5min, part5max),
                                  #RealParameter(part14[j], part5min, part5max),
                                  #RealParameter(part15[j], part5min, part5max)]# we use time_horizon here


    #constraints=[]
                               
 
    #a=[]
    #k=0
    #l=0
    #m=[]
    #n=1
    #for k in range(len(dams)): 
        #for i in range(12):
            #a=Constraint("Constraint_"+str(n),parameter_names=[dams[k][0][i],dams[k][1][i]],function=lambda x,y:max(0,-x+y))
            #m.append(a)
            #n=n+1
            #a=Constraint("Constraint_"+str(n),parameter_names=[dams[k][1][i],dams[k][2][i]],function=lambda x,y:max(0,-x+y))
            #m.append(a)
            #n=n+1
     
        
    basin_model.outcomes = [ScalarOutcome('pwsobjective', kind=ScalarOutcome.MINIMIZE),
                            ScalarOutcome('energyobjective', kind=ScalarOutcome.MAXIMIZE),
                            ScalarOutcome('virobjective', kind=ScalarOutcome.MINIMIZE)]


    convergence_metrics = [EpsilonProgress()]

                       
    #with SequentialEvaluator(basin_model) as evaluator:
        #experiments, outcomes=evaluator.perform_experiments(policies=10)#constraints=constraints)
    with MultiprocessingEvaluator(basin_model) as evaluator:
        results, convergence=evaluator.optimize(nfe=100, searchover='levers',epsilons=[0.20, 0.20,0.20], reference=None, convergence=convergence_metrics)#constraints=m)
        #experiments, outcomes=evaluator.perform_experiments(scenarios=100)

