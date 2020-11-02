# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 15:13:15 2020

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
import os
import glob
import distutils.dir_util
import sys
import shutil
from backports import tempfile
import operator
import time
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
    
reshis = HisFile("C:\\Oer0T_2\\OUMRBIA9.Rbn\\WORK\\pws.his")
reshis.read()
segnam = "b'Pws_Marrakech:P2(CR)'"
sysnam = "b'Gross demand incl.lo'"
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

reshis1 = HisFile("C:\\Oer0T_2\\OUMRBIA9.Rbn\\WORK\\reservoi.his")
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

energy1=datatagzirt
energy2=datahassan
energy3=datatyoughza
energy4=datamoulayyoussef
energy5=dataalmassira
energy6=databinelouidane
energy7=dataahmedhansal
energy8=datasididriss
energytotal=(energy1+energy2+energy3+energy4+energy5+energy6+energy7+energy8)
energyobjective=sum(energytotal[:])
    
reshis = HisFile("C:\\Oer0T_2\\OUMRBIA9.Rbn\\WORK\\vir.his")
reshis.read()
segnam="b'Vir_PMHTessaoutAval:'"
sysnam="b'Demand from network '"
rdate,dataPMHT = reshis.gettimeseries(sysnam, segnam)
segnam="b'Vir_OeRavalMassira:I'"
rdate,dataOeRaval = reshis.gettimeseries(sysnam, segnam)
segnam="b'Vir_MoyenOeR:I4     '"
rdate,dataMoyen = reshis.gettimeseries(sysnam, segnam)
segnam= "b'Vir_Haouz:I7        '"
rdate,dataHaouz= reshis.gettimeseries(sysnam, segnam)
segnam= "b'Vir_GH-TessaoutAvald'"
rdate,dataGHTessaout = reshis.gettimeseries(sysnam, segnam)
segnam="b'Vir_PMH:AmontTadla:I'"
rdate,dataPMHAmont= reshis.gettimeseries(sysnam, segnam)
segnam="b'Vir_PMH:I11         '"
rdate,dataPMHI11 = reshis.gettimeseries(sysnam, segnam)
segnam= "b'Vir_GH_Haut_Doukkala'"
rdate,dataGHHaut = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_GH_BeniAmir:I1  '"
rdate,dataBeni = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_GHBeniMoussa:I5A'"
rdate,dataMoussa = reshis.gettimeseries(sysnam, segnam)
segnam= "b'Vir_GH-TessaoutAvaln'"
rdate,dataGHTessain = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_GH-TessaoutAmont'"
rdate,dataGHTessaamont = reshis.gettimeseries(sysnam, segnam)
segnam= "b'Vir_GH_BasDoukkala:I'"
rdate,dataGHBasDoukkala = reshis.gettimeseries(sysnam, segnam)
segnam= "b'Vir_Shtouka         '"
rdate,datashtouka = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_IPAmontTadlaIp24'"
rdate,dataIp24 = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_IPMoyenOeRIp4   '"
rdate,dataIp4moyen= reshis.gettimeseries(sysnam, segnam)
segnam=  'b"Vir_Droitd\'eauLakhda"'
rdate,dataeau = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_IPLakhdaravalSid'"
rdate,dataiplakhdasid = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_IPTadlaOeRIp3   '"
rdate,dataiptadlaoer = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_PMH_TessaoutInte'"
rdate,datapmhtessa = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_Dir:I2          '"
rdate,dataDir = reshis.gettimeseries(sysnam, segnam)
segnam=  "b'Vir_PMH:OERTadla:I3 '"
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

reshis = HisFile("C:\\Oer0T_2\\OUMRBIA9.Rbn\\WORK\\fir.his")
reshis.read()
sysnam=   "b'Demand (m3/s)       '"
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
print(energyobjective)
print(pwsobjective)
print(virobjective)