import pyfits,string
import PlotLibrary as PlotLibrary
import Loggin
from math import *
import numpy as np
from array import *
import ROOT
import string
from Macro import *
import ReadCatalog as ReadCat
import sys as sys

HztoMeV=4.1357e-21

source2FHL = "2FHL J0639.9-1252"
source3FGL = "3FGL J0640.0-1252"
source = "TXS0637-128"

Emin = 50e3 # 50 GeV
Emax = 2e6  # 2 TeV

Cat2FHL = ReadCat.CatalogReader(source2FHL)
params = Cat2FHL.ReadPL2("2FHL")
Spec = PlotLibrary.Spectrum(params,"PowerLaw2",Emin,Emax)
ener, phi, dphi = Spec.MakeSEDAndError(Emin,Emax)
dphi_down = np.log10(phi) - np.log10(phi-dphi)
dphi_up = np.log10(phi+dphi) - np.log10(phi)

Cat3FGL = ReadCat.CatalogReader(source3FGL)
model = Cat3FGL.GetModels()
if (str(model) == "PowerLaw") :
	params = Cat3FGL.ReadPL("3FGL")
if (str(model) == "LogParabola") :
	params = Cat3FGL.ReadLP("3FGL")
Spec = PlotLibrary.Spectrum(params,str(model),100,300e3)
energ, butt = Spec.GetButterfly(100,300e3)
em,ep,flux,dflux =  Cat3FGL.GetDataPoints('3FGL')
ener = np.sqrt(em*ep)
dem = ener-em
dep = ep-ener
c=Cat3FGL.ReadPL('3FGL')[3]
dnde = (-c+1)*flux*np.power(ener,-c+2)/(np.power((ep),-c+1)-np.power((em),-c+1))*1.6e-6
ddnde = dnde*dflux/flux

for i in range(len(ener)) :
  print ener[i], dnde[i], ddnde[i]

data_SED = np.genfromtxt(str(source)+".txt", unpack=True)
logE = []
sed = []
dsed = []
logE_ul = []
sed_ul = []
for k in range(len(data_SED[0])) :
	if (data_SED[3][k] != 0) :
		logE.append(np.log10(10**data_SED[0][k]*HztoMeV))
		sed.append(data_SED[2][k])
		dsed.append(data_SED[3][k])
	else :
		logE_ul.append(np.log10(10**data_SED[0][k]*HztoMeV))
		sed_ul.append(data_SED[2][k])
if (len(sed_ul) == 1) :
	logE_ul = np.array([logE_ul,logE_ul])
	sed_ul = np.array([sed_ul,sed_ul])
logE = np.array(logE)
sed = np.array(sed)
dsed = np.array(dsed)
logE_ul = np.array(logE_ul)
sed_ul = np.array(sed_ul)

data_ref = np.genfromtxt("ApLibrae.txt", unpack=True)
logE_ref = []
sed_ref = []
dsed_ref = []
logE_reful = []
sed_reful = []
for k in range(len(data_ref[0])) :
	if (data_ref[3][k] != 0) :
		logE_ref.append(np.log10(10**data_ref[0][k]*HztoMeV))
		sed_ref.append(data_ref[2][k])
		dsed_ref.append(data_ref[3][k])
	else :
		logE_reful.append(np.log10(10**data_ref[0][k]*HztoMeV))
		sed_reful.append(data_ref[2][k])
logE_ref = np.array(logE_ref)
sed_ref = np.array(sed_ref)
dsed_ref = np.array(dsed_ref)
logE_reful = np.array(logE_reful)
sed_reful = np.array(sed_reful)








