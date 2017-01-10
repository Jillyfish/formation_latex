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
data_2FHL = np.genfromtxt("sourceslist_LBL_2FHLx3FGL_workon.dat", dtype='string', unpack=True)
# data_2FHL = np.genfromtxt("sourceslist_X-HBL_2FHLx3FGL_workon.dat", dtype='string', unpack=True)

source2FHL = "2FHL J2158.8-3013"
source3FGL = "3FGL J2158.8-3013"
source = "PKS2155-304"

Emin = 50e3 # 50 GeV
Emax = 2e6  # 2 TeV

Cat2FHL = ReadCat.CatalogReader(source2FHL)
params = Cat2FHL.ReadPL2("2FHL")
Spec = PlotLibrary.Spectrum(params,"PowerLaw2",Emin,Emax)
ener, phi, dphi = Spec.MakeSEDAndError(Emin,Emax)
dphi_down = np.log10(phi) - np.log10(phi-dphi)
dphi_up = np.log10(phi+dphi) - np.log10(phi)
tgr_2FHL = MakeTGraph(np.log10(ener),np.log10(phi),9,23)

Cat3FGL = ReadCat.CatalogReader(source3FGL)
model = Cat3FGL.GetModels()
if (str(model) == "PowerLaw") :
	params = Cat3FGL.ReadPL("3FGL")
if (str(model) == "LogParabola") :
	params = Cat3FGL.ReadLP("3FGL")
Spec = PlotLibrary.Spectrum(params,str(model),100,300e3)
energ, butt = Spec.GetButterfly(100,300e3)
tgr_3FGL = MakeTGraph(np.log10(energ),np.log10(butt),6)
em,ep,flux,dflux =  Cat3FGL.GetDataPoints('3FGL')
ener = np.sqrt(em*ep)
dem = ener-em
dep = ep-ener
c=Cat3FGL.ReadPL('3FGL')[3]
dnde = (-c+1)*flux*np.power(ener,-c+2)/(np.power((ep),-c+1)-np.power((em),-c+1))*1.6e-6
ddnde = dnde*dflux/flux
tgr_3FGL_p = MakeTGraphAsymmErrors(np.log10(ener),np.log10(dnde),array('f',len(em)*[0.]),array('f',len(em)*[0.]),np.log10(ddnde),np.log10(ddnde),6,8,0.7)

for i in range(len(ener)) :
  print ener[i], dnde[i], ddnde[i]

# data_analysis = np.genfromtxt("Fermianalysis/SED_"+str(source[0:7])+"_LogParabola.dat", unpack=True)  
# params = np.array[]
# Spec2 = PlotLibrary.Spectrum(params,Model="LogParabola",Emin=100,Emax=300e3)
# energy,butterfly = Spec2.GetButterfly(SED=True)
# tgr_dataanalysis = MakeTGraph(np.log10(energy),np.log10(butterfly),97)

data_SED = np.genfromtxt("SEDbuilder_data/"+str(source)+".txt", unpack=True)
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

tgr_sedsource = MakeTGraphErrors(logE,sed,logE*0,dsed)
tgr_sedul = MakeTGraph(logE_ul,sed_ul)

data_HESS = np.genfromtxt("HESS_analysisdata/PKS0537_flux_ul.dat", unpack = True)
energy = data_HESS[0]*1e6 
ergenergy = energy/624151.
flux = data_HESS[1]*1e-4*ergenergy
tgrUL = MakeTGraph(np.log10(energy),np.log10(flux),8,23)

data_ref = np.genfromtxt("SEDbuilder_data/ApLibrae.txt", unpack=True)
points = np.genfromtxt("SEDbuilder_data/ApLib_SED_2014.dat", unpack=True)
dlogE_down = np.log10(points[0]*1e-6) - np.log10(points[0]*1e-6-points[1]*1e-6)
dlogE_up = np.log10(points[0]*1e-6+points[1]*1e-6) - np.log10(points[0]*1e-6)
dsed_down = np.log10(points[2]) - np.log10(points[2]-points[3])
dsed_up = np.log10(points[2]+points[3]) - np.log10(points[02])
tgrpoints = MakeTGraphAsymmErrors(np.log10(points[0]*1e-6),np.log10(points[2]),dlogE_down,dlogE_up,dsed_down,dsed_up,96)


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
tgr_ref = MakeTGraphErrors(logE_ref,sed_ref,logE_ref*0,dsed_ref,color=15)
tgr_reful = MakeTGraph(logE_reful,sed_reful,color=15)


canvas = ROOT.TCanvas("canvas")
fenetre = ROOT.TH2F("",str(source)+", z = 0.116",10,-14,8,10,-17,-8)
fenetre.SetStats(000)
fenetre.SetXTitle("Energy (MeV)")
fenetre.SetYTitle("#nu F_{#nu} (ergs/cm2/s)")
fenetre.Draw()

# tgr_ref.Draw("Pz")
# tgr_reful.Draw("Pz")
# tgrpoints.Draw("Pz")
tgr_sedsource.Draw("Pz")
tgr_sedul.Draw("Pz")
# tgr_3FGL.Draw("L")
# tgr_2FHL.Draw("L")
# tgrUL.Draw("P")
# tgr_dataanalysis.Draw("L")








