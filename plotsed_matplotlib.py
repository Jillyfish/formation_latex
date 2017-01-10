import pyfits,string
import PlotLibrary as PlotLibrary
import Loggin
from math import *
import numpy as np
from array import *
import matplotlib.pyplot as plt
import ReadCatalog as ReadCat

HztoMeV=4.1357e-21

### source names for the 2FHL and 3FGL catalogs
### watch out for the spaces
source2FHL = "2FHL J0639.9-1252"
source3FGL = "3FGL J0640.0-1252"
source = "TXS0637-128"

Emin = 50e3 # 50 GeV
Emax = 2e6  # 2 TeV

do_2fhl = True
do_3fgl = True
do_source = True
do_ref = True

if (do_2fhl) :
	Cat2FHL = ReadCat.CatalogReader(source2FHL)
	params  = Cat2FHL.ReadPL2("2FHL")
	Spec    = PlotLibrary.Spectrum(params,"PowerLaw2",Emin,Emax)
	ener, phi, dphi = Spec.MakeSEDAndError(Emin,Emax)
	dphi_down = np.log10(phi) - np.log10(phi-dphi)
	dphi_up   = np.log10(phi+dphi) - np.log10(phi)
if (do_3fgl) :
	Cat3FGL = ReadCat.CatalogReader(source3FGL)
	model   = Cat3FGL.GetModels()
	if (str(model) == "PowerLaw") :
		params = Cat3FGL.ReadPL("3FGL")
	if (str(model) == "LogParabola") :
		params = Cat3FGL.ReadLP("3FGL")
	Spec = PlotLibrary.Spectrum(params,str(model),100,300e3)
	energ, butt = Spec.GetButterfly(100,300e3)				# 3FGL butterfly
	em,ep,flux,dflux =  Cat3FGL.GetDataPoints('3FGL')		# 3FGL data points
	ener = np.sqrt(em*ep)
	dem = ener-em
	dep = ep-ener
	c = Cat3FGL.ReadPL('3FGL')[3]
	dnde = (-c+1)*flux*np.power(ener,-c+2)/(np.power((ep),-c+1)-np.power((em),-c+1))*1.6e-6
	ddnde = dnde*dflux/flux

if (do_source) :
	# SED Builder data of the source
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
if (do_ref) :
	data_ref = np.genfromtxt("1ES0229+200.txt", unpack=True)
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

plt.errorbar(logE_ref,sed_ref,yerr=dsed_ref,color='grey',mec='grey',marker='o',capsize=0,ls='None',label='1ES 0229+200')
plt.errorbar(logE,sed,yerr=dsed,color='black',marker='o',capsize=0,ls='None',label='TXS 0637-128')

x = logE_ul
y = sed_ul
yerrl = y*0 + 0.5
yerru = y*0
uplims = np.zeros(x.shape)
i = np.arange(len(x))
uplims[i] = True
plt.errorbar(x,y,yerr=[yerrl,yerru],lolims=uplims,marker='o',ls='None',color='black')

plt.ylabel('SED' )
plt.xlabel('log(Energy)')
plt.show()









