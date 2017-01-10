# author David Sanchez david.sanchez@moi-hd.mpg.de
# Library to plot spectrum and SED. Also can correct for EBL


# Import
from math import *
import numpy as np
from array import *
#import ROOT
import string
#from Macro import *
class Spectrum:
  def __init__(self,Parameters,Model="",Emin=100,Emax=3e5,Npt=1000):

    self.Nbp = Npt
    self.Emin  = Emin
    self.Emax  = Emax
    self.Model = Model
    Model_list = ["PowerLaw","PowerLaw2","LogParabola","PLExpCutoff","file"]

    #Check the input model
    if not(self.Model in Model_list):
      self.usage()
      raise RuntimeError("Model in not in the list")
    #Check the number of parameters
    Parameters = self.CheckNumParemeters(Parameters)

    self.Name     = Parameters[0]
    if not(self.Model=="file") :
      #input units : cm-2 s-1 MeV-1 for diff_flux and 
      #              cm-2 s-1 for integral flux
      self.Norm     = Parameters[1] 
      # Index
      self.Gamma    = Parameters[3]

      # Errors
      self.ErrNorm    = Parameters[2]#input units :  cm-2 s-1 MeV-1
      self.ErrGamma = Parameters[4]

      if self.Model=="PowerLaw"  :
        self.Ed       = Parameters[5]  #Decorrelation Energy in MeV 


      if self.Model=="LogParabola" :
        self.Beta  = Parameters[5] #input units :  cm-2 s-1 MeV-1
        self.ErrBeta = Parameters[6]
        self.Ed  = Parameters[7] #Scale energy  in MeV
        self.Cov1  = Parameters[8] #input units :  cm-2 s-1 MeV-1
        self.Cov2 = Parameters[9]
        self.Cov3  = Parameters[10]

      if self.Model=="PLExpCutoff" :
        self.Ecut       = Parameters[5]  #Cutoff energy Energy in MeV
        self.ErrEcut    = Parameters[6]  #Err Cutoff energy Energy in MeV
        self.Ed       = Parameters[7]  #Scale energy Energy in MeV

    else :
      self.filename = Parameters[1]

  def CheckNumParemeters(self, Params):
    if  self.Model=="PowerLaw" and not(len(Params)==6) :
        self.usage()
        raise RuntimeError("Error : Not good numbers of input parameters\n Got "+str(len(Params))+" Need 6")

    if  self.Model=="PowerLaw2" and not(len(Params)==5) :
        self.usage()
        raise RuntimeError("Error : Not good numbers of input parameters\n Got "+str(len(Params))+" Need 5")


    if  self.Model=="file"  and not(len(Params)==2) :
        self.usage()
        raise RuntimeError("Error : Not good numbers of input parameters\n Got "+str(len(Params))+" Need 2")

    if  self.Model=="LogParabola"  and not(len(Params)==11) :
      if  not(len(Params)==8) :
        self.usage()
        raise RuntimeError("Error : Not good numbers of input parameters\n Got "+str(len(Params))+" Need 11")
      else :
        import logging
        logging.warn("Add 3 zeros for the covariance terms")
        Params.append(0)
        Params.append(0)
        Params.append(0)
    if  self.Model=="PLExpCutoff"  and not(len(Params)==8) :
        self.usage()
        raise RuntimeError("Error : Not good numbers of input parameters\n Got "+str(len(Params))+" Need 8")

    return Params

  def usage(self):
    print "not implemented at the moment"

  def GetValueAndError(self,Emin=0,Emax=10e10,SED=True):
    Emin= max(Emin,self.Emin)
    Emax= min(Emax,self.Emax)
    if  self.Model=="file" :
     # print 'Reading file for source', self.Name
      return self.MakeFluxAndError(Emin,Emax)
    #Commpute either the SED or the differential flux
    if SED :
     # print 'Compute the SED for source', self.Name
      return self.MakeSEDAndError(Emin,Emax)
    else :
     # print 'Compute the Differential butterfly for source', self.Name
      return self.MakeFluxAndError(Emin,Emax)

  def ReadFile(self,Emin=0,Emax=10e10) :
    '''Reads files as produced by enrico'''
    Emin= max(Emin,self.Emin)
    Emax= min(Emax,self.Emax)

    lines = open(self.filename,"r").readlines()
    ilen=len(lines)-1

    ene=[]
    Phi=[]
    dPhi=[]

    for j in xrange(ilen):
      words = string.split(lines[j+1])
      if float(words[0]) >= Emin and float(words[0]) <= Emax :
        ene.append(float(words[0]))
        Phi.append(float(words[1]))
        dPhi.append(float(words[2]))

    return np.array(ene),np.array(Phi),np.array(dPhi)

  def MakeFluxAndError(self,Emin=0,Emax=10e10):
    "return a np array with SED and error (in erg*cm-2*s-1)"

    mylog  = lambda x: np.log(x)

    #Compute the energy array using either the parameters of the class 
    # or the input parameters given by the user when calling this function
    Emin= max(Emin,self.Emin)
    Emax= min(Emax,self.Emax)
    ene = np.logspace(log10(Emin),log10(Emax),self.Nbp)

    Phi = array('f',(self.Nbp)*[0])
    dPhi = array('f',(self.Nbp)*[0])

    if self.Model == 'file':
      ene,Phi,dPhi = self.ReadFile(Emin,Emax) 

    if self.Model == "PowerLaw":
      Phi  = self.Norm*np.power(ene/self.Ed,-self.Gamma)
      dPhi = np.sqrt( (self.ErrNorm/self.Norm)**2 + (mylog(ene/self.Ed)*self.ErrGamma)**2 )*Phi

    if self.Model == "PowerLaw2":
      D = (np.power(self.Emax,1-self.Gamma)-np.power(self.Emin,1-self.Gamma))
      Phi  = self.Norm*(1-self.Gamma)*np.power(ene,-self.Gamma)/D
     # dPhi = np.sqrt( (self.ErrNorm/self.Norm)**2 + ((1./(1-self.Gamma)-mylog(ene/self.Emin))*self.ErrGamma)**2 )*Phi
      dPhi = np.sqrt( (self.ErrNorm/self.Norm)**2 + (((1./(1-self.Gamma)-mylog(ene))-(mylog(self.Emin)*np.power(self.Emin,1-self.Gamma)-mylog(self.Emax)*np.power(self.Emax,1-self.Gamma))/D)*self.ErrGamma)**2 )*Phi

    if self.Model == "PLExpCutoff":
      Phi  = self.Norm*np.power(ene/self.Ed,-self.Gamma)*np.exp(-(ene-self.Ed)/self.Ecut)
      dPhi = np.sqrt( (self.ErrNorm/self.Norm)**2 + (mylog(ene/self.Ed)*self.ErrGamma)**2 + ((ene-self.Ed)/self.Ecut**2*self.ErrEcut)**2 )*Phi

    if self.Model == 'LogParabola' :
      x = ene/self.Ed;
      Phi = self.Norm*np.power(x,-self.Gamma-self.Beta*mylog(x))

      dPhiDNorm = 1/self.Norm
      dPhiDalpha = mylog(x)[0]
      dPhiDbeta = mylog(x)[0]**2

      dPhi = Phi*sqrt((dPhiDNorm*self.ErrNorm)**2+(dPhiDalpha*self.ErrGamma)**2+(dPhiDbeta*self.ErrBeta)**2)

    return ene,Phi,dPhi

  def MakeSEDAndError(self,Emin=0,Emax=10e10):
    ene,Phi,dPhi = self.MakeFluxAndError(Emin=Emin,Emax=Emax)
    Phi *= ene**2*1.602e-6
    dPhi *= ene**2*1.602e-6
    return ene,Phi,dPhi

  def GetModel(self,Emin=0,Emax=10e10,SED=True):
    ene,Phi,_ = self.GetValueAndError(Emin=0,Emax=10e10,SED=SED)
    return ene,Phi

  def GetFlux(self,z=-1,Emin=0,Emax=10e10):
    ene,Phi,dPhi = self.GetValueAndError(Emin=0,Emax=10e10,SED=False)
    if z>0 :
      Phi,dPhi = Absorbtion(z,ene,Phi,dPhi,ModelFile="tau.dat",alpha=1.)
    flux = np.trapz(Phi,ene)
    dflux = np.trapz(dPhi,ene)  
    return flux,dflux

  def GetButterfly(self,Emin=0,Emax=10e10,SED=True):
    ene,Phi,dPhi = self.GetValueAndError(Emin=Emin,Emax=Emax,SED=SED)

    N = 2*len(ene)+1
    but = np.array(N*[1.])
    e_but = np.array(N*[1.])

    #First loop for the butterfly
    for i in xrange(len(ene)):
      e_but[i] = ene[i]
      but[i] = Phi[i]*exp(dPhi[i]/Phi[i])

    #second loop for the butterfly
    for i in xrange(len(ene)):
      e_but[len(ene)+i] = ene[len(ene)-1-i]
      but[len(ene)+i] = Phi[len(ene)-1-i]*exp(-dPhi[len(ene)-1-i]/Phi[len(ene)-1-i])

    #Close the Butterfly
    e_but[-1] = e_but[0]
    but[-1] = but[0]

    return e_but,but


def MakeTGraph(x,y,xtitle="",ytitle="",linecolor=2):
  tgr = ROOT.TGraph(len(x),array('f',x),array('f',y))
  tgr.SetLineColor(linecolor)
  tgr.GetHistogram().SetYTitle(xtitle)
  tgr.GetHistogram().SetXTitle(ytitle)
  return tgr

def GetTau(z_ref,filename='tau.dat'):
  # read a fille with z and tau tabulated, see franchescini
  lines=open(filename,"r").readlines()
  ilen=len(lines)-1

  E = array("f",ilen*[0.])
  tau = array("f",ilen*[0.])
  z = array("f",9*[0.])

  z[0]=0.01
  z[1]=0.03
  z[2]=0.1
  z[3]=0.3
  z[4]=0.5
  z[5]=1.0
  z[6]=1.5
  z[7]=2.0
  z[8]=3.0

  i=0
  ind = find_z_ind(z,z_ref)
  for line in lines[1:] :   
    words = string.split(line)
    taup = map(float,words[1:])

    E[i] = float(words[0])
    tau[i] = np.interp([z_ref],[z[ind],z[ind+1]],[taup[ind],taup[ind+1]])
    i+=1

  if z_ref<0.01 :
    return E,array("f",ilen*[0.])    
  return np.array(E), np.array(tau)

def find_z_ind(ztab,z):
    znp = np.array(ztab)
    return znp.searchsorted(z)-1
   # i=0
   # while ztab[i]<z:
   #   i+=1      
   # return i-1


def Absorbtion(z,ener,flux,dflux=None,ModelFile="tau.dat",alpha=1.):
 # z : redshift of the source
 # flux and ener are the tabulated flux (SED) of the source and the correspondind enegrgy (in MeV)

  E,tau = GetTau(z,ModelFile)
  N=len(flux)
  AbsFlux = np.zeros(N)
  AbsdFlux = np.zeros(N)
  for  i in xrange(N):
    depth = np.interp([ener[i]*1e-6],E,tau,0,0)
    AbsFlux[i]=flux[i]*exp(-alpha*depth[0])
    if not(dflux==None):
      AbsdFlux[i]=dflux[i]*exp(-alpha*depth[0])
  return AbsFlux,AbsdFlux

def DeAbsorbtion(z,ener,flux,ModelFile="tau.dat",alpha=1.):
 # z : redshift of the source
 # flux and ener are the tabulated flux (SED) of the source and the correspondind enegrgy (in MeV)

  E,tau = GetTau(z,ModelFile)
  N=len(flux)
  Result = np.zeros(N)
  for  i in xrange(N):
    detph = np.interp([ener[i]*1e-6],E,tau)
    Result[i]=flux[i]*exp(alpha*detph[0])
  return Result


