# %load imports.py
import numpy as np

import aplpy
from astropy.coordinates import SkyCoord
import astropy.constants as aconst
import scipy as sc
from scipy import ndimage, misc
from scipy import signal

import fdust as fd
import utils as ut 
import const as const
from astropy.io import fits as f

import matplotlib.pyplot as plt 
import matplotlib.patches as mpatche
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.colors import colorConverter

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

import matplotlib as mpl
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
# extras 

import sys 
import pandas as pd
import os 
import math
import numpy.ma as ma

import scipy as sc
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from scipy import ndimage, misc

import cPickle as pickle

from copy import deepcopy
from value import * 

        
class source_info:

    def __init__(self, z, name, position, LIR, cosmo=const.riechers2011):
        self.z = z  # valiur + uncertainity 
        self.position = position 
        self.name = name 
        self.scale = ut.CC_dl.main(z.value,const.riechers2011, wantscale=True)
        self.DL = ut.CC_dl.main(z.value,cosmo) * const.Mpc_to_kpc * const.kpc_to_m * const.m_to_cm 
        self.LIR = LIR 
        self.LIR_SFR = 1e-10 ## conversion factor from the LIR to the SFR - this is assuming a Chabrier IMF; for a Salpeter IMF, change to 1.7e-10 
  
  
  
           
class source: 
        
    '''
    A class for each target, defining the source redshift, position, spectral cubes 
    and their rest frequencies
    '''
    
    def __init__(self, z, name, position, LIR, cosmo=const.riechers2011, fee=value(0.0,1e-10), sil=0.0, fagn=None):
        
        self.add_basics(z, name, position, LIR, cosmo=const.riechers2011)
        self.cubes=[]
        self.z = z  # valiur + uncertainity 
        self.position = position 
        self.name = name 
        self.scale = ut.CC_dl.main(z.value,const.riechers2011, wantscale=True)
        self.DL = ut.CC_dl.main(z.value,cosmo) * const.Mpc_to_kpc * const.kpc_to_m * const.m_to_cm 
        self.LIR = LIR *(1.0-fee.value)
        self.LIR_to_SFR = 1e-10
        self.fee = fee
        self.sil = sil
        self.fagn = fagn
        self.LIR_SFR = self.LIR * (1.0 - self.fagn.value)
        self.LIR_AGN = self.LIR * (self.fagn.value)
        
    def add_basics(self, z, name, position, LIR, cosmo=const.riechers2011): 
        self.info = source_info(z, name, position, LIR, cosmo=const.riechers2011)
        
    def add_cube(self, newcube): 
        flag=0
        for cube in self.cubes: 
            if cube.linename == newcube.linename: 
                print 'Cubename already exists. Pick another'
                flag=1
                break
        if flag<1: 
            self.cubes.append(newcube)
        
    def pick_cube(self, newcube): 
        flag=0
        for cube in self.cubes: 
            if cube.linename == newcube: 
                self.cube = cube 
                flag=1
                break
        if flag<1: 
            print 'Cube not found.'
    
    def dust_props(self, nu_obs, flux, semiMaj, semiMin, Td=25, listem=False):  # all the dust observations should be from 
        
        " The function expects the continuum flux, the size of the source, the dust temperature : alternative functionality will include fitting an image cube with the flux"
        
        self.FluxDust = flux
        self.semiMaj_arcsec = semiMaj
        self.semiMin_arcsec = semiMin
        
        self.semiMaj_kpc = deepcopy(semiMaj)  
        self.semiMin_kpc = deepcopy(semiMin)

        self.semiMaj_kpc.mult(self.scale, unit='kpc')
        self.semiMin_kpc.mult(self.scale, unit='kpc')
        
        self.Td = Td # this is the most recent calculation Td - the one that was used to calculate the dust mass 
        self.Mism = value()
        self.Mism.value = fd.Mism(nu_obs=nu_obs, Snu=self.FluxDust.value, Td=self.Td, z=self.z.value) ## need to propagate error here.
        self.Mism.error = fd.Mism(nu_obs=nu_obs, Snu=self.FluxDust.error, Td=self.Td, z=self.z.value) ## need to propagate error here.
        print('Mism-error not propagated yet - dont trust it.')

        self.area = value()
        self.area.value  = np.pi * self.semiMaj_kpc.value * self.semiMin_kpc.value / (4 *np.log10(2))
        self.area.error  = self.area.value * (self.semiMaj_kpc.frac**2  +  self.semiMin_kpc.frac**2)**0.5

        self.radDust = value()
        self.radDust.value = np.sqrt(self.area.value/np.pi)
        self.radDust.error = self.area.error/ (2 * np.pi * self.radDust.value) 

        self.sf_density = value()
        self.sf_density.value = self.LIR * self.LIR_to_SFR /self.area.value 
        self.sf_density.error = self.sf_density.value * (self.area.error/self.area.value)
        self.sf_density.unit = 'Msun/yr/kpc^2'  
        
        if listem: 
            print "Flux density (mJy): ", np.round(self.FluxDust.value,2), '+/-',np.round(self.FluxDust.error, 2)
            print "size (arcs): ", np.round(self.semiMaj_arcsec.value,2), '+/-', np.round(self.semiMaj_arcsec.error,2), np.round(self.semiMin_arcsec.value, 2), '+/-' ,np.round(self.semiMin_arcsec.error, 2)
            print "size (kpc) : ", np.round(self.semiMaj_kpc.value,2), '+/-',np.round(self.semiMaj_kpc.error,2), np.round(self.semiMin_kpc.value,2), '+/-', np.round(self.semiMin_kpc.error, 2)
            #print "Radius (kpc) : ", np.round(self.radDust.value,2), '+/-',np.round(self.radDust.error, 2)
            print "Td (K): ", np.round(self.Td,2)
            print "Area (kpc^2): ", np.round(self.area.value,2), '+/-',np.round(self.area.error, 2) 
            print "Mism (Msun): ", np.round(self.Mism.value,2),'+/-',np.round(self.Mism.error,2)
            print "sf density (Msun/yr/kpc^{2}): ", np.round(self.sf_density.value,2),'+/-',np.round(self.sf_density.error, 2)
            
            