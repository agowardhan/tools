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


class linecube: 
    
    def __init__(self): 
        self.cubetype = 'line'
        pass 
        
    def cube_props(self, fitsfile, linename, restfreq, axes=[0,1,2], freqaxis=4, z=None): 
        
        self.fitsfile = fitsfile
        self.restfreq = restfreq 
        self.linename = linename 
        self.cubedat = np.squeeze(f.open(self.fitsfile)[0].data).transpose(axes)
        self.header = f.open(self.fitsfile)[0].header
        self.freqs = np.linspace(self.header['CRVAL'+str(freqaxis)], self.header['CRVAL'+str(freqaxis)] + self.header['CDELT'+str(freqaxis)]*self.header['NAXIS'+str(freqaxis)],self.header['NAXIS'+str(freqaxis)])
        self.z = z
        if z!=None: 
            self.vels = ut.Freq2Vel(x0=self.freqs/1e9, z = self.z, restfreq=self.restfreq)

    def extract_Spectra_square(self, restfreq, position, aperture, nbin = 1): 
        ' restfreq, position, aperture, nbin = 1'
        
        self.fitsfile = f.open(self.fitsfile)[0]
        self.w = wcs.WCS(self.fitsfile.header, naxis=2)
        self.position = position 
        self.nbin = nbin 
        self.pixs= self.w.all_world2pix([[self.position.ra.deg, self.position.dec.deg]], 0)
        self.icent, self.jcent = int(self.pixs[0][0]), int(self.pixs[0][1])
        self.spec = 1000*np.mean(self.cube[self.icent-self.nbin:self.icent+self.nbin,self.jcent-self.nbin:self.jcent+self.nbin,:], axis=(0,1)) ## returns in mJy   
        self.lineprops = ['square',aperture, position]

    def extract_Spectra_circ(self, source=True, fitsfile=None, pos=None, radius=0.1, axes= [0,1,2], smth=False, freqsmth = 0, verbose=False, restfreq=115.27120, z = 0.01): 
        'source=True, fitsfile=None, pos=None, radius=0.1, axes= [0,1,2], smth=False, freqsmth = 0, verbose=False, restfreq=115.27120, z = 0.01'
        if source: 
            print "Using cube:", self.fitsfile, "for the line", self.linename 
            fitsfile=self.fitsfile 
            restfreq = self.restfreq
            cube = self.cubedat
            header = self.header
            freqs = self.freqs
            z = self.z
            
        else: 
            header = f.open(fitsfile)[0].header
            cube = np.squeeze(f.open(fitsfile)[0].data).transpose(axes)
            freqs = np.linspace(header['CRVAL'+str(freqaxis)], header['CRVAL'+str(freqaxis)] + header['CDELT'+str(freqaxis)]*header['NAXIS'+str(freqaxis)],header['NAXIS'+str(freqaxis)])

        df, dx, dy = np.shape(cube)
        apertures = SkyCircularAperture(pos, r=radius* aunits.arcsec)
        w = WCS(fitsfile).celestial
        pix_aperture = apertures.to_pixel(w)
        
        if dx!= dy: 
            print('warning : image is not square, may need a transpose')

        mask = pix_aperture.to_mask()[0]
        image = mask.to_image(shape=((dx, dy))).astype('bool')
        vels = ut.Freq2Vel(x0=freqs/1e9, z = z , restfreq=restfreq)

        if smth: 
            cube = ndimage.gaussian_filter(cube, sigma=[freqsmth,0,0], mode='reflect')[::freqsmth, :, :]
            freqs = freqs[::freqsmth]
            vels=vels[::freqsmth]

        spec = np.mean(cube[:,image], axis=(1))

        if verbose: # need to add more messages 
            print('assuming NAXIS 4 is the frequency one')
        
        self.vels = vels 
        self.freqs = freqs
        self.spec = spec 
        self.lineprops = ['circle', radius, pos] # sets the line properties so that we can look them up later 
        return vels, freqs, spec
    
    def fitGaussImg(self, fitsfile, data_mom0=None, p0=None, plot=True): 
        'fitsfile, data_mom0=None, p0=None, plot=True'
        self.TDGauss_popt, self.intflux = GaussFlux(self.fitsfile, data_mom0=data_mom0, p0=p0, plot=plot)

    def fitGaussSpec(self,  p0=[0.0, 1.0, 10.0], frame='VEL'): 
        'p0=[0.0, 1.0, 10.0], frame=VEL'
        self.p0 = p0 
        if frame is 'VEL':
            self.popt_vel, self.pcov_vel = curve_fit(f=ut.Mgauss,xdata=self.vels, ydata=self.spec, p0=self.p0)
            self.popt_vel_err = np.diag(np.sqrt(self.pcov_vel))
        if frame is 'FREQ': 
            self.popt_freq, self.pcov_freq = curve_fit(f=ut.Mgauss,xdata=self.freqs, ydata=self.spec, p0=self.p0)
            self.popt_freq_err = np.diag(np.sqrt(self.pcov_freq))

    def plot(self, frame='VEL'): 
        'frame = VEL or FREQ'
        if frame is 'VEL':
            plt.figure()
            plt.axhline(0, color='grey')
            xnew = np.linspace(min(self.vels), max(self.vels), 1000)
            plt.plot(xnew, ut.Mgauss(xnew, *self.popt_vel), color='black', linestyle='--', linewidth=1.0)
            plt.plot(self.vels, self.spec)
        if frame is 'FREQ':
            plt.figure()
            xnew = np.linspace(min(self.freqs), max(self.freqs), 1000)
            plt.axhline(0, color='grey')
            plt.plot(xnew, ut.Mgauss(xnew, *self.popt_freq), color='black', linestyle='--', linewidth=1.0)
            plt.plot(self.freqs, self.spec)
        
    def calc_Mdyn(self, vel, radius, vel_circ=False):  # can pass it either the FWHM velocity or the circular velocity, radius should be kpc
        ' vel_fwhm, radius'
        self.vcirc = vel
        if not vel_circ:   
            self.vcirc.mult(1/(2 * (np.log(2))**0.5))
            
        self.Mdyn = value() 
        self.Mdyn.value = (self.vcirc.value *1e3)**2  * radius.value * const.kpc_to_m / const.G_si /const.Msun_to_kg
        self.Mdyn.error = self.Mdyn.value * ((2 * self.vcirc.frac)**2 + radius.frac**2)**0.5
        self.Mdyn.unit = 'Msun'
        
    def create_mom0(self, vel1, vel2, mask=None): 
        'vel1, vel2, mask=None'
        idx1, idx2  = ut.find_nearest(self.vels, vel1)[1], ut.find_nearest(self.vels, vel2)[1]
        if mask:
            self.mom0 = np.sum(self.cubedat[idx1:idx2,mask], axis=0)
        else: 
            self.mom0 = np.sum(self.cubedat[idx1:idx2,:,:], axis=0)
        self.mom_props = [vel1, vel2, mask]
        
        
