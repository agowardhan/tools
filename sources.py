import utils as u 
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.io import fits
from astropy import units as units

# file_W  = u.read_casa_fits('../data/light/ir10038_dustW_fit.txt', fits = '../data/light/ir10038_cont.fits')
# file_E  = u.read_casa_fits('../data/light/ir10038_dustE_fit.txt', fits = '../data/light/ir10038_cont.fits')
 
class source:
    
    '''
    A class for each target, defining the source redshift, position, spectral cubes 
    and their rest frequencies
    '''
    
    def __init__(self, z, name, LIR, nu_obs):
        self.z  = value() # valiur + uncertainity 
        self.name = name 
        self.scale = u.CC_dl.main(self.z.value,u.riechers2011, wantscale=True)
        self.DL = u.CC_dl.main(self.z.value,u.riechers2011) * u.Mpc_to_kpc * u.kpc_to_m * u.m_to_cm 
        self.LIR = LIR 
        self.complist = [] # save final fluxes and coordinates and velocity ranges here 
 
    def dust_props(self, nu_obs, flux, semiMaj, semiMin, Td=25): 
        
        self.FluxDust = flux
        self.semiMaj = semiMaj
        self.semiMin = semiMin

        self.Rdust_arcs = (self.semiMaj+self.semiMin)/4.0 * scale
        self.Rdust_kpc =  self.Rdust_arcs * scale # in kpc 
        
        self.Td = Td # this is the most recent calculation Td - the one that was used to calculate the dust mass 
        self.Mism = fd.Mism(nu_obs=nu_obs, Snu=self.FluxDust.value, Td=self.Td, z=self.z) ## need to propagate error here.

        self.area = value()
        self.area.value  = np.pi * self.semiMaj * self.semiMin * self.scale**2
        self.area.err  = np.pi * self.semiMaj.err * self.semiMin *self.scale**2 +  np.pi * self.semiMaj * self.semiMin.err *scale**2

        self.sf_density = value()
        self.sf_density.value = self.LIR * 1e-10/self.area.value 
        self.sf_density.error = self.sf_density.value * self.area.error/self.area.value
    
    def add_coords(self, name, coords, frame, name, flux): 
        "Add component coordinates i.e. blueshifted coordinates, redshifted coordinates"
        temp = value(flux)
        temp_name = 
        
    def add_line(self, name, cube=None)
    
        self.xcsize = -1
        self.ycsize = -1

    sourcefilename_base = 'ir10038'
    filetype='cont'
    openfile = fits.open(u.DATPATH+'/fits/'+ sourcename +'_' + filetype + '.fits')
    header = openfile[0].header
    naxis = header['NAXIS']
    for i in range(naxis):
        if 'FREQ' in header['CTYPE'+str(i+1)]: 
            nu = header['CRVAL'+str(i+1)] 
    
   # to convert between Mpc and centimetres  

    ### NEED TO PUT A BASIC FUNCTION FOR THIS 
    
    flux_dust_W = value(4.75e-3, 0.22e-3, units.mJy)    
    semiMaj = value(0.484, 0.020, units.arcs)
    semiMin = value(0.394, 0.016, units.arcs)

    # for E 
    


#############################
#Dynamical mass measurements# 
#############################

# For this, it is essential to ensure that we are using comparative radii and velocities. 
# so we extract the HCN and HCO+ spectra from an aperture of 0.2

hcn_fwhm_vel = 151.2
hcop_fwhm_vel = 138 # second component may b be broad because of contamination by HCN-vib

#Mdyn3 = (hcn_circ_vel/1.665  *1e3)**2 * radius_hcn * u.kpc_to_m / u.G_si /u.Msun_to_kg # v needs to be in m/s

vfwhm = 65.4002088 * 2.335 # km/s from fitting to the broad line within a  0.23'' diameter
vcirc = vfwhm/(2 * (np.log(2))**0.5)

vfwhm2 = 174 
vcirc2 = vfwhm2/(2 * (np.log(2))**0.5)

Mdyn1 = (vcirc *1e3)**2 * rad_dust_W * u.kpc_to_m / u.G_si /u.Msun_to_kg # v needs to be in m/s
Mdyn2 = (vcirc2 *1e3)**2 * 0.125 * scale_ir10038 * u.kpc_to_m / u.G_si /u.Msun_to_kg # v needs to be in m/s

mass_sf_density1 = Mdyn1 /area_W
mass_sf_density2 = Mdyn2 /area_W
mass_sf_density3 = mism_ir10038_W /area_W




class FitsFile:
    '''
    To describe a fits cube  
    '''
    
    def __init__(self, cube):
        
        self.flip() 
        ## flips it the right way up - 
            
    def flip():
        """ 
        Adjusts a FITS cube to have axes in the order of RA, DEC, and frequency. 
        """
        
    def ReadHeader(self):
     
    
    def spectrum(self, xcoord, ycoord, pix=False): 
        '''
        Returns the spectrum from a fixed position/aperture from a cube 
        '''
        return -1 
    
    def moments(self): 
        return -1 
    
    def SpecMap(self):
        return -1 

    def pv(self): 
        return -1 

    
class value(self, val, err, unit=None):
    '''
     Structure to save a value, error, and unit. 
     '''
    self.val = val 
    self.err = err
    self.unit = unit 
    self.physics = physics()

               
class lines:
    '''
    A class to deal with fits files 
    
    '''
    
    def __init__(self, cube, restfreq, name):
         
            
    def info(self):
        print things 
        return -1 
    
    def spectrum(self): 
        '''
        Returns the spectrum from a fixed position/aperture from a cube 
        '''
        return -1 
    
    def moments(self): 
        return -1 
    
    def SpecMap(self):
        return -1 

    def pv(self): 
        return -1 
    
    def dust(self): 
        
        # dust continuum looks like a single nucleus, unresolved. 
        # RA : 04:46:49.529 -048.33.30.05531
        # pix center: 251.45, 248.23 
        # deconvoled image size : 137+/-13, 119+/-14 
        # beam: 0.15x0.13, 6.95 degrees 
        # flux : integrated, 17.38 \pm 0.74 mJy 
        # peak : 9.58+-/0.28 mJy 

class continuum:
    '''
    To describe a fits cube 
    '''
    
    def __init__(self, cube, restfreq, name):
         
            
    def info(self):
        print things 
        return -1 
    
    def spectrum(self): 
        '''
        Returns the spectrum from a fixed position/aperture from a cube 
        '''
        return -1 
    
    def moments(self): 
        return -1 
    
    def SpecMap(self):
        return -1 

    def pv(self): 
        return -1 
    
    def dust(self): 
        

        
        
        # dust continuum looks like a single nucleus, unresolved. 
        # RA : 04:46:49.529 -048.33.30.05531
        # pix center: 251.45, 248.23 
        # deconvoled image size : 137+/-13, 119+/-14 
        # beam: 0.15x0.13, 6.95 degrees 
        # flux : integrated, 17.38 \pm 0.74 mJy 
        # peak : 9.58+-/0.28 mJy 

#### For the HCN 

# dust continuum looks like a single nucleus, unresolved. 
# RA : 04:46:49.529 -048.33.30.05531
# pix center: 251.45, 248.23 
# deconvoled image size : 137+/-13, 119+/-14 
# beam: 0.15x0.13, 6.95 degrees 
# flux : integrated, 17.38 \pm 0.74 
# peak : 9.58+-/0.28 

### for the HCO+ -- forHCN and HCO+ binning over the velocity -300 to 300 km.s


# asuming it is all SFR 
# the mass estimates are accurate enogh for a first pass -- all basic fits 

# file_W  = u.read_casa_fits('../data/light/ir10038_dustW_fit.txt', fits = '../data/light/ir10038_cont.fits')
# file_E  = u.read_casa_fits('../data/light/ir10038_dustE_fit.txt', fits = '../data/light/ir10038_cont.fits')


