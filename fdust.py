import numpy as np
import CC_dl 
import utils as u 
import math
import const as const 

def GammaRJ(Td = 25, nu_obs=352.6970, z=0.0): 
    '''
    Function to calculate the rayleigh-jeans correcting factor GammaRJ for the dust continuum emission 
    
        Td = dust temperature in Kelvin 
    
        nu = observed frequency of emission in GHz 
    
        z = source redshift 
    
    Returns : dimensionless RJ correction factor 
    
    '''
    xx = 0.672 * nu_obs*(1.0+z)/350 /(Td/25.)
    return xx/(math.exp(xx) - 1.0)
    
    
def Lnu_850(Snu, nu_obs=352.6970, z=0.0, cosmo=const.riechers2011):
    '''
    Function to calculate 850 micron luminosity from other observed dust continuum emission -- hopefully on the RJ tail 
    
        Snu = flux density in Jy
            
        nu = observed frequency of emission in GHz 
    
        z = source redshift 
        
        dL = luminosity distance in Mpc 
    
    Returns : specific luminosity in ergs/s/Hz 
    
    '''
    dL = CC_dl.main(z,cosmo)# * 3.086e+24 # to convert between Mpc and centimetres  
    snu_850 = Snu_850(Snu, nu_obs, z)
    return 1.19e27  * snu_850 * dL**2 /(1+z) * GammaRJ(25, u.Wave2Freq(850),0.0 )/GammaRJ(25, nu_obs ,z )
   
    
def Snu_850(Snu, nu_obs=352.6970, z=0.0):
    '''
    Function to calculate 850 micron luminosity from other observed dust continuum emission -- hopefully on the RJ tail 
    
        Snu = flux density in Jy
            
        nu = observed frequency of emission in GHz 
    
        z = source redshift 
        
        dL = luminosity distance in Mpc 
    
    Returns : specific luminosity in ergs/s/Hz 
    
    '''
    return Snu * (u.Wave2Freq(850)/(nu_obs*(1+z)))**3.8 

    
def Mism(Snu, Td =25, nu_obs = 350.0, z=0.0, cosmo=const.riechers2011):
    
    '''
    Formula to calculate total ISM mass, from Scoville+14 paper on Arp 220. 
    The value 0.83 uses the alpha_870 value of 1.0 \pm 0.23 e20 ergs/s/Hz/Msun, check Scoville+14.    
    '''
    
    dL =  CC_dl.main(z,cosmo)
    return 0.83e10 * Snu *(dL/1e3)**2 /(1.0+z)**4.8  /(nu_obs/350.)**3.8 /(Td/25) /GammaRJ(Td=Td, nu_obs=nu_obs ,z=z )
    
    
def Mism2(Snu, nu_obs = 350.0, z=0.0, alp=6.7e19, Td = 25, cosmo=const.riechers2011):
    '''
    Formula to calculate total ISM mass, from Scoville+16 paper erratum. Use this with caution,  there is a Td dependence in alpha-850,
    which is not reflected in the function Mism2, which assumes a constant alpha. 
    
    '''
    dL =  CC_dl.main(z,cosmo)
    return 1.78e10 * Snu*(dL/1e3)**2 /(1.0+z)**4.8 * (u.Wave2Freq(850)/nu_obs)**3.8 *(6.7e19/alp) * 0.71/GammaRJ(Td=Td, nu_obs=nu_obs, z=z)