import numpy as np
import CC_dl 
from scipy.special import erf
from astropy  import constants as aconst
from astropy  import units as aunits
from astropy.io import fits as f 
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture
from astropy.wcs import WCS
from scipy import ndimage, misc


def extractSpectra(fitsfile, pos, radius, axes= [0,1,2], smth=False, freqsmth = 0, verbose=False, restfreq=115.27120,z = 0.01): 
    
    apertures = SkyCircularAperture(pos, r=radius* aunits.arcsec)
    w = WCS(fitsfile).celestial
    pix_aperture = apertures.to_pixel(w)
    header = f.open(fitsfile)[0].header
    
    cube = np.squeeze(f.open(fitsfile)[0].data).transpose(axes)
    df, dx, dy = np.shape(cube)
    if dx!= dy: 
        print('warning : image is not square, may need a transpose')
    
    mask = pix_aperture.to_mask()[0]
    image = mask.to_image(shape=((dx, dy))).astype('bool')
    freqs = np.linspace(header['CRVAL4'], header['CRVAL4'] + header['CDELT4']*header['NAXIS4'],header['NAXIS4'])
    vels = Freq2Vel(x0=freqs/1e9, z = z , restfreq=restfreq)

    if smth: 
        cube = ndimage.gaussian_filter(cube, sigma=[freqsmth,0,0], mode='reflect')[::freqsmth, :, :]
        freqs = freqs[::freqsmth]
        vels=vels[::freqsmth]

    spec = np.mean(cube[:,image], axis=(1))

    if verbose:
        print('assuming NAXIS 4 is the frequency one')

    return vels, freqs, spec


def unique(a):
    """ return the list with duplicate elements removed """
    return list(set(a))

def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))


def Lprime(z, sco, nu_rest, cosmo): 
    '''
    Given a redshift,velocity integrated line flux, rest frame frequncym and a cosmology, calculates 
    L' based on formula in Solomon+05. 
    
    Inputs: 
    
    z: float 
         redshift 
    
    sco: float 
        velocity integrated line flux in Jy km/s 
    
    nu_rest: float 
        rest frame frequency in GHz 
    
    cosmo: [float, float, float]
        selected cosmology in the format eg. riechers2011 = [71., 0.27,0.73]
        
    Output: float 
        L' value in K km/s pc^-2  
        
    '''
    
    dl = CC_dl.main(z,cosmo,verbose=-1)
    temp = 3.25*10.**7. *dl**2. *sco
    temp2 = ((1.+z) *nu_rest**(2.))
    return temp/temp2

# defining all needed functions 

def rebin(x,y,bins): 
    n = len(x)
    xnew = np.zeros((n/bins))
    ynew = np.zeros((n/bins))
    for i in range(1,n/bins): 
        xnew[i] = np.mean(x[bins*(i-1): bins*i])
        ynew[i] = np.sum(y[bins*(i-1): bins*i])/bins
    
    return xnew, ynew

def smooth(x,window_len=10,window='hanning'):
    
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]


def sed_inplace(filename, pattern, repl):
    '''
    Perform the pure-Python equivalent of in-place `sed` substitution: e.g.,
    `sed -i -e 's/'${pattern}'/'${repl}' "${filename}"`.
    '''
    # For efficiency, precompile the passed regular expression.
    pattern_compiled = re.compile(pattern)

    # For portability, NamedTemporaryFile() defaults to mode "w+b" (i.e., binary
    # writing with updating). This is usually a good thing. In this case,
    # however, binary writing imposes non-trivial encoding constraints trivially
    # resolved by switching to text writing. Let's do that.
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
        with open(filename) as src_file:
            for line in src_file:
                tmp_file.write(pattern_compiled.sub(repl, line))

    # Overwrite the original file with the munged temporary file in a
    # manner preserving file attributes (e.g., permissions).
    shutil.copystat(filename, tmp_file.name)
    shutil.move(tmp_file.name, filename)

# Do it for Johnny.
#sed_inplace('file', r'^\# deb', 'deb')

def find_nearest(array,value):
    
    '''
    Given a 1-D array and a float value, searched for the closest value in the array and returns the index. 
    
    Inputs: 
    
    array: array[float]
        array in which to search 
    
    value: float 
        number for which to search 
    
    Output: int 
        index of closest value 
        
    '''
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx

def sed_integral(att,lamb1,lamb2, verbose=0): 
    [lambs, flux] = att 
    integ = 0. 
    indx1 = find_nearest(lambs,lamb1)[1]
    indx2 = find_nearest(lambs,lamb2)[1]
    if verbose: 
        print 'Found indices:', indx1, indx2
        print 'Corresponding to', lambs[indx1], lambs[indx2]
    for i in range(indx1,indx2): 
        delL = np.abs(lambs[i] - lambs[i-1])
        integ = integ + flux[i]*delL
    if verbose: 
        print 'Integrated SED is', integ, 'in L_{\odot}'
    return integ 

def spec_integral(att,lamb1,lamb2, verbose=0): 
    [lambs, flux] = att 
    integ = 0. 
    indx1 = find_nearest(lambs,lamb1)[1]
    indx2 = find_nearest(lambs,lamb2)[1]
    if verbose: 
        print 'Found indices:', indx1, indx2
        print 'Corresponding to', lambs[indx1], lambs[indx2]
    for i in range(min(indx1,indx2), max(indx1,indx2)+1): 
        delL = abs(lambs[i] - lambs[i-1])
        integ = integ + flux[i]*delL
    if verbose: 
        print 'Integrated spectrum', integ, 'in L_{\odot}'
    return integ 

def sed_integral_lsuns(att,lamb1,lamb2, verbose=0): 
    [lambs, flux] = att 
    integ = 0. 
    indx1 = find_nearest(lambs,lamb1)
    indx2 = find_nearest(lambs,lamb2)
    if verbose: 
        print 'Found indices:', indx1, indx2
        print 'Corresponding to', lambs[indx1], lambs[indx2]
    for i in range(indx1,indx2): 
        integ = integ + flux[i]
    if verbose: 
        print 'Integrated SED is', integ, 'in L_{\odot}'
    return integ 

# function to calculate the mass outflow rate , Maiolino2012 

def Mgauss(x, *params):

    '''
    Given a 1-D array and a set of initial parameters, fit M gaussians.  
    
    Inputs: 
    
    x: array[float]
        list of x-axis values 
    
    params: array[float], 
        len(params) must be multiple of 3. Each set of 3 consequetive values are the center, amplitude, width, for a gaussian. 
    
    Output: array
        returns sum of all gaussians 
        
    '''
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp*1./(wid*np.sqrt(2*np.pi)) *np.exp(-1*(x - ctr)**2/(2*wid*wid))
    return y

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()



def phi(x0, mu, sigma, n): 
    '''
    Finds the area under the a Gaussian curve from -inf to x = x0, given the normalization, mean, and sigma. 
    '''
    return (n/2) * (1 + erf((x0 - mu) / (sigma * np.sqrt(2))))


def Freq2Vel(x0, restfreq, z=0.0): 
    '''
    Given a redshift, a rest frequency, and an observed spectrum, converts between the observed frequency to the rest frame v
    velocity 
    Reference: https://www.iram.fr/IRAMFR/ARN/may95/node4.html
    Assumes that the x0 and restfreq are in the same frequency units. 
    
    '''
    
    ele2 =  299792.458*(restfreq - x0*(1.0 + z))/restfreq
    return ele2 


def Wave2Vel(x0, restwave, z=0.0): 
    '''
    Given a redshift, a rest frequency, and an observed spectrum, converts between the observed wavelength to the rest frame v
    velocity 
    
    Reference: https://www.iram.fr/IRAMFR/ARN/may95/node4.html

    
    '''
    ele2 =  299792.458*(x0/(1.0 + z) - restwave)/ (x0/(1.0 + z))
    return ele2 

def Vel2Freq(v0, restfreq, z=0.0):
    '''
    Given a redshift, a rest frequency, and an observed spectrum, converts the velocity spectrum to the observed spectrum 
    velocity 
    '''
    fs = restfreq * (1.0 - v0/299792.458)/(1.0 + z)
    return fs 

def Wave2Freq(lamb):
    '''
    Given a wavelength in microns, converts to Frequncy in GHz.  
    '''
    fs = 299792.458/lamb 
    return fs 


def Vesc(M, r):
    '''
    Given mass in Msun, and radius in kpc, returns escape velocity in km/s. 
    '''
    return np.sqrt(2* G_si * M * Msun_to_kg/(r * kpc_to_m)) * 1e-3 # m should be in solar, r should be in kpc 

def Mout(m, v, r):
    '''
    Mass in solar masses 
    v in km/s 
    r in kpc 
    Returns outflow properties of mass outflow rate (in Msun/yr), outflow momentum flux (in cgs), and outflow energy flux (in cgs)
    '''
    mout = m*v * 3.24078e-17 /r * 31556926
    pout = mout * v  *2e33 * 1e5 # converting to cgs 
    eout = 0.5 * mout*2e33 * v**2 *1e10 # converting to cgs 
    return mout, pout, eout

def MoutRate(Mout, R, v0): 
    '''
    Given outflowing gas, with mass Mout in outflow, radius of outflow (in kpc)
    and maximal outflow velocity v(in km/s), calculate the mass outflow rate, kinetic energy of outflow, 
    and momentum in outflow \citep{Maiolino+12}
    '''
    v = 1.02201e-9*v0
    t1 = v*Mout/R         # in the case of a spherical outflow # make sure mass is the H2 gas mass. 
    t2 = 0.5*v0*v0*t1*6.3e35   # in units of ergs /s 
    t3 = v0*t1*6.3e30 # in CGS units 
    return t1, t2 ,t3 

def smooth_spec(arr, window_len=3):
    
    xx = arr[:,0][::window_len]
    yy = smooth(arr[:,1], window_len=window_len, window='hanning')[::window_len]
    nearr = np.zeros((len(xx), 2))
    nearr[:,0] = xx 
    nearr[:,1] = yy 
    return nearr 

def ImgDefaults(fig, sb =3, ls = 30 , xfm = 'hh:mm:ss', yfm =('dd:mm:ss'), corner='bottom right'): 
    fig.ticks.set_linewidth(2)
    fig.ticks.set_length(8)
    fig.axis_labels.set_font(size=ls)
    fig.tick_labels.set_font(size=ls)
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')
    fig.axis_labels.set_ypad(-30)
    fig.ticks.set_color('black')
    fig.add_scalebar(sb/3600.)
    fig.scalebar.set_corner(corner)
    fig.scalebar.set_linestyle('solid')
    fig.scalebar.set_color('black')
    fig.scalebar.set_linewidth(4)
    fig.scalebar.set_font_size(ls)
    fig.scalebar.set_label(str(int(sb)) + ' arcsec')
    return fig 




def sSersic2D((x,y), amplitude, r_eff, n, x_0, y_0,ellip, theta): 
    mod = Sersic2D(amplitude, r_eff, n, x_0, y_0, ellip, theta)
    img = mod(x, y)
    return img.ravel()

def printpopt(popt, pcov): 
    "A quick function to print results from scipy.optimize.curvefit"

    for i in range(len(popt)): 
        print(popt[i]),(np.diag(pcov)**0.5)[i]
        
def radial_profile(data, center):
    
    'Gets a radial profile centered on a given pixel coordinate for an image'
    
    x, y = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 


def GaussFlux(fitsfile, data_mom0=None, p0=None, plot=True): 
    
    if data_mom0 is None: 
        data_mom0 = np.squeeze(f.open(fitsfile)[0].data)
    
    xlen, ylen = np.shape(data_mom0)
    
    if p0 is None: 
        p0 = [1e-3, xlen/2, ylen/2, 1, 1, 0, 0]
        
    header = f.open(fitsfile)[0].header
    x = np.arange(xlen)
    y = np.arange(ylen)
    (x,y) = np.meshgrid(x,y)
    popt, pcov = sc.optimize.curve_fit(twoD_Gaussian, (x, y), data_mom0.flatten(), p0=p0)
    
    amp = value()
    x0 = value()
    y0 = value()
    sigx = value()
    sigy = value()
    pa = value()
    offset = value()
    
    amp.value, x0.value, y0.value, sigx.value, sigy.value, pa.value, offset.value = popt 
    amp.error, x0.error, y0.error, sigx.error, sigy.error, pa.error, offset.error = np.diag(np.sqrt(pcov)) 

    arcsec_per_pix = header['CDELT1']*3600
    beam_ma =  header['BMAJ'] / header['CDELT1']
    beam_mi = header['BMIN']/ header['CDELT1'] # in pixels 
    
    fact_corr = (beam_mi*beam_ma)**0.5
    popt_ret = [amp, x0, y0, sigx, sigy, pa, offset] 

    amp.calcFrac()
    x0.calcFrac()
    y0.calcFrac()
    sigx.calcFrac()
    sigy.calcFrac()
    pa.calcFrac()
    offset.calcFrac()
    
    nbeams = value()
    nbeams.value = 2.355**2 * sigx.value * sigy.value / (beam_ma*beam_mi)
    nbeams.error = nbeams.value * ((sigx.frac)**2  + (sigy.frac)**2)**0.5
    nbeams.calcFrac()

    int_flux = value()
    int_flux.value = amp.value*nbeams.value
    int_flux.error = int_flux.value * (amp.frac**2 + nbeams.frac**2)**0.5

    ### should return position in RA, DEC 
    ## and size in pixels 
    
    if plot: 
        data_fitted = twoD_Gaussian((x, y), *popt)
        fig, ax = plt.subplots(1, 1)
        ax.imshow(np.squeeze(data_mom0), cmap=plt.cm.jet, origin='bottom',
            extent=(x.min(), x.max(), y.min(), y.max()))
        ax.contour(x, y, data_fitted.reshape(xlen,ylen), 8, colors='w')
        plt.show()
    
    return popt_ret, int_flux, arcsec_per_pix, fact_corr


def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset, beamx=1,beamy=1):
    
    '(x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset, beamx=1,beamy=1'
    
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    
    g = ndimage.gaussian_filter(g, sigma=[beamx,beamy], mode='reflect')

    return g.ravel()

