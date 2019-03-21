import aplpy
from astropy.coordinates import SkyCoord
import numpy as np 
import const as const
from astropy.io import fits as f

import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.colors import colorConverter

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

from matplotlib import rcParams
import matplotlib

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


def makemom0(name,  vmin, vmax, centx, centy, width, height, stretch='arcsinh', unit='mJy beam$^{-1}$km s$^{-1}$', spacing=0.5*2.7e-4, grid=False, typename=None, typecolor='white',noborders=False, colorbarticks=None):
    
    f1= plt.figure(figsize=(5,5))
    fig = aplpy.FITSFigure(const.PATH+'data/fits/'+name+'.fits', figure=f1)

    if grid: 
        fig.add_grid()
        fig.grid.set_xspacing(spacing)  # degrees
        fig.grid.set_yspacing(spacing)  # degrees
        fig.grid.set_color('grey')
        fig.grid.set_alpha(0.5)

    fig.show_colorscale(stretch=stretch,cmap='Spectral_r', vmin=vmin, vmax=vmax)#,vmax=None,pmax=98)#96#
    fig.recenter(centx, centy, width=width, height=height)
    fig.set_nan_color('black')
    
    
    if typename: 
        fig.add_label(0.1, 0.9, typename, relative=True, size=15, color=typecolor,horizontalalignment='left')

    if noborders: 
        fig.axis_labels.hide_y()
        fig.axis_labels.hide_x()
        fig.tick_labels.hide_x()
        fig.tick_labels.hide_y()
    else: 
        fig.axis_labels.set_font(size=13)
        fig.tick_labels.set_font(size=13)
        fig.tick_labels.set_xformat('mm:ss.s')
        fig.tick_labels.set_yformat('mm:ss.s')
        fig.axis_labels.set_ypad(-15)

        fig.add_colorbar()
        fig.colorbar.set_width(width=0.1)
        fig.colorbar.set_font(size=10)

        if unit: 
            fig.colorbar.set_axis_label_text('mJy beam$^{-1}$km s$^{-1}$')
            fig.colorbar.set_axis_label_font(size=10)
        if colorbarticks: 
            fig.colorbar.set_ticks(ticks=colorbarticks)

    fig.save(const.SPATH+name+'.png', adjust_bbox='tight', format='png', dpi=200)
    
def makemom1(name, vmin, vmax, centx, centy, width,height, spacing=0.5*2.7e-4, grid=False, typename=None, noborders=False, colorbarticks=None):
    
    f1= plt.figure(figsize=(5,5))
    fig = aplpy.FITSFigure(const.PATH+'data/fits/'+name+'.fits', figure=f1)
    
    if grid: 
        fig.add_grid()
        fig.grid.set_xspacing(spacing)  # degrees
        fig.grid.set_yspacing(spacing)  # degrees
        fig.grid.set_color('grey')
        fig.grid.set_alpha(0.5)
    fig.show_colorscale(stretch='linear',cmap='seismic', vmin=vmin, vmax=vmax)#,vmax=None,pmax=98)#96#
    fig.recenter(centx, centy, width=width, height=height)
    fig.set_nan_color('black')
    
    
    if typename: 
        fig.add_label(0.1, 0.9, typename, relative=True, size=15, color='white', horizontalalignment='left')

    if noborders: 
        fig.axis_labels.hide_y()
        fig.axis_labels.hide_x()
        fig.tick_labels.hide_x()
        fig.tick_labels.hide_y()

    else: 
        fig.axis_labels.set_font(size=13)
        fig.tick_labels.set_font(size=13)
        fig.tick_labels.set_xformat('mm:ss.s')
        fig.tick_labels.set_yformat('mm:ss.s')
        fig.axis_labels.set_ypad(-15)
        
        fig.add_colorbar()
        fig.colorbar.set_width(width=0.1)
        fig.colorbar.set_axis_label_text('Velocity (km s$^{-1}$)')
        fig.colorbar.set_axis_label_font(size=10)
        fig.colorbar.set_font(size=10)
        #fig.colorbar.set_ticks(ticks=[0.1,0.2])

    fig.save(const.SPATH+name+'.png', adjust_bbox='tight', format='png', dpi=200)

def makeHST(outfile, width=0.004,height=0.004, vmin=1e-2, vmax=10.0,centx=const.x_ir04454, centy=const.y_ir04454,scalebar=5, labels=False, labelcoords=None, 
            labelnames=None, labelcolor='white', box=None, wbox=None, hbox=None, name=None, returnfig=False): 

    fig= plt.figure(figsize=(5,5))
    f1 = aplpy.FITSFigure(outfile,  figure=fig)#, north=True)
    f1.show_colorscale(stretch='arcsinh',cmap='inferno', vmax=vmax, vmin =vmin)#,vmax=None,pmax=98)#96#
    f1.recenter(centx,centy , width=width, height=height)

    f1.axis_labels.set_font(size=15)
    f1.tick_labels.set_font(size=15)
    f1.add_scalebar(scalebar/3600., color='black')
    f1.scalebar.set_corner('top right')
    f1.scalebar.set_label(str(scalebar) + ' arcsec')
    f1.scalebar.set_linestyle('solid')
    f1.scalebar.set_color('white')
    f1.scalebar.set_linewidth(4)
    f1.scalebar.set_font_size(25)
    f1.tick_labels.set_xformat('mm:ss.s')
    f1.tick_labels.set_yformat('mm:ss.s')
    f1.ticks.set_color('white')
    f1.axis_labels.set_ypad(-20)

    if labels: 
        for i in range(len(labelnames)): 
            f1.add_label(labelcoords[2*i], labelcoords[2*i+1], labelnames[i], relative=True, size=15, color=labelcolor)

    if box: 
        f1.show_rectangles(centx, centy, width=wbox/3600, height = hbox/3600, edgecolor='lime', linewidth=2.0)

    if returnfig: 
        return f1
    else: 
        fig.savefig(const.SPATH+str(name)+'_hst.png', format='png', bbox_inches = 'tight' , dpi=200)
        

def makeCont(outfile, vmin=1e-2, vmax=10.0, width=3.0/3600,height=3.0/3600, centx=const.x_ir04454, centy=const.y_ir04454, scalebar=1, 
             levels=np.array([-5,-3, 3, 5, 9, 15, 30, 60]), rms=0.002, spacing=0.2*2.7e-4, beamx=0.1, beamy = 0.1, pa=10.0, name=None, 
             circles=None, circleparams=None, returnfig=False): 
    
    f1= plt.figure(figsize=(6,6))
    fig = aplpy.FITSFigure(outfile,  figure=f1)#, north=True)

    fig.show_colorscale(stretch='arcsinh',cmap='Spectral_r', vmin=vmin, vmax=vmax) #,vmax=None,pmax=98)#96#
    fig.set_nan_color('black')
    fig.recenter(centx,centy, width=width, height=height)

    fig.show_contour(outfile,
                     levels=np.multiply(levels, rms),
                     linewidths=2.0, overlap=True, colors='black')


    fig.add_grid()
    fig.grid.set_xspacing(spacing)  # degrees
    fig.grid.set_yspacing(spacing)  # degrees
    fig.grid.set_color('grey')
    fig.grid.set_alpha(0.5)

    fig.show_beam(beamx/3600.,beamy/3600.,pa,fill=True,color='grey')

    fig.beam.set_edgecolor('blue')
    fig.beam.set_frame(True)
    fig.beam.set_pad(0.5)
    fig.beam.set_alpha(0.9)
    fig.beam.set_hatch('/')

    fig.ticks.set_color('black')
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=15)
    fig.tick_labels.set_xformat('mm:ss.s')
    fig.tick_labels.set_yformat('mm:ss.s')
    fig.axis_labels.set_ypad(-5)

    fig.add_colorbar()
    fig.colorbar.set_width(width=0.1)
    fig.colorbar.set_axis_label_text('Jy beam$^{-1}$')
    fig.colorbar.set_axis_label_font(size=13)
    fig.colorbar.set_font(size=13)

    fig.add_label(0.1, 0.9,  "continuum",  relative=True, size=15, color='black', horizontalalignment='left')
    
    if circles: 
        for i in range(0, len(circleparams),3): 
            fig.show_circles(circleparams[i], circleparams[i+1], radius=circleparams[i+2]/3600, color='white', linestyle='--', linewidth=2.0, zorder=1)
    fig.save(const.SPATH+str(name)+'_cont.png', format='png', adjust_bbox='tight', dpi=200)
    