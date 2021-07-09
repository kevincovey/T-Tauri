# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 15:35:14 2021

@author: Table
"""
import pathlib

from master import *

def Wavelength_Flux_Plot(plate,mjd,fiber):
    
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

    if len(str(fiber)) == 2:
        fiber = '0' + str(fiber)
    elif len(str(fiber)) == 1:
        fiber = '00' + str(fiber)
    else:
        fiber = str(fiber)

# sets up reading the csv files in emitters folder only for me (EK)
    #serverpath = '/Users/ellio/Desktop/Research/DR15/Spectra Files/Emitters/'
    #serverpath = '/Users/khilfek/Covey/APOGEE_Spectra/Elliott/DR15/Spectra Files/Emitters/'    

#### NOTE!  KEVIN CHANGED THE PATH!
    serverpath = '/Users/coveyk/Dropbox/python/T-Tauri/Summer_2021/dataFiles/DR15/Spectra Files/Emitters/'    
    
    filepath = serverpath + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
    openfile = pd.read_csv(filepath)

# I think this is making individual lists for each variable
    wave = openfile['Wavelength']
    flux = openfile['Flux']
    a = np.where(flux>-100)[0]
    flux=np.array(flux[a])    #Kevin is not used to pandas, so he converted the flux and wave columns to numpy arrays.
    wave=np.array(wave[a])    #having these elements as numpy arrays allowed Kevin to use np.where and np.percentile below
    
    #define the blue, green, and red chips
    blue = np.where( (wave > 15140) & (wave < 15810) )    #having the wave column as a numpy array made selecting the wavelengths
    green = np.where( (wave > 15855) & (wave < 16435) )   #in each chip much easier.
    red = np.where( (wave > 16470) & (wave < 16955) )

    #find the values that should be used to set the ylims for each panel
    #do this using the np.percentile command to exclude the few highest and lowest pixels in each array,
    #as they often ruin the y scaling
    blue_max = 1.1*np.percentile(flux[blue], 99.5)
    blue_min = 0.9*np.percentile(flux[blue], 0.5)
    
    green_max = 1.1*np.percentile(flux[green], 99.5)
    green_min = 0.9*np.percentile(flux[green], 0.5)

    red_max = 1.1*np.percentile(flux[red], 99.5)
    red_min = 0.9*np.percentile(flux[red], 0.5)

# This makes the plots of flux vs emitter 
    
    fig, (ax, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Flux vs. Wavelength of Emitter ' + str(plate) + '-' + str(mjd) + '-' + str(fiber))
    fig.tight_layout()
    #kevin set the plot size to be large, so it would print well.
    fig.set_size_inches(11, 8)
    
    ax.plot(wave[blue], flux[blue], color='k', linewidth=.4)
    ax2.plot(wave[green], flux[green], color='k', linewidth=.4)
    ax3.plot(wave[red], flux[red], color='k', linewidth=.4)

    ax.set_xlim(15140,15810)
    ax2.set_xlim(15855,16435)
    ax3.set_xlim(16470,16955)
    
    ax.set_ylim(blue_min,blue_max)
    ax2.set_ylim(green_min,green_max)
    ax3.set_ylim(red_min,red_max)
    
    #plt.autoscale(enable=True, axis='y' ,tight='True')
    
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_major_formatter('{x:.0f}')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    
    ax2.xaxis.set_major_locator(MultipleLocator(100))
    ax2.xaxis.set_major_formatter('{x:.0f}')
    ax2.xaxis.set_minor_locator(MultipleLocator(10))
    
    ax3.xaxis.set_major_locator(MultipleLocator(100))
    ax3.xaxis.set_major_formatter('{x:.0f}')
    ax3.xaxis.set_minor_locator(MultipleLocator(10))
    
    #ax.yaxis.set_major_locator(MultipleLocator(100))
    #ax.yaxis.set_major_formatter('{x:.0f}')
    #ax.yaxis.set_minor_locator(MultipleLocator(10))
    
    #ax2.yaxis.set_major_locator(MultipleLocator(100))
    #ax2.yaxis.set_major_formatter('{x:.0f}')
    #ax2.yaxis.set_minor_locator(MultipleLocator(10))
    
   # ax3.yaxis.set_major_locator(MultipleLocator(100))
   # ax3.yaxis.set_major_formatter('{x:.0f}')
   #ax3.yaxis.set_minor_locator(MultipleLocator(10))
  
    plt.xlabel('Wavelength (\u00C5)')
    plt.ylabel('Flux')

    #save a digital copy of the figure as a file; do this in the emitters folder where the csv lives.
    plt.savefig(serverpath+str(plate) + '-' + str(mjd) + '-' + str(fiber)+'.png', dpi=100)

    #close the plot because failing to do so would consume a lot of memory.
    plt.close()
    
    #ax2.xlabel('Wavelength (\u00C5)')
    #ax2.ylabel('Flux (10^-17 erg/s/cm^2/\u00C5)')
    #ax3.xlabel('Wavelength (\u00C5)')
    #ax3.ylabel('Flux (10^-17 erg/s/cm^2/\u00C5)')
    
    #plt.show()


#modify the main script so that it makes plots for everything in the emitters list csv file.

#import pandas so we can open the csv as a panda
import pandas as pd

#read in the emitters list
openfile = pd.read_csv('/Users/coveyk/Dropbox/python/T-Tauri/Summer_2021/dataFiles/DR15/Emitters_List_DR15.csv')

#pull the plate, mjd and fiber columns out of the csv file
plate = openfile['Plate']
mjd = openfile['MJD']
fiber = openfile['Fiber']

#loop through every object in the input list, and use the Wavelength_Flux_Plot function to make a plot of it.
for x in range(len(plate)):
    Wavelength_Flux_Plot(plate[x],mjd[x],fiber[x])
    

 

######## COPY PASTED FROM EQUIVALENT_WIDTHS ########

# for k in range(10):
 #        line = 11 + k
 #        rest = Rest_Emission(line)
 #        elrest = Find_Nearest(wave,rest)
        


 #        L1 = elrest - 241
 #        L2 = elrest - 120
 #        R1 = elrest + 120
 #        R2 = elrest + 241

   

 #        ''' Calculating local continuum via median '''
 #        # Setting up windows
 #        left_win = np.asarray(flux[L1:L2])
 #        right_win = np.asarray(flux[R1:R2])
 #        left_win_err = np.asarray(errors[L1:L2])
 #        right_win_err = np.asarray(errors[R1:R2])
 #        left_snr = np.asarray(snr[L1:L2])
 #        right_snr = np.asarray(snr[R1:R2])

 #        # Only import non-junk data
 #        left_win_err = left_win_err[left_win != -9999]
 #        right_win_err = right_win_err[right_win != -9999]
 #        left_snr = left_snr[left_win != -9999]
 #        right_snr = right_snr[right_win != -9999]
 #        left_win = left_win[left_win != -9999]   
 #        right_win = right_win[right_win != -9999]

 #        # Only import data above an SNR of 10
 #        left_win_err = left_win_err[left_snr > 10]
 #        right_win_err = right_win_err[right_snr > 10]
 #        left_win = left_win[left_snr > 10]   
 #        right_win = right_win[right_snr > 10]
 #        left_snr = left_snr[left_snr > 10]
 #        right_snr = right_snr[right_snr > 10]

 #        # Median
 #        local_range = np.concatenate((left_win,right_win))
 #        local = np.median(local_range)     

##########################################################

def Fits_Info(plate,mjd,fiber,version):
#not sure what "version" is?
    from astropy.io import fits
    from PyAstronomy.pyasl import helcorr

    # Now we stitch together the server path to each file. Careful though as there is a folder for the APO telescope and one for LCO that   
    # contain different field. So we will check against LCO as it contains fewer plates

#lco only has these 7 plates
#confused because I couldn't find any of these plates in apo25m... is the purpose for them to NOT have overlapping data?
    lco = [9753,9761,9762,9856,9906,9907,9908]
    apopath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
    lcopath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/lco25m/'

    int_plate = int(plate)

    # Setting path to file
    if int_plate in lco:
        filepath = lcopath + str(plate) + '/' + str(mjd) + '/asVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
    else:
        filepath = apopath + str(plate) + '/' + str(mjd) + '/apVisit-' + str(version) + '-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
    
    # Now we can open the fits file on the server
    try:
        fitsfile = fits.open(filepath)
        header = fitsfile[0].header
    except:
        print('Bad file path!')

    ''' Get: wavegrid, flux, error, VHELIO, RA, DEC, TELESCOP '''

    # Grabbing the data for the wavelength and flux
#what are these numbers 4,2,1? -EK
    wavedata = fitsfile[4].data
    fluxdata = fitsfile[1].data
    ferrordata = fitsfile[2].data

    # Getting location information and the barycentric correction
    ra = header['RA']
    dec = header['DEC']
    jd = header['JD-MID']
    telescope = header['TELESCOP']

    # Sometimes there is no reported BC so we need to get one from the helcorr package
#is this BC as in barycentric correction? -EK
    try:
        vbc = header['BC']

#again, what are these numbers? -EK
    except:
        if telescope == 'apo25m':
            vbc,hjd = helcorr(-105.4913,36.4649,2788,ra,dec,jd)
        else:
            vbc,hjd = helcorr(-70.413336,-29.05256,2380,ra,dec,jd)

    fitsfile.close()
    
    return(wavedata,fluxdata,ferrordata,vbc,ra,dec)
