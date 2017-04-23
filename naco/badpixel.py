import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from multiprocessing import Pool
from .misc import read_fits
from .params import *

def make_masterflat_bpm(flatdir,fnmflat,fnbpm,sigma=3):
    '''Using median of images to sort in increasing or decreasing intensity. Then clip the
    Values which are of. Using for loops. So slow. Returns a matrix which is 1 whre a bp is found.'''
    print('Making BPM from files in {}'.format(flatdir))
    fns,flats,header = read_fits(flatdir,only_last=True)
    findices = []
    for iim,head in enumerate(header):
        if 'FLAT' in head[hobstype].upper():
            findices.append(iim)
    flats = flats[findices,:,:]
    sortord = np.argsort(np.median(flats,axis=(1,2)))
    bpm = np.full(flats.shape[1:],0,dtype=np.int64)
    gradmap = np.full(flats.shape[1:],np.nan,dtype=np.float64)
    nflats = len(sortord)
    #sort flats by sum of their entries (e.g. darkest first)
#    sums = []
#    for ii in range(nflats):
#        sums.append(np.sum(flats[ii,::]))
#    sumsorted = np.argsort(sums)
    flats = flats[sortord,::]
    iflats = np.linspace(0,nflats-1,nflats).astype(int)
    for yy in range(flats.shape[1]):
        if yy%100 ==0: print('At row %d of %d'%(yy,flats.shape[1]))
        for xx in range(flats.shape[2]):
            values = sigma_clip(flats[:,yy,xx],sigma=5) #filter cosmic rays
            gradmap[yy,xx] = np.polyfit(iflats[~values.mask],\
                                        values[~values.mask],deg=1)[1]
    bpm[sigma_clip(gradmap,sigma=sigma).mask] = 1
    masterflat = gradmap/np.median(gradmap[bpm == 1])
    bpm[masterflat <= 0] = 1
    masterflat = gradmap/np.median(gradmap[bpm == 1])
    fits.writeto(fnmflat,masterflat,overwrite=True)
    fits.writeto(fnbpm,bpm,overwrite=True)
    return masterflat,bpm
    
    
