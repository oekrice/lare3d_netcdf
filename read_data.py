#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file

data_directory = './Data/'



nx = 96#144#int(paras[1])
ny = 96#144#int(paras[2])
nz = 96#144#int(paras[3])

xs = np.linspace(-200,200, nx+1)
ys = np.linspace(-200,200, ny+1)
zs = np.linspace(0,400, nz+1)

xc = 0.5*(xs[1:] + xs[:-1])
yc = 0.5*(ys[1:] + ys[:-1])
zc = 0.5*(zs[1:] + zs[:-1])

dx = xs[1] - xs[0]
dy = ys[1] - ys[0]
dz = zs[1] - zs[0]

xc = 0.5*(xs[1:] + xs[:-1])
yc = 0.5*(ys[1:] + ys[:-1])
zc = 0.5*(zs[1:] + zs[:-1])

dx = xs[1] - xs[0]
dy = ys[1] - ys[0]
dz = zs[1] - zs[0]


class Grid():
    def __init__(self):
        self.x0 = xs[0]; self.x1 = xs[-1]
        self.y0 = ys[0]; self.y1 = ys[-1]
        self.z0 = zs[0]; self.z1 = zs[-1]
        self.nx = nx ; self.ny = ny; self.nz = nz

i = 0
run_number = 0
nsnaps = 200

for plot_num in range(0,nsnaps,1):

    fname = '%s%04d.nc' % (data_directory, plot_num)

    slice_index = ny//2
    i = plot_num
    wait = 0
    fname = '%s%04d.nc' % (data_directory, i)
    print('Loading file name', fname)
    print('Making plot', i, 'from run number', run_number)

    while not os.path.exists(fname):
        time.sleep(0.5)
    try:
        data = netcdf_file(fname, 'r', mmap=False)
        print('File', fname, 'found')

    except:
        print('File', fname, 'not found')
        continue

    bx = np.zeros((nx+1,ny+2,nz+2))
    by = np.zeros((nx+2,ny+1,nz+2))
    bz = np.zeros((nx+2,ny+2,nz+1))

    bx[:,1:-1,1:-1] = np.swapaxes(data.variables['bx'][:],0,2)
    by[1:-1,:,1:-1] = np.swapaxes(data.variables['by'][:],0,2)
    bz[1:-1,1:-1,:] = np.swapaxes(data.variables['bz'][:],0,2)

    en = np.zeros((nx+2,ny+2,nz+2))
    rho = np.zeros((nx+2,ny+2,nz+2))

    en[1:-1,1:-1,1:-1] = np.swapaxes(data.variables['en'][:],0,2)
    rho[1:-1,1:-1,1:-1] = np.swapaxes(data.variables['rho'][:],0,2)

    pr = rho*en*2/3

    def magfield(bx, by, bz):
        bx1 = 0.5*(bx[1:,slice_index,1:-1] + bx[:-1,slice_index,1:-1])
        by1 = 0.5*(by[1:-1,slice_index,1:-1] + by[1:-1,slice_index,1:-1])
        bz1 = 0.5*(bz[1:-1,slice_index,1:] + bz[1:-1,slice_index,:-1])
        return 0.5*(bx1**2 + by1**2+ bz1**2)

    data.close()

    beta = pr[1:-1,slice_index,1:-1].T/magfield(bx,by,bz).T

    if True:
        fig, axs = plt.subplots(2,3, figsize = (10,4))

        im = axs[0,0].pcolormesh(xc,zs,bx[1:-1,slice_index,1:-1].T)
        plt.colorbar(im, ax=axs[0,0])
        axs[0,0].set_title('Bx')
        
        im = axs[0,1].pcolormesh(xs,zs,by[1:-1,slice_index,1:-1].T)
        plt.colorbar(im, ax=axs[0,1])
        axs[0,1].set_title('By')
        
        im = axs[0,2].pcolormesh(xs,zc,bz[1:-1,slice_index,1:-1].T)
        plt.colorbar(im, ax=axs[0,2])
        axs[0,2].set_title('Bz')
        
        im = axs[1,0].pcolormesh(xs,zs,en[1:-1,slice_index,1:-1].T)
        plt.colorbar(im, ax=axs[1,0])
        axs[1,0].set_title('Energy')

        im = axs[1,1].pcolormesh(xs,zs,np.log(rho[1:-1,slice_index,1:-1].T))
        plt.colorbar(im, ax=axs[1,1])
        axs[1,1].set_title('Density (log)')


        im = axs[1,2].pcolormesh(xs,zs,np.log(pr[1:-1,slice_index,1:-1].T/magfield(bx,by,bz).T))
        plt.colorbar(im, ax=axs[1,2])
        axs[1,2].set_title('Plasma Beta (log)')

        plt.tight_layout()
        #plt.show()
        plt.savefig('plots/a%04d' % plot_num)
        plt.close()
