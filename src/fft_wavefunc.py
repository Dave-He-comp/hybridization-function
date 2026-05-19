#!/usr/bin/env python
# coding: utf-8


import numpy as np
import copy
import os
import matplotlib.pyplot as plt
import time
import scipy.interpolate as interpolate
import Fortran_read



# Band energy window. If the chosen_sys(env)_bands file is not available, 
# the program will read the bands within [band_width_for_***_i, band_width_for_***_f].
band_width_for_sys_i = -0.2
band_width_for_tip_i = -0.2
band_width_for_sub_i = -0.2
#
band_width_for_sys_f = 0.1
band_width_for_tip_f = 0.1
band_width_for_sub_f = 0.1


# Read the necessary VASP calculation parameters from WAVECAR_***,
# including the number of atoms, number of bands, eigenvalues,
# lattice vectors, reciprocal lattice vectors, etc.
# in_title_wavecar: input file name
# in_title: output file name
def read_wave_para(in_title_wavecar,energy_windows_i=-0.4,energy_windows_f=0.1,in_title='wavecar_output'):
    print('read VASP parameters from the file', in_title_wavecar,':\n')
    Fortran_read.read_wavecar(in_title_wavecar,energy_windows_i,energy_windows_f,1)
    #
    self_para = np.zeros(4,dtype=np.int32)
    #
    self_para[0] = int(np.loadtxt(in_title,skiprows=1,max_rows=1))
    self_para[1] = int(np.loadtxt(in_title,skiprows=3,max_rows=1))
    self_para[2] = int(np.loadtxt(in_title,skiprows=5,max_rows=1))
    self_para[3] = int(np.loadtxt(in_title,skiprows=7,max_rows=1))
    #
    print(' number of spin:\n', self_para[0],'\n',\
    'number of k points:\n', self_para[1],'\n',\
    'number of bands:\n', self_para[2],'\n',\
    'number of plane wave:\n', self_para[3],'\n')
    #
    self_a = np.zeros((3,3),dtype=np.float64)
    self_a[0] = np.loadtxt(in_title,skiprows=12,max_rows=1)
    self_a[1] = np.loadtxt(in_title,skiprows=14,max_rows=1)
    self_a[2] = np.loadtxt(in_title,skiprows=16,max_rows=1)
    #
    print('real lattice vectors:\n',\
    'ax:', self_a[0],'\n',\
    'ay:', self_a[1],'\n',\
    'az:', self_a[2],'\n')
    #
    self_b = np.zeros((3,3),dtype=np.float64)
    self_b[0] = np.loadtxt(in_title,skiprows=19,max_rows=1)
    self_b[1] = np.loadtxt(in_title,skiprows=21,max_rows=1)
    self_b[2] = np.loadtxt(in_title,skiprows=23,max_rows=1)
    #
    print('recip. lattice vectors:\n',\
    'bx:', self_b[0],'\n',\
    'by:', self_b[1],'\n',\
    'bz:', self_b[2],'\n')
    #
    self_k = np.loadtxt(in_title,skiprows=25,max_rows=1)
    #
    print('k point value:\n',self_k, '\n')
    #
    self_g_mesh = np.loadtxt('G_MAX')
    print('G MAX mesh:\n',self_g_mesh, '\n')
    return self_para, self_a, self_b, self_k, self_g_mesh


# Read the wavefunction parameters from WAVECAR_***
def read_wavecar(in_title):
    self_nplane = int(np.loadtxt(in_title,skiprows=7,max_rows=1))
    #
    self_igall = np.zeros((self_nplane,3),dtype=np.int32)
    self_igall = np.loadtxt(in_title,skiprows=31,usecols=(0,1,2))
    #
    self_coeff = np.zeros(self_nplane,dtype=np.complex128)
    self_coeff = np.loadtxt(in_title,skiprows=31,usecols=(3))\
    + 1.j * np.loadtxt(in_title,skiprows=31,usecols=(4))
    return self_igall, self_coeff


os.system('grep fermi OUTCAR_SYS > fermi_sys')
sys_para, sys_a, sys_b, sys_k, sys_g_mesh = read_wave_para('WAVECAR_SYS')
#
sys_nspin = sys_para[0]
sys_nwk = sys_para[1]
sys_nband = sys_para[2]
sys_nplane = sys_para[3]
np.savetxt('sys_para',sys_para[2:], fmt='%d', delimiter=' ')
#
del sys_para


i = int(np.loadtxt('chosen_sys_bands'))
print(i)
title_real = 'wavefunc-real-spin-1-wk-1-band-' + str(i) +'.vasp'
title_imag = 'wavefunc-imag-spin-1-wk-1-band-' + str(i) +'.vasp'
title_fft = 'wavefunc-fft-spin-1-wk-1-band-' + str(i)
mesh = np.loadtxt(title_real,max_rows=1)
tmp = int(mesh[0]*mesh[1]*mesh[2])
print('mesh:\n',mesh)

Nx = int(mesh[0])
Ny = int(mesh[1])
Nz = int(mesh[2])

lx = sys_a[0,0]
ly = sys_a[1,1]
lz = sys_a[2,2]


wavecar_real_x = np.fft.fftshift(np.fft.fftfreq(Nx,1.0e0/Nx))
wavecar_real_y = np.fft.fftshift(np.fft.fftfreq(Ny,1.0e0/Ny))
wavecar_real_z = np.fft.fftshift(np.fft.fftfreq(Nx,1.0e0/Nz))


# Convert the x, y, and z dimensions of LOCPOT to the corresponding dimensions in WAVECAR
NGx = int(sys_g_mesh[0]) * 2 -1
NGy = int(sys_g_mesh[1]) * 2 -1
NGz = int(sys_g_mesh[2]) * 2 -1


wavecar_fft_x = np.fft.fftshift(np.fft.fftfreq(NGx,1.0e0/NGx))
wavecar_fft_y = np.fft.fftshift(np.fft.fftfreq(NGy,1.0e0/NGy))
wavecar_fft_z = np.fft.fftshift(np.fft.fftfreq(NGz,1.0e0/NGz))

sample_freqx = (NGx-1) / lx
sample_freqy = (NGy-1) / ly
sample_freqz = (NGz-1) / lz


np.savetxt('freqx',wavecar_fft_x,['%.24E'])
np.savetxt('freqy',wavecar_fft_y,['%.24E'])
np.savetxt('freqz',wavecar_fft_z,['%.24E'])


#
print('Band index ',i, 'spin  1')
a = np.loadtxt(title_real,skiprows=1,usecols=0) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=0)
read_coloum = len(a)
print(read_coloum)
if (tmp%5 == 1):
    read_coloum = len(a)-1
    print(tmp%5,read_coloum)
b = np.loadtxt(title_real,skiprows=1,usecols=1,max_rows=read_coloum) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=1,max_rows=read_coloum)
#read_coloum = len(a)
print(read_coloum)
if (tmp%5 == 2 and read_coloum == len(a)):
    read_coloum = len(a)-1
    print(tmp%5,read_coloum)
c = np.loadtxt(title_real,skiprows=1,usecols=2,max_rows=read_coloum) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=2,max_rows=read_coloum) 
#read_coloum = len(a)
print(read_coloum)
if (tmp%5 == 3 and read_coloum == len(a)):
    read_coloum = len(a)-1
    print(tmp%5,read_coloum)
d = np.loadtxt(title_real,skiprows=1,usecols=3,max_rows=read_coloum) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=3,max_rows=read_coloum) 
#read_coloum = len(a)
print(read_coloum)
if (tmp%5 == 4 and read_coloum == len(a)):
    read_coloum = len(a)-1
    print(tmp%5,read_coloum)
e = np.loadtxt(title_real,skiprows=1,usecols=4,max_rows=read_coloum) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=4,max_rows=read_coloum) 
#
print('reading finished')
max_col = max(a.size,b.size,c.size,d.size,e.size)
print('size of a', a.size,'\n',\
      ' size of b', b.size,'\n',\
      ' size of c', c.size,'\n',\
      ' size of d', d.size,'\n',\
      ' size of e', e.size)
print('maxium coloumn is ', max_col)
#
ncount = 0
if (a.size == (max_col-1)):
    a = np.append(a,np.array([0]))
    ncount =  ncount + 1
if (b.size == (max_col-1)):
    b = np.append(b,np.array([0]))
    ncount =  ncount + 1
if (c.size == (max_col-1)):
    c = np.append(c,np.array([0]))
    ncount =  ncount + 1
if (d.size == (max_col-1)):
    d = np.append(d,np.array([0]))
    ncount =  ncount + 1
if (e.size == (max_col-1)):
    e = np.append(e,np.array([0]))
    ncount =  ncount + 1
print('ncount=', ncount)
#
wavefunc = np.hstack((a.reshape(a.size,1), b.reshape(b.size,1), c.reshape(c.size,1), d.reshape(d.size,1), e.reshape(e.size,1)))
wavefunc = wavefunc.reshape(wavefunc.size)
if (ncount != 0):
    wavefunc = wavefunc[0:-ncount]
wavefunc = wavefunc.reshape(int(mesh[2]),int(mesh[1]),int(mesh[0]))
#
print('fft begin')
fft_shifted_wavefunc = np.fft.fftshift(np.fft.fftn(wavefunc))/(sample_freqx*sample_freqy*sample_freqz)
#
np.savetxt(title_fft,fft_shifted_wavefunc.reshape(int(mesh[0]*mesh[1]*mesh[2])),['%.8E    %.8E'])
print('fft done')


title_real = 'wavefunc-real-spin-2-wk-1-band-' + str(i) +'.vasp'
title_imag = 'wavefunc-imag-spin-2-wk-1-band-' + str(i) +'.vasp'
title_fft = 'wavefunc-fft-spin-2-wk-1-band-' + str(i)
print('Band index ',i, 'spin  2')
a = np.loadtxt(title_real,skiprows=1,usecols=0) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=0)
read_coloum = len(a)
print(read_coloum)
if (tmp%5 == 1):
    read_coloum = len(a)-1
    print(tmp%5,read_coloum)
b = np.loadtxt(title_real,skiprows=1,usecols=1,max_rows=read_coloum) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=1,max_rows=read_coloum)
#read_coloum = len(a)
print(read_coloum)
if (tmp%5 == 2 and read_coloum == len(a)):
    read_coloum = len(a)-1
    print(tmp%5,read_coloum)
c = np.loadtxt(title_real,skiprows=1,usecols=2,max_rows=read_coloum) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=2,max_rows=read_coloum) 
#read_coloum = len(a)
print(read_coloum)
if (tmp%5 == 3 and read_coloum == len(a)):
    read_coloum = len(a)-1
    print(tmp%5,read_coloum)
d = np.loadtxt(title_real,skiprows=1,usecols=3,max_rows=read_coloum) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=3,max_rows=read_coloum) 
#read_coloum = len(a)
print(read_coloum)
if (tmp%5 == 4 and read_coloum == len(a)):
    read_coloum = len(a)-1
    print(tmp%5,read_coloum)
e = np.loadtxt(title_real,skiprows=1,usecols=4,max_rows=read_coloum) +\
    1.j * np.loadtxt(title_imag,skiprows=1,usecols=4,max_rows=read_coloum) 
#
print('reading finished')
max_col = max(a.size,b.size,c.size,d.size,e.size)
print('size of a', a.size,'\n',\
      ' size of b', b.size,'\n',\
      ' size of c', c.size,'\n',\
      ' size of d', d.size,'\n',\
      ' size of e', e.size)
print('maxium coloumn is ', max_col)
#
ncount = 0
if (a.size == (max_col-1)):
    a = np.append(a,np.array([0]))
    ncount =  ncount + 1
if (b.size == (max_col-1)):
    b = np.append(b,np.array([0]))
    ncount =  ncount + 1
if (c.size == (max_col-1)):
    c = np.append(c,np.array([0]))
    ncount =  ncount + 1
if (d.size == (max_col-1)):
    d = np.append(d,np.array([0]))
    ncount =  ncount + 1
if (e.size == (max_col-1)):
    e = np.append(e,np.array([0]))
    ncount =  ncount + 1
print('ncount=', ncount)
#
wavefunc = np.hstack((a.reshape(a.size,1), b.reshape(b.size,1), c.reshape(c.size,1), d.reshape(d.size,1), e.reshape(e.size,1)))
wavefunc = wavefunc.reshape(wavefunc.size)
if (ncount != 0):
    wavefunc = wavefunc[0:-ncount]
wavefunc = wavefunc.reshape(int(mesh[2]),int(mesh[1]),int(mesh[0]))
#
print('fft begin')
fft_shifted_wavefunc = np.fft.fftshift(np.fft.fftn(wavefunc))/(sample_freqx*sample_freqy*sample_freqz)
#
np.savetxt(title_fft,fft_shifted_wavefunc.reshape(int(mesh[0]*mesh[1]*mesh[2])),['%.8E    %.8E'])
print('fft done')