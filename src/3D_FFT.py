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
band_width_for_sys_i = -2.0
band_width_for_tip_i = -0.2
band_width_for_sub_i = -2.0
#
band_width_for_sys_f = 2.0
band_width_for_tip_f = 0.1
band_width_for_sub_f = 2.0


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
    print(' number of spin:\n'     ,self_para[0],'\n',\
          'number of k points:\n'  ,self_para[1],'\n',\
          'number of bands:\n'     ,self_para[2],'\n',\
          'number of plane wave:\n',self_para[3],'\n')
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


# Read the dimensions of LOCPOT along the x, y, and z directions
def read_fft_para(in_title_locpot='LOCPOT'):
    Fortran_read.read_locpot(in_title_locpot,1)
    #
    self_fft_mesh = np.loadtxt('FFT_MESH')
    print('FFT mesh:\n',self_fft_mesh)
    return self_fft_mesh


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


os.system('grep fermi OUTCAR_SUB > fermi_sub')
sub_para, sub_a, sub_b, sub_k, sub_g_mesh = read_wave_para('WAVECAR_SUB')
#
sub_nspin = sub_para[0]
sub_nwk = sub_para[1]
sub_nband = sub_para[2]
sub_nplane = sub_para[3]
np.savetxt('sub_para',sub_para[2:], fmt='%d', delimiter=' ')
#
del sub_para


# Check whether the parameters are consistent among different WAVECAR files
print('checking the consistence of different parameters \n')
if (sub_g_mesh.all() == sys_g_mesh.all()):
    print('the G mesh is OK')
    g_mesh = copy.deepcopy(sub_g_mesh)
    del sys_g_mesh, sub_g_mesh
#
if (sub_a.all() == sys_a.all()):
    print('the recip. lattice is OK')
    a = copy.deepcopy(sub_a)
    del sys_a, sub_a
#
if (sub_b.all() == sys_b.all()):
    print('the real lattice is OK')
    b = copy.deepcopy(sub_b)
    del sys_b, sub_b
#
if (sub_k.all() == sys_k.all()):
    print('the k points is OK')
    k = copy.deepcopy(sub_k)
    del sys_k, sub_k
#
if (sub_nwk == sys_nwk):
    print('the number of k points is OK')
    nwk = copy.deepcopy(sub_nwk)
    del sys_nwk, sub_nwk
#
if (sub_nspin == sys_nspin):
    print('the nspin is OK')
    nspin = copy.deepcopy(sub_nspin)
    del sys_nspin, sub_nspin

if (sub_nplane == sys_nplane):
   print('the nplane is OK')


# Read the wavefunctions from WAVECAR_SYS
print('reading WAVECAR of system \n')
Fortran_read.read_wavecar('WAVECAR_SYS',band_width_for_sys_i,band_width_for_sys_f,0)
os.system('mkdir SYS')
os.system('mv spin-* ./SYS/')


# Read the wavefunctions from WAVECAR_SUB
print('reading WAVECAR of substrate \n')
Fortran_read.read_wavecar('WAVECAR_SUB',band_width_for_sub_i,band_width_for_sub_f,0)
os.system('mkdir SUB')
os.system('mv spin-* ./SUB/')