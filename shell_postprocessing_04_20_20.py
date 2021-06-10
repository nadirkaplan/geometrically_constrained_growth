#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 16:18:53 2020

@author: nadir
"""

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from IPython import get_ipython


import sys, math
import numpy as np


shape='coral'
kk=6
d3 = '21_05_23'

theta=-2.*np.pi/float(kk)

rotmatrix=[np.matrix([[np.cos(theta*ii), -np.sin(theta*ii)], [np.sin(theta*ii), np.cos(theta*ii)]]) for ii in range(1, int(kk/2))]

filenamemesh=shape+'_mesh_kk_'+str(kk)+'_'+d3+'.txt'
filenameX=shape+'_Xdata_kk_'+str(kk)+'_'+d3+'.txt'
filenameY=shape+'_Ydata_kk_'+str(kk)+'_'+d3+'.txt'
filenameZ=shape+'_Zdata_kk_'+str(kk)+'_'+d3+'.txt'
filenametime=shape+'_timearray_kk_'+str(kk)+'_'+d3+'.txt'

filenamekg=shape+'_kappag_kk_'+str(kk)+'_'+d3+'.txt'
filenameU=shape+'_U_kk_'+str(kk)+'_'+d3+'.txt'
filenamekN2=shape+'_kappaN2_kk_'+str(kk)+'_'+d3+'.txt'
filenamegg=shape+'_gg_kk_'+str(kk)+'_'+d3+'.txt'
filenamekN1=shape+'_kappaN1_kk_'+str(kk)+'_'+d3+'.txt'
filenametaug=shape+'_taug_kk_'+str(kk)+'_'+d3+'.txt'


coordinates=np.loadtxt(filenamemesh)
X1array=np.loadtxt(filenameX, delimiter='\t')
X2array=np.loadtxt(filenameY, delimiter='\t')
X3array=np.loadtxt(filenameZ, delimiter='\t')
timearray=np.loadtxt(filenametime, delimiter='\t')


kappag_array=np.loadtxt(filenamekg, delimiter='\t')
U_array=np.loadtxt(filenameU, delimiter='\t')
kN2_array=np.loadtxt(filenamekN2, delimiter='\t')
gg_array=np.loadtxt(filenamegg, delimiter='\t')
kN1_array=np.loadtxt(filenamekN1, delimiter='\t')
taug_array=np.loadtxt(filenametaug, delimiter='\t')

ntime=np.max(np.nonzero(timearray))+1
nmesh=coordinates.shape[0];

aatemp=[[np.matmul(rotmatrix[ii], [X1array[:, jj], X3array[:, jj]]) for ii in range(int(kk/2)-1)] for jj in range(ntime)]
aatemp2=np.asarray(aatemp)

X1arrayF=np.zeros((nmesh*int(kk/2), ntime))
X2arrayF=np.zeros((nmesh*int(kk/2), ntime))
X3arrayF=np.zeros((nmesh*int(kk/2), ntime))
#X1arrayF=[np.append(X1arrayF, [X1array[:, jj], [aatemp[jj][ii][0, :] for ii in range(int(kk/2)-1)]]) for jj in range(ntime)]
#X1arrayF=np.append(X1arrayF, [[X1array[:, jj], [aatemp[jj][ii][0, :] for ii in range(int(kk/2)-1)]] for jj in range(ntime)])
#X1arrayF=[np.append(X1arrayF, [aatemp[jj][ii][0, :] for ii in range(int(kk/2)-1)]) for jj in range(ntime)]


X1arrayF[0:nmesh, :]=X1array[0:nmesh, :]
X2arrayF[0:nmesh, :]=X2array[0:nmesh, :]
X3arrayF[0:nmesh, :]=X3array[0:nmesh, :]

for ii in range(1, int(kk/2)):
    print(ii)
    X1arrayF[ii*nmesh:(ii+1)*nmesh, :]=aatemp2[:, ii-1, 0, :].transpose()
    X2arrayF[ii*nmesh:(ii+1)*nmesh, :]=X2array[0:nmesh, :]
    X3arrayF[ii*nmesh:(ii+1)*nmesh, :]=aatemp2[:, ii-1, 1, :].transpose()
    
Hmean=np.zeros((nmesh, ntime))
    
#mean curvature
for jj in range(ntime):
    
    Hmean[:, jj]=np.sum([1./2.*kN1_array[:,jj], 1./2.*kN2_array[:, jj]], axis=0)
    #Hmean[:,jj]=1./2.*Hmean[:, jj]
    
    
#X1array=np.append(X1array, [aatemp2[:, 1-1, 0, :].transpose()])
# [ (X1arrayF[ii*(nmesh+1):nmesh*(ii+1)+ii, :]=aatemp2[:, ii-1, 0, :].transpose()) for ii in range(1, int(kk/2))]
#=aatemp2[:, ii-1, 0, :].transpose()
#=aatemp[:][ii-1][0, :]
#X1arrayF[1*(nmesh+1):nmesh*(1+1)+1, :]=aatemp2[:, 1-1, 0, :].transpose()

#get_ipython().run_line_magic('matplotlib', 'qt')


#ntime=50

#fig=plt.figure()
plt.figure(0)
plt.ion()
ax= plt.axes(projection='3d')
[ax.plot3D(X1arrayF[:, jj], X2arrayF[:, jj], X3arrayF[:, jj]) for jj in range(ntime)]
plt.show()


plt.figure(1)
plt.ion()
[plt.plot(X1array[:, jj], X2array[:, jj]) for jj in range(ntime)]
plt.axes().set_aspect(aspect=1)
plt.show()


plt.figure(2)
plt.ion()
[plt.plot(X1array[:, jj], X3array[:, jj]) for jj in range(ntime)]
plt.axes().set_aspect(aspect=1)
plt.show()


#ntime=1

#get_ipython().run_line_magic('matplotlib', 'inline')

plt.figure(3)
plt.ioff()

plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\kappa_g$')

[plt.plot(coordinates, kappag_array[:, jj]) for jj in range(ntime)]
plt.show()

plt.figure(4)

plt.xlabel(r'$\sigma$')
plt.ylabel(r'$U$')
[plt.plot(coordinates, U_array[:, jj]) for jj in range(ntime)]
plt.show()

plt.figure(5)

plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\kappa_{N,2}$')
[plt.plot(coordinates, kN2_array[:, jj]) for jj in range(ntime)]
plt.show()

plt.figure(6)

plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\sqrt{g}$')
[plt.plot(coordinates, gg_array[:, jj]) for jj in range(ntime)]
plt.show()

plt.figure(7)

plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\kappa_{N}$')
[plt.plot(coordinates, kN1_array[:, jj]) for jj in range(ntime)]
plt.show()

plt.figure(8)

plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\tau_g$')
[plt.plot(coordinates[:], taug_array[:, jj]) for jj in range(ntime)]
plt.show()
#nmesh-10:nmesh

plt.figure(9)

plt.xlabel(r'$\sigma$')
plt.ylabel('mean curvature')
[plt.plot(coordinates[:], Hmean[:, jj]) for jj in range(ntime)]
plt.show()








#[plt.plot(mesh.coordinates(), uarray[jj, 1, :]) for jj in range(ntime)]
#plt.show()

#[plt.plot(qarray[jj, 0, :], qarray[jj, 2, :]) for jj in range(ntime)]
#plt.show()

