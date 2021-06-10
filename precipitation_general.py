#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 14:28:39 2020

@author: nadir
"""

#from fenics import *
from dolfin import *
#import matplotlib.pyplot as plt
#plt.ion()
import sys, math
import numpy as np
import time

#input parameters in a class format
        
class Bend: #eevolution equation of KN2
    def __init__(self, chitemp = None):
        if chitemp is None:
            self.gamma=Constant(2.)
            self.zeta=Constant(4.)
        else:
            self.chi=Constant(chitemp)
        
class Domain:
        
    def __init__(self, kk = None):
        if kk is None:
            self.min=Constant(1e-8)
            self.max=Constant(1./4.)
        else:
            self.min=Constant(0.)
            self.max=Constant(1./float(kk))

       

def geomseq(x, smin, smax, tau):
    return smin + (smax-smin)*((x-smin)/(smax-smin))**tau
    #return ((smax-smin)*(np.exp(x/tau)-np.exp(smin/tau))
    #        /(np.exp(smax/tau)-np.exp(smin/tau))+smin)



#mean curvature, Gaussian curvature, and speed functions

def H(kn1i, kn2i): #mean curvature
    return (kn1i+kn2i)/2.
def kappaG(kn1i, kn2i, taugi): #Gaussian curvature
    return kn1i*kn2i-taugi**2
def Utemp(kgi, kn1i, kn2i, taugi, param, fmean): #param is a class type Rate
    #param=Rate()
    return (-param.alpha1*kgi+param.alpha2*kgi**2+param.alpha3*kgi**3+param.eta1*fmean(kn1i, kn2i)**2
            +param.eta21*kn1i*kn2i+param.eta22*taugi**2+param.eta3*kgi*fmean(kn1i, kn2i)**2-param.eta41*kgi*kn1i*kn2i
            +param.eta42*kgi*taugi**2+param.eta5*taugi*fmean(kn1i, kn2i)+param.eta6*kgi*taugi*fmean(kn1i, kn2i))

#class definitions for boundary labels

class BoundaryX0(SubDomain):
    def __init__(self, kk=None):
        SubDomain.__init__(self)
        self.kk=kk
        
    def inside(self, x, on_boundary):
        return on_boundary and near(x, float(Domain(self.kk).min), DOLFIN_EPS)
    
    
class BoundaryX1(SubDomain):
    def __init__(self, kk=None):
        SubDomain.__init__(self)
        self.kk=kk
        
    def inside(self, x, on_boundary):
        return on_boundary and near(x, float(Domain(self.kk).max), DOLFIN_EPS)
    
#flux and source functions

def Q2(Ui, ggi):
    return grad(Ui)/ggi

def f2(Ui, ggi, kgi, kn1i, kn2i, taugi, fGaussian):
    return ggi*fGaussian(kn1i, kn2i, taugi)*Ui

'''
Science paper:
def Q3(ggi, kn2i):
    param=Bend()
    return param.gamma*grad(kn2i)/ggi


def f3(Ui, ggi, kgi, kn1i, kn2i, qb):
    param=Bend()
    return param.zeta*ggi*kgi*kn1i*Ui*(kn2i-qb)-ggi*kgi*kn2i*Ui
'''
#alternative closure; source term:
    
def f3(chitemp, Ui, kn1i, kn2i, fmean):
    param=Bend(chitemp)
    return param.chi*Ui*fmean(kn1i, kn2i)


def f4(Ui, ggi, kgi):
    return -ggi*kgi*Ui

def f5(Ui, ggi, kgi, kn1i, kn2i, taugi):
    return Ui*taugi.dx(0)+2*taugi*Ui.dx(0)-ggi*kgi*kn2i*Ui

def f6(Ui, ggi, kgi, kn1i, kn2i, taugi):
    return Ui*kn2i.dx(0)+(kn2i-kn1i)*Ui.dx(0)+ggi*kgi*taugi*Ui

#Backward differentiation formula variables & parameters
    
def BDFparameters(nn, ccn, ccn_1, ccn_2, ccn_3, ccn_4, ccr):
    if nn==0:
        ccn.assign(Constant(-1.)); ccr.assign(Constant(1.)); ccn_1.assign(Constant(0.))
        ccn_2.assign(Constant(0.)); ccn_3.assign(Constant(0.)); ccn_4.assign(Constant(0.))                
    elif nn==1:
        ccn.assign(Constant(-4./3.)); ccr.assign(Constant(2./3.)); ccn_1.assign(Constant(1./3.))
        ccn_2.assign(Constant(0.)); ccn_3.assign(Constant(0.)); ccn_4.assign(Constant(0.))
    elif nn==2:
        ccn.assign(Constant(-18./11.)); ccr.assign(Constant(6./11.)); ccn_1.assign(Constant(9./11.))
        ccn_2.assign(Constant(-2./11.)); ccn_3.assign(Constant(0.)); ccn_4.assign(Constant(0.))
    elif nn==3:
        ccn.assign(Constant(-48./25.)); ccr.assign(Constant(12./25.)); ccn_1.assign(Constant(36./25.))
        ccn_2.assign(Constant(-16./25.)); ccn_3.assign(Constant(3./25.)); ccn_4.assign(Constant(0.))
    else:
        ccn.assign(Constant(-300./137.)); ccr.assign(Constant(60./137.)); ccn_1.assign(Constant(300./137.))
        ccn_2.assign(Constant(-200./137.)); ccn_3.assign(Constant(75./137.)); ccn_4.assign(Constant(-12./137.))

def BDFvariables(nn, u, un, un_1, un_2, un_3, un_4):
    if nn==0:
        un.assign(u)    
    elif nn==1:
        un_1.assign(un); un.assign(u)  
    elif nn==2:
        un_2.assign(un_1); un_1.assign(un); un.assign(u)  
    elif nn==3:
        un_3.assign(un_2); un_2.assign(un_1); un_1.assign(un); un.assign(u)  
    else:
        un_4.assign(un_3); un_3.assign(un_2); un_2.assign(un_1); un_1.assign(un); un.assign(u)  
        