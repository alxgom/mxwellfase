# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 21:43:26 2016

@author: Alexis
"""

import numpy as np
from scipy.integrate import odeint

def intmbfase2(yinit,time,k,mu,Dphi0,d,g,D0,m,wf):

    def mb(y, t):
        """ y[0],y[1] campo electrico en x. y[2],y[3] campo electrico en y, y[4],y[5]  polarizacion en x, y[6],y[7]  polarizacion en y, y[8] poblacion. """
        dfxr=-k*y[0]+mu*y[4]
        dfxi=-k*y[1]+mu*y[5]
        dfyr=-k*y[2]+mu*y[6]-y[3]*(Dphi0+m*np.cos(wf*t))
        dfyi=-k*y[3]+mu*y[7]+y[2]*(Dphi0+m*np.cos(wf*t))
        drxr=-(1*y[4]-d*y[5])+y[0]*y[8]
        drxi=-(1*y[5]+d*y[4])+y[1]*y[8]
        dryr=-(1*y[6]-d*y[7])+y[2]*y[8]
        dryi=-(1*y[7]+d*y[6])+y[3]*y[8]
        ddelta=-g*(y[8]-D0+ 0.5*(y[0]*y[4]+y[1]*y[5]+y[2]*y[6]+y[3]*y[7]))
        return [dfxr,dfxi,dfyr,dfyi,drxr,drxi,dryr,dryi,ddelta]

    
    y = odeint(mb, yinit, time)
    return y, time
