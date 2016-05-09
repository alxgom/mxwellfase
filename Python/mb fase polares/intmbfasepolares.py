# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 23:29:37 2015

@author: Alexis
"""

"""
Created on Tue Dec 01 19:56:47 2015

@author: Alexis
"""

import numpy as np
from scipy.integrate import odeint

def intmbfase(yinit,time,k,g,dw,D0,gama):

    pi=np.pi

    def mbpolares(y, t):
        """ y[0],y[1] campo electrico en x. y[2],y[3] campo electrico en y, y[4],y[5]  polarizacion en x, y[6],y[7]  polarizacion en y, y[8] poblacion. """
        dr=-k*y[0]+g*y[0]*y[2]/((1+dw**2))
        dphi=-g*dw*y[2]/((1+dw**2))
        dde=-gama*(y[2]-D0+y[0]*y[2]/(1+dw**2))
        return [dr,dphi,dde]
        
        

    y = odeint(mbpolares, yinit, time)
    return y, time
