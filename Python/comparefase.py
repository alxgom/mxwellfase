# -*- coding: utf-8 -*-
"""
Created on Fri Dec 04 16:54:23 2015

@author: Alexis
"""
import numpy as np
import scipy as sc
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import linspace
import math
from matplotlib.collections import LineCollection
from intmbfase import intmbfase as int


def comparembfase(y,yinit,time,k,mu,Dphi0,d,g,D0,m,wf):

    plt.ion()
    
    pi=np.pi #defino pi
    gperp=10**8. #gamma perpendicular
    scale=1*10.**6/gperp #scale to micro seconds
    wscale=1000/scale#scale frequency to khz
    w_res=np.sqrt(k*g*((D0*mu/k)-1.))*wscale #resonance frequency
    a=D0*mu/k
    w=np.sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale #Relaxation oscilations frequency
    wf_real=wf*wscale
    #wf2=np.sqrt(k*g*((D0*mu/k)-1.))
    a=D0*mu/k
    
    fig=plt.figure()
    fig.suptitle('E_y vs tiempo', fontsize=12, fontweight='bold')
    ax1 = fig.add_subplot(111)
    
    plt.plot(time*scale,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2),'r',label='|E| w/modulation')
    plt.plot(time*scale,np.sqrt(y[:,2]**2+y[:,3]**2),'--r',dashes=[100,3,15,3], label='|$E_y$| w/modulation (dashed)')
    
    '''integration with m=0, Dphi0=0'''
    y, time=int(yinit,time,k,mu,0.,d,g,D0,0.,wf)
         
    plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2),'b',label='|E| No modulation')
    plt.plot(time,np.sqrt(y[:,2]**2+y[:,3]**2),'--b',dashes=[100,3,15,3], label='|$E_y$| no modulation (dashed)')
    ax1.set_xlabel('time($\mu s$)')
    ax1.set_ylabel('$E$')
    plt.xlim(min(time*scale), max(time*scale))
    plt.text(-0.1,-.33, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
    plt.subplots_adjust(bottom=0.22)
    plt.legend(fontsize = 'small')
    fig.set_size_inches(8, 6)


    
    