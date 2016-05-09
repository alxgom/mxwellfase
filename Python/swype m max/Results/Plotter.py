# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 06:23:41 2016

@author: Alexis
"""

''' Code to Replot previous run'''


import numpy as np
import matplotlib.pyplot as plt
#from scipy.integrate import odeint
from scipy.signal import argrelextrema
from time import localtime
import time as timing
import os


#%%


m_maxs=''
max_peaks=''


##Read parameters from 'Params.txt', set.
#txtfile=open("Params.txt")
#tempvar=txtfile.readlines()
#
#'''parameters for normalization'''
#a= float(tempvar[1].split()[1])
#gperp= float(tempvar[2].split()[1])  #gamma perpendicular, loss rate
#scale=1*(10.**6)/gperp #scale to micro seconds
#wscale=1000*gperp/(10.**6)#scale frequency to khz
#
#'''parameters for the equation'''
#kr= float(tempvar[3].split()[1])
#k=kr/gperp #normalized loss rate
#mu= float(tempvar[4].split()[1])
#Dphi0= float(tempvar[5].split()[1]) #phase shift [-pi,pi]
#d= float(tempvar[6].split()[1])#detuning
#gr= float(tempvar[7].split()[1])
#g=gr/gperp #*((2*pi)**2) #sigma parallel, normalized loss rate
#D0=a*k/mu #Poblation
#m= float(tempvar[8].split()[1]) #modulation amplitud [0,1]
#wf= float(tempvar[9].split()[1])
#print 'Params: ', a,gperp,kr,mu,Dphi0,d,gr,m,wf
#
#'''parameters to compare with the results'''
#w_res=np.sqrt(k*g*((D0*mu/k)-1.))*wscale #resonance frequency
#a=D0*mu/k
#w=np.sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale #Relaxation oscilations frequency
#wf_real=wf*wscale
#
#txtfile.close()#deberia cerrarlo aca?
#



with open('Status.txt','w') as f:
        f.write('Status: 0') #ends the run, sets status back to 0.


'''plot swype'''
m_peks=np.fromfile(m_maxs,dtype=np.float64)
max_peak=np.fromfile(max_peaks,dtype=np.float64)
print np.shape(m_peks), np.shape(max_peak)
fig5=plt.figure()
fig5.suptitle('Local Max of intensity vs. m', fontsize=12, fontweight='bold')
ax1 = fig5.add_subplot(111)
ax1.plot(m_peks,max_peak,',b')
ax1.set_xlabel('m')
ax1.set_ylabel('Local Max. Intensity |E|')
ax1.set_xlim(np.min(m_peks), np.max(m_peks))
ax1.set_ylim(0, 2.5)
#plt.text(-0.1,-.32, "\n Parameters: $wf= $ %s , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (wf*wscale, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
plt.subplots_adjust(bottom=0.22)
fig5.set_size_inches(7,7)
#if save==True: 
#    fname='%d_%d_%d-%d.%d.%d-Max_vs_m (w=%s, m=(%s - %s)).png' % tuple(list(localtime()[0:6])+[wf*wscale, m_peks[0], m_peks[-1]])
#    fig5.savefig('Results/'+fname, dpi = 120)# when saving, specify the DPI
#plt.show