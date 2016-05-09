# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 19:11:29 2016

@author: Alexis
"""

import numpy as np
#import scipy as sc
#import matplotlib
import matplotlib.pyplot as plt
#from scipy.integrate import odeint
#from numpy import linspace
#import math
#from matplotlib.collections import LineCollection
from intmbfase import intmbfase as integ
from intfaserev import intfaserev as intrev

#import datetime as dt
from time import localtime

'''parameters for normalization'''
a=2
gperp=10**8. #gamma perpendicular, loss rate
scale=1*(10.**6)/gperp #scale to micro seconds
wscale=1000*gperp/(10.**6)#scale frequency to khz

'''parameters for the equation'''
k=0.9*10.**7/gperp #normalized loss rate
mu=.25*10**(-4) #g
Dphi0=0.0 #phase shift [-pi,pi]
d=1.0 #detuning
g=2.5*10.**4/gperp #*((2*pi)**2) #sigma parallel, normalized loss rate
D0=a*k/mu #Poblation
m=.02 #modulation amplitud [0,1]
wf=wf=0.00370    #0.0026089#0.000474340.001*2*pi #modulation frequency

'''parameters to compare with the results'''
w_res=np.sqrt(k*g*((D0*mu/k)-1.))*wscale #resonance frequency
a=D0*mu/k
w=np.sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale #Relaxation oscilations frequency
wf_real=wf*wscale

'''initial conditions'''
intime=500.*7*10**(-6)*gperp #integration time

def  initial():
    initial=raw_input('Initial contitions: \n  Write "new" to use new defined initial conditions.  \n write "l" to use last result as initial condition. \n yinit: ')
    if initial=='new':
        '''User defined initial condition'''
        timeinit = np.arange(0., intime, 0.05)
        dfxinit=[-80., 80.] 
        dfyinit=[-80.,  -8.9]  
        drxinit=[-1.,   1.]
        dryinit=[-1.34720367e+02 ,  1.31585718e+02]
        ddeltainit=[0.65973518e+03]
        yinit=np.array(dfxinit+dfyinit+drxinit+dryinit+ddeltainit)
    if initial=='l':
        '''initial condition from last simulation'''
        timeinit = np.arange(time[-1] ,intime+time[-1] , 0.5)
        yinit=y[-1]
    if initial=='':
        '''reuse last initial condition '''
    return yinit, timeinit

yinit, time=initial()
yinit1=np.copy(yinit)
yinit1[0]=yinit1[0]+yinit1[0]*10*np.finfo(float).eps
print '%.28f' %yinit[0], '%.28f' %yinit1[0]

'''integration'''
y, time=intrev(yinit,time,k,mu,Dphi0,d,g,D0,m,wf)

print 'Last result to use as new initial condition: \n', 'dfxinit=[%e, %e] \n dfyinit=[%e,  %e] \n drxinit=[%e, %e] \n dryinit=[%e,  %e] \n  ddeltainit=[%e] ' %(y[-1][0], y[-1][1], y[-1][2], y[-1][3], y[-1][4], y[-1][5], y[-1][6], y[-1][7], y[-1][8])

'''intensitys'''
intensity_ex=np.sqrt(y[:,0]**2+y[:,1]**2)
intensity_ey=np.sqrt(y[:,2]**2+y[:,3]**2)
intensity=np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2)

'''plots'''
save=False #set True if i want to save files automatically


fig0=plt.figure()
fig0.suptitle('Temporal evolution of |E|, |P| and the Population', fontsize=12, fontweight='bold')
ax1 = fig0.add_subplot(3, 1, 1)
ax1.plot(time*scale,intensity)
ax1.set_ylabel('|E| ')
ax1.set_xlim(min(time*scale), max(time*scale))
ax2 = fig0.add_subplot(3, 1, 2)
ax2.set_xlim(min(time*scale), max(time*scale))
ax2.plot(time*scale,intensity)
ax2.set_ylabel('|P|')
ax3 = fig0.add_subplot(3, 1, 3)
ax3.plot(time*scale,y[:,8])
plt.ylim(min(y[:,8]), max(y[:,8]))
ax3.set_ylabel('Population')
ax3.set_xlim(min(time*scale), max(time*scale))
plt.xlabel('Time($\mu s$)')
plt.text(-0.09,-1.05, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)
fig0.set_size_inches(11, 7)

if save==True: 
    fname='%d_%d_%d-%d.%d.%d-Time_series.png' % localtime()[0:6]
    fig0.savefig(fname) 

#mpld3.display(fig0)


