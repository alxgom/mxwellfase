# -*- coding: utf-8 -*-
"""
Created on Fri Dec 04 17:14:22 2015

@author: Alexis
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 17:05:30 2015

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

plt.ion()
pi=np.pi #defino pi

'''parameters'''
a=2
gperp=1000000000. #gamma perpendicular
scale=1*1000000./gperp #scale to micro seconds
wscale=1000*gperp#scale frequency to khz

k=0.01
mu=20000. 
Dphi0=0. #phase shift
d=0.000001 #detuning
g=0.00001 #*((2*pi)**2) #sigma parallel
D0=k*a/mu #Poblation0
m=1. #modulation amplitud [-pi,pi]
wf=16 #modulation frequency

w_res=np.sqrt(k*g*((D0*mu/k)-1.)) #resonance frequency
a=D0*mu/k
w=np.sqrt(k*g*(a-1.)-(g*g*a*a)/4) #Relaxation oscilations frequency
wf_real=wf

'''equation'''
def mb(y, t):
    """ y[0],y[1] Electric field x. y[2],y[3] Electric field y, y[4],y[5]  polarization x, y[6],y[7]  polarization y, y[8] population. """
    dfxr=-k*y[0]+mu*y[4]
    dfxi=-k*y[1]+mu*y[5]
    dfyr=-k*y[2]+mu*y[6]-y[3]*(Dphi0+m*np.cos(wf*t))
    dfyi=-k*y[3]+mu*y[7]+y[2]*(Dphi0+m*np.cos(wf*t))
    drxr=-(1*y[4]-d*y[5])+y[0]*y[8]
    drxi=-(1*y[5]+d*y[4])+y[1]*y[8]
    dryr=-(1*y[6]-d*y[7])+y[2]*y[8]
    dryi=-(1*y[7]+d*y[6])+y[3]*y[8]
    ddelta=-g*(y[8]-D0+(y[0]*y[4]+y[1]*y[5]+y[2]*y[6]+y[3]*y[7]))
    return [dfxr,dfxi,dfyr,dfyi,drxr,drxi,dryr,dryi,ddelta]

'''initial conditions'''
time = np.arange(0.0001, 15, 0.001)
dfxinit=[.1, .1] 
dfyinit=[.1, .1]  
drxinit=[.1, .1]
dryinit=[.01, .01]
ddeltainit=[.04]
yinit=np.array(dfxinit+dfyinit+drxinit+dryinit+ddeltainit)

print mb(yinit, time[1])
y = odeint(mb, yinit, time)

#fig0=plt.figure()
#plt.text(0,0,"\n Parameters: $m= $ %s , $w_{forced}$= %s , $\Delta \phi_0=$ %s , \n $k$=%s , $\mu'=$ %s , $\delta= $ %s , $\gamma_{parallel}=$ %s , $D_0=$ %s" % (m,wf, Dphi0,k,mu, d, g, D0), fontsize=11)

'''plots'''

fig0=plt.figure()
fig0.suptitle('Temporal evolution of |E|, |P| and the Population', fontsize=12, fontweight='bold')
ax1 = fig0.add_subplot(3, 1, 1)
f1=plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2))
ax1.set_ylabel('|E| ')
ax2 = fig0.add_subplot(3, 1, 2)
f2=plt.plot(time,np.sqrt(y[:,4]**2+y[:,5]**2+y[:,6]**2+y[:,7]**2))
ax2.set_ylabel('|P|')
ax3 = fig0.add_subplot(3, 1, 3)
f3=plt.plot(time,y[:,8])
plt.ylim(min(y[:,8]), max(y[:,8]))
ax3.set_ylabel('Population')
plt.xlabel('Time')
plt.text(-0.09,-1.05, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f  , $\Delta \phi_0=$ %s ,  $\Omega=$ %f  \n $k$=%.2f , $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)
fig0.savefig('maxt_reales_temporales_fasemod.png')


fig4=plt.figure()
fig4.suptitle('Comparison between |E| and the Modulation', fontsize=12, fontweight='bold')
ax1 = fig4.add_subplot(3, 1, 1)
f1=plt.plot(time,np.cos(wf*time))
plt.xlim(min(time), min(time)+(10*2*pi/wf))
ax1.set_ylabel('Modulation ')
ax2 = fig4.add_subplot(3, 1, 2)
f2=plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2))
ax2.set_ylabel('|E|')
plt.xlim(min(time), min(time)+10*2*pi/wf)
ax3 = fig4.add_subplot(3, 1, 3)
plt.xlim([max(time)-10*2*pi/wf, max(time)])
f3=plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2))
ax3.set_ylabel('|E|')
plt.text(-0.09,-1.05,  "\n Parameters: $m= $ %s , $w_{mod}$= %.2f  , $\Delta \phi_0=$ %s ,  $\Omega=$ %f  \n $k$=%.2f , $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)
plt.xlabel('Time')


fig4=plt.figure()
fig4.suptitle('Comparison between Re(E_y) and the Modulation', fontsize=12, fontweight='bold')
ax1 = fig4.add_subplot(3, 1, 1)
f1=plt.plot(time,np.cos(wf*time))
plt.xlim(min(time), min(time)+(15*2*pi/wf))
ax1.set_ylabel('Modulation ')
ax2 = fig4.add_subplot(3, 1, 2)
f2=plt.plot(time,y[:,2])
ax2.set_ylabel('$Re(E_y)$')
plt.xlim(min(time), min(time)+15*2*pi/wf)
ax3 = fig4.add_subplot(3, 1, 3)
plt.xlim([max(time)-10*2*pi/wf, max(time)])
f3=plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2))
ax3.set_ylabel('|$E_y$|')
plt.text(-0.09,-1.05,  "\n Parameters: $m= $ %s , $w_{mod}$= %.2f  , $\Delta \phi_0=$ %s ,  $\Omega=$ %f  \n $k$=%.2f , $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)
plt.xlabel('Time')
fig4.savefig('modulo_comparadion_fasemod.png')


fig1=plt.figure()
fig1.suptitle('|E| vs tiempo', fontsize=12, fontweight='bold')
ax1 = fig1.add_subplot(111)
f1=plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2),'b')
f2=plt.plot(time,y[:,0],'g')
f3=plt.plot(time,y[:,2],'r')
ax1.set_xlabel('time ')
ax1.set_ylabel('|E| ')
plt.xlim(min(time), max(time))
plt.text(-0.1,-1.1, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f  , $\Delta \phi_0=$ %s ,  $\Omega=$ %f  \n $k$=%.2f , $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f \n |E|(blue), Re($E_x$) (Green), Re($E_y$) (Red)" % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)
plt.legend()
fig1.savefig('moduloE_fasemodulada.png')


fig5=plt.figure()
fig5.suptitle('E_y vs tiempo', fontsize=12, fontweight='bold')
ax1 = fig5.add_subplot(111)
f1=plt.plot(time,y[:,3],'b')
f1=plt.plot(time,y[:,2],'g')
f1=plt.plot(time,np.sqrt(y[:,2]**2+y[:,3]**2),'k')
ax1.set_xlabel('time')
ax1.set_ylabel('$E_y$')
plt.xlim(min(time), max(time))
plt.text(-0.1,-1.04,  "\n Parameters: $m= $ %s , $w_{mod}$= %.2f  , $\Delta \phi_0=$ %s ,  $\Omega=$ %f  \n $k$=%.2f , $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f  \n $|E_y|$ (black dashed), $Re(E_y)$(green) $Im(E_y)$(blue)" % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)
fig1.savefig('moduloE_fasemodulada.png')


fig2=plt.figure()
fig2.suptitle('|E| vs population', fontsize=12, fontweight='bold')
ax2 = fig2.add_subplot(111)
plt.plot(np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2),y[:,8])
ax2.set_xlabel('|E|')
ax2.set_ylabel('Population')
plt.text(-0.09,-1.05,  "\n Parameters: $m= $ %s , $w_{mod}$= %.2f  , $\Delta \phi_0=$ %s ,  $\Omega=$ %f  \n $k$=%.2f , $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)


fig2=plt.figure()
fig2.suptitle('|P| vs population', fontsize=12, fontweight='bold')
ax2 = fig2.add_subplot(111)
plt.plot(np.sqrt(y[:,4]**2+y[:,5]**2+y[:,6]**2+y[:,7]**2),y[:,8])
ax2.set_xlabel('|P| ')
ax2.set_ylabel('Population')
plt.ylim(min(y[:,8]), max(y[:,8]))
plt.text(-0.09,-1.05,  "\n Parameters: $m= $ %s , $w_{mod}$= %.2f  , $\Delta \phi_0=$ %s ,  $\Omega=$ %f  \n $k$=%.2f , $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)


#fig5=plt.figure()
#ax1 = fig5.add_subplot(111)
#f1=plt.plot(time,y[:,0])
#f1=plt.plot(time,y[:,2])
#fig5.suptitle('E_x  & E_y vs tiempo', fontsize=12, fontweight='bold')
#ax1.set_xlabel('time(ms) ')
#ax1.set_ylabel('E ')
#plt.xlim(min(time), max(time))
#plt.text(-0.1,-1.04, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
#plt.subplots_adjust(bottom=0.22)
#fig1.savefig('moduloE_fasemodulada.png')

print 'wf=', wf ,'\n a=', a, '\n w=', w , '\n w_res=', w_res
