# -*- coding: utf-8 -*-
"""
Created on Tue Dec 01 20:00:54 2015

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
import datetime as dt
from time import localtime
import matplotlib.animation as animation

'''---- |2>
    ||
    || R21
    ||
   ---- |1>'''

plt.ion()
pi=np.pi #defino pi

'''parameters for normalization'''
a=2
gperp=10**8. #gamma perpendicular
scale=1*(10.**6)/gperp #scale to micro seconds
wscale=1000*gperp/(10.**6)#scale frequency to khz

'''parameters for the equation'''
k=0.9*10.**7/gperp
mu=.25*10**(-4) 
Dphi0=0.0 #phase shift [-pi,pi]
d=0.0 #detuning
g=2.5*10.**4/gperp #*((2*pi)**2) #sigma parallel
D0=a*k/mu #Poblation
m=.02 #modulation amplitud [0,1]
wf=0.0047434#0.001*2*pi #modulation frequency

'''parameters to compare with the results'''
w_res=np.sqrt(k*g*((D0*mu/k)-1.))*wscale #resonance frequency
a=D0*mu/k
w=np.sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale #Relaxation oscilations frequency
wf_real=wf*wscale

'''initial conditions'''
time = np.arange(0.0001, 5*500*10**(-6)*gperp, 0.5)
dfxinit=[.75, .1] 
dfyinit=[.75, .1]  
drxinit=[2600., .5]
dryinit=[2600., .5]
ddeltainit=[3600.]
yinit=np.array(dfxinit+dfyinit+drxinit+dryinit+ddeltainit)


for n in np.arange(100,1000,100):
    time = np.arange(0.0001, 50*500*10**(-6)*gperp, 0.5)
    y, time=int(yinit,time,k,mu,Dphi0,d,g,a*k/mu,m,wf/n)
    print wf/n
    
    '''intensitys'''
    intensity_ex=np.sqrt(y[:,0]**2+y[:,1]**2)
    intensity_ey=np.sqrt(y[:,2]**2+y[:,3]**2)
    intensity=np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2)    
        
    D0=a*k/mu
    
    fig4=plt.figure()
    periods=5
    fig4.suptitle('Comparison between |E| and the Modulation', fontsize=12, fontweight='bold')
    ax1 = fig4.add_subplot(3, 1, 1)
    ax1.plot(time*scale,np.cos(wf*time), label='modulation')
    ax1.set_xlim(min(time*scale), min(time*scale)+(periods*2*pi/wf)*scale)
    ax1.set_ylabel('Modulation ')
    ax2 = fig4.add_subplot(3, 1, 2)
    ax2.plot(time*scale,intensity,label='First %i periods' %periods)
    ax2.set_ylabel('|E|')
    ax2.set_xlim(min(time*scale), min(time*scale)+(periods*2*pi/wf)*scale)
    ax2.legend(fontsize = 'small')
    ax3 = fig4.add_subplot(3, 1, 3)
    ax3.plot(time*scale,intensity,label='Last %i periods' %periods)
    ax3.set_xlim(max(time*scale)-(periods*2*pi/wf)*scale, max(time*scale))
    ax3.set_ylim(intensity[len(intensity)-1]-8*(max(intensity[len(intensity)-90000:])-intensity[len(intensity)-1]),intensity[len(intensity)-1]+8*(max(intensity[len(intensity)-90000:])-intensity[len(intensity)-1]))
    ax3.set_ylabel('|E|')
    ax3.legend(fontsize = 'small')
    plt.text(-0.09,-1.06, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf*wscale/n, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
    plt.subplots_adjust(bottom=0.22)
    plt.xlabel('Time ($\mu s$)')
    
#    
#def animate(i):
#    time = np.arange(0.0001, 50*500*10**(-6)*gperp, 0.5)
#    y, time=int(yinit,time,k,mu,Dphi0,d,g,a*k/mu,m,wf/n)
#        
#    '''intensitys'''
#    intensity_ex=np.sqrt(y[:,0]**2+y[:,1]**2)
#    intensity_ey=np.sqrt(y[:,2]**2+y[:,3]**2)
#    intensity=np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2)    
#        
#    D0=a*k/mu
#    
#    fig4=plt.figure()
#    periods=5
#    fig4.suptitle('Comparison between |E| and the Modulation', fontsize=12, fontweight='bold')
#    ax1 = fig4.add_subplot(3, 1, 1)
#    ax1.plot(time*scale,np.cos(wf*time), label='$w_{mod}=$' %(wf*wscale/n))
#    ax1.set_xlim(min(time*scale), min(time*scale)+(periods*2*pi/wf)*scale)
#    ax1.set_ylabel('Modulation ')
#    ax2 = fig4.add_subplot(3, 1, 2)
#    ax2.plot(time*scale,intensity,label='First %i periods' %periods)
#    ax2.set_ylabel('|E|')
#    ax2.set_xlim(min(time*scale), min(time*scale)+(periods*2*pi/wf)*scale)
#    ax2.legend(fontsize = 'small')
#    ax3 = fig4.add_subplot(3, 1, 3)
#    ax3.plot(time*scale,intensity,label='Last %i periods' %periods)
#    ax3.set_xlim(max(time*scale)-(periods*2*pi/wf)*scale, max(time*scale))
#    ax3.set_ylim(intensity[len(intensity)-1]-8*(max(intensity[len(intensity)-90000:])-intensity[len(intensity)-1]),intensity[len(intensity)-1]+8*(max(intensity[len(intensity)-90000:])-intensity[len(intensity)-1]))
#    ax3.set_ylabel('|E|')
#    ax3.legend(fontsize = 'small')
#    plt.text(-0.09,-1.06, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf*wscale/n, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
#    plt.subplots_adjust(bottom=0.22)
#    plt.xlabel('Time ($\mu s$)') # update the data
#    return line,
#
#
## Init only required for blitting to give a clean slate.
#def init():
#    line.set_ydata(np.ma.array(x, mask=True))
#    return line,
#
#ani = animation.FuncAnimation(fig, animate, np.arange(1, 200), init_func=init,
#                              interval=25, blit=True)
#plt.show()    
    
#    fig0=plt.figure()
#    ax1 = fig0.add_subplot(3, 1, 1)
#    f1=plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2))
#    ax1.set_ylabel('|E|')
#    ax2 = fig0.add_subplot(3, 1, 2)
#    f2=plt.plot(time,np.sqrt(y[:,4]**2+y[:,5]**2+y[:,6]**2+y[:,7]**2))
#    ax2.set_ylabel('|P|')
#    ax3 = fig0.add_subplot(3, 1, 3)
#    f3=plt.plot(time,y[:,8])
#    plt.ylim(min(y[:,8]), max(y[:,8]))
#    ax3.set_ylabel('Population')
#    fig0.suptitle('Temporal evolution of |E|, |P| and the Population', fontsize=12, fontweight='bold')
#    plt.xlabel('Time(ms)')
#    plt.text(-0.09,-1.05, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf*wscale/n, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
#    plt.subplots_adjust(bottom=0.22)  
#    #plt.tight_layout()
#    plt.subplots_adjust(bottom=0.22)
#    plt.show()
    #fig0.savefig('maxt_reales_temporales_fasemod.png')
#
#fig4=plt.figure()
#fig4.suptitle('Comparison between |E| and the Modulation', fontsize=12, fontweight='bold')
#ax1 = fig4.add_subplot(3, 1, 1)
#f1=plt.plot(time,m*np.cos(wf*time))
#plt.xlim(min(time), min(time)+10*2*pi/wf)
##real del campo vs real polarizacion
#ax1.set_ylabel('Modulation ')
#ax2 = fig4.add_subplot(3, 1, 2)
#f2=plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2))
#ax2.set_ylabel('|E|')
#plt.xlim(min(time), min(time)+10*2*pi/wf)
#ax3 = fig4.add_subplot(3, 1, 3)
#plt.xlim([max(time)-10*2*pi/wf, max(time)])
#f3=plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2))
#ax3.set_ylabel('|E|')
#plt.text(-0.09,-1.05, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $w_{rel}=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf, Dphi0, w ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
#plt.subplots_adjust(bottom=0.22)
##plt.axes('tight')
#plt.xlabel('Time(ms)')
#fig4.savefig('modulo_comparadion_fasemod.png')
#
#
#fig1=plt.figure()
#ax1 = fig1.add_subplot(111)
#f1=plt.plot(time,np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2))
#fig1.suptitle('modulo E vs tiempo', fontsize=12, fontweight='bold')
#ax1.set_xlabel('time(ms) ')
#ax1.set_ylabel('|E| ')
#plt.xlim(min(time), max(time))
#plt.text(-0.1,-0.3, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $w_{rel}=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf, Dphi0, w ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
#plt.subplots_adjust(bottom=0.22)
#fig1.savefig('moduloE_fasemodulada.png')
#
#
#fig2=plt.figure()
#ax2 = fig2.add_subplot(111)
#plt.plot(np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2),y[:,8])#real del campo vs real poblacion
#fig2.suptitle('|E| vs population', fontsize=12, fontweight='bold')
#ax2.set_xlabel('|E|')
#ax2.set_ylabel('Population')
#plt.text(-0.09,-1.05, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $w_{rel}=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf, Dphi0, w ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
#plt.subplots_adjust(bottom=0.22)


print 'wf=', wf ,'\n a=', a, '\n w=', w #, '\n wf2=', wf2



