# -*- coding: utf-8 -*-
"""
Created on Fri Dec 04 23:29:00 2015

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
from comparefase import comparembfase
from intmbfase import intmbfase as int

'''---- |2>
    ||
    || R21
    ||
   ---- |1>'''

plt.ion()
pi=np.pi #defino pi

'''parameters'''
a=2
gperp=10**8. #gamma perpendicular
scale=1*(10.**6)/gperp #scale to micro seconds
wscale=1000*gperp/(10.**6)#scale frequency to khz

k=0.9*10.**7/gperp
mu=.25*10**(-4) 
Dphi0=0.0 #phase shift [-pi,pi]
d=0.0 #detuning
g=2.5*10.**4/gperp #*((2*pi)**2) #sigma parallel
D0=a*k/mu #Poblation
m=.02 #modulation amplitud [0,1]
wf=0.001*2*pi #modulation frequency

w_res=np.sqrt(k*g*((D0*mu/k)-1.))*wscale #resonance frequency
a=D0*mu/k
w=np.sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale #Relaxation oscilations frequency
wf_real=wf*wscale

'''initial conditions'''
time = np.arange(0.0001, 5000*10**(-6)*gperp, 0.01)
dfxinit=[.75, .1] 
dfyinit=[.75, .1]  
drxinit=[2600., .5]
dryinit=[2600., .5]
ddeltainit=[3600.]
yinit=np.array(dfxinit+dfyinit+drxinit+dryinit+ddeltainit)


'''integration'''
y, time=int(yinit,time,k,mu,Dphi0,d,g,D0,m,wf)

'''comparison between the solution with and without modulation'''
def comp():
    comparembfase(yinit,time,k,mu,Dphi0,d,g,D0,m,wf)

'''intensitys'''
intensity_ex=np.sqrt(y[:,0]**2+y[:,1]**2)
intensity_ey=np.sqrt(y[:,2]**2+y[:,3]**2)
intensity=np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2)
#fig0=plt.figure()
#plt.text(0,0,"\n Parameters: $m= $ %s , $w_{forced}$= %s , $\Delta \phi_0=$ %s , \n $k$=%s , $\mu'=$ %s , $\delta= $ %s , $\gamma_{parallel}=$ %s , $D_0=$ %s" % (m,wf, Dphi0,k,mu, d, g, D0), fontsize=11)

'''plots'''

fig0=plt.figure()
fig0.suptitle('Temporal evolution of |E|, |P| and the Population', fontsize=12, fontweight='bold')
ax1 = fig0.add_subplot(3, 1, 1)
f1=plt.plot(time*scale,intensity)
ax1.set_ylabel('|E| ')
ax2 = fig0.add_subplot(3, 1, 2)
f2=plt.plot(time*scale,intensity)
ax2.set_ylabel('|P|')
ax3 = fig0.add_subplot(3, 1, 3)
f3=plt.plot(time*scale,y[:,8])
plt.ylim(min(y[:,8]), max(y[:,8]))
ax3.set_ylabel('Population')
plt.xlabel('Time($\mu s$)')
plt.text(-0.09,-1.05, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)



fig5=plt.figure()
fig5.suptitle('Comparison between |E| and modulation', fontsize=12, fontweight='bold')
ax1 = fig5.add_subplot(111)
f1=plt.plot(time*scale,intensity,'b', label='|E|')
ax2=ax1.twinx()
ax2=plt.plot(time*scale,((max(intensity)-min(intensity))/2)*np.cos(wf*time),'g', label='modulation')
ax1.set_xlabel('time($\mu s$)')
ax1.set_ylabel('|E|')
plt.xlim(min(time*scale), max(time*scale))
plt.text(-0.1,-0.33, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f" % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
plt.subplots_adjust(bottom=0.22)
plt.legend()


fig4=plt.figure()
fig4.suptitle('Comparison between |E| and the Modulation', fontsize=12, fontweight='bold')
ax1 = fig4.add_subplot(3, 1, 1)
f1=plt.plot(time,np.cos(wf*time))
plt.xlim(min(time), min(time)+(10*2*pi/wf))
ax1.set_ylabel('Modulation ')
ax2 = fig4.add_subplot(3, 1, 2)
f2=plt.plot(time,intensity)
ax2.set_ylabel('|E|')
plt.xlim(min(time), min(time)+10*2*pi/wf)
ax3 = fig4.add_subplot(3, 1, 3)
plt.xlim([max(time)-10*2*pi/wf, max(time)])
f3=plt.plot(time,intensity)
ax3.set_ylabel('|E|')
plt.text(-0.09,-1.05, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
plt.subplots_adjust(bottom=0.22)
plt.xlabel('Time renormalized')


#fig4=plt.figure()
#fig4.suptitle('Comparison between Re(E_y) and the Modulation', fontsize=12, fontweight='bold')
#ax1 = fig4.add_subplot(3, 1, 1)
#f1=plt.plot(time,np.cos(wf*time))
#plt.xlim(min(time), min(time)+(15*2*pi/wf))
#ax1.set_ylabel('Modulation ')
#ax2 = fig4.add_subplot(3, 1, 2)
#f2=plt.plot(time,y[:,2])
#ax2.set_ylabel('$Re(E_y)$')
#plt.xlim(min(time), min(time)+15*2*pi/wf)
#ax3 = fig4.add_subplot(3, 1, 3)
#plt.xlim([max(time)-10*2*pi/wf, max(time)])
#f3=plt.plot(time,intensity_ey)
#ax3.set_ylabel('|$E_y$|')
#plt.text(-0.09,-1.05, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
#plt.subplots_adjust(bottom=0.22)
#plt.xlabel('Time renormalized')

fig1=plt.figure()
fig1.suptitle('|E| vs time', fontsize=12, fontweight='bold')
ax1 = fig1.add_subplot(111)
f1=plt.plot(time*scale,intensity,'b',label='$|E|$')
f2=plt.plot(time*scale,y[:,0],'g', label='$Re(E_x)$')
f3=plt.plot(time*scale,y[:,2],'r', label='$Re(E_y)$')
ax1.set_xlabel('time($\mu s$)')
ax1.set_ylabel('Electric field')
plt.xlim(min(time*scale), max(time*scale))
plt.text(-0.1,-.33, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f \n |E|(blue), Re($E_x$) (Green), Re($E_y$) (Red)" % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
plt.subplots_adjust(bottom=0.22)
plt.legend(fontsize = 'medium')

fig2=plt.figure()
fig2.suptitle('|E| vs time', fontsize=12, fontweight='bold')
ax1 = fig2.add_subplot(111)
f1=plt.plot(time*scale,intensity,'b', label='$|E|$')
f2=plt.plot(time*scale,y[:,1],'g', label='$Im(E_x)$')
f3=plt.plot(time*scale,y[:,3],'r', label='$Im(E_y)$')
ax1.set_xlabel('time($\mu s$)')
ax1.set_ylabel('Electric field')
plt.xlim(min(time*scale), max(time*scale))
plt.text(-0.1,-.33, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f \n |E|(blue), Im($E_x$) (Green), Im($E_y$) (Red)" % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
plt.subplots_adjust(bottom=0.22)
plt.legend(fontsize = 'medium')


fig3=plt.figure()
fig3.suptitle('|E| vs time', fontsize=12, fontweight='bold')
ax1 = fig3.add_subplot(111)
f1=plt.plot(time*scale,intensity,'b', label='$|E|$')
f2=plt.plot(time*scale,intensity_ex,'g', label='$|E_x|$')
f3=plt.plot(time*scale,intensity_ey,'r', label='$|E_y|$')
ax1.set_xlabel('time($\mu s$)')
ax1.set_ylabel('Electric field ')
plt.xlim(min(time*scale), max(time*scale))
plt.text(-0.1,-.33, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f \n |E|(blue), |$E_x$| (Green), |$E_y$| (Red)" % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
plt.subplots_adjust(bottom=0.22)
plt.legend(fontsize = 'medium')



#fig5=plt.figure()
#fig5.suptitle('E_y vs tiempo', fontsize=12, fontweight='bold')
#ax1 = fig5.add_subplot(111)
#f1=plt.plot(time,y[:,3],'b')
#f1=plt.plot(time,y[:,2],'g')
#f1=plt.plot(time,np.sqrt(y[:,2]**2+y[:,3]**2),'k')
#ax1.set_xlabel('time(ms)')
#ax1.set_ylabel('$E_y$')
#plt.xlim(min(time), max(time))
#plt.text(-0.1,-1.04, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f \n $|E_y|$ (black dashed), $Re(E_y)$(green) $Im(E_y)$(blue)" % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
#plt.subplots_adjust(bottom=0.22)
#fig1.savefig('moduloE_fasemodulada.png')


fig9=plt.figure()
fig9.suptitle('|E| vs population', fontsize=12, fontweight='bold')
ax2 = fig9.add_subplot(111)
plt.plot(y[:,8], intensity)
ax2.set_xlabel('Population')
ax2.set_ylabel('|E|')
plt.text(-0.1,-.3, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax2.transAxes)   
plt.subplots_adjust(bottom=0.22)


#fig2=plt.figure()
#fig2.suptitle('|P| vs population', fontsize=12, fontweight='bold')
#ax2 = fig2.add_subplot(111)
#plt.plot(np.sqrt(y[:,4]**2+y[:,5]**2+y[:,6]**2+y[:,7]**2),y[:,8])
#ax2.set_xlabel('|P| ')
#ax2.set_ylabel('Population')
#plt.ylim(min(y[:,8]), max(y[:,8]))
#plt.text(-0.09,-1.05, "\n Parameters: $m= $ %s , $w_{mod}$= %.2f khz , $\Delta \phi_0=$ %s ,  $\Omega=$ %f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m,wf_real, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax3.transAxes)   
#plt.subplots_adjust(bottom=0.22)


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

#
#fig2=plt.figure()
#ax2 = fig2.add_subplot(111)
#plt.plot(np.sqrt(y[:,4]**2+y[:,5]**2+y[:,6]**2+y[:,7]**2),np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2))#real del campo vs real poblacion
#fig2.suptitle('modulo P vs poblacion', fontsize=12, fontweight='bold')
#ax2.set_xlabel('|p| ')
#ax2.set_ylabel('|E|')
#

#fig2=plt.figure()
#ax2 = fig2.add_subplot(111)
#plt.plot(y[:,0],y[:,2])#real del campo vs real poblacion
#fig2.suptitle('$E_x$ vs $E_y$', fontsize=12, fontweight='bold')
#ax2.set_xlabel('$Re(E_x)$ ')
#ax2.set_ylabel('$Re(E_y)$')

#points = np.array([y[:,0], y[:,4]]).T.reshape(-1, 1, 2)
#segments = np.concatenate([points[:-1], points[1:]], axis=1)
#lc = LineCollection(segments, cmap=plt.get_cmap('rainbow'))
#lc.set_array(time)#si en lugar de time pongo otro intervalo puedo visualizar que tarn rapido se mueve la curva.
#lc.set_linewidth(1)
#
#fig3 = plt.figure()
#plt.gca().add_collection(lc)
#plt.xlim(min(y[:,0]), max(y[:,0]))
#plt.ylim(min(y[:,4]), max(y[:,4]))
#plt.show()
#
#points = np.array([y[:,0], y[:,2]]).T.reshape(-1, 1, 2)
#segments = np.concatenate([points[:-1], points[1:]], axis=1)
#lc = LineCollection(segments, cmap=plt.get_cmap('rainbow'))
#lc.set_array(time)#si en lugar de time pongo otro intervalo puedo visualizar que tarn rapido se mueve la curva.
#lc.set_linewidth(1)
#
#fig3 = plt.figure()
#plt.gca().add_collection(lc)
#plt.xlim(min(y[:,0]), max(y[:,0]))
#plt.ylim(min(y[:,2]), max(y[:,2]))
#plt.show()


