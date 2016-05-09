# -*- coding: utf-8 -*-
"""
Integrates the maxwell-block equations with a modulation on the polarization phase.
Swypes the value of the amplitud modulation and output files needed to plot an
stroboscopic map.


Created on Thu Apr 21 16:20:03 2016

@author: Alexis
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy.integrate import odeint
from intmbfase import intmbfase as integ
from scipy.signal import argrelextrema
from time import localtime
import time as timing
import os


#%%
#Read parameters from 'Params.txt', set.
txtfile=open("Params.txt")
tempvar=txtfile.readlines()

'''parameters for normalization'''
a= float(tempvar[1].split()[1])
gperp= float(tempvar[2].split()[1])  #gamma perpendicular, loss rate
scale=1*(10.**6)/gperp #scale to micro seconds
wscale=1000*gperp/(10.**6)#scale frequency to khz

'''parameters for the equation'''
kr= float(tempvar[3].split()[1])
k=kr/gperp #normalized loss rate
mu= float(tempvar[4].split()[1])
Dphi0= float(tempvar[5].split()[1]) #phase shift [-pi,pi]
d= float(tempvar[6].split()[1])#detuning
gr= float(tempvar[7].split()[1])
g=gr/gperp #*((2*pi)**2) #sigma parallel, normalized loss rate
D0=a*k/mu #Poblation
m= float(tempvar[8].split()[1]) #modulation amplitud [0,1]
wf= float(tempvar[9].split()[1])
print 'Params: ', a,gperp,kr,mu,Dphi0,d,gr,m,wf

'''parameters to compare with the results'''
w_res=np.sqrt(k*g*((D0*mu/k)-1.))*wscale #resonance frequency
a=D0*mu/k
w=np.sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale #Relaxation oscilations frequency
wf_real=wf*wscale

txtfile.close()#deberia cerrarlo aca?

txtfile=open("Status.txt",'r')
tempvar=txtfile.readlines()
stat= int(tempvar[0].split()[1])
print 'Status=', stat
txtfile.close()#deberia cerrarlo aca?

#%%

def editline(file,n_line,text):
    with open(file) as infile:
        lines = infile.readlines()
    lines[n_line] = text+' \n'
    with open(file, 'w') as outfile:
        outfile.writelines(lines)        

#%%        
def swipe(m,k,mu,Dphi0,d,g,D0,wf):
    '''Swype parameters'''
    m_max=0.0260
    m_min=0.023
    h=0.0000025
    intime=500.*15*10**(-6)*gperp #integration time FOR TRANSITORY    
    print 'total number of steps: %s' %((m_max-m_min)/h), 'Integration time per step: %s' %intime
    
    if stat==0:    
        '''User defined initial condition'''
        dfxinit=[1., 1.] 
        dfyinit=[1., -1.9]  
        drxinit=[1., 1.]
        dryinit=[1., -1.9] 
        ddeltainit=[6.65973518e+03]
        yinit=np.array(dfxinit+dfyinit+drxinit+dryinit+ddeltainit)
        time=np.array([0])
        m=m_min
        editline('Params.txt',9,'wf: %f' %wf)#update next wf to file

        binwrite=open('m_max.in','wb')#clear w_max and strobo
        binwrite.close()
        binwrite=open('max_peak.in','wb')
        binwrite.close()
        binwrite=open('benchmark.in','wb')# to save runtimes
        np.array(intime).tofile(binwrite)
        binwrite.close()

    elif stat==1:  
        yinit=np.fromfile('yinitials.in',dtype=np.float64)
        time=np.fromfile('tinitials.in',dtype=np.float64)
          
    with open('Status.txt','w') as f:
        f.write('Stat: 1')
   
    '''swipe''' 
    while m<m_max:
        #error print
        tin=timing.time()

        if m>m+h: print 'man, m esta disminuyendo, se hace infinito el loop.' 

        print 'integrating m= %f' %m, 
        '''Transitory integration'''
        timeinit = np.arange(time[-1] ,intime+time[-1] , 15.)
        y, time=integ(yinit,timeinit,k,mu,Dphi0,d,g,D0,m,wf)
        '''integration'''
        timeinit = np.arange(time[-1] ,intime*5/15+time[-1] , 15.)
        y, time=integ(y[-1],timeinit,k,mu,Dphi0,d,g,D0,m,wf)
        '''set intensitys'''
#        intensity_ex=np.sqrt(y[:,0]**2+y[:,1]**2)
#        intensity_ey=np.sqrt(y[:,2]**2+y[:,3]**2)
        intensity=np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2)
        '''map'''
        
        peak_max=intensity[argrelextrema(intensity, np.greater)[0]]#intensity  puedo scar el set
        m_peaks=m*np.ones_like(peak_max)#vector or m, the same lenght as peak_max        
        m=m+h   #ste next wf 

        '''save'''
        editline('Params.txt',8,'m: %f' %m)#update next wf to file
        
        binwrite=open('m_max.in','ab')
        m_peaks.tofile(binwrite)
        binwrite.close()
        binwrite=open('max_peak.in','ab')
        peak_max.tofile(binwrite)
        binwrite.close()

        binwrite=open('yinitals.in','wb')
        np.ndarray.tofile(y[-1],'yinitials.in')
        binwrite.close()
        
        binwrite=open('tinitals.in','wb')
        np.ndarray.tofile(np.array(time[-1]),'tinitials.in')
        binwrite.close() 
        
#       tout=timing.time()  
        runtime=timing.time()-tin
        binwrite=open('benchmark.in','ab')
        np.array([m-h,runtime]).tofile(binwrite)
        binwrite.close() 
        print 'runtime: %.2f s' %(runtime)
        
    return #strobo_map  

#%%    
strobo_map=swipe(m,k,mu,Dphi0,d,g,D0,wf)    
#%%
with open('Status.txt','w') as f:
        f.write('Status: 0') #ends the run, sets status back to 0.

#%%
'''write the results to a subfolder ../Results  (so i can plot them again if needed)'''
filename = "Results/m_max.in"
if not os.path.exists(os.path.dirname(filename)):
    try:
        os.makedirs(os.path.dirname(filename))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
binwrite=open('Results/'+'%d_%d_%d-%d.%d.%d-m_max.in'  %localtime()[0:6],'wb')
m_peks=np.fromfile('m_max.in',dtype=np.float64)        
m_peks.tofile(binwrite)
binwrite.close()        
binwrite=open('Results/'+'%d_%d_%d-%d.%d.%d-max_peak.in' %localtime()[0:6],'wb') 
peak_max=np.fromfile('max_peak.in',dtype=np.float64)
peak_max.tofile(binwrite)
binwrite.close()
#%%    
'''Plots!''' 
save=True #set True if i want to save files automatically
pi=np.pi

#%%    
'''plot swype'''
m_peks=np.fromfile('m_max.in',dtype=np.float64)
max_peak=np.fromfile('max_peak.in',dtype=np.float64)
print np.shape(m_peks), np.shape(max_peak)
fig5=plt.figure()
fig5.suptitle('Local Max of intensity vs. m', fontsize=12, fontweight='bold')
ax1 = fig5.add_subplot(111)
ax1.plot(m_peks,max_peak,',b')
ax1.set_xlabel('m')
ax1.set_ylabel('Local Max. Intensity |E|')
ax1.set_xlim(np.min(m_peks), np.max(m_peks))
ax1.set_ylim(0, 2.5)
plt.text(-0.1,-.32, "\n Parameters: $wf= $ %s , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (wf*wscale, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
plt.subplots_adjust(bottom=0.22)
fig5.set_size_inches(7,7)
if save==True: 
    fname='%d_%d_%d-%d.%d.%d-Max_vs_m (w=%s, m=(%s - %s)).png' % tuple(list(localtime()[0:6])+[wf*wscale, m_peks[0], m_peks[-1]])
    fig5.savefig('Results/'+fname, dpi = 180)# when saving, specify the DPI
plt.show
#%%
'''plot swype benchmark'''
bench=np.fromfile('benchmark.in',dtype=np.float64)
fig5=plt.figure()
fig5.suptitle('Run time vs. m  (integration time= %s ) ' %bench[0] , fontsize=12, fontweight='bold')
ax1 = fig5.add_subplot(111)
for i in arange(1,np.shape(bench[1:])[0],2):
    ax1.plot(bench[i],bench[i+1],'.b')
ax1.set_xlabel('m')
ax1.set_ylabel('Runtime [s]')
plt.text(-0.1,-.32, "\n Parameters: $wf= $ %s  , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (wf*wscale, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
plt.subplots_adjust(bottom=0.22)
fig5.set_size_inches(7,7)
if save==True: 
    fname='%d_%d_%d-%d.%d.%d-benchmark_vs_m (w=%s, m=(%s - %s)).png' % tuple(list(localtime()[0:6])+[wf*wscale, m_peks[0], m_peks[-1]])
    fig5.savefig('Results/'+fname, dpi = 100)# when saving, specify the DPI
plt.show