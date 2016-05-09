# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 01:50:41 2015

@author: Alexis
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 06 17:09:46 2015

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
#j=complex(0,1) #defino el complejo

tita=0.
k=1.
ro=2.
d=-0.18
g=0.01
c=4.
m=5.
a=4.

print 'Parametros:', 'tita= ',  tita, 'k=', k, 'ro=', ro, '....' 


l=1. #(longitud)
n=20. #cant de divisiones
h=l/n #estencil
dt=0.5 #paso temporal.
r=dt/(h**2)


U[i]=[y0=np.array(), y1=np.array(), y2=np.array(), y3=np.array(), y4=np.array()]
Rta=U


def mb(y, t):
    """ y[0],y[1] campo electrico. y[2],y[3] polarizacion, y[4] poblacion. """
    dfr=-k*(((1+m*np.cos(tita*t))*y[0]+(d+a*(laplaciano/4+1-ro**2))*y[1])-y[2])
    dfi=-k*(((1+m*np.cos(tita*t))*y[1]-(d+a*(laplaciano/4+1-ro**2))*y[0])-y[3])
    drr=-(1*y[2]-d*y[3])+y[0]*y[4]
    dri=-(1*y[3]+d*y[2])+y[1]*y[4]
    ddelta=-g*(y[4]-c+(y[0]*y[2]+y[1]*y[3]))
    return [dfr, dfi, drr, dri, ddelta]

time = np.arange(0.01, 70 * pi, 0.001)
dfinit=[30., 30.]   
drinit=[1., 2.]
ddeltainit=[.5]
yinit=np.array(dfinit+drinit+ddeltainit)

print mb(yinit, 0)
y = odeint(mb, yinit, time)

fig0=plt.figure()
ax1 = fig0.add_subplot(3, 1, 1)
f1=plt.plot(time,y[:,0])
plt.ylim(min(y[:,0]), max(y[:,0]))
ax1.set_ylabel('Campo F real')
ax2 = fig0.add_subplot(3, 1, 2)
f2=plt.plot(time,y[:,2])
plt.ylim(min(y[:,2]), max(y[:,2]))
ax2.set_ylabel('Campo R real')
ax3 = fig0.add_subplot(3, 1, 3)
f3=plt.plot(time,y[:,4])
plt.ylim(min(y[:,4]), max(y[:,4]))
ax3.set_ylabel('Poblacion')

fig0.suptitle('Evolucion temporal de Fr, Rr y Poblacion', fontsize=12, fontweight='bold')
plt.xlabel('Tiempo')
fig0.savefig('maxt_reales_temporales.png')


fig1=plt.figure()
ax1 = fig1.add_subplot(111)
f1=plt.plot(y[:,0],y[:,2])
#real del campo vs real polarizacion
fig1.suptitle('Campo F real vs Campo R real', fontsize=12, fontweight='bold')
ax1.set_xlabel('Campo F real')
ax1.set_ylabel('Campo R real')
plt.xlim(min(y[:,0]), max(y[:,0]))
plt.ylim(min(y[:,2]), max(y[:,2]))
fig1.savefig('testreales.png')

fig1=plt.figure()
ax1 = fig1.add_subplot(111)
f1=plt.plot(np.sqrt(y[:,0]**2+y[:,1]**2),np.sqrt(y[:,2]**2+y[:,3]**2))
#real del campo vs real polarizacion
fig1.suptitle('modulo F vs modulo R', fontsize=12, fontweight='bold')
ax1.set_xlabel('Modulo F ')
ax1.set_ylabel('modulo R ')
plt.xlim(min(y[:,0]), max(y[:,0]))
plt.ylim(min(y[:,2]), max(y[:,2]))
fig1.savefig('testreales.png')


fig2=plt.figure()
ax2 = fig2.add_subplot(111)
plt.plot(y[:,0],y[:,4])#real del campo vs real poblacion
fig2.suptitle('Campo F real vs poblacion', fontsize=12, fontweight='bold')
ax2.set_xlabel('Campo F real')
ax2.set_ylabel('Poblacion')
plt.xlim(min(y[:,0]), max(y[:,0]))
plt.ylim(min(y[:,4]), max(y[:,4]))

fig2=plt.figure()
ax2 = fig2.add_subplot(111)
plt.plot(np.sqrt(y[:,0]**2+y[:,1]**2),y[:,4])#real del campo vs real poblacion
fig2.suptitle('modulo F vs poblacion', fontsize=12, fontweight='bold')
ax2.set_xlabel('modulo F ')
ax2.set_ylabel('Poblacion')
plt.xlim(min(y[:,0]), max(y[:,0]))
plt.ylim(min(y[:,4]), max(y[:,4]))



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
