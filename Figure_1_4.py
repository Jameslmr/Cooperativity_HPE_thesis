# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:44:04 2024

@author: marti119
"""

import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import matplotlib.ticker as ticker
import scipy as sp
from scipy.interpolate import splrep, BSpline

plt.rcParams.update({'font.size':7})
plt.rcParams.update({'figure.autolayout': 1})#forcesa tight lauout
plt.rcParams['legend.fontsize'] =6
plt.rcParams["figure.figsize"] = (12/2.54,8/2.54)
plt.rcParams['legend.title_fontsize'] = 'small'
plt.rcParams['lines.markersize'] = 2
plt.rcParams['lines.linewidth'] = 1
#cycles through colours and line style for the plots
custom_cycler =(cycler(color= ['#000000', '#0f94d1','#f7ad05','#009e74', '#f02495', '#f2e41d', '#0071b2', '#d55c00']) + cycler(linestyle=['solid', 'dotted', 'dashed', 'dashdot',(0, (5, 10)),(0, (3, 10, 1, 10)),(0, (5, 1)),(0, (1, 1))]))
plt.rc('axes',prop_cycle=custom_cycler)

mcount=0
def markers(): 
    global mcount 
    mcount+=1
    for i, j in zip(range(mcount), marker_cycle()):
        a=list(j.values())
    return a[0]
    

marker_cycle=(cycler(marker=["o", "v", "*","s","P","d","p","+",]))



#%%

data=np.loadtxt('PMAA_titration.txt',delimiter=',')

fig,ax=plt.subplots()
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)
x=data[:,0]
y=data[:,1]

tck = splrep(x, y, s=0)

x_dim=np.linspace(0.0,1,101)



ax.scatter(x,y,marker=markers())
ax.plot(x_dim, BSpline(*tck,extrapolate=0)(x_dim))

data=np.loadtxt('PMAA_viscosity.txt',delimiter=',')

x=data[:,0]
y=data[:,1]

secax=ax.twinx()

tck = splrep(x, y, s=0)






secax.scatter(x,y,marker=markers(),color='#f7ad05')
secax.plot(x_dim, BSpline(*tck,extrapolate=0)(x_dim),color='#f7ad05')
secax.set_ylim(-0.009,0.0346)
secax.set_yticks((0,0.005,0.01,0.015,0.02,0.025,0.03))


data=np.loadtxt('PAA_titration.txt',delimiter=',')

x=data[:,0]
y=data[:,1]

tck = splrep(x, y, s=0)




ax.scatter(x,y,marker=markers())
ax.plot(x_dim, BSpline(*tck,extrapolate=0)(x_dim))

ax.set_xlim(0,1)

ax.set_xlabel(r'$\theta$')
ax.set_ylabel(r'$pK$')
secax.set_ylabel(r'$\eta$')


