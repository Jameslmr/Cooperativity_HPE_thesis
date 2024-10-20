# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:27:22 2023

@author: marti119
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 10:57:06 2022

@author: james
"""

#simple script to predict some transitions and illustrate the dependence of the sharpness on the value of M

import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from cycler import cycler
import matplotlib.ticker as ticker
import scipy as sp

plt.rcParams.update({'font.size':7})
plt.rcParams.update({'figure.autolayout': 1})#forcesa tight lauout
plt.rcParams['legend.fontsize'] =6
plt.rcParams["figure.figsize"] = (12/2.54,5.5/2.54)
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


def theta(p,K):
    return np.array(((p*K/(1+p*K))),dtype=float)

#defining our p dimension
p_range = np.linspace(0, 100, 1001)
        




K=np.array([1E-2,5E-2,1E-1,5E-1])



fig,(ax1,ax2)=plt.subplots(1,2,num=1)
ax1.text(-0.25, 0.95, '(a)', transform=ax1.transAxes, size=8)
ax2.text(-0.25, 0.95, '(b)', transform=ax2.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)





ax1.plot(p_range,theta(p_range,K[0]),label='{:.1e}'.format(K[0]))
ax1.plot(p_range,theta(p_range,K[1]),label='{:.1e}'.format(K[1]))
ax1.plot(p_range,theta(p_range,K[2]),label='{:.1e}'.format(K[2]))
ax1.plot(p_range,theta(p_range,K[3]),label='{:.1e}'.format(K[3]))
ax1.set_xlabel(r'$\frac{p}{p_0}$')
ax1.set_ylabel(r'$\theta$')
ax1.legend(frameon=False,title='$K_b$')







p_range2 = np.logspace(-2, 4, 1001)


ax2.plot(np.log10(p_range2),theta(p_range2,K[0]),label='{:.1e}'.format(K[0]))
ax2.plot(np.log10(p_range2),theta(p_range2,K[1]),label='{:.1e}'.format(K[1]))
ax2.plot(np.log10(p_range2),theta(p_range2,K[2]),label='{:.1e}'.format(K[2]))
ax2.plot(np.log10(p_range2),theta(p_range2,K[3]),label='{:.1e}'.format(K[3]))
ax2.set_xlabel(r'$\log_{10}(\frac{p}{p_0})$')
ax2.set_ylabel(r'$\theta$')
ax2.legend(frameon=False,title='$K_b$')


