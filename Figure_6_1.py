# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 13:11:19 2022

@author: marti119
"""

import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from scipy.optimize import curve_fit
import time
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from cycler import cycler



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

colours= ['#000000', '#0f94d1','#f7ad05','#009e74', '#f02495', '#f2e41d', '#0071b2', '#d55c00']




pKa=4.56
def fH(x,g,M): #fraction in hydrophobic state - no of hydrophobic groups = no of ionic groups
    return (1 + np.exp(-M*g)*(1 + 10**(x - pKa))**M)**-1     #x is gonnabe pH
                                                                    #g hydrophobic penalty per hydrophobic monomer
x=np.linspace(5,9,1000)
fig,ax=plt.subplots()
for n in [2,10]:
#for n in [2,10]:
    plt.plot(x, fH(x,5.6,n))
#plt.title('Influence of M')
plt.xlabel('$pH$')
#plt.ylabel('f_H')
plt.ylabel('$f_H$')

y=np.linspace(-1,2,1000)

def straight_line(y,a):
    
    return a+0*y

plt.ylim(-0.05,1.05)
plt.plot(straight_line(y,6.7),y,linestyle='dashed')
plt.plot(straight_line(y,7.2),y,linestyle='dashed')

ax.xaxis.set_ticks(np.arange(5, 9.5, 0.5))

plt.annotate('Tumor cell $pH$',(5.25,0.1),color='#f7ad05')
plt.annotate('Healthy cell $pH$',(7.35,0.5),color='#009e74')
plt.annotate('$M=2$',(5.9,0.84),color='#000000')
plt.annotate('$M=10$',(8.1,0.9),color='#0f94d1')
plt.arrow(8,0.92,-1.1,-0.1,color='#0f94d1',length_includes_head=True,
          head_width=0.03, head_length=0.1)




#less cluttered 2nd version

fig,ax=plt.subplots()
plt.plot(x, fH(x,5.3,1),label='1')
plt.plot(x, fH(x,5.3,10),label='10',linestyle='dashdot')
#plt.title('Influence of M')
plt.xlabel('$pH$')
#plt.ylabel('f_H')
plt.ylabel('$f_H$')
plt.legend(frameon=0,title=r'$M_{eff}$')
y=np.linspace(-1,2,1000)

def straight_line(y,a):
    
    return a+0*y

plt.ylim(-0.05,1.05)
plt.plot(straight_line(y,6.7),y,linestyle='dashed')
plt.plot(straight_line(y,7.2),y,linestyle='dashed')

ax.xaxis.set_ticks(np.arange(5, 9.5, 0.5))

plt.annotate('Tumor cell $pH$',(5.25,0.3),color='#f7ad05')
plt.annotate('Healthy cell $pH$',(7.35,0.5),color='#009e74')
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)


