# -*- coding: utf-8 -*-
"""
Created on Tue May 17 11:49:57 2022

@author: james
"""
#simple script to predict some transitions and illustrate the dependence of the sharpness on the value of M

import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from cycler import cycler
import matplotlib.ticker as ticker
import scipy as sp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import gridspec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

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


#In this case we will set the hydrophobic partition function to 1.
def Ph():
    return 1



#hydrophobic partition function, instead I will add an effective M parameter,e 

def Pa(pH,Gh2,pka2,M,r,e):
    return mp.power((mp.exp(-Gh2*(r))*(1+(mp.power(10,(pH-pka2))))),e*M) #for an acidic hpe




#def Pd(pH,Gh3,pka3,pka4,M,f,r):
   # return np.vectorize(Pd2(pH,Gh3,pka3,pka4,M,f,r))




#total partition function

def Pt(pH,Gh2,pka2,M,r,e):
    return Ph()+Pa(pH,Gh2,pka2,M,r,e)


#defining the fraction of polymer in each state

def fh(pH,Gh2,pka2,M,r,e):
    return np.array(Ph()/Pt(pH,Gh2,pka2,M,r,e),dtype=float)
    
def fa(pH,Gh2,pka2,M,r,e):
    return np.array(Pa(pH,Gh2,pka2,M,r,e)/Pt(pH,Gh2,pka2,M,r,e),dtype=float)

def theta(pH,Gh2,pka2,M,r,e):
    return np.array(((mp.power(10,(pH-pka2))/(1+mp.power(10,(pH-pka2))))*(1-fh(pH,Gh2,pka2,M,r,e))),dtype=float)

def theta_1p(pH,Gh2,pka2,M,r,e):
    return np.array(((mp.power(10,(pH-pka2))/(1+mp.power(10,(pH-pka2))))),dtype=float)


#defining our pH dimension
pH_range = np.linspace(2, 8, 1401)
        

def fh2(pH,M,Gh2):
    e=1
    r=1
    pka2=4.56
    return fh(pH,Gh2,pka2,M,r,e)


fh3=np.vectorize(fh2)


def thet2(pH,M,Gh2):
    e=1
    r=1
    pka2=4.56
    return theta(pH,Gh2,pka2,M,r,e)


thet3=np.vectorize(thet2)

def thet_1p2(pH,M,Gh2):
    e=1
    r=1
    pka2=4.56
    return theta_1p(pH,Gh2,pka2,M,r,e)


thet_1p3=np.vectorize(thet_1p2)

Gh2=3

fig = plt.figure(figsize=(12/2.54,13.5/2.54))

gs = gridspec.GridSpec(3, 2,height_ratios=[0.7272727272727,0.7272727272727,1])
ax = plt.subplot(gs[0:2, 0:])
ax3 = plt.subplot(gs[2, 0])
ax4 = plt.subplot(gs[2, 1])
ax3.text(-0.25, 0.95, '(b)', transform=ax3.transAxes, size=8)
ax4.text(-0.25, 0.95, '(c)', transform=ax4.transAxes, size=8)
ax.text(-0.109, 0.98, '(a)', transform=ax.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)


ax.plot(pH_range,fh3(pH_range,1,Gh2),label='$1$')
ax.plot(pH_range,fh3(pH_range,2,Gh2),label='$2$')
ax.plot(pH_range,fh3(pH_range,5,Gh2),label='$5$')
ax.plot(pH_range,fh3(pH_range,20,Gh2),label='$20$')
ax.plot(pH_range,fh3(pH_range,40,Gh2),label='$40$')
ax.set_xlim(4.5,8)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.set_xlabel('$pH$')
ax.set_ylabel(r'$f_H$')
h, l = ax.get_legend_handles_labels()
empty = mpatches.Patch(color='white')
ax.legend(frameon=0,handles=[empty,h[0],h[1],h[2],h[3],h[4]],labels=['$M_{eff}$',l[0],l[1],l[2],l[3],l[4]],loc='center left')





#it would also be intersting to look at the width of the transition that comes out of the model, it's definitely not linear which is ideal.

def bounds(M,Gh2):
    interp_data=sp.interpolate.interp1d(fh3(pH_range,M,Gh2),pH_range)
    upper=interp_data(0.9)
    lower=interp_data(0.1)
    return lower-upper

i=0
bounds_array=np.zeros(100)
for M in np.arange(1,101):
    bounds_array[i]=bounds(M,Gh2)
    i+=1
    



ax2=inset_axes(ax, width='55%', height='50%', loc='center right',bbox_to_anchor=(0.14,0.3,0.85,0.85), bbox_transform=ax.transAxes)
ax2.scatter(np.arange(1,101),bounds_array)
ax2.set_xlabel('$M_{eff}$')
ax2.set_ylabel(r'$\Delta pH_{90-10}$')






#graph for theta
Gh2=3


ax3.plot(pH_range,fh3(pH_range,20,Gh2),label='$f_{H}$')
ax3.plot(pH_range,1-fh3(pH_range,20,Gh2),label='$f_{aq}$')
ax3.plot(pH_range,thet3(pH_range,20,Gh2),label=r'$\theta$')
ax3.plot(pH_range,thet_1p3(pH_range,20,Gh2),label=r'$\theta_{HH}$')
ax3.set_xlabel('$pH$')
ax3.set_ylabel(r'$f$ , $\theta$')
ax3.legend(frameon=False)
ax3.set_xlim(2,8)
ax3.xaxis.set_major_locator(ticker.MultipleLocator(1))






#same for a basic polylectrolyte



#hydrophobic partition function, instead I will add an effective M parameter,e 

def Pa_b(pH,Gh2,pka2,M,r,e):
    return mp.power((mp.exp(-Gh2*(r))*(1+(mp.power(10,(pka2-pH))))),e*M) #for an acidic hpe




#def Pd(pH,Gh3,pka3,pka4,M,f,r):
   # return np.vectorize(Pd2(pH,Gh3,pka3,pka4,M,f,r))




#total partition function

def Pt_b(pH,Gh2,pka2,M,r,e):
    return Ph()+Pa_b(pH,Gh2,pka2,M,r,e)


#defining the fraction of polymer in each state

def fh_b(pH,Gh2,pka2,M,r,e):
    return np.array(Ph()/Pt_b(pH,Gh2,pka2,M,r,e),dtype=float)
    
def fa_b(pH,Gh2,pka2,M,r,e):
    return np.array(Pa_b(pH,Gh2,pka2,M,r,e)/Pt_b(pH,Gh2,pka2,M,r,e),dtype=float)

def theta_b(pH,Gh2,pka2,M,r,e):
    return np.array(((mp.power(10,(pka2-pH))/(1+mp.power(10,(pka2-pH))))*(1-fh_b(pH,Gh2,pka2,M,r,e))),dtype=float)

def theta_1p_b(pH,Gh2,pka2,M,r,e):
    return np.array(((mp.power(10,(pka2-pH))/(1+mp.power(10,(pka2-pH))))),dtype=float)


#defining our pH dimension
pH_range = np.linspace(2, 8, 1401)
        

def fh2_b(pH,M,Gh2):
    e=1
    r=1
    pka2=10.1
    return fh_b(pH,Gh2,pka2,M,r,e)


fh3_b=np.vectorize(fh2_b)


def thet2_b(pH,M,Gh2):
    e=1
    r=1
    pka2=10.1
    return theta_b(pH,Gh2,pka2,M,r,e)


thet3_b=np.vectorize(thet2_b)

def thet_1p2_b(pH,M,Gh2):
    e=1
    r=1
    pka2=10.1
    return theta_1p_b(pH,Gh2,pka2,M,r,e)

thet_1p3_b=np.vectorize(thet_1p2_b)

#graph for theta
Gh2=3
pH_range = np.linspace(5, 14, 1401)
ax4.plot(pH_range,fh3_b(pH_range,20,Gh2),label='$f_{H}$')
ax4.plot(pH_range,1-fh3_b(pH_range,20,Gh2),label='$f_{aq}$')
ax4.plot(pH_range,thet3_b(pH_range,20,Gh2),label=r'$\theta$')
ax4.plot(pH_range,thet_1p3_b(pH_range,20,Gh2),label=r'$\theta_{HH}$')
ax4.set_xlabel('$pH$')
ax4.legend(frameon=False)
ax4.set_xlim(7,13)
ax4.xaxis.set_major_locator(ticker.MultipleLocator(1))












