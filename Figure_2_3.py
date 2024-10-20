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
from scipy import optimize
import matplotlib.patches as mpatches
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


#In this case we will set the hydrophobic partition function to 1.
def Ph():
    return 1



#hydrophobic partition function, instead I will add an effective M parameter,e 

def Pa(pH,Gh2,pka2,M,r,e):
    return mp.power((mp.exp(-Gh2*(r-1))*(1+(mp.power(10,(pH-pka2))))),e*M) #for an acidic hpe




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


#defining our pH dimension
pH_range = np.linspace(0, 14, 1401)
        

def fh2(pH,Gh2):
    e=1
    r=2
    M=40
    pka2=4.56
    return fh(pH,Gh2,pka2,M,r,e)


fh3=np.vectorize(fh2)

def theta2(pH,Gh2):
    e=1
    r=2
    M=40
    pka2=4.56
    return theta(pH,Gh2,pka2,M,r,e)

theta3=np.vectorize(theta2)

fig,(ax1,ax2)=plt.subplots(1,2,num=0)
ax1.text(-0.25, 0.95, '(a)', transform=ax1.transAxes, size=8)
ax2.text(-0.25, 0.95, '(b)', transform=ax2.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)
ax1.plot(pH_range,1-fh3(pH_range,0.2),label='$f_{aq}$')
ax1.plot(pH_range,theta3(pH_range,0.2),label=r'$\theta $')
ax1.set_xlabel('$pH$')
ax1.set_ylabel(r'$ f_{aq},\theta $')
ax1.legend(frameon=0)
ax1.set_xlim(2.5,6.5)
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.5))

ax2.plot(pH_range,1-fh3(pH_range,3),label='3',linestyle='-')
ax2.plot(pH_range,1-fh3(pH_range,1.5),label='1.5',linestyle='-')
ax2.plot(pH_range,1-fh3(pH_range,0.5),label='0.5',linestyle='-')
ax2.plot(pH_range,1-fh3(pH_range,0.2),label='0.2',linestyle='-')
plt.gca().set_prop_cycle(None)
plt.rc('axes',prop_cycle=custom_cycler)
ax2.plot(pH_range,theta3(pH_range,3),linestyle='dashed',label=' ')
ax2.plot(pH_range,theta3(pH_range,1.5),linestyle='dashed',label=' ')
ax2.plot(pH_range,theta3(pH_range,0.5),linestyle='dashed',label=' ')
ax2.plot(pH_range,theta3(pH_range,0.2),linestyle='dashed',label=' ')
ax2.set_xlim(3,6.5)
ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax2.set_xlabel('$pH$')
#ax2.set_ylabel(r'$f_{aq}, \theta$')
#ax2.legend(frameon=False,title='$g_{H}$')


handles_, labels_ = ax2.get_legend_handles_labels()
empty = mpatches.Patch(color='white')
ax2.legend(frameon=0,handles=[handles_[0],handles_[4],handles_[1],handles_[5],handles_[2],handles_[6],handles_[3],handles_[7]],labels=[labels_[0],labels_[4],labels_[1],labels_[5],labels_[2],labels_[6],labels_[3],labels_[7]],title='$g_{H}$')


#%%
fig,ax=plt.subplots(num=2)
plt.plot(pH_range,theta3(pH_range,0),label='HH, $g_H=0$',linestyle='-')
plt.xlabel('$pH$')
plt.ylabel(r'$\theta$')
plt.xlim(3,6.5)
plt.legend()

# def fh_fit2(pH,M,Gh2):
#     e=1
#     r=2
#     pka2=2
#     return fh(pH,Gh2,pka2,M,r,e)

# fh_fit=np.vectorize(fh_fit2)

# bounds=((0,0),(50,20))

# p0=[5,3]

# gh_array=np.linspace(0.1,5,50)
# M_fit_array=np.zeros(50)
# for i in np.arange(0,50):
#     params_a, params_covariance_a = optimize.curve_fit(fh_fit,pH_range,np.array(fh3(pH_range,gh_array[i]),dtype=float),bounds=bounds, p0=p0)
#     M_fit_array[i]=params_a[0]

# plt.figure()
# plt.scatter(gh_array,M_fit_array)    
# plt.xlabel(r'$g_{H}$')
# plt.ylabel(r'$M_{fit}$')













