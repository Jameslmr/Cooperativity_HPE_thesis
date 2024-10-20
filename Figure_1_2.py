# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:44:04 2024

@author: marti119
"""

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
from scipy.special import factorial
from scipy import optimize
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

data=np.loadtxt('Imai_data.txt',delimiter=',')

x=data[:,0]
y=data[:,1]

plt.figure(-1)
plt.scatter(x,y)


thet= (10**y)/(1+10**y)

plt.figure(0)
plt.scatter(10**x,thet)

plt.figure(1)
plt.scatter(x,thet)

#%%
#fitting

#tense parittion function
def Pt(p,n,et):

    return (factorial(4)/(factorial(n)*factorial(4-n)))*(p**n)*np.exp(-et*n)




#Relaxed partition function

def Pr(p,n,er,Gh):

    return (factorial(4)/(factorial(n)*factorial(4-n)))*(p**n)*np.exp(-er*n)*np.exp(-Gh)



def Ptot(p,er,et,Gh):
    return Pt(p,0,et)+Pt(p,1,et)+Pt(p,2,et)+Pt(p,3,et)+Pt(p,4,et)+Pr(p,0,er,Gh)+Pr(p,1,er,Gh)+Pr(p,2,er,Gh)+Pr(p,3,er,Gh)+Pr(p,4,er,Gh)

def ft(p,n,er,et,Gh):
    return Pt(p,n,et)/Ptot(p,er,et,Gh)

def fr(p,n,er,et,Gh):
    return Pr(p,n,er,Gh)/Ptot(p,er,et,Gh)

def ft2(p,er,et,Gh):
    return (Pt(p,0,et)+Pt(p,1,et)+Pt(p,2,et)+Pt(p,3,et)+Pt(p,4,et))/Ptot(p,er,et,Gh)

def fr2(p,er,et,Gh):
    return (Pr(p,0,er,Gh)+Pr(p,1,er,Gh)+Pr(p,2,er,Gh)+Pr(p,3,er,Gh)+Pr(p,4,er,Gh))/Ptot(p,er,et,Gh)

def theta(p,er,et,Gh):

    return (1/Ptot(p,er,et,Gh))*((p*np.exp(-et)*(1+p*np.exp(-et))**3)+(np.exp(-Gh)*p*np.exp(-er)*(1+p*np.exp(-er))**3))

def theta2(p,er,Gh):
    et=-0.9985 #found from a discussion in https://doi.org/10.1371/journal.pone.0182871, dont fully grasp the units tbh
    return (1/Ptot(p,er,et,Gh))*((p*np.exp(-et)*(1+p*np.exp(-et))**3)+(np.exp(-Gh)*p*np.exp(-er)*(1+p*np.exp(-er))**3))

p_axis=np.linspace(0,500,1000)

bounds=((-100,-1000),(10,10))

p0=[-4,10]

params, params_covariance = optimize.curve_fit(theta2, 10**x,thet,bounds=bounds, p0=p0)

#%%

fig = plt.figure(figsize=(12/2.54,13.5/2.54))

gs = gridspec.GridSpec(3, 2,height_ratios=[0.7272727272727,0.7272727272727,1])
ax = plt.subplot(gs[0:2, 0:])
ax3 = plt.subplot(gs[2, 0])
ax4 = plt.subplot(gs[2, 1])
ax3.text(-0.25, 0.95, '(b)', transform=ax3.transAxes, size=8)
ax4.text(-0.25, 0.95, '(c)', transform=ax4.transAxes, size=8)
ax.text(-0.109, 0.98, '(a)', transform=ax.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)

ax.scatter(x,thet,label='Exp. data')

def theta_lang(p,K):
    return np.array(((p*K/(1+p*K))),dtype=float)
ax.plot(np.log10(p_axis),theta_lang(p_axis,0.0956),label=r"$\theta_{lang}$")


ax.plot(np.log10(p_axis),theta2(p_axis,*params),label=r"$\theta_{MWC}$")
ax.plot(np.log10(p_axis),ft2(p_axis,params[0],4.85,params[1]),label=r"$f_T$")
ax.plot(np.log10(p_axis),fr2(p_axis,params[0],4.85,params[1]),label=r"$f_R$")
ax.set_xlabel(r'log$_{10}(p_{\mathrm{O}_{2}})$')
ax.set_ylabel(r'$\theta$ , $f$')
ax.legend(frameon=0)



ax3.plot(np.log10(p_axis),ft(p_axis,0,params[0],4.85,params[1]),label=r"$0$")
ax3.plot(np.log10(p_axis),ft(p_axis,1,params[0],4.85,params[1]),label=r"$1$")
ax3.plot(np.log10(p_axis),ft(p_axis,2,params[0],4.85,params[1]),label=r"$2$")
ax3.plot(np.log10(p_axis),ft(p_axis,3,params[0],4.85,params[1]),label=r"$3$")
ax3.plot(np.log10(p_axis),ft(p_axis,4,params[0],4.85,params[1]),label=r"$4$")

ax3.set_xlabel(r'log$_{10}(p_{\mathrm{O}_{2}})$')
ax3.set_ylabel(r'$f$')
ax3.legend(frameon=0,title='$T$')



ax4.plot(np.log10(p_axis),fr(p_axis,0,params[0],4.85,params[1]),label=r"$0$")
ax4.plot(np.log10(p_axis),fr(p_axis,1,params[0],4.85,params[1]),label=r"$1$")
ax4.plot(np.log10(p_axis),fr(p_axis,2,params[0],4.85,params[1]),label=r"$2$")
ax4.plot(np.log10(p_axis),fr(p_axis,3,params[0],4.85,params[1]),label=r"$3$")
ax4.plot(np.log10(p_axis),fr(p_axis,4,params[0],4.85,params[1]),label=r"$4$")

ax4.set_xlabel(r'log$_{10}(p_{\mathrm{O}_{2}})$')
ax4.set_ylabel(r'$f$')
ax4.legend(frameon=0,title='$R$')

