# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 10:57:06 2022

@author: james
"""

#note this is all for the gao data, so for a basic system not an acidic one
#script to try and understand how theta works. I think the3 absolute values of the pka's do matter in this case, so it's really quite hard to fit stuff 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import mpmath as mp
from cycler import cycler
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

# also no disk state
#total/ionizable groups. so 2 is for a 1 to 1 monomer. set so the hydrophobic parameter is a group specific parameter and not dependent on composition. Easier to see if physical values come out.


#hydrophobic partition function

#f in this case (assuming we set pka3=pka means that f is the fraction of groups that can still ionize in the collapsed conformation)
def Ph():
    return 1




#aqueous partition function

def Pa(pH,Gh,pka,M,r):
    return mp.power((mp.exp(-Gh*(r))),M)*mp.power((1+(mp.power(10,(pH-pka)))),(M))

def Pa3(pH,Gh,pka,M,r):
    return np.array(mp.power((mp.exp(-Gh*(r))),M)*mp.power((1+(mp.power(10,(pH-pka)))),(M)),dtype=float)

Pa2=np.vectorize(Pa3)



#total partition function

def Pt(pH,Gh,pka,M,r):
    return Ph()+Pa(pH,Gh,pka,M,r)


#defining the fraction of polymer in each state

def fh2(pH,Gh,pka,M,r):
    return np.array(Ph()/Pt(pH,Gh,pka,M,r),dtype=float)

fh=np.vectorize(fh2)
    
def fa2(pH,Gh,pka,M,r):
    return np.array(Pa(pH,Gh,pka,M,r)/Pt(pH,Gh,pka,M,r),dtype=float)

fa=np.vectorize(fa2)


def fhv2(pH,Gh,pka,M,r,vol_frac):
    return np.array((1/(1+((Pa(pH,Gh,pka,M,r)/Ph())*vol_frac))),dtype=float)

fh_v=np.vectorize(fhv2)

#defining our pH dimension
pH_range = np.linspace(3, 8, 1401)


#polymer parameters

Gh=3
pka=4.56
M=20
r=1
fig,(ax,ax2)=plt.subplots(1,2)

ax.plot(pH_range,Pa2(pH_range,Gh,pka,M,r),label='$M=20$')
ax.plot(pH_range,Pa2(pH_range,Gh,pka,10,r),label='$M=10$')
ax.plot(pH_range,Pa2(pH_range,Gh,pka,5,r),label='$M=5$')
ax.plot(pH_range,Pa2(pH_range,Gh,pka,2,r),label='$M=2$')
ax.plot(pH_range,Pa2(pH_range,Gh,pka,1,r),label='$M=1$')

ax.set_yscale('log')
ax.set_xlim(3,8)
ax.set_ylim(1E-20,1E+8)
ax.set_xticks((3,4,5,6,7,8))
ax.grid(alpha=0.5)
ax.set_xlabel('$pH$')
ax.set_ylabel('$K_{H-aq}$')
ax.legend(frameon=0)
vol_frac=100
ax2.plot(pH_range,1-fh_v(pH_range,Gh,pka,M,r,vol_frac),label='$M=20$')
    
ax2.plot(pH_range,1-fh_v(pH_range,Gh,pka,5,r,vol_frac),label='$M=5$')

ax2.plot(pH_range,1-fh_v(pH_range,Gh,pka,1,r,vol_frac),label='$M=1$')

vol_frac=1

ax2.plot(pH_range,1-fh_v(pH_range,Gh,pka,M,r,vol_frac),label='$M=20$')
    
ax2.plot(pH_range,1-fh_v(pH_range,Gh,pka,5,r,vol_frac),label='$M=5$')

ax2.plot(pH_range,1-fh_v(pH_range,Gh,pka,1,r,vol_frac),label='$M=1$')
ax2.set_xticks((3,4,5,6,7,8))
ax2.set_xlim(3,8)
ax2.set_xlabel('$pH$')
ax2.set_ylabel('$f_{aq}$')

ax.text(-0.25, 0.95, '(a)', transform=ax.transAxes, size=8,zorder=100)

ax2.text(-0.25, 0.95, '(b)', transform=ax2.transAxes, size=8,zorder=100)

fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)



h,l=ax2.get_legend_handles_labels()
empty = mpatches.Patch(color='white')

ax2.legend(frameon=0,handles=[empty,h[0],h[1],h[2],empty,empty,h[3],h[4],h[5]],labels=['$V_r=100$' , l[0],l[1],l[2],'','$V_r=1$',l[3],l[4],l[5]],loc='center left')





