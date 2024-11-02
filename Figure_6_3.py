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

# also no disk state
#total/ionizable groups. so 2 is for a 1 to 1 monomer. set so the hydrophobic parameter is a group specific parameter and not dependent on composition. Easier to see if physical values come out.
plt.rcParams.update({'font.size':7})
plt.rcParams.update({'figure.autolayout': 1})
plt.rcParams['legend.fontsize'] =6
plt.rcParams["figure.figsize"] = (12/2.54,5.5/2.54)
plt.rcParams['legend.title_fontsize'] = 'small'
plt.rcParams['lines.markersize'] = 2
plt.rcParams['lines.linewidth'] = 1
from cycler import cycler

#cycles through colours and line style for the plots
custom_cycler =(cycler(color= ['#000000', '#0f94d1','#f7ad05','#009e74', '#f02495', '#f2e41d', '#0071b2', '#d55c00']) + cycler(linestyle=['solid', 'dotted', 'dashed', 'dashdot',(0, (5, 10)),(0, (3, 10, 1, 10)),(0, (5, 1)),(0, (1, 1))]))
plt.rc('axes',prop_cycle=custom_cycler)
#hydrophobic partition function

#r=1 for a homopolymer
def Ph():
    return 1




#aqueous partition function

def Pa(pH,Gh,pka1,pka2,M1,M2):
    return mp.exp(-Gh*(M1+M2))*mp.power(((1+(mp.power(10,(pH-pka1))))),M1)*mp.power(((1+(mp.power(10,(pka2-pH))))),M2)





#total partition function

def Pt(pH,Gh,pka1,pka2,M1,M2):
    return Ph()+Pa(pH,Gh,pka1,pka2,M1,M2)


#defining the fraction of polymer in each state

def fh2(pH,Gh,pka1,pka2,M1,M2):
    return np.array(Ph()/Pt(pH,Gh,pka1,pka2,M1,M2),dtype=float)

fh=np.vectorize(fh2)
    
def fa2(pH,Gh,pka1,pka2,M1,M2):
    return np.array(Pa(pH,Gh,pka1,pka2,M1,M2)/Pt(pH,Gh,pka1,pka2,M1,M2),dtype=float)

fa=np.vectorize(fa2)

def theta1_2(pH,Gh,pka1,pka2,M1,M2):
    return (1/Pt(pH,Gh,pka1,pka2,M1,M2))*mp.power(10,(pH-pka1))* mp.exp(-Gh*(M1+M2))*mp.power(((1+(mp.power(10,(pka2-pH))))),M2)*mp.power(((1+(mp.power(10,(pH-pka1))))),M1-1)

thet1=np.vectorize(theta1_2)

def theta2_2(pH,Gh,pka1,pka2,M1,M2):
    return (1/Pt(pH,Gh,pka1,pka2,M1,M2))*mp.power(10,(pka2-pH))* mp.exp(-Gh*(M1+M2))*mp.power(((1+(mp.power(10,(pka2-pH))))),M2-1)*mp.power(((1+(mp.power(10,(pH-pka1))))),M1)

thet2=np.vectorize(theta2_2)

#we know need a graph where we can vary the quantities for pka's and Gh's and see the effect on the fractions. Slider would be ideal. Copied from the slider scripts

#defining our pH dimension
pH_range = np.linspace(1, 13, 1001)



fig,(ax1,ax2)=plt.subplots(1,2)

M=20
Gh=6.9
pka1=4.69
pka2=10.1


ax1.plot(pH_range,fh(pH_range,Gh,pka1,pka2,M,M))
ax1.plot(pH_range,fa(pH_range,Gh,pka1,pka2,M,M))
ax1.plot(pH_range,thet1(pH_range,Gh,pka1,pka2,M,M))
ax1.plot(pH_range,thet2(pH_range,Gh,pka1,pka2,M,M))
ax1.set_xlabel('$pH$')
ax1.set_ylabel(r'$f$ , $\theta$')
ax1.set_xticks([1,3,5,7,9,11,13])
ax1.set_xlim(1,13)


M=20
Gh=2
pka1=4.69
pka2=6

ax2.plot(pH_range,fh(pH_range,Gh,pka1,pka2,M,M),label='$f_{H}$')
ax2.plot(pH_range,fa(pH_range,Gh,pka1,pka2,M,M),label='$f_{aq}$')
ax2.plot(pH_range,thet1(pH_range,Gh,pka1,pka2,M,M),label=r'$\theta_{acid}$')
ax2.plot(pH_range,thet2(pH_range,Gh,pka1,pka2,M,M),label=r'$\theta_{basic}$')
ax2.legend(frameon=0)
ax2.set_xlabel('$pH$')
#ax2.set_ylabel(r'$f$ , $\theta$')
ax2.set_xticks([1,3,5,7,9,11,13])
ax2.set_xlim(1,13)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)

ax1.text(-0.25, 0.95, '(a)', transform=ax1.transAxes, size=8,zorder=100)
ax2.text(-0.25, 0.95, '(b)', transform=ax2.transAxes, size=8,zorder=100)
