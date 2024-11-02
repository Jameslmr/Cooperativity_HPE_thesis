# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:22:05 2023

@author: marti119
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import mpmath as mp
from matplotlib import ticker
#%% ####plotting stuff##

from cycler import cycler

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

colours= ['#000000', '#0f94d1','#f7ad05','#009e74', '#f02495', '#f2e41d', '#0071b2', '#d55c00']
#########
#exactly the same as the other 3 state system but we have explicitely linked M and the Gh which makes the change in M a bit more clear. So we have an r value to represent the hydrophilic to hydrophobic ratio.

#The main thing in this script is to unify the description for the partition function of the polymer in all enviroments under one general formula, M will be inmutable and we will instead change pka and Gh. Will this work? I think the relative change in pka between states will yield the same effect as an artificially constrained effective M

#In this script we will look at 3 coexisting states the "oily", the "aqueous", and an intermediate state, described here as a disk state.

#the parameters that are important n this case are pH, pka, Gh. We assume acidic ionizable groups.

#we will define the patitition functions as P and the sum of all of them Ptot. Derivation comes from assuming a template model for the polymer base (or proton holes) and a enviroment dependen Gh and pka

#the Gh parameters are in KT! So no need for units, all is essentially unitless

#hydrophobic partition function

#total/ionizable groups. so 2 is for a 1 to 1 monomer

# use of mp math allows for arbitrary floating point precision, in this case we use to have massive values of M not produce overflows.
def Ph2(pH,Gh1,pka1,M,r):
    return mp.power((mp.exp(-Gh1*(r))*(1+(mp.power(10,(pH-pka1))))),M)

Ph=np.vectorize(Ph2)


#hydrophobic partition function

def Pa2(pH,Gh2,pka2,M,r):
    return mp.power((mp.exp(-Gh2*(r))*(1+(mp.power(10,(pH-pka2))))),M)

Pa=np.vectorize(Pa2)


#intermediate, disk, partition function

def Pd2(pH,Gh3,pka3,pka4,M,f,r):
    return mp.power((mp.exp(-Gh3*(r))),f*M)*mp.power((1+(mp.power(10,(pH-pka3)))),(f*M))*mp.power(((1+(mp.power(10,(pH-pka4))))),((1-f)*M))

#def Pd(pH,Gh3,pka3,pka4,M,f,r):
   # return np.vectorize(Pd2(pH,Gh3,pka3,pka4,M,f,r))

Pd=np.vectorize(Pd2)


#total partition function

def Pt(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r):
    return Ph(pH,Gh1,pka1,M,r)+Pa(pH,Gh2,pka2,M,r)+Pd(pH,Gh3,pka3,pka4,M,f,r)


#defining the fraction of polymer in each state

def fh(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r):
    return Ph(pH,Gh1,pka1,M,r)/Pt(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r)
    
def fa(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r):
    return Pa(pH,Gh2,pka2,M,r)/Pt(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r)

def fd(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r):
    return Pd(pH,Gh3,pka3,pka4,M,f,r)/Pt(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r)

#the sum of the aqueous fractions and hydrophobic fractions will be the intensity measured in an experiment
def I(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r):
     return fh(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r)+fa(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r)
 
 
def theta(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r):
    return (mp.power(10,(pH))/(M*np.array(Pt(pH,Gh1,pka1,Gh2,pka2,Gh3,pka3,pka4,M,f,r)))*(((f*M)*mp.power(10,(-pka3))*mp.power(mp.exp(-Gh3*(r)),f*M)*mp.power((1+(mp.power(10,(pH-pka3)))),((f*M)-1)))+((M)*mp.power(10,(-pka2))*(mp.power(mp.exp(-Gh2*(r)),M))*mp.power((1+(mp.power(10,(pH-pka2)))),((M)-1)))+((M)*mp.power(10,(-pka1))*(mp.power(mp.exp(-Gh1*(r)),M))*mp.power((1+(mp.power(10,(pH-pka1)))),((M)-1)))))

theta3=np.vectorize(theta)
    


#we know need a graph where we can vary the quantities for pka's and Gh's and see the effect on the fractions. Slider would be ideal. Copied from the slider scripts

#defining our pH dimension
pH_range = np.linspace(3, 14, 1401)

# Define initial parameters
Gh1 = 0
Gh2 = 3
Gh3 = 1.7
pka1 = 14
pka2 = 4.56
pka3 = 4.56
pka4 = 14
M = 20
f=0.9
r=1


# Create the figure
fig, (ax1,ax2) = plt.subplots(1,2)

ax1.text(-0.25, 0.95, '(a)', transform=ax1.transAxes, size=8)
ax2.text(-0.25, 0.95, '(b)', transform=ax2.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)

ax1.plot(pH_range, fh(pH_range,  Gh1, pka1,  Gh2, pka2,  Gh3, pka3,  pka4, M, f, r),label='$f_{H}$')
ax1.plot(pH_range, fa(pH_range,  Gh1, pka1,  Gh2, pka2,  Gh3, pka3,  pka4, M, f, r), label='$f_{aq}$')
ax1.plot(pH_range, fd(pH_range,  Gh1, pka1,  Gh2, pka2,  Gh3, pka3,  pka4, M, f, r), label='$f_{D}$')
ax1.plot(pH_range, theta3(pH_range,  Gh1, pka1,  Gh2, pka2,  Gh3, pka3,  pka4, M, f, r), label=r'$\theta$')
ax1.set_ylabel(r'$f$ , $\theta$')
ax1.set_xlabel('$pH$')
ax1.set_xlim(3,13.5)
ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
#ax1.legend(frameon=0)


Gh1 = 0
Gh2 = 3
Gh3 = 1
pka1 = 14
pka2 = 4.56
pka3 = 4.56
pka4 = 14
M = 20
f=0.5
r=1


ax2.plot(pH_range, fh(pH_range,  Gh1, pka1,  Gh2, pka2,  Gh3, pka3,  pka4, M, f, r),label='$f_{H}$')
ax2.plot(pH_range, fa(pH_range,  Gh1, pka1,  Gh2, pka2,  Gh3, pka3,  pka4, M, f, r), label='$f_{aq}$')
ax2.plot(pH_range, fd(pH_range,  Gh1, pka1,  Gh2, pka2,  Gh3, pka3,  pka4, M, f, r), label='$f_{D}$')
ax2.plot(pH_range, theta3(pH_range,  Gh1, pka1,  Gh2, pka2,  Gh3, pka3,  pka4, M, f, r), label=r'$\theta$')
#ax2.set_ylabel(r'$f$ , $\theta$')
ax2.set_xlabel('$pH$')
ax2.set_xlim(3,9)
ax2.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax2.legend(frameon=0)
