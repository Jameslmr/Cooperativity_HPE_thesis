# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 17:46:06 2024

@author: marti119
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import mpmath as mp
from cycler import cycler
from scipy import optimize

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
#exactly the same as the other 3 state system but we have explicitely linked M and the Gh which makes the change in M a bit more clear. So we have an r value to represent the hydrophilic to hydrophobic ratio.

#The main thing in this script is to unify the description for the partition function of the polymer in all enviroments under one general formula, M will be inmutable and we will instead change pka and Gh. Will this work? I think the relative change in pka between states will yield the same effect as an artificially constrained effective M

#In this script we will look at 3 coexisting states the "oily", the "aqueous", and an intermediate state, described here as a disk state.

#the parameters that are important n this case are pH, pka, Gh. We assume acidic ionizable groups.

#we will define the patitition functions as P and the sum of all of them Ptot. Derivation comes from assuming a template model for the polymer base (or proton holes) and a enviroment dependen Gh and pka

#the Gh parameters are in KT! So no need for units, all is essentially unitless

#hydrophobic partition function

#total/ionizable groups. so 2 is for a 1 to 1 monomer

# use of mp math allows for arbitrary floating point precision, in this case we use to have massive values of M not produce overflows.
def Ph2():
    return 1

Ph=np.vectorize(Ph2)


#hydrophobic partition function

def Pa2(pH,Gh2,pka2,M,r):
    return mp.exp(-Gh2*(r-1)*M)*mp.power((1+(mp.power(10,(pH-pka2)))),M)

Pa=np.vectorize(Pa2)


#intermediate, disk, partition function

def Pd2(pH,Gh3,pka3,M,f,r):
    return mp.exp(-Gh3*(r-1)*M)*mp.power((1+(mp.power(10,(pH-pka3)))),(f*M))

#def Pd(pH,Gh3,pka3,pka4,M,f,r):
   # return np.vectorize(Pd2(pH,Gh3,pka3,pka4,M,f,r))

Pd=np.vectorize(Pd2)


#total partition function

def Pt(pH,Gh2,pka2,Gh3,pka3,M,f,r):
    return Ph()+Pa(pH,Gh2,pka2,M,r)+Pd(pH,Gh3,pka3,M,f,r)


#defining the fraction of polymer in each state

def fh(pH,Gh2,pka2,Gh3,pka3,M,f,r):
    return Ph()/Pt(pH,Gh2,pka2,Gh3,pka3,M,f,r)
    
def fa(pH,Gh2,pka2,Gh3,pka3,M,f,r):
    return Pa(pH,Gh2,pka2,M,r)/Pt(pH,Gh2,pka2,Gh3,pka3,M,f,r)

def fd(pH,Gh2,pka2,Gh3,pka3,M,f,r):
    return Pd(pH,Gh3,pka3,M,f,r)/Pt(pH,Gh2,pka2,Gh3,pka3,M,f,r)

#the sum of the aqueous fractions and disk fractions will be the intensity measured in an experiment for a coil to globule transition for example assuming they both scatter equally.

 
def I(pH,Gh2,pka2,Gh3,pka3,M,f,r):
     return (fh(pH,Gh2,pka2,Gh3,pka3,M,f,r)+fd(pH,Gh2,pka2,Gh3,pka3,M,f,r))  
 
def theta(pH,Gh2,pka2,Gh3,pka3,M,f,r):
    return (mp.power(10,(pH))/(M*np.array(Pt(pH,Gh2,pka2,Gh3,pka3,M,f,r)))*((mp.exp(-Gh2*(r-1)*M)*M*mp.power((1+(mp.power(10,(pH-pka2)))),M-1)*mp.power(10,(-pka2)))+(mp.exp(-Gh3*(r-1)*M)*(f*M)*mp.power((1+(mp.power(10,(pH-pka3)))),((f*M)-1))*mp.power(10,(-pka3)))))

theta3=np.vectorize(theta)
    


#we know need a graph where we can vary the quantities for pka's and Gh's and see the effect on the fractions. Slider would be ideal. Copied from the slider scripts

#defining our pH dimension
pH_range = np.linspace(2, 10, 401)



#naive fitting of what we have here

def Ph_n():
    return 1




#hydrophobic partition function

def Pa_n(pH,Gh2,pka2,M,r):
    return mp.power((mp.exp(-Gh2*(r-1))*(1+(mp.power(10,(pH-pka2))))),M)


def Pt_n(pH,Gh2,pka2,M,r):
    return Ph_n()+Pa_n(pH,Gh2,pka2,M,r)


#defining the fraction of polymer in each state

def fh_n(pH,Gh2,pka2,M,r):
    return Ph_n()/Pt_n(pH,Gh2,pka2,M,r)
    
def fa_n(pH,Gh2,pka2,M,r):
    return Pa_n(pH,Gh2,pka2,M,r)/Pt_n(pH,Gh2,pka2,M,r)


def fh_fit2(pH,Gh2,M):
    pka2=1.5
    r=2
    return np.array(fh_n(pH,Gh2,pka2,M,r),dtype=float)


fh_fit=np.vectorize(fh_fit2)

bounds=((0,0),(30,100))

p0=[5,40]

def fit_line(pH,Gh2,pka2,Gh3,pka3,M,f,r):
    params_a, params_covariance_a = optimize.curve_fit(fh_fit, pH_range,I(pH,Gh2,pka2,Gh3,pka3,M,f,r),bounds=bounds, p0=p0)
    return params_a[0],params_a[1]





fig,(ax1,ax2)=plt.subplots(1,2,num=0)
ax1.text(-0.25, 0.95, '(a)', transform=ax1.transAxes, size=8)
ax2.text(-0.25, 0.95, '(b)', transform=ax2.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)


# Define parameters

Gh2 = 4.5
Gh3 = 2
pka2 = 4.5
pka3 = 4.5
M = 40
f=0.5
r=2

ax1.plot(pH_range, fd(pH_range, Gh2,pka2, Gh3,pka3,M,f,r),label='$f_{int}$')
ax1.plot(pH_range, fh(pH_range, Gh2,pka2, Gh3,pka3,M,f,r), label='$f_{H}$')
ax1.plot(pH_range, I(pH_range, Gh2,pka2, Gh3,pka3,M,f,r),  label='$f_{exp}$')
ax1.plot(pH_range, theta3(pH_range, Gh2,pka2, Gh3,pka3,M,f,r), label=r'$\theta$')

params=fit_line(pH_range, Gh2,pka2, Gh3,pka3,M,f,r)

ax1.plot(pH_range,fh_fit(pH_range,params[0],params[1]), label="$M_{eff}=$"+str(np.round(params[1],decimals=1)))

ax1.set_xlim(5.8,7.5)
ax1.set_xlabel('$pH$')
ax1.set_ylabel(r'$f$')
ax1.legend(frameon=0,title='$M_{aq}-M_{int}=20$')



# Define parameters

Gh2 = 4.5
Gh3 = 4
pka2 = 4.5
pka3 = 4.5
M = 40
f=0.9
r=2

ax2.plot(pH_range, fd(pH_range, Gh2,pka2, Gh3,pka3,M,f,r),label='$f_{int}$')
ax2.plot(pH_range, fh(pH_range, Gh2,pka2, Gh3,pka3,M,f,r), label='$f_{H}$')
ax2.plot(pH_range, I(pH_range, Gh2,pka2, Gh3,pka3,M,f,r),  label='$f_{exp}$')
ax2.plot(pH_range, theta3(pH_range, Gh2,pka2, Gh3,pka3,M,f,r), label=r'$\theta$')

params2=fit_line(pH_range, Gh2,pka2, Gh3,pka3,M,f,r)

ax2.plot(pH_range,fh_fit(pH_range,params2[0],params2[1]), label="$M_{eff}=$"+str(np.round(params2[1],decimals=1)))

ax2.set_xlim(5.8,7.5)
ax2.set_xlabel('$pH$')
ax2.set_ylabel(r'$f$')
ax2.legend(frameon=0,title='$M_{aq}-M_{int}=4$')





