# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 13:23:27 2023

@author: marti119
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
from scipy.optimize import curve_fit
import mpmath as mp
import matplotlib.patches as mpatches


####plotting stuff###

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
    

marker_cycle=(cycler(marker=["o", "v", "*","s","P","d","p"]))

colours= ['#000000', '#0f94d1','#f7ad05','#009e74', '#f02495', '#f2e41d', '#0071b2', '#d55c00']
#########

#functions used

r=2
pka=4.56
def fit_function(pH,gh,M):
    return 1/(1+((m.exp(-gh*(r-1)))*(1+10**(pH-pka)))**M)

p0a=[10,1]
bounds=((-10,0.2),(14,20))
pH_range=np.linspace(1,14,651)





data2=pd.read_csv('Acetic acid.csv') 
data2=pd.DataFrame.to_numpy(data2)



#scaling the high pH points
data_scaled=(data2[:,3]-data2[:,3][-1])

fig,ax=plt.subplots(num=6)
ax.scatter(data2[:,0],data_scaled,label='Acetic',marker=markers())
pka=4.56
params_a, params_covariance_a = curve_fit(fit_function,data2[:,0],data_scaled, bounds=bounds,p0=p0a)  
ax.plot(pH_range, fit_function(pH_range, params_a[0],params_a[1]),label='$M_{eff}$:'+'%.2f'%(round(params_a[1],2))+'  $g_H$:'+'%.2f'%(round(params_a[0],2))+' $pK_a$:'+'%.2f'%(pka))



data2=pd.read_csv('Salicylic acid.csv') 
data2=pd.DataFrame.to_numpy(data2)


#scaling the high pH points
data_scaled=(data2[:,3]-data2[:,3][-1])


plt.figure(6)
ax.scatter(data2[:,0],data_scaled,label='Salicylic',marker=markers())

pka=2.80
params_a, params_covariance_a = curve_fit(fit_function,data2[:,0],data_scaled, bounds=bounds,p0=p0a)  
ax.plot(pH_range, fit_function(pH_range, params_a[0],params_a[1]),label='$M_{eff}$:'+'%.2f'%(round(params_a[1],2))+'  $g_H$:'+'%.2f'%(round(params_a[0],2))+'  $pK_a$:'+'%.2f'%(pka))


data2=pd.read_csv('Naphtoic acid.csv') 
data2=pd.DataFrame.to_numpy(data2)


#scaling the high pH points
data_scaled=(data2[:,2]-data2[:,2][-1])
plt.figure(6)
ax.scatter(data2[:,0],data_scaled,label='Napthoic',marker=markers())
pka=3.7
params_a, params_covariance_a = curve_fit(fit_function,data2[:,0],data_scaled, bounds=bounds,p0=p0a)  
ax.plot(pH_range, fit_function(pH_range, params_a[0],params_a[1]),label='$M_{eff}$:'+'%.2f'%(round(params_a[1],2))+'  $g_H$:'+'%.2f'%(round(params_a[0],2))+'  $pK_a$:'+'%.2f'%(pka))


data2=pd.read_csv('Naphtoic acid_0.5MNacl.csv') 
data2=pd.DataFrame.to_numpy(data2)


#scaling the high pH points
data_scaled=(data2[:,2]-data2[:,2][-1])
plt.figure(6)
ax.scatter(data2[:,0],data_scaled,label='Napthoic (0.5M NaCl)',marker=markers())
pka=3.7
params_a, params_covariance_a = curve_fit(fit_function,data2[:,0],data_scaled, bounds=bounds,p0=p0a)  
ax.plot(pH_range, fit_function(pH_range, params_a[0],params_a[1]),label='$M_{eff}$:'+'%.2f'%(round(params_a[1],2))+'  $g_H$:'+'%.2f'%(round(params_a[0],2))+'  $pK_a$:'+'%.2f'%(pka))

handles_, labels_ = ax.get_legend_handles_labels()
empty = mpatches.Patch(color='white')
ax.legend(frameon=0,handles=[handles_[0],handles_[1],empty,handles_[2],handles_[3],empty,handles_[4],handles_[5],empty,handles_[6],handles_[7]],labels=[labels_[0],labels_[1],' ',labels_[2],labels_[3],' ',labels_[4],labels_[5],' ',labels_[6],labels_[7]])
plt.xlabel('$pH$')
ax.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14])
ax.set_xlim(1,14)
plt.ylabel('$f_H$')
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)

