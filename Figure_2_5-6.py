# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 16:53:11 2023

@author: marti119
"""
#this version of the chemical dispersity/polydispersity code might is mostly suited for homopolymers where the chemical dispersity is 0. The ratio of M to Mt is in essence fixed. M=Mt and any scaling can be done via the gH parameter.
import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from scipy.optimize import curve_fit
import time
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# Grab Currrent Time Before Running the Code
start = time.time()

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
    

marker_cycle=(cycler(marker=["o", "v", "*","s","P","d","p","+",]))

colours= ['#000000', '#0f94d1','#f7ad05','#009e74', '#f02495', '#f2e41d', '#0071b2', '#d55c00']
#########

#want to remake the random copolymer with a distribution in both Mt and M in a way tat it is possible to fit the width of the M distrbution and it is compatible with mpmath.
#this has already been done in the allostery_prediction_binomialMpolydisperseMt_v5 file but it is frankly a mess.

#it might be worth to change the mpmath functions into numpy functions that apply the mpmath functions but which then allow for matrices as inputs directly. Again i've done it before but it is messy.

exp_f=np.frompyfunc(mp.exp,1,1)
power_f= np.frompyfunc(mp.power,2,1)

#want to define a function for the distribution of the weights of the length


#In this case Mt is the number of repeating units and PDI is the polydispersity. PDI is not a nice measure of polydispersity but I have an intuitive understanding of it now and it is used EVERYWHERE in the literature
def sig_Mt(av_Mt,PDI):
    return av_Mt*((PDI-1)**(1/2))

#then we need to go from u and s to a discrete distribution
    
def Mt_dist_weights(Mt,av_Mt,PDI):

    #changing the real variables into the lognormal ones, this ensures that the average and standard deviation of the lognormal distribution matches the "real" ones. 
    
    mean=av_Mt
    std=sig_Mt(av_Mt,PDI)
    var=std**2
    mudist=np.log((mean**2)/((var+(mean**2))**(1/2)))
    stdist=(np.log((var/(mean**2))+1))**(1/2)
    
    return (1/(Mt*stdist*(2*np.pi)**(1/2)))*exp_f((-(np.log(Mt)-mudist)**2)/(2*(stdist**2))) 

#################
#optional gaussian Mt instead log-normal. in practive really doesn't make much of a difference
# def Mt_dist_weights2(Mt,av_Mt,PDI):

    
#     return (1/(sig_Mt(av_Mt,PDI)*((2*np.pi)**(1/2))))*exp_f((-(Mt-av_Mt)**2)/(2*(sig_Mt(av_Mt,PDI)**2)))
###################

def f_H(pH,pka,gH,Mt):
    return power_f(1+((exp_f(-gH*(Mt)))*np.transpose(power_f((1+(power_f(10,(pH-pka)))),(np.transpose(Mt))))),(-1))
 #to make sure it was possible to broadcast everything together a fake thid dimension was given ysin np.newaxis to the (Mt,M) matrix to get (Mt,M,1) and then the (pH,M) matrix was transposed into (M,pH) so on multiplication it would give us a (Mt,M,pH) matrix. REMEMBER that numpy broadcasts from the furthest right axis to left.

def theta(pH,pka,gH,Mt):
    return (np.transpose((power_f(10,pH-pka)/(1+power_f(10,pH-pka)))*np.transpose(Mt))*(1-power_f(1+((exp_f(-gH*(Mt)))*np.transpose(power_f((1+(power_f(10,(pH-pka)))),(np.transpose(Mt))))),(-1)))) #extra multiplication by M is to get the total number of ionized groups 


#the following functions are the "naive" fitting functions
#hydrophobic partition function
def Ph():
    return 1



#aqueous partition function, r is the value of MH/M
def Pa(pH,gH,pka,M):
    return power_f((exp_f(-gH)*(1+(power_f(10,(pH-pka))))),M)


#total partition function

def Pt(pH,gH,pka,M):
    return Ph()+Pa(pH,gH,pka,M)


#defining the fraction of polymer in each state

def fh(pH,gH,pka,M):
    return Ph()/Pt(pH,gH,pka,M)

#will set some parameters before the fit
def fh_fit(pH,gH,M):
    pka2=pka
    return np.array(fh(pH,gH,pka2,M),dtype=float)

#polymer inout values
av_Mt=20
PDI=1.05
gH=6
pka=4.56
#in general we will keep Mt on the rows and M on the columns.
Mt=np.transpose(np.array([np.linspace(0,av_Mt*30,av_Mt*30)])+1) #this dimension must be very large when dealing with the most disperse polymers
pH=np.transpose(np.array([np.linspace(5,11,1001)]))

#much more useful to have this in the form of a function

def calculation (av_Mt,PDI,gH,pka):
#weight factors for Mt and M, I think all these normalizations are kinda redundant but it's nice to be sure.

    #factors for Mt
    dist_Mt=np.nan_to_num(np.array(Mt_dist_weights(Mt,av_Mt,PDI),dtype=float)) 
    dist_Mt=dist_Mt/np.sum(dist_Mt)#normalizing for when it becomes not well behaved
#check for well-behavedness
    print(np.sum(dist_Mt))


#this gives us the fraction of the polymer f_h at each pH for each M and Mt combination.
    fraction_vals=f_H(pH,pka,gH,Mt)

# #we then multiply the fraction by the weight of each Mt M polymer combination.
    weighted_fraction_vals2=np.transpose(np.transpose(fraction_vals)*np.transpose(dist_Mt)) #weird transpose stuff is to keep the Mt and M axis at the same

#Finally we must collapse this matrix by summing over all of the fractions at the same pH.
    fraction2=np.sum(weighted_fraction_vals2,axis=(0))

#finding the average value of Mt
    av_Mt_calc=np.sum(dist_Mt*Mt)
    print(av_Mt)
    
#fiting the calculation fraction with the naive equation fh_fit

    bounds=((0,0),(20,100))
    init_parameters=[5,1]

    fit_params,fit_cov_params=curve_fit(fh_fit,np.ndarray.flatten(pH),fraction2,bounds=bounds,p0=init_parameters) #need to flatten pH as it it needs 1D arrays in the form (n,) to work for some reason. 

#ionization fraction section so it can be easily commented out
    theta_vals=theta(pH,pka,gH,Mt)/(av_Mt_calc) #alternative is to find the av_Mt separately. Ideally we want this value to be as close as possible to av_Mt but we cannot have infintely large dimensions
    weighted_theta_vals2=np.transpose(np.transpose(theta_vals)*np.transpose(dist_Mt))
    theta2=np.sum(weighted_theta_vals2,axis=(0))



    return dist_Mt, fraction2, theta2, fit_params 


fig,(ax1,ax2) = plt.subplots(2,1,figsize=(12/2.54,13.5/2.54))
ax1.text(-0.16, 0.95, '(a)', transform=ax1.transAxes, size=8)
ax2.text(-0.16, 0.95, '(b)', transform=ax2.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)

output=calculation(av_Mt, 1.01, 6, pka)
dist_Mt=output[0]
fraction_vals=output[1]
theta_vals=output[2]
fit_params=output[3]

PDI=1.01

ax1.plot(Mt,dist_Mt,label='%.2f'%(PDI))
ax1.legend(frameon=0)
ax1.set_xlabel('$M$')
ax1.set_ylabel(r'$\omega (M)$')
ax1.set_xlim(0,60)




ax2.scatter(pH,fraction_vals,label='\u0110$=$ '+'%.2f'%(PDI),marker=markers())
#plt.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)
ax2.set_xlim(6.8,7.6)
ax2.set_xlabel('$pH$')
ax2.set_ylabel('$f_H$')
ax2.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)

fig2,ax=plt.subplots(num=2)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)
ax.scatter(pH,fraction_vals,label='$f_H$: ' +'\u0110$=$ '+'%.2f'%(PDI),marker=markers())
#plt.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax.set_xlim(6.8,7.6)
ax.set_xlabel('$pH$')
ax.set_ylabel('$f_H$')


ax.scatter(pH,1-theta_vals,label=r'$1- \theta$: ' +'\u0110$=$ '+'%.2f'%(PDI),marker=markers())
plt.ylabel(r'$f_H$ , $1-\theta$')

PDI=1.05
output=calculation(av_Mt, 1.05, 6, pka)
dist_Mt=output[0]
fraction_vals=output[1]
theta_vals=output[2]
fit_params=output[3]


ax1.plot(Mt,dist_Mt,label='%.2f'%(PDI))
ax1.legend(frameon=0)
ax1.set_xlabel('$M$')
ax1.set_ylabel(r'$\omega (M)$')
ax1.set_xlim(0,60)




ax2.scatter(pH,fraction_vals,label='\u0110$=$ '+'%.2f'%(PDI),marker=markers())
#plt.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)
ax2.set_xlim(6.8,7.6)
ax2.set_xlabel('$pH$')
ax2.set_ylabel('$f_H$')
ax2.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)

PDI=1.1
output=calculation(av_Mt, 1.1, 6, pka)
dist_Mt=output[0]
fraction_vals=output[1]
theta_vals=output[2]
fit_params=output[3]


ax1.plot(Mt,dist_Mt,label='%.2f'%(PDI))
ax1.legend(frameon=0)
ax1.set_xlabel('$M$')
ax1.set_ylabel(r'$\omega (M)$')
ax1.set_xlim(0,60)




ax2.scatter(pH,fraction_vals,label='\u0110$=$ '+'%.2f'%(PDI),marker=markers())
#plt.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)
ax2.set_xlim(6.8,7.6)
ax2.set_xlabel('$pH$')
ax2.set_ylabel('$f_H$')
ax2.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)

PDI=1.5
output=calculation(av_Mt, 1.5, 6, pka)
dist_Mt=output[0]
fraction_vals=output[1]
theta_vals=output[2]
fit_params=output[3]


ax1.plot(Mt,dist_Mt,label='%.2f'%(PDI))
ax1.legend(frameon=0)
ax1.set_xlabel('$M$')
ax1.set_ylabel(r'$\omega (M)$')
ax1.set_xlim(0,60)




ax2.scatter(pH,fraction_vals,label='\u0110$=$ '+'%.2f'%(PDI),marker=markers())
#plt.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)
ax2.set_xlim(6.8,7.6)
ax2.set_xlabel('$pH$')
ax2.set_ylabel('$f_H$')
ax2.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)


PDI=2
output=calculation(av_Mt, 2, 6, pka)
dist_Mt=output[0]
fraction_vals=output[1]
theta_vals=output[2]
fit_params=output[3]


ax1.plot(Mt,dist_Mt,label='%.2f'%(PDI))
ax1.legend(frameon=0,title='\u0110')
ax1.set_xlabel('$M$')
ax1.set_ylabel(r'$\omega (M)$')
ax1.set_xlim(0,50)




ax2.scatter(pH,fraction_vals,label='\u0110$=$ '+'%.2f'%(PDI),marker=markers())
#plt.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)
ax2.set_xlim(6.9,7.5)
ax2.set_xlabel('$pH$')
ax2.set_ylabel('$f_H$')
ax2.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax2.legend(frameon=0)



#%% 

output=calculation(av_Mt, 1.5, 6.1, pka)
dist_Mt=output[0]
fraction_vals=output[1]
theta_vals=output[2]
fit_params=output[3]

# plt.figure(0)
# plt.plot(Mt,dist_Mt,label='\u0110$=$ '+'%.2f'%(PDI))
# plt.legend(frameon=0)
# plt.xlabel('$M$')
# plt.ylabel(r'$\omega (M)$')
# plt.xlim(0,60)



# plt.figure(1)
# plt.scatter(pH,fraction_vals,label'\u0110$=$ '+'%.2f'%(PDI),marker=markers())
# #plt.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
# plt.xlim(6.8,7.6)
# plt.xlabel('$pH$')
# plt.ylabel('$f_H$')
# plt.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
# plt.legend(frameon=0)


ax.scatter(pH,fraction_vals,label='$f_H$: ' +'\u0110$=$ '+'%.2f'%(PDI),marker=markers())
#plt.plot(pH,fh_fit(pH,*fit_params),label=' $M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
ax.set_xlim(6.8,7.6)
ax.set_xlabel('$pH$')



ax.scatter(pH,1-theta_vals,label=r'$1- \theta$: ' +'\u0110$=$ '+'%.2f'%(PDI),marker=markers())
ax.set_ylabel(r'$f_H$ , $1-\theta$')
empty = mpatches.Patch(color='white')
marker1 = mlines.Line2D([], [], color='#000000', marker="o", linestyle='None')
marker2 = mlines.Line2D([], [], color='#0f94d1', marker="v", linestyle='None' )
marker3 = mlines.Line2D([], [], color='#f7ad05', marker="*", linestyle='None')
marker4 = mlines.Line2D([], [], color='#009e74', marker="s", linestyle='None' )
purple_triangle = mlines.Line2D([], [], color='purple', marker='^', linestyle='None')
ax.legend(frameon=0,handles=[empty,marker1,marker2,empty,empty,marker3,marker4],labels=['\u0110$=1.05$' , '$f_H$',r'$1-\theta$',' ','\u0110$=1.5$','$f_H$',r'$1-\theta$'])













# Grab Currrent Time After Running the Code
end = time.time()

#Subtract Start Time from The End Time
total_time = end - start
print("\n"+ str(total_time))





