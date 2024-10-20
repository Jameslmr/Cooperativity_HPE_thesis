# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:45:28 2022

@author: marti119
"""

#want to see what the effect of dispersity in the value of Gh is on the transition. In this case M=Mt as all the groups have ionizable groups on them!
#I suspect it is more than we might expect.


#will use the gao data from 2014 as an example. SO the polymer we look at is basic

#the distribution of the values of will be binomial



import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from scipy.optimize import curve_fit
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Grab Currrent Time Before Running the Code
start = time.time()

####plotting stuff###

from cycler import cycler

plt.rcParams.update({'font.size':7})
plt.rcParams.update({'figure.autolayout': 1})#forcesa tight lauout
plt.rcParams['legend.fontsize'] =6
plt.rcParams["figure.figsize"] = (12/2.54,13.5/2.54)
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

#this is a function for the distribution of the values of M for each value of Mt spit out. 

def N1_dist_weights(N1,Mt,N1_Mt_ratio,sig_N1):

    #in this case we simply want 
    
    
    return (1/(sig_N1*((2*np.pi)**(1/2))))*exp_f((-(N1-(Mt*N1_Mt_ratio))**2)/(2*(sig_N1**2)))

def f_H(pH,pka,gH_1,gH_2,N1,Mt):
    return power_f(1+(np.transpose(exp_f(-(gH_1*(N1)+gH_2*(Mt-N1))))[:,:,np.newaxis]*np.transpose(power_f((1+(power_f(10,(pH-pka)))),(np.transpose(Mt))))),(-1)) #to make sure it was possible to broadcast everything together a fake thid dimension was given ysin np.newaxis to the (Mt,M) matrix to get (Mt,M,1) and then the (pH,M) matrix was transposed into (M,pH) so on multiplication it would give us a (Mt,M,pH) matrix. REMEMBER that numpy broadcasts from the furthest right axis to left.

def theta(pH,pka,gH_1,gH_2,N1,Mt):
    return (np.transpose((power_f(10,pH-pka)/(1+power_f(10,pH-pka)))*np.transpose(Mt))*(1-power_f(1+(np.transpose(exp_f(-(gH_1*(N1)+gH_2*(Mt-N1))))[:,:,np.newaxis]*np.transpose(power_f((1+(power_f(10,(pH-pka)))),(np.transpose(Mt))))),(-1))))


#will now fit this equation with a naive approach that assumes an effective value of M for the polymer

#hydrophobic partition function
def Ph():
    return 1



#aqueous partition function, r is the value of MH/M
def Pa(pH,gH,pka,M,r):
    return power_f((exp_f(-gH*(r))*(1+(power_f(10,(pH-pka))))),M)


#total partition function

def Pt(pH,gH,pka,M,r):
    return Ph()+Pa(pH,gH,pka,M,r)


#defining the fraction of polymer in each state

def fh(pH,gH,pka,M,r):
    return Ph()/Pt(pH,gH,pka,M,r)

#will set some parameters before the fit
def fh_fit(pH,gH,M):
    pka=4.56
    r=1 #as all groups have a hydrophobic group in essence
    return np.array(fh(pH,gH,pka,M,r),dtype=float)

#polymer input values
av_Mt=20
N1_Mt_ratio=1/2
PDI=1.05
sig_N1= (av_Mt*N1_Mt_ratio**2)**(1/2) #for reference for the binomial distribution the value of sigma will be (av_Mt*M_Mt_ratio**2)**(1/2)
gH_1=6 #for the hydrophobic groups of type 1
gH_2=6 #for the hydrophobic groups of type 2. 
pka=4.56

#in general we will keep Mt on the rows and M on the columns.
Mt=np.transpose(np.array([np.linspace(0,av_Mt*2,av_Mt*2)])+1)
N1=np.array([np.linspace(0,int(np.ceil(av_Mt*N1_Mt_ratio*2)),int(np.ceil(av_Mt*N1_Mt_ratio*2))+1)])
pH=np.transpose(np.array([np.linspace(5,11,1001)]))


#weight factors for Mt and M, I think all these normalizations are kinda redundant but it's nice to be sure.

def calculation(av_Mt,N1_Mt_ratio,PDI,sig_N1,gH_1,gH_2,pka):
    #factors for Mt
    dist_Mt=np.nan_to_num(np.array(Mt_dist_weights(Mt,av_Mt,PDI),dtype=float)) #makes sure that the nan crated for the Mt=0 value doesn't break everything. dodgy
    dist_Mt=dist_Mt/np.sum(dist_Mt)#normalizing for when it becomes not well behaved
    #check
    print(np.sum(dist_Mt))
    
    #factors for M
    dist_N1=N1_dist_weights(N1,Mt,N1_Mt_ratio,sig_N1)
    
    #we know want to remove any values of M>Mt
    dist_N1_2=dist_N1
    for j in range(0,np.size(dist_N1,1)):
        for i in range (0,j):
            dist_N1_2[i,j]=0
    
    #test=np.array(dist_M2,dtype=float)
    
    dist_N1_3=np.transpose(np.transpose(dist_N1_2)/np.array([np.sum(dist_N1_2,axis=1)]))#normalizing
    #test2=np.array(dist_M3,dtype=float)
    #check 
    #print(np.sum(dist_N1_3,axis=1))
    #factors for Mt and M multiply out them together
    dist_Mt_N1=dist_Mt*dist_N1_3
    
    
    
    #check normalization before
    print(np.sum(dist_Mt_N1))
    #normalizing even though it's not really needed as the error is tiny for low values of polydispersity. But actually does become useful for larger values
    dist_Mt_N1_2=dist_Mt_N1/np.sum(dist_Mt_N1)
    print(np.sum(dist_Mt_N1_2))
    

    

    
    #now we need to calculate the fraction of the polymer f_h at each pH for each M and Mt combination.
    
    #defining our pH dimension
    
    
    
    
    #this gives us the fraction of the polymer f_h at each pH for each M and Mt combination.
    fraction_vals=f_H(pH,pka,gH_1,gH_2,N1,Mt)
    
    
    # #we then multiply the fraction by the weight of each Mt M polymer combination.
    # weighted_fraction_vals=np.transpose(np.transpose(fraction_vals)*np.transpose(dist_Mt_M)) #weird transpose stuff is to keep the Mt and M axis at the same location as all the way through. 
    
    weighted_fraction_vals2=np.transpose(np.transpose(fraction_vals)*(dist_Mt_N1_2)) #weird transpose stuff is to keep the Mt and M axis at the same
    
    
    #Finally we must collapse this matrix by summing over all of the fractions at the same pH.
    
    # fraction=np.sum(weighted_fraction_vals,axis=(0,1))
    fraction2=np.sum(weighted_fraction_vals2,axis=(0,1))
    
    #finding the average value of Mt
    av_Mt_calc=np.sum(dist_Mt*Mt)
    
    #ionization fraction section so it can be easily commented out
    #theta_vals=theta(pH,pka,gH,Mt)/(av_Mt) # using the av_Mt requires us to have the Mt dimension be very large as to keep this av_Mt accurat enough for this calculation
    theta_vals=theta(pH,pka,gH_1,gH_2,N1,Mt)/(av_Mt_calc) #alternative is to find the av_Mt separately. Ideally we want this value to be as close as possible to av_Mt but we cannot have infintely large dimensions
    weighted_theta_vals2=np.transpose(np.transpose(theta_vals)*(dist_Mt_N1_2))
    theta2=np.sum(weighted_theta_vals2,axis=(0,1))

    
    
    
    
    
        
    #fiting the calculation fraction with the naive equation fh_fit
    
    bounds=((0,0),(20,20))
    init_parameters=[5,1]
    
    fit_params,fit_cov_params=curve_fit(fh_fit,np.ndarray.flatten(pH),fraction2,bounds=bounds,p0=init_parameters) #need to flatten pH as it it needs 1D arrays in the form (n,) to work for some reason. 
    
    return dist_Mt_N1_2, fraction2,theta2,fit_params



def plotting(av_Mt, N1_Mt_ratio, PDI, sig_N1, gH_1, gH_2, pka):
    
    global counter
    
    output=calculation(av_Mt, N1_Mt_ratio, PDI, sig_N1, gH_1, gH_2, pka)
    dist_Mt_N1_2=output[0]
    fraction_vals=output[1]
    theta_vals=output[2]
    fit_params=output[3]
    
    
    
    if counter==0:
    
        #plotting the Mt and M factors 
        X,Y=np.meshgrid(N1,Mt) #again M are the columns so the x axis in this case.
    
    
        im=ax1.pcolormesh(X,Y,np.array(dist_Mt_N1_2,dtype=float),cmap='Reds')
        fig.colorbar(im, cax=cax, orientation='vertical')
        ax1.set_xlabel('$N_1$')
        ax1.set_ylabel('$M_t$')
        ax1.set_xticks([0,2,4,6,8,10,12,14,16,18,20])
        
        
    
    #plotting the calculated fraction

    #plt.scatter(pH,fraction,label='$M_{t,av}=$ '+str(av_Mt)+' $PDI=$ '+str(PDI)+' $M_{av}=$ '+str(sig_M)+' $\sigma _M =$ '+str(sig_M)+' $g _H =$ '+str(gH))
    #plt.scatter(pH,fraction2,label='$M_{t,av}=$ '+str(av_Mt)+' $PDI=$ '+str(PDI)+'\n'+' $M_{av}=$ '+str(np.round(av_Mt*N1_Mt_ratio,2))+' $\sigma _{N_1} =$ '+str(np.round(sig_N1,2))+' $g _{H,1} =$ '+str(gH_1)+' $g _{H,2} =$ '+str(gH_2))
    ax2.scatter(pH,fraction_vals,label=r'$\Delta g _{H} =$ '+'%.2f'%(gH_2-gH_1),marker=markers())
    ax2.legend()
    # plt.scatter(pH,1-theta_vals,label=r'$\Delta g _{H} =$ '+'%.2f'%(gH_2-gH_1),marker=markers())
    # plt.legend()
    
    
    plt.figure(1)
    #plt.plot(pH,fh_fit(pH,*fit_params),label='$M_{eff}:$ '+'%.2f'%(np.round(fit_params[1],2))+' $g_{H}=$ '+'%.2f'%(np.round(fit_params[0],2)))
    ax2.plot(pH,fh_fit(pH,*fit_params),label='$M_{eff}:$ '+'%.2f'%(np.round(fit_params[1],2)))
    ax2.legend(frameon=0)
    ax2.set_xlim((6.5,8.0))
    ax2.set_xlabel('$pH$')
    ax2.set_ylabel('$f_H$')

    counter+=1

fig,(ax1,ax2) = plt.subplots(2,1)
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cax.set_title(r'$\omega$')
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)
ax1.text(-0.16, 0.95, '(a)', transform=ax1.transAxes, size=8)
ax2.text(-0.16, 0.95, '(b)', transform=ax2.transAxes, size=8)

counter=0

plotting(av_Mt, N1_Mt_ratio, PDI, sig_N1, gH_1, gH_2, pka)
plotting(av_Mt, N1_Mt_ratio, PDI, sig_N1, 5.5, 6.5, pka)
plotting(av_Mt, N1_Mt_ratio, PDI, sig_N1, 5, 7, pka)
plotting(av_Mt, N1_Mt_ratio, PDI, sig_N1, 4.5, 7.5, pka)
plotting(av_Mt, N1_Mt_ratio, PDI, sig_N1, 4, 8, pka)



# Grab Currrent Time After Running the Code
end = time.time()

#Subtract Start Time from The End Time
total_time = end - start
print("\n"+ str(total_time))








