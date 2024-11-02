# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:11:04 2023

@author: marti119
"""

#In this script we want to be able to predict what the transition curve that would be that comes of a fractionated polymer

import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from scipy.optimize import curve_fit
import matplotlib.image as mpimg
from matplotlib import gridspec
import string
#%% ####plotting stuff##

from cycler import cycler

plt.rcParams.update({'font.size':7})
plt.rcParams.update({'figure.autolayout': 1})#forcesa tight lauout
plt.rcParams['legend.fontsize'] =6
plt.rcParams["figure.figsize"] = (12/2.54,9/2.54)
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


#%% defining functions



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

def M_dist_weights(M,Mt,M_Mt_ratio,sig_M):

    #in this case we simply want 
    
    
    return (1/(sig_M*((2*np.pi)**(1/2))))*exp_f((-(M-(Mt*M_Mt_ratio))**2)/(2*(sig_M**2)))

def f_H(pH,pka,gH_1,gH_2,M,Mt):
    return power_f(1+((exp_f(-(gH_1*(M)+gH_2*(Mt-M))))[:,:,np.newaxis]*np.transpose(power_f((1+(power_f(10,(pH-pka)))),(M)))),(-1)) #to make sure it was possible to broadcast everything together a fake thid dimension was given ysin np.newaxis to the (Mt,M) matrix to get (Mt,M,1) and then the (pH,M) matrix was transposed into (M,pH) so on multiplication it would give us a (Mt,M,pH) matrix. REMEMBER that numpy broadcasts from the furthest right axis to left.

def theta(pH,pka,gH_1,gH_2,M,Mt):
    return (np.transpose(power_f(10,pH-pka)/(1+power_f(10,pH-pka)))*np.transpose(M)*(1-power_f(1+((exp_f(-(gH_1*(M)+gH_2*(Mt-M))))[:,:,np.newaxis]*np.transpose(power_f((1+(power_f(10,(pH-pka)))),(M)))),(-1)))) #extra multiplication by M is to get the total number of ionized groups  
#this gives us the fraction of the polymer f_h at each pH for each M and Mt combination.



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
    r=1
    return np.array(fh(pH,gH,pka,M,r),dtype=float)

def plot_heatmap(data):
    X,Y=np.meshgrid(M,Mt) #again M are the columns so the x axis in this case.

    plt.pcolormesh(X,Y,np.array(data,dtype=float),cmap='Reds')
    plt.colorbar()
    plt.xlabel('M')
    plt.ylabel('$M_t$')


#polydisperse and binomial transition calculation
def calculation(pH,av_Mt,M_Mt_ratio,PDI,sig_M,gH_1,gH_2,pka):

    #weight factors for Mt and M, I think all these normalizations are kinda redundant but it's nice to be sure.
    
    #factors for Mt
    dist_Mt=np.nan_to_num(np.array(Mt_dist_weights(Mt,av_Mt,PDI),dtype=float)) #makes sure that the nan crated for the Mt=0 value doesn't break everything
    dist_Mt=dist_Mt/np.sum(dist_Mt)#normalizing for when it becomes not well behaved
    
    #check
    print(np.sum(dist_Mt))
    
    #factors for M
    dist_M=M_dist_weights(M,Mt,M_Mt_ratio,sig_M)
    
    #we know want to remove any values of M>Mt
    dist_M2=dist_M
    for j in range(0,np.size(dist_M,1)):
        for i in range (0,j):
            dist_M2[i,j]=0
    

    
    dist_M3=np.transpose(np.transpose(dist_M2)/np.array([np.sum(dist_M2,axis=1)]))#normalizing

    #factors for Mt and M multiply out them together
    dist_Mt_M=dist_Mt*dist_M3
    
    
    
    #check normalization before
    print(np.sum(dist_Mt_M))
    #normalizing even though it's not really needed as the error is tiny for low values of polydispersity. But actually does become useful for larger values
    dist_Mt_M2=dist_Mt_M/np.sum(dist_Mt_M)
    print(np.sum(dist_Mt_M2))
    
    
     
    #now we need to calculate the fraction of the polymer f_h at each pH for each M and Mt combination.
    
    fraction_vals=f_H(pH,pka,gH_1,gH_2,M,Mt)
    
    # #we then multiply the fraction by the weight of each Mt M polymer combination.
    weighted_fraction_vals2=np.transpose(np.transpose(fraction_vals)*np.transpose(dist_Mt_M2)) #weird transpose stuff is to keep the Mt and M axis at the same
    
    
    #Finally we must collapse this matrix by summing over all of the fractions at the same pH.
    fraction2=np.sum(weighted_fraction_vals2,axis=(0,1))
    
    # #ionization fraction section so it can be easily xommented out
    # theta_vals=theta(pH,pka,gH_1,gH_2,M,Mt)/(av_Mt*M_Mt_ratio)
    # weighted_theta_vals2=np.transpose(np.transpose(theta_vals)*np.transpose(dist_Mt_M2))
    # theta2=np.sum(weighted_theta_vals2,axis=(0,1))

        
    #fiting the calculation fraction with the naive equation fh_fit
    
    bounds=((0,0),(20,20))
    init_parameters=[5,1]
    
    fit_params,fit_cov_params=curve_fit(fh_fit,np.ndarray.flatten(pH),fraction2,bounds=bounds,p0=init_parameters) #need to flatten pH as it it needs 1D arrays in the form (n,) to work for some reason. 
    
    return dist_Mt_M2, fraction2, fit_params,fraction_vals


def calculation2(pH,dist_Mt_M,gH_1,gH_2,pka):

     
    #now we need to calculate the fraction of the polymer f_h at each pH for each M and Mt combination.
    
    #defining our pH dimension
    
    
    dist_Mt_M2=dist_Mt_M/np.sum(dist_Mt_M)
    print(np.sum(dist_Mt_M2))
    
    
    fraction_vals=f_H(pH,pka,gH_1,gH_2,M,Mt)
    
    # #we then multiply the fraction by the weight of each Mt M polymer combination.
    weighted_fraction_vals2=np.transpose(np.transpose(fraction_vals)*np.transpose(dist_Mt_M2)) #weird transpose stuff is to keep the Mt and M axis at the same
    
    #Finally we must collapse this matrix by summing over all of the fractions at the same pH.
    
    fraction2=np.sum(weighted_fraction_vals2,axis=(0,1))
    
    #fiting the calculation fraction with the naive equation fh_fit
    
    bounds=((0,0),(20,20))
    init_parameters=[5,1]
    
    fit_params,fit_cov_params=curve_fit(fh_fit,np.ndarray.flatten(pH),fraction2,bounds=bounds,p0=init_parameters) #need to flatten pH as it it needs 1D arrays in the form (n,) to work for some reason. 
    
    return dist_Mt_M2, fraction2, fit_params


        
def calculation3(pH,av_Mt,M_Mt_ratio,PDI,sig_M,gH_1,gH_2,pka):

    #weight factors for Mt and M, I think all these normalizations are kinda redundant but it's nice to be sure.
    
    #factors for Mt
    dist_Mt=np.nan_to_num(np.array(Mt_dist_weights(Mt,av_Mt,PDI),dtype=float)) #makes sure that the nan crated for the Mt=0 value doesn't break everything
    dist_Mt=dist_Mt/np.sum(dist_Mt)#normalizing for when it becomes not well behaved
    #check
    print(np.sum(dist_Mt))
    
    #factors for M
    dist_M=M_dist_weights(M,Mt,M_Mt_ratio,sig_M)
    
    #we know want to remove any values of M>Mt
    dist_M2=dist_M
    for j in range(0,np.size(dist_M,1)):
        for i in range (0,j):
            dist_M2[i,j]=0
    
    #test=np.array(dist_M2,dtype=float)
    
    dist_M3=np.transpose(np.transpose(dist_M2)/np.array([np.sum(dist_M2,axis=1)]))#normalizing

    dist_Mt_M=dist_Mt*dist_M3
    
    
    
    #check normalization before
    print(np.sum(dist_Mt_M))
    #normalizing even though it's not really needed as the error is tiny for low values of polydispersity. But actually does become useful for larger values
    dist_Mt_M2=dist_Mt_M/np.sum(dist_Mt_M)
    print(np.sum(dist_Mt_M2))
    
    
     
    #now we need to calculate the fraction of the polymer f_h at each pH for each M and Mt combination.
    
    #defining our pH dimension
    
    
    fraction_vals=f_H(pH,pka,gH_1,gH_2,M,Mt)
    
    # #we then multiply the fraction by the weight of each Mt M polymer combination.


    weighted_fraction_vals2=np.transpose(np.transpose(fraction_vals)*np.transpose(dist_Mt_M2)) #weird transpose stuff is to keep the Mt and M axis at the same
    
    
    #Finally we must collapse this matrix by summing over all of the fractions at the same pH.
    
    # fraction=np.sum(weighted_fraction_vals,axis=(0,1))
    fraction2=np.sum(weighted_fraction_vals2,axis=(0,1))
    
    


    return  fraction2

#this is the fit function that we will use to fit the data with

def fit_function_disperse2(pH,gH_1,gH_2):
    av_Mt2=av_Mt
    M_Mt_ratio2=M_Mt_ratio
    PDI2=PDI
    pka2=pka
    sig_M=sig_M=(av_Mt*M_Mt_ratio**2)**(1/2)
    return np.array(calculation3(pH,av_Mt2,M_Mt_ratio2,PDI2,sig_M,gH_1,gH_2,pka2),dtype=float)   


fit_function_disperse=np.vectorize(fit_function_disperse2)


def slicing_func(pH_array,av_Mt, M_Mt_ratio, PDI, sig_M, gH_1, gH_2, pka, pH1,pH2):

    pH=np.array([pH1])
    calculation_vals1=calculation(pH,av_Mt, M_Mt_ratio, PDI, sig_M, gH_1, gH_2, pka)
    fraction_vals_slice_1=calculation_vals1[3]
    dist_Mt_M2_vals=calculation_vals1[0]

#in the current fractionation procedure where we end wth the polymer in the oil phase we start by using the polymer in the water phase, therefore 1-fraction_vals is what we want

    slice_1_dist=np.squeeze(1-fraction_vals_slice_1)*dist_Mt_M2_vals

    slice_1_dist_n=slice_1_dist/np.sum(slice_1_dist)
    print(np.sum(slice_1_dist_n))

#plotting the Mt and M factors
    plt.figure() 
    plot_heatmap(slice_1_dist_n)
    
    #slice 2 calculations

    pH=np.array([pH2])
    calculation_vals2=calculation(pH,av_Mt, M_Mt_ratio, PDI, sig_M, gH_1, gH_2, pka)
    fraction_vals_slice_2=calculation_vals2[3]
    
    #in the current fractionation procedure where we end wth the polymer in the oil phase we start by using the polymer in the water phase, therefore fraction_vals is what we want for the second separation
    
    #for the atom efficiency calculation
    slice_2_dist_eff=np.squeeze(fraction_vals_slice_2)*slice_1_dist
    atom_eff=np.array(np.sum(slice_2_dist_eff),dtype=float)
     
    slice_2_dist=np.squeeze(fraction_vals_slice_2)*slice_1_dist_n
    
    slice_2_dist_n=slice_2_dist/np.sum(slice_2_dist)
    print(np.sum(slice_2_dist_n))
    
    #plotting the Mt and M factors
    plt.figure()
    plot_heatmap(slice_2_dist_n)
    
    calculation_vals_sliced=calculation2(pH_array,slice_2_dist_n, gH_1, gH_2, pka)
      
    fraction_vals=calculation_vals_sliced[1]
    fit_params=calculation_vals_sliced[2]
    
    return fraction_vals, fit_params, atom_eff


#%% making the original data to be fractionated
#polymer input values for the polymer to be fractionated
av_Mt=20
M_Mt_ratio=1/2
PDI=1.04
sig_M=(av_Mt*M_Mt_ratio**2)**(1/2) #for reference for the binomial distribution the value of sigma will be (av_Mt*M_Mt_ratio**2)**(1/2)
gH_1=1.94 # this is the value of gH associated with the M directly.
gH_2=2.48 # this is the value of gH associated with the hydrophobic group without an ionizable group.
pka=4.56

#in general we will keep Mt on the rows and M on the columns.
Mt=np.transpose(np.array([np.linspace(0,av_Mt*2,av_Mt*2)+1])) 
M=np.array([np.linspace(0,int(np.ceil(av_Mt*M_Mt_ratio*2)),int(np.ceil(av_Mt*M_Mt_ratio*2))+1)])
pH=np.transpose(np.array([np.linspace(3,14,101)])) 
pH_range=np.linspace(1,14,521)


#what does the transition look like originally
original_transition=calculation(pH,av_Mt,M_Mt_ratio,PDI,sig_M,gH_1,gH_2,pka)
#plotting it out for interest

dist_Mt_M2=original_transition[0]
fraction_vals=original_transition[1]
fit_params=original_transition[2]

plt.figure(-2)
plot_heatmap(dist_Mt_M2)
plt.xticks([0,2,4,6,8,10,12,14,16,18,20])


#plotting the calculated fraction
plt.figure(-1)
plt.scatter(pH,fraction_vals,label='$\sigma _M =$ '+'%.2f'%(np.round(sig_M,2)),marker=markers())
plt.plot(pH,fh_fit(pH,*fit_params),label='$M_{eff}=$ '+'%.2f'%(np.round(fit_params[1],2)))
plt.legend(frameon=0)
plt.xlabel('$pH$')
plt.ylabel('$f_H$')
plt.xlim((4,9))
 

#%%calculating the fractionation stuff
#so in  way becuase we have the fraction of states in the oil phase in the fraction_vals matrix the fractionantion has been done for us. The approach is to update the weights arrray distMt_Mt wth the slices, then recalculate the transition. 
#so i need to run the calculation with a pH array that is only one value. 

pH1_1=7
pH2_1=6.5

fractionation_result1=slicing_func(pH,av_Mt, M_Mt_ratio, PDI, sig_M, gH_1, gH_2, pka, pH1_1,pH2_1)


pH1_2=6.9
pH2_2=6

fractionation_result2=slicing_func(pH,av_Mt, M_Mt_ratio, PDI, sig_M, gH_1, gH_2, pka, pH1_2,pH2_2)


#%%setting up the figure and plotting it
fig = plt.figure(figsize=(12/2.54,16.5/2.54))

gs = gridspec.GridSpec(3, 2)
ax1 = plt.subplot(gs[:2, :2])
ax3 = plt.subplot(gs[2, 0])
ax2 = plt.subplot(gs[2, 1])

mcount=0
ax2.scatter(pH,fractionation_result1[0],label='$pH_{1,2}:$'+'(%.2f'%pH1_1 + ' , '+'%.2f'%pH2_1+ ')',marker=markers())
ax2.plot(pH_range,fh_fit(pH_range,*fractionation_result1[1]),label='$M_{eff}=$ '+'%.2f'%(np.round(fractionation_result1[1][0],2)))


ax2.scatter(pH,fractionation_result2[0],label='$pH_{1,2}:$'+'(%.2f'%pH1_2 + ' , '+'%.2f'%pH2_2+ ')',marker=markers())
ax2.plot(pH_range,fh_fit(pH_range,*fractionation_result2[1]),label='$M_{eff}=$ '+'%.2f'%(np.round(fractionation_result2[1][0],2)))

#plotting the corresponding lines on the ax3 for one of the fractionation procedures as two gets a bit messy.

ax3.axvline(x=pH1_2,color='#0f94d1',linestyle='dashed')
ax3.annotate('$pH_1$',(pH1_2+0.1,1),color='#0f94d1',fontsize=6)
ax3.axvline(x=pH2_2,color='#0f94d1',linestyle='dashed')
ax3.annotate('$pH_2$',(pH2_2+0.1,1),color='#0f94d1',fontsize=6)

ax2.legend(frameon=0)
ax2.set_xlabel('$pH$')
ax2.set_ylabel('$f_H$')

ax2.set_xlim(5.5,9)
ax2.set_xticks([5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0])

#%% will now add the data from the pbaa titration to the grah

pH_plot=np.array([2.734184513999999844e+00,
2.734184513999999844e+00,
4.229242175000000437e+00,
4.229242175000000437e+00,
4.657825371000000381e+00,
4.657825371000000381e+00,
5.814003295000000016e+00,
5.814003295000000016e+00,
6.272487643999999918e+00,
6.272487643999999918e+00,
6.621334432000000270e+00,
6.621334432000000270e+00,
6.910378912999999734e+00,
6.910378912999999734e+00,
7.099752882999999848e+00,
7.099752882999999848e+00,
7.438632618999999835e+00,
7.438632618999999835e+00,
8.176194398999999891e+00,
8.176194398999999891e+00,
9.003459638000000709e+00,
9.003459638000000709e+00,
9.342339373999999808e+00,
9.342339373999999808e+00,
9.621416804000000766e+00,
9.621416804000000766e+00,
9.920428336000000513e+00,
9.920428336000000513e+00,
1.136565073999999953e+01,
1.136565073999999953e+01
])

f_h_plot=np.array([9.597192432051734334e-01,
9.588388769728207750e-01,
1.004761980984089309e+00,
9.952380190159106910e-01,
1.008987738899382203e+00,
1.035190639305951343e+00,
9.698994781829242129e-01,
9.563738515222335312e-01,
6.880222172423728910e-01,
6.908073758683612731e-01,
3.556279412235490711e-01,
3.591013861766495507e-01,
2.231088132663188173e-01,
2.177785958958927404e-01,
1.723997182828056474e-01,
1.801629477862790873e-01,
1.049316515670518785e-01,
1.008659602394596372e-01,
5.904056087332330144e-02,
5.702372186829721828e-02,
3.469443288407977194e-02,
3.511060601210103105e-02,
1.218186765694529161e-01,
1.250040016646925267e-01,
1.564650894772223744e-02,
1.702308160194637709e-02,
1.250920382879274276e-02,
1.334155008483525404e-02,
-5.362230687966184186e-04,
5.362230687966084439e-04
])

mcount=0
ax3.scatter(pH_plot,f_h_plot,marker=markers())

#want to fit with a disperse model,this section is really slow

p0d=[6,2.5]

params_d, params_covariance_d = curve_fit(fit_function_disperse,pH_plot,np.ndarray.flatten(f_h_plot), bounds=((0,0),(14,5)),p0=p0d)  
    
ax3.plot(pH_range, fit_function_disperse(pH_range, params_d[0],params_d[1]),label='$g_{H,1}$:'+'%.2f'%(round(params_d[0],2))+' $g_{H,2}$:'+'%.2f'%(round(params_d[1],2)))

ax3.set_xlabel('$pH$')
ax3.set_ylabel('$f_H$')
ax3.set_xlim((4,9))
ax3.set_xticks([4,5.0,6.0,7.0,8.0,9])    

#%%plotting the image
img=mpimg.imread('Fractionation_figure_top.png')
ax1.imshow(img, aspect='equal',zorder=100)#,aspect='auto')
ax1.set_xlim(0,2148)
ax1.set_ylim(1548,-600)
ax1.set_xlabel('pH',color='white')

#weird change of colour instead f simply removing the axes so that it doesn't change the easpect ratio of the image. 
ax1.spines['left'].set_color('white')
ax1.spines['right'].set_color('white')
ax1.spines['bottom'].set_color('white')
ax1.spines['top'].set_color('white')
ax1.tick_params(axis='x', colors='white')
ax1.tick_params(axis='y', colors='white')
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)

for n, ax in enumerate((ax3,ax2)):
    
    ax.text(-0.25, 0.95, '('+string.ascii_lowercase[n+1]+')', transform=ax.transAxes, 
            size=8,zorder=100)

ax1.text(-0.105, 0.67, '('+string.ascii_lowercase[0]+')', transform=ax1.transAxes, size=8,zorder=100)

plt.savefig('fractionation_figure_v6.pdf',format='pdf',dpi=1200)