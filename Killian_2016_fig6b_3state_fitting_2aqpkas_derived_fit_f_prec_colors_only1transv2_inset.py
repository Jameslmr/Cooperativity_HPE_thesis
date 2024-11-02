# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 17:15:06 2021

@author: marti119
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import mpmath as mp
from cycler import cycler
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as ticker

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
                

#this can be inserter to make open symbols
#facecolors='none'
#,edgecolors=colours()

#plt.rc("marker",prop_cycle= marker_cycler)
#here we are deriving the parameters from the midpoints plot and assuming their description of x to be correct. We then fit f as this will be dependent on the sharpness of the transition. f affects the transition point so will have to make delta free again and then try to use that f with the midpoints again to get a final graph.

#exactly the same as the other 3 state system but we have explicitely linked M and the Gh which makes the change in M a bit more clear. So we have an r value to represent the hydrophilic to hydrophobic ratio.

#The main thing in this script is to unify the description for the partition function of the polymer in all enviroments under one general formula, M will be inmutable and we will instead change pka and Gh. Will this work? I think the relative change in pka between states will yield the same effect as an artificially constrained effective M

#In this script we will look at 3 coexisting states the "oily", the "aqueous", and an intermediate state, described here as a disk state.

#the parameters that are important n this case are pH, pka, Gh. We assume acidic ionizable groups.

#we will define the patitition functions as P and the sum of all of them Ptot. Derivation comes from assuming a template model for the polymer base (or proton holes) and a enviroment dependen Gh and pka

#the Gh parameters are in KT! So no need for units, all is essentially unitless

#hydrophobic partition function

#total/ionizable groups. so 2 is for a 1 to 1 monomer. set so the hydrophobic parameter is a group specific parameter and not dependent on composition. Easier to see if physical values come out.


#total/ionizable groups. so 2 is for a 1 to 1 monomer. set so the hydrophobic parameter is a group specific parameter and not dependent on composition. Easier to see if physical values come out.








def Ph():
    return 1

#intermediate, disk, partition function

def Pd(pH,Gh3,pka3,M,r):
    return mp.power((mp.exp(-Gh3*(r-1))),(M))*mp.power((1+(mp.power(10,(pH-pka3)))),(M))

#def Pd(pH,Gh3,pka3,pka4,M,f,r):
   # return np.vectorize(Pd2(pH,Gh3,pka3,pka4,M,f,r))




#total partition function

def Pt(pH,Gh3,pka3,M,r):
    return Ph()+Pd(pH,Gh3,pka3,M,r)


#defining the fraction of polymer in each state

def fh(pH,Gh3,pka3,M,r):
    return Ph()/Pt(pH,Gh3,pka3,M,r)
    


def fd(pH,Gh3,pka3,M,r):
    return Pd(pH,Gh3,pka3,M,r)/Pt(pH,Gh3,pka3,M,r)



#these values are derived from the polymer characterisation provided in the paper.
r_array=np.array([2.4,3,4,5])
#r_array=np.array([2,2,2,2])
M_array=np.array([19,31,23,20])
#M_array=np.array([20,30,30,30])

#Actual data
    
#SMA 1.4:1
a=np.array([[4.0193, 0.1577],
[4.9868, 0.115],
[5.2582, 0.1503],
[5.5083, 0.3697],
[5.7431, 0.4262],
[6.0003, 0.5724],
[7.0046, 0.5556],
[8.0241, 0.6804],
[8.9698, 0.6496]])


#SMA 2:1
b=np.array([[4.0193, 0.9808],
[4.9973, 0.9808   ],
[5.2466, 0.983     ],
[5.4996, 4.0318e-3 ],
[6.0055, 6.2715e-3 ],
[6.5187, 3.7914e-3 ],
[6.9952, 8.397e-3  ],
[8.0147, 0.1049    ],
[8.9982, 0.4137    ]])


#SMA 3:1
c=np.array([[4.0193, 0.9808],
[4.9826, 0.9666],
[5.261, 0.9264 ],
[5.5323, 0.924 ],
[5.7596, 0.9122],
[6.013, 0.0393 ],
[6.5188, 0.0274],
[7.0027, 0.032 ],
[7.9933, 0.3007],
[9.0064, 0.6449]])


#SMA 4:1
d=np.array([[4.0193, 0.9808],
[4.9973, 0.9808],
[5.7598, 0.9806],
[5.9868, 0.9192],
[6.5, 0.9002   ],
[7.0112, 0.3575],
[8.0167, 0.6522],
[8.9703, 0.7957]])


    


#SMA 2:1
b1=np.array([[4.0193, 0.9808],
[4.9973, 0.9808   ],
[5.2466, 0.983     ],
[5.4996, 4.0318e-3 ],
[6.0055, 6.2715e-3 ],
[6.5187, 3.7914e-3 ],
[6.9952, 8.397e-3  ]])

b2=np.array([
[8.0147, 0.1049    ],
[8.9982, 0.4137    ]])
            
            


#SMA 3:1
c1=np.array([[4.0193, 0.9808],
[4.9826, 0.9666],
[5.261, 0.9264 ],
[5.5323, 0.924 ],
[5.7596, 0.9122],
[6.013, 0.0393 ],
[6.5188, 0.0274],
[7.0027, 0.032 ]])


c2=np.array([
[7.9933, 0.3007],
[9.0064, 0.6449]])

#SMA 4:1
d1=np.array([[4.0193, 0.9808],
[4.9973, 0.9808],
[5.7598, 0.9806],
[5.9868, 0.9192],
[6.5, 0.9002   ]])




d2=np.array([
[7.0112, 0.3575],
[8.0167, 0.6522],
[8.9703, 0.7957]])


plt.figure(0)

plt.scatter(b1[:,0],b1[:,1],marker=markers(),label='SMA 2:1')
plt.scatter(c1[:,0],c1[:,1],marker=markers(),label='SMA 3:1')
plt.scatter(d1[:,0],d1[:,1],marker=markers(),label='SMA 4:1')
mcount=0
plt.scatter(b2[:,0],b2[:,1],marker=markers(),facecolors='none'
,edgecolors=colours[0])
plt.scatter(c2[:,0],c2[:,1],marker=markers(),facecolors='none'
,edgecolors=colours[1])
plt.scatter(d2[:,0],d2[:,1],marker=markers(),facecolors='none'
,edgecolors=colours[2])
plt.scatter(a[:,0],a[:,1],marker=markers(),facecolors='none'
,edgecolors=colours[3],label='SMA 1.4:1')
plt.legend()



trans1=[5.410000000000000142e+00,
5.800000000000000711e+00,
6.549999999999999822e+00] #this comes from a full fit of the data inclduing the other side of the transition

#note the change in defiinition of r in this case from Mt/M to Mh/M. M+Mh=1 therefore mh+m +1 =r . therefore r'=r-1
plt.figure(1)
plt.scatter([2,3,4],trans1,label='Oil-disk transition',s=100)




#will not fit the 1.4 curve because it is obivously in a different regime to the others and does not even come close to fitting the trend. the questions is of course why?
#so can fit straight lines to both of them.







#note the change in defiinition of r in this case from Mt/M to Mh/M. M+Mh=1 therefore mh+m +1 =r . therefore r'=r-1

def fit_function1(r,Gh3,intr):
    return intr+np.log10((np.exp(Gh3*(r)))-1)

def fit_function2(r,Gh3,intr):
    return intr+0.4343*(Gh3*(r))




p0=[0.01,3]

params_1,covariance_1=optimize.curve_fit(fit_function1,[2,3,4],trans1,p0=p0,bounds=((0),(100)))
params_2,covariance_2=optimize.curve_fit(fit_function2,[2,3,4],trans1,p0=p0,bounds=((0),(100)))




Gh3=params_1[0]
pka3=params_1[1] #worth mentioning that differs subtanstially from pka2 calculated from the solubilization data

x_vals=np.linspace(1,6,61)

plt.figure(1)
#plt.plot(x_vals,fit_function1(x_vals,*params_1),label='Full calculation')
plt.plot(x_vals,fit_function2(x_vals,*params_2),label='Eq. 9')

#plt.title('Midpoint values')
plt.xlabel(r'$\frac{M_H}{M}$ , Styrene to Maleic Acid ratio')
plt.ylabel('$pH(f_H=0.5)$')
plt.legend(frameon=False)
#plt.savefig("SMA_memb_2.pdf")

#now we can make a predicition instead of a fit

fig,ax=plt.subplots(num=5)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)


mcount=0
# plt.scatter(b2[:,0],b2[:,1],marker=markers(),facecolors='none'
# ,edgecolors=colours[0],s=s)
# plt.scatter(c2[:,0],c2[:,1],marker=markers(),facecolors='none'
# ,edgecolors=colours[1],s=s)
# plt.scatter(d2[:,0],d2[:,1],marker=markers(),facecolors='none'
# ,edgecolors=colours[2],s=s)
# plt.scatter(a[:,0],a[:,1],marker=markers(),facecolors='none'
# ,edgecolors=colours[3],label='1.4:1',s=s)
mcount=0
ax.scatter(b1[:,0],b1[:,1],marker=markers(),label='2')
ax.scatter(c1[:,0],c1[:,1],marker=markers(),label='3')
ax.scatter(d1[:,0],d1[:,1],marker=markers(),label='4')
ax.legend()


def fh2(pH,M):
    r=r_array[i]
    Gh3=params_1[0]
    pka3=params_1[1]
    return np.array(fh(pH,Gh3,pka3,M,r),dtype=float)

fh3=np.vectorize(fh2)    


pH_range = np.linspace(2, 12, 1001)
p0=6
bounds=[1,50]
i=1
params_b, params_covariance_b = optimize.curve_fit(fh3, b1[:,0],b1[:,1],bounds=bounds, p0=p0)

i=2
params_c, params_covariance_c = optimize.curve_fit(fh3, c1[:,0],c1[:,1],bounds=bounds, p0=p0)

M_array=np.array([19,19.5,14.8,20])

def fh4(pH,i):
    M=M_array[i]
    r=r_array[i]
    Gh3=params_1[0]
    pka3=params_1[1]
    return np.array(fh(pH,Gh3,pka3,M,r),dtype=float)

fh5=np.vectorize(fh4)    


ax.plot(pH_range,fh5(pH_range,1))
ax.plot(pH_range,fh5(pH_range,2))
ax.plot(pH_range,fh5(pH_range,3))
#ax.plot(pH_range,fh5(pH_range,0))
ax.set_xlim((3.8,10))
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.legend(title="         $N'_S/M$",frameon=False, loc='center left')
plt.text(4.35,0.64,"SMA",fontsize="small")

#plt.title('Predicted transitions')
ax.set_xlabel('$pH$')
ax.set_ylabel('$f_H$')
#plt.savefig("SMA_memb.pdf")

ax2=inset_axes(ax, width='40%', height='40%', loc='lower right',bbox_to_anchor=(0,0.19,1,1), bbox_transform=ax.transAxes)
ax2.scatter([2,3,4],trans1,label='Oil-disk transition')
ax2.plot(x_vals,fit_function2(x_vals,*params_2),label='Eq. 9')
ax2.tick_params(axis='both', which='major')

#plt.title('Midpoint values')
ax2.set_xlabel(r"$N'_S/M$")
ax2.set_ylabel('pH($f_H$=0.5)')
#ax2.legend(frameon=False)




#attempting to add an image to the figure

im = plt.imread('SMA.png') # insert local path of the image.

newax = fig.add_axes([0.62,0.665,0.3,0.3], anchor='SE', zorder=1)
newax.imshow(im)
newax.axis('off')
plt.show()

#%%
plt.savefig("Killian_membrane_sol.png",dpi=3000)

