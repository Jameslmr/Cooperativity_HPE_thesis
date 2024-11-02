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
    


#here we are deriving the parameters from the midpoints plot and assuming their description of x to be correct. We then fit f as this will be dependent on the sharpness of the transition. f affects the transition point so will have to make delta free again and then try to use that f with the midpoints again to get a final graph.

#exactly the same as the other 3 state system but we have explicitely linked M and the Gh which makes the change in M a bit more clear. So we have an r value to represent the hydrophilic to hydrophobic ratio.

#The main thing in this script is to unify the description for the partition function of the polymer in all enviroments under one general formula, M will be inmutable and we will instead change pka and Gh. Will this work? I think the relative change in pka between states will yield the same effect as an artificially constrained effective M

#In this script we will look at 3 coexisting states the "oily", the "aqueous", and an intermediate state, described here as a disk state.

#the parameters that are important n this case are pH, pka, Gh. We assume acidic ionizable groups.

#we will define the patitition functions as P and the sum of all of them Ptot. Derivation comes from assuming a template model for the polymer base (or proton holes) and a enviroment dependen Gh and pka

#the Gh parameters are in KT! So no need for units, all is essentially unitless

#hydrophobic partition function

#total/ionizable groups. so 2 is for a 1 to 1 monomer. set so the hydrophobic parameter is a group specific parameter and not dependent on composition. Easier to see if physical values come out.
def Ph():
    return 1




#hydrophobic partition function

def Pa(pH,Gh2,pka2,M,r):
    return mp.exp(-Gh2*(r-1))*mp.power(((1+(mp.power(10,(pH-pka2))))),M)




#intermediate, disk, partition function

def Pd(pH,Gh3,pka3,M,Md,r):
    return (mp.exp(-Gh3*(r-1)))*mp.power((1+(mp.power(10,(pH-pka3)))),(Md))

#def Pd(pH,Gh3,pka3,pka4,M,f,r):
   # return np.vectorize(Pd2(pH,Gh3,pka3,pka4,M,f,r))




#total partition function

def Pt(pH,Gh2,pka2,Gh3,pka3,M,Md,r):
    return Ph()+Pa(pH,Gh2,pka2,M,r)+Pd(pH,Gh3,pka3,M,Md,r)


#defining the fraction of polymer in each state

def fh(pH,Gh2,pka2,Gh3,pka3,M,Md,r):
    return Ph()/Pt(pH,Gh2,pka2,Gh3,pka3,M,Md,r)
    
def fa(pH,Gh2,pka2,Gh3,pka3,M,Md,r):
    return Pa(pH,Gh2,pka2,M,r)/Pt(pH,Gh2,pka2,Gh3,pka3,M,Md,r)

def fd(pH,Gh2,pka2,Gh3,pka3,M,Md,r):
    return Pd(pH,Gh3,pka3,M,Md,r)/Pt(pH,Gh2,pka2,Gh3,pka3,M,Md,r)


#actual data

#100
a=np.array([[6.3006, 0.0133],
[6.4085, 0.0158],
[6.6085, 0.4132],
[6.6245, 0.6764],
[6.7804, 0.9154],
[6.8843, 0.991 ],
[6.9523, 0.9935],
[7.4042, 1.0278]]
)

#73
b=np.array([[6.1806, 3.5266e-3],
[6.2646, 5.9982e-3],
[6.3565, 0.0402   ],
[6.3406, 0.279    ],
[6.3925, 0.5423   ],
[6.4885, 0.5911   ],
[6.4605, 0.7032   ],
[6.6284, 0.8593   ],
[6.5845, 0.986    ],
[6.6564, 0.9885   ],
[6.7084, 0.9983   ],
[6.8484, 1.0056   ],
[6.9723, 0.9935   ]]
)

#58
c=np.array([[5.7048, 0.0106 ],
[5.8247, 3.382e-3],
[5.9047, 0.0132  ],
[5.9847, 0.0522  ],
[6.0526, 0.1497  ],
[6.0886, 0.3813  ],
[6.1406, 0.6421  ],
[6.2086, 0.708   ],
[6.3086, 0.9346  ],
[6.5085, 1.0055  ],
[6.6085, 0.9885  ],
[6.8124, 0.9983  ]])

#49
d=np.array([[5.5088, -1.6214e-3],
[5.6728, 0.0204    ],
[5.7727, 0.1837    ],
[5.8247, 0.1813    ],
[6.0207, 0.7518    ],
[6.0606, 0.8103    ],
[6.1646, 0.9444    ],
[6.2526, 1.0029    ],
[6.4005, 0.9933    ],
[6.4365, 0.9934  ]]
)

#0
e=np.array([[5.5288, 0.988 ],
[5.6488, 1.01  ],
[5.9967, 1.0297],
[5.8607, 0.9929],
[6.9603, 0.9959]])



plt.figure(0)
plt.scatter(a[:,0],a[:,1],label='100',marker=markers())
plt.scatter(b[:,0],b[:,1],label='73',marker=markers())
plt.scatter(c[:,0],c[:,1],label='58',marker=markers())
plt.scatter(d[:,0],d[:,1],label='49',marker=markers())
plt.scatter(e[:,0],e[:,1],label='0',marker=markers())
plt.legend()






#defining our pH dimension
pH_range = np.linspace(4, 7.5, 351)


#need to reduce the parameter space for fitting, we want to fix all pkas and M. For this iteration (number i don' even want to think about it) we will fit delta and f to find an appropriate number for f considering a presuposed value of M. f affects the sharpness.

#setting gh3 to 0 means making it the reference state which makes sense
def II(pH,Gh2,dM):
    Gh3=0
    pka2=4.56
    pka3=4.56
    r=2
    Md=60
    M=Md+dM
    return np.array(fa(pH,Gh2,pka2,Gh3,pka3,M,Md,r),dtype=float)

#note the absolute value pf Md is not relevant here just dM

I2=np.vectorize(II)

bounds=((0,0),(50,40))

p0=[0.008,2]
params_a, params_covariance_a = optimize.curve_fit(I2, a[:,0],a[:,1],bounds=bounds, p0=p0)

plt.figure(0)
plt.plot(pH_range,I2(pH_range,*params_a))


p0=[0.008,2]
params_b, params_covariance_b = optimize.curve_fit(I2, b[:,0],b[:,1],bounds=bounds, p0=p0)

plt.figure(0)
plt.plot(pH_range,I2(pH_range,*params_b))


p0=[0.4,2]
params_c, params_covariance_c = optimize.curve_fit(I2, c[:,0],c[:,1],bounds=bounds, p0=p0)

plt.figure(0)
plt.plot(pH_range,I2(pH_range,*params_c))



p0=[0.007,2]
params_d, params_covariance_d = optimize.curve_fit(I2, d[:,0],d[:,1],bounds=bounds, p0=p0)

plt.figure(0)
plt.plot(pH_range,I2(pH_range,*params_d))



p0=[0.006,2]
params_e, params_covariance_e = optimize.curve_fit(I2, e[:,0],e[:,1],bounds=bounds, p0=p0)

plt.figure(0)
plt.plot(pH_range,I2(pH_range,*params_e))
plt.xlim((4,7.75))



plt.figure(1)
names=['deltaG','dM']
for i in range(0,2):
    plt.plot([100,73,58,49,0],[params_a[i],params_b[i],params_c[i],params_d[i],params_e[i]],'--o',label=names[i])

plt.legend()




#plotting midpoints


#will do the weird shifting a matrix over itself should work quite well. Kinda redundant process in my opinion but good practice nonetheless
#the thing is we know what the general shape for the function is, first transition is always descending, second is ascending
def midpoint_finder(*params):
    bool_ar=(I2(pH_range,*params)>0.5).astype(int)
    val_ar=(np.diff(bool_ar)!=0).astype(int)#problem with this is that we then have an indexing array that is one index smaller. a bit of a problem
    vals_index=np.where(val_ar==1)#gives us the index where the sign change is occuring. Can now linear interpolate between that value and an adjacent one but which direction to take the adjacent one in depends on the sign of the gradient. Should we include some logic to know what direction the gradient it is going in. From the bool array we know that if the value of the array at the transition is 1 then it is descending and if it is a 0 it is descending. I don't think this really matters though. What's important is that we separate both transitions. might do it manually for just 4 curves
    trans_pH=pH_range[vals_index]
    return trans_pH



#transition from aq to disk, no transition for series e 
trans=np.array([midpoint_finder(*params_a)[0],midpoint_finder(*params_b)[0],midpoint_finder(*params_c)[0],midpoint_finder(*params_d)[0]])


plt.figure(4)
plt.scatter([1,0.73,0.58,0.49],trans,label='Disk-aqueous transition')

#plt.legend()



pka=4.56
r=2


ge=(trans[0]-pka)/0.4343

def Gh(x,gm):
    return x*ge+(1-x)*gm

def fit_function(x,gm):
    return pka+np.log10((np.exp(Gh(x,gm)))-1) #why have i not fit the pka???




p0=[0.005]

params_1,covariance_1=optimize.curve_fit(fit_function,[1,0.73,0.58,0.49],trans,p0=p0,bounds=((0),(ge)))

gm=params_1[0]

x_vals=np.linspace(0,1,101)

plt.figure(4)
plt.plot(x_vals,fit_function(x_vals,*params_1),label="Eq.10")
#plt.title('Midpoint values')
plt.legend(frameon=False)
plt.xlabel('$x$, EAA mole fraction ')
plt.ylabel('pH($f_{aq}=0.5$)')
#plt.savefig('Tirrell_memb_2.pdf')








#defining our pH dimension
pH_range = np.linspace(4, 10, 601)


#need to reduce the parameter space for fitting, we want to fix all pkas and M. These leaves the hydrophobic energies and f to be fitted.
#will redefine the function then


def III(pH,x,dM):
    geaa=ge
    gmaa=gm
    Gh2=(x*geaa+(1-x)*gmaa)*dM
    Gh3=0
    pka2=4.56
    pka3=4.56
    Md=50
    M=Md+dM
    r=2
    return fa(pH,Gh2,pka2,Gh3,pka3,M,Md,r)


I3=np.vectorize(III)

#just a reminder about what parameters the equations have (Gh1,Gh2,Gh3)

mcount=0
fig,ax=plt.subplots(num=5)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)
ax.scatter(a[:,0],a[:,1],label='100',marker=markers())
ax.scatter(b[:,0],b[:,1],label='73',marker=markers())
ax.scatter(c[:,0],c[:,1],label='58',marker=markers())
ax.scatter(d[:,0],d[:,1],label='49',marker=markers())
ax.scatter(e[:,0],e[:,1],label='0',marker=markers())


dM=6.5

ax.plot(pH_range,I3(pH_range,1,dM))
ax.plot(pH_range,I3(pH_range,0.73,dM))
ax.plot(pH_range,I3(pH_range,0.58,dM))
ax.plot(pH_range,I3(pH_range,0.49,dM))
ax.plot(pH_range,I3(pH_range,0,dM)) 
ax.set_xlim((4.4,9.4))
leg=plt.legend(frameon=False,title=r'$\%_{EAA}$',loc='center left')

plt.text(4.45,0.7,'P(EAA$_i$-$r$-MAA$_j$)',fontsize="small")

#plt.title('Predicted transitions')
ax.set_xlabel('$pH$')
ax.set_ylabel('$f_{aq}$')
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
#plt.savefig('Tirrell_memb.pdf')


ax2=inset_axes(ax, width='40%', height='40%', loc='lower right',bbox_to_anchor=(0,0.12,1,1), bbox_transform=ax.transAxes)
ax2.scatter([1,0.73,0.58,0.49],trans,label='Disk-aqueous transition')
ax2.plot(x_vals,fit_function(x_vals,*params_1),label="Eq.10")
ax2.tick_params(axis='both', which='major')
#plt.title('Midpoint values')
#plt.legend(frameon=False)
ax2.set_xlabel('$x_{EAA}$')
ax2.set_ylabel('$pH(f_{aq}=0.5)$')
#plt.savefig('Tirrell_memb_2.pdf')


#attempting to add an image to the figure

im = plt.imread('PEAA-PMAA_2.png') # insert local path of the image.

newax = fig.add_axes([0.68,0.61,0.27,0.27], anchor='SE', zorder=1)
newax.imshow(im)
newax.axis('off')
plt.show()

#%%


#plt.savefig('Tirrell_membrane_sol.png',dpi=3000)
