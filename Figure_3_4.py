_# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 17:15:06 2021

@author: marti119
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import stats
import mpmath as mp
from cycler import cycler
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FormatStrFormatter
plt.rcParams.update({'font.size':7})
plt.rcParams.update({'figure.autolayout': 1})
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







# also no disk state
#total/ionizable groups. so 2 is for a 1 to 1 monomer. set so the hydrophobic parameter is a group specific parameter and not dependent on composition. Easier to see if physical values come out.



#In this case we will set the hydrophobic partition function to 1.
def Ph():
    return 1



#hydrophobic partition function, instead I will add an effective M parameter,e 
#redefined r as MH/M to keep stuff consistent in the paper.

def Pa(pH,Gh2,pka2,M,r,e):
    return mp.power((mp.exp(-Gh2*(r))*(1+(mp.power(10,(pH-pka2))))),e*M)




#def Pd(pH,Gh3,pka3,pka4,M,f,r):
   # return np.vectorize(Pd2(pH,Gh3,pka3,pka4,M,f,r))




#total partition function

def Pt(pH,Gh2,pka2,M,r,e):
    return Ph()+Pa(pH,Gh2,pka2,M,r,e)


#defining the fraction of polymer in each state

def fh(pH,Gh2,pka2,M,r,e):
    return np.array(Ph()/Pt(pH,Gh2,pka2,M,r,e),dtype=float)
    
def fa(pH,Gh2,pka2,M,r,e):
    return np.array(Pa(pH,Gh2,pka2,M,r,e)/Pt(pH,Gh2,pka2,M,r,e),dtype=float)



    
#actual data from https://doi.org/10.1016/j.bpj.2016.09.025

#4:1

#4:1

d_raw=np.array(
[[4.0057, 0.7158   ],
[4.2478, 1.1383   ],
[4.4967, 1.2221   ],
[4.7456, 1.3502   ],
[4.9872, 1.4224   ],
[5.5101, 1.4126   ],
[6.0034, 0.8306   ],
[7.0084, 0.0135   ],
[8.0108, 6.7793e-3],
[9.0024, 6.6435e-3],
[5.2541, 1.4192   ],
[5.7625, 1.3649   ],
[6.2441, 0.0761   ],
[6.5001, 0.076    ],
[6.7525, 0.048    ]]
)




#3:1

c_raw=np.array(
[[4.0208, 1.3372   ],
[4.2769, 1.3651   ] ,
[4.4969, 1.3897   ] ,
[4.7419, 1.2516   ] ,
[5.0045, 0.7305   ] ,
[5.5337, 5.4748e-3] ,
[5.9952, 5.4116e-3] ,
[7.0084, 0.0135   ] ,
[8.0108, 6.7793e-3] ,
[9.0024, 6.6435e-3] ,
[5.2596, 0.0203   ] ,
[5.7536, 0.0153   ]]
)

#2:1

b_raw=np.array(
[[3.9992, 1.4424   ],
[4.2469, 0.4577   ],
[4.5061, 0.0418   ],
[4.7765, 8.8662e-3],
[4.9964, 7.1923e-3],
[5.5337, 5.4748e-3],
[5.9952, 5.4116e-3],
[7.0084, 0.0135   ],
[8.0108, 6.7793e-3],
[9.0024, 6.6435e-3],
[5.2596, 0.0203   ],
[5.7536, 0.0153   ]]
)


#1.4:1

a_raw=np.array(
[[3.9832, 0.0106   ],
[4.2392, 7.296e-3 ],
[4.4988, 2.3289e-3],
[4.7765, 8.8662e-3],
[4.9964, 7.1923e-3],
[5.5337, 5.4748e-3],
[5.9952, 5.4116e-3],
[7.0084, 0.0135   ],
[8.0108, 6.7793e-3],
[9.0024, 6.6435e-3]]
)


mcount=0 
plt.figure(0)
plt.scatter(a_raw[:,0],a_raw[:,1],label='SMA 1.4:1',marker=markers())
plt.scatter(b_raw[:,0],b_raw[:,1],label='SMA 2:1',marker=markers())
plt.scatter(c_raw[:,0],c_raw[:,1],label='SMA 3:1',marker=markers())
plt.scatter(d_raw[:,0],d_raw[:,1],label='SMA 4:1',marker=markers())
plt.xlim((3.5,9.5))
plt.legend()



b2=(b_raw[:,1]-np.min(b_raw[:,1]))/(np.max(b_raw[:,1]-np.min(b_raw[:,1])))
c2=(c_raw[:,1]-np.min(c_raw[:,1]))/(np.max(c_raw[:,1]-np.min(c_raw[:,1])))
d2=(d_raw[:,1]-np.min(d_raw[:,1]))/(np.max(d_raw[:,1]-np.min(d_raw[:,1])))








a=a_raw
b=np.transpose(np.vstack((b_raw[:,0],b2)))
c=np.transpose(np.vstack((c_raw[:,0],c2)))
d=np.transpose(np.vstack((d_raw[:,0],d2)))



mcount=0 
plt.figure(1)
plt.scatter(a[:,0],a[:,1],label='1.4:1',marker=markers())
plt.scatter(b[:,0],b[:,1],label='2:1',marker=markers())
plt.scatter(c[:,0],c[:,1],label='3:1',marker=markers())
plt.scatter(d[:,0],d[:,1],label='4:1',marker=markers())
plt.xlim((3.5,9.5))
plt.legend()





r_array=np.array([1.4,2,3,4])
#r_array=np.array([2,2,2,2])
M_array=np.array([19,31,23,20]) #where did this come from?
#M_array=np.array([20,30,30,30])


    
#defining our pH dimension
pH_range = np.linspace(3, 10, 1001)




def fh2(pH,Gh2,e):
    pka2=2.9
    M=M_array[i]
    r=r_array[i]
    return fh(pH,Gh2,pka2,M,r,e)


fh3=np.vectorize(fh2)

bounds=((0,0),(13,1))

p0=[1,0.8]

i=0
params_a, params_covariance_a = optimize.curve_fit(fh3, a[:,0],a[:,1],bounds=bounds, p0=p0)

i=1
params_b, params_covariance_b = optimize.curve_fit(fh3, b[:,0],b[:,1],bounds=bounds, p0=p0)

i=2
params_c, params_covariance_c = optimize.curve_fit(fh3, c[:,0],c[:,1],bounds=bounds, p0=p0)

i=3
params_d, params_covariance_d = optimize.curve_fit(fh3, d[:,0],d[:,1],bounds=bounds, p0=p0)

mcount=0 
plt.figure(1)
i=0
plt.plot(pH_range,fh3(pH_range,*params_a))
i=1
plt.plot(pH_range,fh3(pH_range,*params_b))
i=2
plt.plot(pH_range,fh3(pH_range,*params_c))
i=3
plt.plot(pH_range,fh3(pH_range,*params_d))
plt.xlim((1,11))



e_array=np.array([params_a[1],params_b[1],params_c[1],params_d[1]])

plt.figure(2)
names=['Gh2','e']
for i in range(0,2):
    plt.plot([1.4,2,3,4],[params_a[i],params_b[i],params_c[i],params_d[i]],'--o',label=names[i])

plt.legend()




#finding midpoints. interpolation was not working very well. Will just extract it directly from the function we've already fit. Seems a bit dodge but actually probably the best way to do it.  Will just ifnd the values of pH that cross 0.5. I've done this before somewhere


#will do the weird shifting a matrix over itself should work quite well. Kinda redundant process in my opinion but good practice nonetheless
#the thing is we know what the general shape for the function is, first transition is always descending, second is ascending
def midpoint_finder(ii,*params):
    global i
    i=ii
    bool_ar=(fh3(pH_range,*params)>0.5).astype(int)
    val_ar=(np.diff(bool_ar)!=0).astype(int)#problem with this is that we then have an indexing array that is one index smaller. a bit of a problem
    vals_index=np.where(val_ar==1)#gives us the index where the sign change is occuring. Can now linear interpolate between that value and an adjacent one but which direction to take the adjacent one in depends on the sign of the gradient. Should we include some logic to know what direction the gradient it is going in. From the bool array we know that if the value of the array at the transition is 1 then it is descending and if it is a 0 it is descending. I don't think this really matters though. What's important is that we separate both transitions. might do it manually for just 4 curves
    trans_pH=pH_range[vals_index]
    return trans_pH





#transition from oil to aq, no transition for series a 
trans1=np.array([midpoint_finder(1,*params_b)[0],midpoint_finder(2,*params_c)[0],midpoint_finder(3,*params_d)[0]])

fig,(ax3,ax1)=plt.subplots(1,2,num=4)

ax3.text(-0.25, 0.95, '(a)', transform=ax3.transAxes, size=8)
ax1.text(-0.25, 0.95, '(b)', transform=ax1.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)

ax1.scatter([2,3,4],trans1,label='Hydrophobic-aqueous transition',marker=markers()) #redefined r as MH/M to keep stuff consistent in the paper.


#plt.legend()



#will not fit the 1.4 curve because it is obivously in a different regime to the others and does not even come close to fitting the trend. the questions is of course why?
#so can fit straight lines to both of them.


e=(params_b[1]+params_c[1]+params_d[1])/3 #this f is the one we extract from the first fit
#we need to define the equations that desribe the transition point



def fit_function(r,g,intr):
    return intr+np.log10((np.exp(g*(r))-1))

def fit_function2(r,g,intr):
    return intr+0.4343*(g*(r)) #redefined r as MH/M to keep stuff consistent in the paper.

p0=[0.5,4.5]

params_1,covariance_1=optimize.curve_fit(fit_function,[2,3,4],trans1,p0=p0,bounds=((0,0),(100,10)))

params_2,covariance_2=optimize.curve_fit(fit_function2,[2,3,4],trans1,p0=p0,bounds=((0,0),(100,10)))

Gh2=params_1[0]



x_vals=np.linspace(1,6,61)



#ax.plot(x_vals,fit_function(x_vals,*params_1),label='Full calculation')
ax1.plot(x_vals,fit_function2(x_vals,*params_2),label='Eq. 9')
#plt.legend()
#plt.title('Midpoint values')
ax1.set_xlabel("$N'_S/M$")
ax1.set_ylabel('pH($f_H=0.5$)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
#plt.legend(frameon=False)
#now we can make a predicition instead of a fit
newax = fig.add_axes([0.66,0.61,0.20,0.34], anchor='SW', zorder=1)
im = plt.imread('Figure_3_4_structure.png') # insert local path of the image.
newax.imshow(im)
newax.axis('off')
plt.show()

# mcount=0 
# fig,ax=plt.subplots(num=4)

# ax.scatter(a[:,0],a[:,1],label='1.4:1',marker=markers())
# ax.scatter(b[:,0],b[:,1],label='2:1',marker=markers())
# ax.scatter(c[:,0],c[:,1],label='3:1',marker=markers())
# ax.scatter(d[:,0],d[:,1],label='4:1',marker=markers())


# def fh4(pH,i):
#     M=M_array[i]
#     r=r_array[i]
#     e=e_array[i]
#     Gh2=params_1[0]
#     pka2=params_1[1]
#     return fh(pH,Gh2,pka2,M,r,e)


# fh5=np.vectorize(fh4)







# ax.plot(pH_range,fh5(pH_range,0))
# ax.plot(pH_range,fh5(pH_range,1))
# ax.plot(pH_range,fh5(pH_range,2))
# ax.plot(pH_range,fh5(pH_range,3))


# #plt.title('Predicted transitions')
# ax.set_xlabel('pH')
# ax.set_ylabel('$f_H$')






ta=np.array([[3.2572, 1.7492  ],
[3.3411, 1.7379  ],
[3.425, 1.7334   ],
[3.5367, 1.7266  ],
[3.6709, 1.713   ],
[3.8218, 1.6947  ],
[3.9616, 1.6676  ],
[4.1348, 1.6381  ],
[4.2634, 1.5995  ],
[4.3753, 1.5587  ],
[4.3977, 1.5088  ],
[4.4872, 1.4657  ],
[4.5655, 1.4248  ],
[4.6493, 1.3795  ],
[4.7612, 1.3341  ],
[4.8115, 1.2887  ],
[4.8451, 1.2411  ],
[5.0352, 1.2025  ],
[5.0632, 1.1775  ],
[5.2085, 1.1526  ],
[5.3091, 1.1367  ],
[5.4098, 1.1185  ],
[5.5718, 1.0981  ],
[5.7339, 1.08    ],
[5.8849, 1.0618  ],
[6.0525, 1.0414  ],
[6.2091, 1.0233  ],
[6.3822, 1.0074  ],
[6.5388, 0.987   ],
[6.6618, 0.9665  ],
[6.8127, 0.9484  ],
[6.9636, 0.928   ],
[7.0587, 0.9144  ],
[7.176, 0.8961   ],
[7.3046, 0.869   ],
[7.4163, 0.8531  ],
[7.5282, 0.8349  ],
[7.6232, 0.8168  ],
[7.7295, 0.7941  ],
[7.8245, 0.776   ],
[7.9195, 0.7555  ],
[8.0369, 0.7306  ],
[8.1655, 0.7056  ],
[8.2605, 0.6829  ],
[8.3556, 0.6602  ],
[8.4506, 0.6398  ],
[8.5455, 0.6194  ],
[8.6574, 0.5967  ],
[8.7525, 0.5763  ],
[8.8363, 0.5535  ],
[8.9369, 0.5309  ],
[9.0376, 0.5105  ],
[9.155, 0.4833   ],
[9.2444, 0.4651  ],
[9.3394, 0.4402  ],
[9.4401, 0.422   ],
[9.5295, 0.4016  ],
[9.6078, 0.3834  ],
[9.6972, 0.3676  ],
[9.7978, 0.3494  ],
[9.8873, 0.329   ],
[9.9767, 0.3131  ],
[10.0941, 0.2881 ],
[10.2339, 0.2655 ],
[10.3624, 0.245  ],
[10.4575, 0.2292 ],
[10.5356, 0.2133 ],
[10.6196, 0.2019 ],
[10.7425, 0.1883 ],
[10.9382, 0.1407 ]])

tb=np.array([[3.3633, 1.8355 ],
[3.4695, 1.8264 ],
[3.6037, 1.8218 ],
[3.7825, 1.8195 ],
[3.9502, 1.806  ],
[4.2073, 1.7788 ],
[4.4029, 1.7492 ],
[4.6265, 1.7153 ],
[4.8054, 1.6812 ],
[4.8893, 1.6449 ],
[5.1129, 1.604  ],
[5.2359, 1.5655 ],
[5.3868, 1.5292 ],
[5.4819, 1.4861 ],
[5.6328, 1.443  ],
[5.7, 1.409     ],
[5.8117, 1.3681 ],
[5.9292, 1.3273 ],
[6.041, 1.2887  ],
[6.1808, 1.2479 ],
[6.2926, 1.2093 ],
[6.4268, 1.1707 ],
[6.4939, 1.1412 ],
[6.5051, 1.0868 ],
[6.6896, 1.055  ],
[6.8462, 1.0187 ],
[7.0083, 0.9733 ],
[7.1089, 0.9348 ],
[7.2711, 0.8961 ],
[7.3717, 0.8599 ],
[7.4667, 0.8327 ],
[7.5674, 0.8032 ],
[7.6512, 0.7737 ],
[7.7463, 0.7442 ],
[7.8469, 0.7147 ],
[7.942, 0.6852  ],
[8.0258, 0.658  ],
[8.1097, 0.6284 ],
[8.1992, 0.5944 ],
[8.2998, 0.5695 ],
[8.4003, 0.5445 ],
[8.5066, 0.5196 ],
[8.6296, 0.4901 ],
[8.7247, 0.4651 ],
[8.8141, 0.4424 ],
[8.9315, 0.4175 ],
[9.0433, 0.3948 ],
[9.1495, 0.3698 ],
[9.2278, 0.3517 ],
[9.3117, 0.3313 ],
[9.4011, 0.3108 ],
[9.4849, 0.2927 ],
[9.58, 0.2768   ],
[9.6918, 0.2609 ],
[9.7812, 0.2405 ],
[9.8707, 0.2269 ],
[9.9601, 0.2065 ],
[10.044, 0.1883 ],
[10.1166, 0.1634],
[10.2396, 0.1475],
[10.3626, 0.1248],
[10.4576, 0.1089],
[10.5471, 0.0976],
[10.6421, 0.0817],
[10.7315, 0.0749],
[10.9383, 0.0408]]
)

tc=np.array([[3.4751, 1.8514 ],
[3.5756, 1.8378 ],
[3.7545, 1.8332 ],
[4.0116, 1.8218 ],
[4.2576, 1.7947 ],
[4.5537, 1.7629 ],
[4.8389, 1.7243 ],
[5.0737, 1.6835 ],
[5.3084, 1.6381 ],
[5.5041, 1.595  ],
[5.6942, 1.5542 ],
[5.8841, 1.5043 ],
[5.9961, 1.4543 ],
[6.147, 1.4158  ],
[6.3538, 1.3727 ],
[6.4825, 1.3296 ],
[6.6167, 1.2796 ],
[6.7675, 1.2365 ],
[6.7955, 1.1865 ],
[6.9577, 1.1435 ],
[7.1031, 1.1027 ],
[7.2373, 1.0528 ],
[7.3827, 1.0096 ],
[7.4889, 0.9665 ],
[7.5448, 0.9257 ],
[7.5673, 0.8917 ],
[7.6679, 0.8576 ],
[7.7629, 0.8304 ],
[7.8468, 0.8054 ],
[7.9363, 0.776  ],
[8.0201, 0.7419 ],
[8.1096, 0.7102 ],
[8.1935, 0.6806 ],
[8.2717, 0.6489 ],
[8.3333, 0.6126 ],
[8.406, 0.5854  ],
[8.4731, 0.5558 ],
[8.5457, 0.5286 ],
[8.6129, 0.499  ],
[8.6911, 0.4742 ],
[8.7638, 0.4424 ],
[8.87, 0.422    ],
[8.9539, 0.3925 ],
[9.0378, 0.3721 ],
[9.1384, 0.3471 ],
[9.2278, 0.3244 ],
[9.3173, 0.2995 ],
[9.4123, 0.2791 ],
[9.5018, 0.2564 ],
[9.58, 0.2337   ],
[9.6862, 0.2201 ],
[9.7757, 0.1929 ],
[9.8428, 0.177  ],
[9.9434, 0.1656 ],
[10.0328, 0.152 ],
[10.1055, 0.1316],
[10.2005, 0.1179],
[10.3738, 0.0953],
[10.2564, 0.1112],
[10.4687, 0.0839],
[10.5638, 0.0726],
[10.6309, 0.0567],
[10.726, 0.0454 ],
[10.8154, 0.0363],
[10.9048, 0.0204]]
)
       
# ax2=ax.twinx()
# mcount=0
# ax2.scatter(ta[:,0],(ta[:,1]-2)*-1,marker=markers(),facecolors='none'
# ,edgecolors=colours[0])
# ax2.scatter(tb[:,0],(tb[:,1]-2)*-1,marker=markers(),facecolors='none'
# ,edgecolors=colours[1])
# ax2.scatter(tc[:,0],(tc[:,1]-2)*-1,marker=markers(),facecolors='none'
# ,edgecolors=colours[2])
# #plt.xlabel('$pH$')
# ax2.set_ylabel(r'$\theta$')
# #plt.title('Titration of SMA variants')
# fig.legend(loc="upper right", bbox_to_anchor=(1,0.615), bbox_transform=ax.transAxes)

mcount=0 


ax3.scatter(a[:,0],1-a[:,1],label='1.4',marker=markers())
ax3.scatter(b[:,0],1-b[:,1],label='2',marker=markers())
ax3.scatter(c[:,0],1-c[:,1],label='3',marker=markers())
ax3.scatter(d[:,0],1-d[:,1],label='4',marker=markers())


def fh4(pH,i):
    M=M_array[i]
    r=r_array[i]
    e=e_array[i]
    Gh2=params_1[0]
    pka2=params_1[1]
    return fh(pH,Gh2,pka2,M,r,e)


fh5=np.vectorize(fh4)







ax3.plot(pH_range,1-fh5(pH_range,0))
ax3.plot(pH_range,1-fh5(pH_range,1))
ax3.plot(pH_range,1-fh5(pH_range,2))
ax3.plot(pH_range,1-fh5(pH_range,3))


#plt.title('Predicted transitions')
ax3.set_xlabel('$pH$')
ax3.set_ylabel('$1-f_H$')
ax3.set_ylim((-0.05,1.05))
ax3.xaxis.set_major_locator(ticker.MultipleLocator(1))

ax2=ax3.twinx()
mcount=0
ax2.scatter(ta[:,0],(ta[:,1]-2)*-1,marker=markers(),facecolors='none'
,edgecolors=colours[0])
ax2.scatter(tb[:,0],(tb[:,1]-2)*-1,marker=markers(),facecolors='none'
,edgecolors=colours[1])
ax2.scatter(tc[:,0],(tc[:,1]-2)*-1,marker=markers(),facecolors='none'
,edgecolors=colours[2])
#plt.xlabel('$pH$')
ax2.set_ylabel(r'$\theta$')
ax2.set_ylim((-0.1,2.1))
#plt.title('Titration of SMA variants')
h1,l1=ax3.get_legend_handles_labels()
h2,l2=ax2.get_legend_handles_labels()
h=h1+h2
l=l1+l2
fig.legend(loc="upper right", bbox_to_anchor=(1,0.615), bbox_transform=ax3.transAxes,frameon=False,title="     $N'_S/M$",handles=h,labels=l)
plt.text(9.8,1.24,"SMA",fontsize="small")
#plt.savefig('SMA_prec.pdf')

plt.savefig('SMA_precip.png', dpi=3000)
