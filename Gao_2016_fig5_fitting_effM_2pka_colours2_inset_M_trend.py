# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 17:23:42 2022

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
# also no disk state
#total/ionizable groups. so 2 is for a 1 to 1 monomer. set so the hydrophobic parameter is a group specific parameter and not dependent on composition. Easier to see if physical values come out.
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

def Pa(pH,Gh2,pka2,M,r,e):
    return mp.power((mp.exp(-Gh2*(r-1))*(1+(mp.power(10,(-pH+pka2))))),e*M)




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

def theta(pH,Gh2,pka2,M,r,e):
    return np.array(((mp.power(10,(-pH+pka2))/(1+mp.power(10,(-pH+pka2))))*(1-fh(pH,Gh2,pka2,M,r,e))),dtype=float)




eff_pka_array=np.array([[5,	6.663],
[10,	6.414],
[20,	6.335],
[60,	6.256],
[100,	6.255]]
)
    
#actual data



#5:


a=np.array([[0.615 ,	0.96429],
[0.389	,0.92857   ],
[0.327	,0.89286   ],
[0.298	,0.85714   ],
[0.269	,0.82143   ],
[0.252	,0.78571   ],
[0.223	,0.75      ],
[0.198	,0.71429   ],
[0.155	,0.67857   ],
[0.128	,0.64286   ],
[0.097	,0.60714   ],
[0.066	,0.57143   ],
[0.035	,0.53571   ],
[0	   , 0.5       ],
[-0.037,	0.46429],
[-0.077,	0.42857],
[-0.123,	0.39286],
[-0.174,	0.35714],
[-0.225,	0.32143],
[-0.28	,0.28571   ],
[-0.346,	0.25   ],
[-0.416,	0.21429],
[-0.489,	0.17857],
[-0.572,	0.14286],
[-0.667,	0.10714],
[-0.765,	0.07143],
[-0.882,	0.03971],
[-0.981,	0.00799]]
)


#10:

b=np.array([[0.535	,0.98148],
[0.265	,0.96296],
[0.179	,0.92593],
[0.151	,0.88889],
[0.142	,0.85185],
[0.133	,0.81481],
[0.12	,0.77778],
[0.105	,0.74074],
[0.092	,0.7037 ],
[0.083	,0.66667],
[0.069	,0.62963],
[0.058	,0.59259],
[0.04	,0.55556],
[0.023	,0.51852],
[0	    ,0.48148],
[-0.022	,0.44444],
[-0.053	,0.40741],
[-0.088	,0.37037],
[-0.123	,0.33333],
[-0.162	,0.2963 ],
[-0.217	,0.25926],
[-0.279	,0.22222],
[-0.358	,0.18519],
[-0.444	,0.14815],
[-0.553	,0.11111],
[-0.688	,0.07407],
[-0.873	,0.03704]])

#20:

c=np.array([[0.762	,0.986  ],
[0.344	,0.97265],
[0.222	,0.93568],
[0.149	,0.8988 ],
[0.106	,0.86185],
[0.064	,0.82   ],
[0.044	,0.786  ],
[0.029	,0.74074],
[0.023	,0.7037 ],
[0.016	,0.66667],
[0.01	,0.62963],
[0.006	,0.59259],
[0	    ,0.55556],
[-0.003	,0.51852],
[-0.011	,0.48148],
[-0.018	,0.44444],
[-0.03	,0.40741],
[-0.043	,0.37037],
[-0.063	,0.33333],
[-0.088	,0.2963 ],
[-0.119	,0.25926],
[-0.154	,0.22222],
[-0.204	,0.18519],
[-0.265	,0.14815],
[-0.343	,0.11111],
[-0.448	,0.07407],
[-0.573	,0.03704],
[-0.791	,0.01852]])

#60:

d=np.array([[0.85	,0.985   ],
[0.626	,0.98305 ],
[0.291	,0.9661  ],
[0.153	,0.94915 ],
[0.103	,0.9322  ],
[0.079	,0.91525 ],
[0.063	,0.88136 ],
[0.051	,0.84746 ],
[0.043	,0.81356 ],
[0.038	,0.77966 ],
[0.034	,0.74576 ],
[0.029	,0.71186 ],
[0.024	,0.67797 ],
[0.019	,0.64407 ],
[0.015	,0.61017 ],
[0.008	,0.57627 ],
[0	    ,0.54237 ],
[-0.009	,0.50847 ],
[-0.017	,0.47458 ],
[-0.023	,0.44068 ],
[-0.029	,0.40678 ],
[-0.037	,0.37288 ],
[-0.045	,0.33898 ],
[-0.057	,0.30508 ],
[-0.073	,0.27119 ],
[-0.086	,0.23729 ],
[-0.107	,0.20339 ],
[-0.124	,0.16949 ],
[-0.171	,0.13559 ],
[-0.22	,0.10169 ],
[-0.309	,0.0678  ],
[-0.512	,0.0339  ],
[-0.95	,0.0085  ]])

#100:

e=np.array([[0.93	  ,0.982    ],
[0.225	  ,0.97143  ],
[0.049	  ,0.94286  ],
[0.028	  ,0.91429  ],
[0.024	  ,0.88571  ],
[0.019	  ,0.85714  ],
[0.015	  ,0.82857  ],
[0.013	  ,0.8      ],
[0.012	  ,0.77143  ],
[0.011	  ,0.74286  ],
[0.01	  ,0.71429  ],
[0.009	  ,0.68571  ],
[0.008	  ,0.65714  ],
[0.007	  ,0.62857  ],
[0.006	  ,0.6      ],
[0.005	  ,0.57143  ],
[0.004	  ,0.54286  ],
[0.003	  ,0.51429  ],
[0.002	  ,0.48571  ],
[1.00E-03	,0.45714],
[0	      ,0.42857  ],
[-1.00E-03,	0.4     ],
[-0.002	  ,0.37143  ],
[-0.003	  ,0.34286  ],
[-0.004	  ,0.31429  ],
[-0.005	  ,0.28571  ],
[-0.008	  ,0.25714  ],
[-0.011	  ,0.22857  ],
[-0.017	  ,0.2      ],
[-0.029	  ,0.17143  ],
[-0.037	  ,0.14286  ],
[-0.07	  ,0.11429  ],
[-0.11	  ,0.08571  ],
[-0.217	  ,0.051    ],
[-0.433	  ,0.02857  ],
[-0.779	  ,0.01     ]])

a=np.transpose(np.vstack((eff_pka_array[0,1]-a[:,0],a[:,1])))
b=np.transpose(np.vstack((eff_pka_array[1,1]-b[:,0],b[:,1])))
c=np.transpose(np.vstack((eff_pka_array[2,1]-c[:,0],c[:,1])))
d=np.transpose(np.vstack((eff_pka_array[3,1]-d[:,0],d[:,1])))
e=np.transpose(np.vstack((eff_pka_array[4,1]-e[:,0],e[:,1])))

plt.figure(1)
plt.scatter(a[:,0],a[:,1],label='5',marker=markers())
plt.scatter(b[:,0],b[:,1],label='10',marker=markers())
plt.scatter(c[:,0],c[:,1],label='20',marker=markers())
plt.scatter(d[:,0],d[:,1],label='60',marker=markers())
plt.scatter(e[:,0],e[:,1],label='100',marker=markers())
plt.xlim((4,8))
plt.legend()




M_array=np.array([5,10,20,60,100])


#will attempt a trial fit, although don't know any parameters, maybe the assymetry is due to pka3 not being equal to pka2.worth a try


#defining our pH dimension
pH_range = np.linspace(3, 10, 1001)




def theta3(pH,Gh2,pka2,e):
    M=M_array[i]
    r=2
    return theta(pH,Gh2,pka2,M,r,e)



theta2=np.vectorize(theta3)






bounds=((0,4,0),(15,14,1))

p0=[8,10,0.5]

i=0
params_a, params_covariance_a = optimize.curve_fit(theta2, a[:,0],a[:,1],bounds=bounds, p0=p0)



i=1
params_b, params_covariance_b = optimize.curve_fit(theta2, b[:,0],b[:,1],bounds=bounds, p0=p0)




i=2
params_c, params_covariance_c = optimize.curve_fit(theta2, c[:,0],c[:,1],bounds=bounds, p0=p0)



i=3
params_d, params_covariance_d = optimize.curve_fit(theta2, d[:,0],d[:,1],bounds=bounds, p0=p0)




i=4
params_e, params_covariance_e = optimize.curve_fit(theta2, e[:,0],e[:,1],bounds=bounds, p0=p0)

plt.figure(1)
i=0
plt.plot(pH_range,theta2(pH_range,*params_a))
i=1
plt.plot(pH_range,theta2(pH_range,*params_b))
i=2
plt.plot(pH_range,theta2(pH_range,*params_c))
i=3
plt.plot(pH_range,theta2(pH_range,*params_d))
i=4
plt.plot(pH_range,theta2(pH_range,*params_e))
plt.xlim((1,11))



e_array=np.array([params_a[2],params_b[2],params_c[2],params_d[2],params_e[2]])

plt.figure(2)
names=['Gh2','pka2','e']
for i in range(0,3):
    plt.plot([5,10,20,60,100],[params_a[i],params_b[i],params_c[i],params_d[i],params_e[i]],'--o',label=names[i])

plt.legend()



#finding midpoints. interpolation was not working very well. Will just extract it directly from the function we've already fit. Seems a bit dodge but actually probably the best way to do it.  Will just ifnd the values of pH that cross 0.5. I've done this before somewhere


#will do the weird shifting a matrix over itself should work quite well. Kinda redundant process in my opinion but good practice nonetheless
#the thing is we know what the general shape for the function is, first transition is always descending, second is ascending
def midpoint_finder(ii,*params):
    global i
    i=ii
    bool_ar=(theta2(pH_range,*params)>0.5).astype(int)
    val_ar=(np.diff(bool_ar)!=0).astype(int)#problem with this is that we then have an indexing array that is one index smaller. a bit of a problem
    vals_index=np.where(val_ar==1)#gives us the index where the sign change is occuring. Can now linear interpolate between that value and an adjacent one but which direction to take the adjacent one in depends on the sign of the gradient. Should we include some logic to know what direction the gradient it is going in. From the bool array we know that if the value of the array at the transition is 1 then it is descending and if it is a 0 it is descending. I don't think this really matters though. What's important is that we separate both transitions. might do it manually for just 4 curves
    trans_pH=pH_range[vals_index]
    return trans_pH





#transition from oil to aq, no transition for series a 
trans1=np.array([midpoint_finder(0,*params_a)[0],midpoint_finder(1,*params_b)[0],midpoint_finder(2,*params_c)[0],midpoint_finder(3,*params_d)[0],midpoint_finder(4,*params_e)[0]])

mcount=0

#plt.legend()

#have to set the plka becuase if not stuff breaks

pka=10.1


def Gh2_func(M,a,b):
    return a*(1-((b)/M))
def fit_function(M,a,b):
    return pka-np.log10(np.exp(Gh2_func(M,a,b))-1)






bounds=((0,0),(20,10))
p0=[4,2]

params_1, params_covariance_1 = optimize.curve_fit(fit_function, M_array,trans1,bounds=bounds,p0=p0)    


def fit_function2(M,a,b,pk):
    return pk-np.log10(np.exp(Gh2_func(M,a,b))-1)






bounds=((0,0,4),(20,10,12))
p0=[4,2,10]

params_2, params_covariance_2 = optimize.curve_fit(fit_function2, M_array,trans1,bounds=bounds,p0=p0)    


M_range=np.linspace(0,150,151)

fig,(ax3,ax1)=plt.subplots(1,2,num=3)

ax3.text(-0.25, 0.95, '(a)', transform=ax3.transAxes, size=8)
ax1.text(-0.25, 0.95, '(b)', transform=ax1.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)


ax1.scatter([5,10,20,60,100],trans1,label='Hydrophobic-aqueous transition',marker=markers())

ax1.plot(M_range,fit_function(M_range,*params_1),label="Eq. 8")
#plt.plot(M_range,fit_function2(M_range,*params_2))
ax1.set_xlabel('$i$')
ax1.set_ylabel(r'$pH(\theta=0.5)$')
#ax.legend()
im = plt.imread('PEG-pdpa.png') # insert local path of the image.

newax = fig.add_axes([0.1,0.4,0.3,0.34], anchor='SW', zorder=1)
newax.imshow(im)
newax.axis('off')
plt.show()


#plt.savefig('Gao_2_v2.pdf')

def theta4(pH,i):
    a=params_1[0]
    b=params_1[1]
    pka2=pka
    M=M_array[i]
    Gh2=Gh2_func(M,a,b)
    e=e_array[i]
    r=2
    return theta(pH,Gh2,pka2,M,r,e)

theta5=np.vectorize(theta4)


mcount=0


ax3.scatter(a[:,0],a[:,1],label='5',marker=markers())
ax3.scatter(b[:,0],b[:,1],label='10',marker=markers())
ax3.scatter(c[:,0],c[:,1],label='20',marker=markers())
ax3.scatter(d[:,0],d[:,1],label='60',marker=markers())
ax3.scatter(e[:,0],e[:,1],label='100',marker=markers())
ax3.set_ylabel(r'$\theta$')
ax3.set_xlabel('$pH$')


ax3.plot(pH_range,theta5(pH_range,0))
ax3.plot(pH_range,theta5(pH_range,1))
ax3.plot(pH_range,theta5(pH_range,2))
ax3.plot(pH_range,theta5(pH_range,3))
ax3.plot(pH_range,theta5(pH_range,4))

ax3.set_xlim((3.8,7.88))
ax3.legend(bbox_to_anchor=(0.7,0.4),frameon=False,title="$i$", bbox_transform=ax3.transAxes)
ax3.text(6.55,0.94,"PEO-$b$-P(DPA)$_i$",fontsize="small")



ax2=inset_axes(ax1, width='70%', height='65%', loc='upper right',bbox_to_anchor=(0.2,0.1,0.78,0.9), bbox_transform=ax1.transAxes)
ax2.scatter(M_array,e_array*M_array)
ax2.set_xlabel('$i$')
ax2.set_ylabel('$M_{eff}$')
#%%

plt.savefig('Gao2016_l_depen_fig.png', dpi=3000)