# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 13:39:56 2022

@author: james
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





    
#actual data from https://doi.org/10.1021/ja5053158



#PDPA

a=np.array([[7.4342, -6.276e-5],
[7.1329, 2.0675e-3],
[6.942, 1.0364e-4 ],
[6.6577, 2.2282e-3],
[6.4412, 8.3865e-3],
[6.2926, 0.1241   ],
[6.0923, 0.9477   ],
[5.9055, 0.9923   ],
[5.6806, 0.9924   ],
[5.4812, 0.9884   ],
[5.2945, 0.9885   ]])



#P(DPA60-DBA20)

b=np.array([[6.9844, -1.9391e-3],
[6.7722, -1.8674e-3],
[6.3479, 6.3897e-3 ],
[6.1781, 0.0247    ],
[5.9786, 0.1181    ],
[5.825, 0.9518     ],
[5.5788, 0.9884    ],
[5.3497, 0.9864    ],
[4.8679, 1.0018    ]]
)




#P(DPA40-DBA40)

c=np.array([[4.6803, 0.9943    ],
[5.0102, 0.9759    ],
[5.3529, 0.9966    ],
[5.6139, 0.9341    ],
[5.7993, 0.0588    ],
[5.9713, 0.039     ],
[6.1824, 0.0141    ],
[6.3511, 9.4312e-3 ],
[6.7711, 6.6848e-4 ],
[4.8691, 1.0008    ]])



#P(DPA20-DBA60)

d=np.array([[4.6282, 0.9948  ],
[4.7811, 1.0089  ],
[4.9529, 1.0058  ],
[5.2287, 0.9987  ],
[5.4926, 0.2836  ],
[5.62, 0.0766    ],
[5.7602, 0.0269  ],
[5.932, 0.0187   ],
[6.1421, 0.0156  ],
[6.3458, 9.433e-3]]
)


#PDBA

e=np.array([[4.484, 0.9898      ],
[4.6261, 0.9796     ],
[4.7832, 0.9866     ],
[4.9487, 0.9906     ],
[5.2504, 0.6162     ],
[5.3101, 0.2704     ],
[5.4504, 0.0381     ],
[5.5735, 5.6373e-3  ],
[5.7602, -4.5678e-3 ],
[5.9618, -6.6644e-3 ],
[6.14, -6.3943e-4   ]])


#might be worth renormalizing them as the data point extraction was a little wavy haha
"""
a2=(ar[:,1]-np.min(ar[:,1]))/(np.max(ar[:,1]-np.min(ar[:,1])))
b2=(br[:,1]-np.min(br[:,1]))/(np.max(br[:,1]-np.min(br[:,1])))
c2=(cr[:,1]-np.min(cr[:,1]))/(np.max(cr[:,1]-np.min(cr[:,1])))
d2=(dr[:,1]-np.min(dr[:,1]))/(np.max(dr[:,1]-np.min(dr[:,1])))
e2=(er[:,1]-np.min(er[:,1]))/(np.max(er[:,1]-np.min(er[:,1])))

a=np.transpose(np.vstack((ar[:,0],a2)))
b=np.transpose(np.vstack((br[:,0],b2)))
c=np.transpose(np.vstack((cr[:,0],c2)))
d=np.transpose(np.vstack((dr[:,0],d2)))
e=np.transpose(np.vstack((er[:,0],e2)))
"""

plt.figure(1)
plt.scatter(a[:,0],a[:,1],label='PDPA',marker=markers())
plt.scatter(b[:,0],b[:,1],label='P(DPA60-DBA20)',marker=markers())
plt.scatter(c[:,0],c[:,1],label='P(DPA40-DBA40)',marker=markers())
plt.scatter(d[:,0],d[:,1],label='P(DPA20-DBA60)',marker=markers())
plt.scatter(e[:,0],e[:,1],label='PDBA',marker=markers())
plt.xlim((4,8))
plt.legend()






#will attempt a trial fit, although don't know any parameters, maybe the assymetry is due to pka3 not being equal to pka2.worth a try


#defining our pH dimension
pH_range = np.linspace(3, 10, 1001)




def fa3(pH,Gh2,e):
    M=80
    r=2
    pka2=10.1
    return fa(pH,Gh2,pka2,M,r,e)



fa2=np.vectorize(fa3)






bounds=((0,0),(15,1))

p0=[10,0.5]

params_a, params_covariance_a = optimize.curve_fit(fa2, a[:,0],a[:,1],bounds=bounds, p0=p0)




params_b, params_covariance_b = optimize.curve_fit(fa2, b[:,0],b[:,1],bounds=bounds, p0=p0)


p0=[12,0.5]


params_c, params_covariance_c = optimize.curve_fit(fa2, c[:,0],c[:,1],bounds=bounds, p0=p0)



params_d, params_covariance_d = optimize.curve_fit(fa2, d[:,0],d[:,1],bounds=bounds, p0=p0)





params_e, params_covariance_e = optimize.curve_fit(fa2, e[:,0],e[:,1],bounds=bounds, p0=p0)

plt.figure(1)
i=0
plt.plot(pH_range,fa2(pH_range,*params_a))
i=1
plt.plot(pH_range,fa2(pH_range,*params_b))
i=2
plt.plot(pH_range,fa2(pH_range,*params_c))
i=3
plt.plot(pH_range,fa2(pH_range,*params_d))
i=4
plt.plot(pH_range,fa2(pH_range,*params_e))
plt.xlim((1,11))



e_array=np.array([params_a[1],params_b[1],params_c[1],params_d[1],params_e[1]])

#will set  pdba as the 100 percent fraction

hf=np.array([0,.25,0.5,0.75,1])

plt.figure(2)
names=['Gh2','e']
for i in range(0,2):
    plt.plot(hf,[params_a[i],params_b[i],params_c[i],params_d[i],params_e[i]],'--o',label=names[i])

plt.legend()



#finding midpoints. interpolation was not working very well. Will just extract it directly from the function we've already fit. Seems a bit dodge but actually probably the best way to do it.  Will just ifnd the values of pH that cross 0.5. I've done this before somewhere


#will do the weird shifting a matrix over itself should work quite well. Kinda redundant process in my opinion but good practice nonetheless
#the thing is we know what the general shape for the function is, first transition is always descending, second is ascending
def midpoint_finder(*params):

    bool_ar=(fa2(pH_range,*params)>0.5).astype(int)
    val_ar=(np.diff(bool_ar)!=0).astype(int)#problem with this is that we then have an indexing array that is one index smaller. a bit of a problem
    vals_index=np.where(val_ar==1)#gives us the index where the sign change is occuring. Can now linear interpolate between that value and an adjacent one but which direction to take the adjacent one in depends on the sign of the gradient. Should we include some logic to know what direction the gradient it is going in. From the bool array we know that if the value of the array at the transition is 1 then it is descending and if it is a 0 it is descending. I don't think this really matters though. What's important is that we separate both transitions. might do it manually for just 4 curves
    trans_pH=pH_range[vals_index]
    return trans_pH





#transition from oil to aq, no transition for series a 
trans1=np.array([midpoint_finder(*params_a)[0],midpoint_finder(*params_b)[0],midpoint_finder(*params_c)[0],midpoint_finder(*params_d)[0],midpoint_finder(*params_e)[0]])

mcount=0


#plt.legend()

#have to set the plka becuase if not stuff breaks

pka=10.1 #vslue from science gao paper form monomer dpa 


def Gh2_func(x,a,b):
    return a*x+(1-x)*b
def fit_function(x,a,b):
    return pka-np.log10(np.exp(Gh2_func(x,a,b))-1)






bounds=((0,0),(20,15))
p0=[14,12]

params_1, params_covariance_1 = optimize.curve_fit(fit_function,hf,trans1,bounds=bounds,p0=p0)
r_squared=stats.linregress(hf,trans1)[2]**2    


# def fit_function2(M,a,b,pk):
#     return pk-np.log10(np.exp(Gh2_func(M,a,b))-1)






# bounds=((0,0,0),(14,14,14))
# p0=[11,9,11]

# params_2, params_covariance_2 = optimize.curve_fit(fit_function2, hf,trans1,bounds=bounds,p0=p0)    


hf_range=np.linspace(0,1,101)

fig,ax= plt.subplots(num=3)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)


#plt.savefig("Figure_comparison_2_v2.pdf")
#plt.savefig("Figure_comparison_2_v2.png")
e_av=np.mean(e_array)


def fa4(pH,i):
    a=params_1[0]
    b=params_1[1]
    pka2=pka
    hfval=hf[i]
    Gh2=Gh2_func(hfval,a,b)
    e=e_av
    r=2
    M=80
    return fa(pH,Gh2,pka2,M,r,e)

fa5=np.vectorize(fa4)


mcount=0

#pdpa to pdba
ax.scatter(a[:,0],a[:,1],label='80:0',marker=markers())
ax.scatter(b[:,0],b[:,1],label='60:20',marker=markers())
ax.scatter(c[:,0],c[:,1],label='40:40',marker=markers())
ax.scatter(d[:,0],d[:,1],label='20:60',marker=markers())
ax.scatter(e[:,0],e[:,1],label='  0:80',marker=markers())
ax.set_ylabel(r'$f_{aq}$')
ax.set_xlabel('$pH$')

ax.plot(pH_range,fa5(pH_range,0))
ax.plot(pH_range,fa5(pH_range,1))
ax.plot(pH_range,fa5(pH_range,2))
ax.plot(pH_range,fa5(pH_range,3))
ax.plot(pH_range,fa5(pH_range,4))

ax.set_xlim((2.9,7.1))
ax.set_ylim((-0.05,1.1))




#overlapping the data from figure 11 in the acs supplementary

#titration
c3=np.array([[0,3.8766],
[0.0274,     4.6936],
[0.0533,     5.3352],
[0.0804,     5.5434],
[0.1074,     5.5764],
[0.1354,     5.5929],
[0.1614,     5.6038],
[0.1885,     5.6148],
[0.2164,     5.6313],
[0.2425,     5.6367],
[0.2695,     5.6477],
[0.2966,     5.6532],
[0.3235,     5.6642],
[0.3506,     5.6696],
[0.3776,     5.6806],
[0.4047,     5.6806],
[0.4317,     5.6861],
[0.4586,     5.6971],
[0.4857,     5.6971],
[0.5137,     5.708 ],
[0.5408,     5.7135],
[0.5678,     5.7135],
[0.5948,     5.73  ],
[0.6209,     5.7354],
[0.6489,     5.7464],
[0.6759,     5.7519],
[0.7029,     5.7683],
[0.7299,     5.7738],
[0.757,      5.7848],
[0.784,      5.8012],
[0.811,      5.8286],
[0.838,      5.8615],
[0.8651,     5.9219],
[0.9191,     6.2288],
[0.9461,     6.4099],
[0.9731,     6.5853],
[1.0,        7.3859]]
)

#fluorescence
c4=np.array([[4.6837, 0.9963   ],
[5.0134, 0.9826   ],
[5.3531, 0.995    ],
[5.6128, 0.9367   ],
[5.7991, 0.0608   ],
[6.1837, 0.0161   ],
[6.351, 9.9256e-3 ],
[6.5609, 4.9628e-3],
[6.7707, 0.       ],
[7.4351, 0.       ]]
)




ax2=ax.twinx()
ax2.scatter(c3[:,1],1-c3[:,0],label='40:40 \nTitration',marker=markers(),color='black')
#plt.scatter(c4[:,0],c4[:,1],label=r'P(DPA40-DBA40) fluor',marker=markers())
ax2.set_ylabel(r'$\theta$')
ax2.set_ylim((-0.05,1.1))

ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
fig.legend(loc="upper right", bbox_to_anchor=(1.,0.96), bbox_transform=ax.transAxes,frameon=False,title="     $i:j$")
plt.text(6.15,1.045,"PEO-$b$-P($\mathregular{DPA_{\mathit{i}}}$-$r$-$\mathregular{DBA_{\mathit{j}}}$)",fontsize="small")
#plt.savefig('Figure_comparison_v2_theta.png')

ax3=inset_axes(ax, width='40%', height='40%', loc='lower left',bbox_to_anchor=(0.12,0.12,1,1), bbox_transform=ax.transAxes)
ax3.scatter(hf,trans1,label='Hydrophobic-aqueous transition',marker=markers())
ax3.plot(hf_range,fit_function(hf_range,*params_1),label="Eq. 8")
#ax.plot(hf_range,fit_function2(hf_range,*params_2))
ax3.set_xlabel('$x$')
ax3.set_ylabel('$pH(f_{aq}=0.5)$')
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax.annotate("$R^2=$"+str(np.round(r_squared,2)),(0.19,5.7))
#ax.legend(loc='lower left')
im = plt.imread('Figure_3_1_structure.png') # insert local path of the image.

newax = fig.add_axes([0.11,0.6,0.39,0.8], anchor='SE', zorder=1)
newax.imshow(im)
newax.axis('off')
plt.show()

plt.savefig('Gao2014_corr_fig.png', dpi=3000)
