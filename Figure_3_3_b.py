# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 18:38:42 2021

@author: james
"""



import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import mpmath as mp
from cycler import cycler


plt.rcParams.update({'font.size':20})
plt.rcParams.update({'figure.autolayout': 1})
plt.rcParams['legend.fontsize'] = 'small'
plt.rcParams["figure.figsize"] = (8,6)


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


#only one data set so can't look at midpoints

# also no disk state
#total/ionizable groups. so 2 is for a 1 to 1 monomer. set so the hydrophobic parameter is a group specific parameter and not dependent on composition. Easier to see if physical values come out.



#in this version of the script we will actually take an effective M, this means leaving Ph as 1
def Ph():
    return 1




#hydrophobic partition function

def Pa(pH,Gh2,pka2,M,r):
    return mp.power((mp.exp(-Gh2*(r-1))*(1+(mp.power(10,(pH-pka2))))),M)




#def Pd(pH,Gh3,pka3,pka4,M,f,r):
   # return np.vectorize(Pd2(pH,Gh3,pka3,pka4,M,f,r))




#total partition function

def Pt(pH,Gh2,pka2,M,r):
    return Ph()+Pa(pH,Gh2,pka2,M,r)


#defining the fraction of polymer in each state

def fh(pH,Gh2,pka2,M,r):
    return np.array(Ph()/Pt(pH,Gh2,pka2,M,r),dtype=float)
    
def fa(pH,Gh2,pka2,M,r):
    return np.array(Pa(pH,Gh2,pka2,M,r)/Pt(pH,Gh2,pka2,M,r),dtype=float)



#need to scale the data with the axis coordinates
#format of axis coordinate matrices is like a clockwise square from the origin.

#data from https://doi.org/10.1021/ma00196a038

a=np.array([[0.0335, 4.791   ],
[0.055, 5.1122   ],
[0.0724, 5.279   ],
[0.0898, 5.4458  ],
[0.1093, 5.5868  ],
[0.1338, 5.7532  ],
[0.1573, 5.8554  ],
[0.1849, 5.9445  ],
[0.2114, 6.0594  ],
[0.242, 6.1356   ],
[0.2706, 6.2504  ],
[0.29, 6.3012    ],
[0.3135, 6.4162  ],
[0.335, 6.5312   ],
[0.3595, 6.6979  ],
[0.385, 6.8128   ],
[0.4086, 6.9923  ],
[0.4341, 7.1587  ],
[0.4536, 7.3512  ],
[0.4705, 7.4601  ],
[0.4858, 7.5948  ],
[0.499, 7.7359   ],
[0.5165, 7.8705  ],
[0.5304, 8.0568  ],
[0.5472, 8.2301  ],
[0.5606, 8.4099  ],
[0.5764, 8.5381  ],
[0.5913, 8.7308  ],
[0.6082, 8.9104  ],
[0.6215, 9.0774  ],
[0.6374, 9.2185  ],
[0.6512, 9.3982  ],
[0.6644, 9.4944  ],
[0.6835, 9.687   ],
[0.7024, 9.8537  ],
[0.7264, 10.0202 ],
[0.7473, 10.1547 ],
[0.773, 10.3341  ],
[0.797, 10.4491  ],
[0.8159, 10.5899 ],
[0.8322, 10.7118 ],
[0.8496, 10.8013 ],
[0.8639, 10.8652 ],
[0.8833, 10.9353 ],
[0.8997, 11.0119 ]])

b=np.array([[6.4007, 12.3018],
[6.0998, 24.3453 ],
[5.9529, 47.0753 ],
[5.8995, 47.6124 ],
[5.7457, 54.5699 ],
[5.4985, 74.8984 ],
[5.204, 73.3075  ],
[6.9026, 11.7448 ],
[7.1971, 11.999  ],
[7.6991, 11.7093 ]])

c=np.array([[7.7964, 6.0935 ],
[6.8993, 5.8695 ],
[6.505, 6.2957  ],
[6.1497, 5.8845 ],
[5.8936, 4.4631 ],
[5.6524, 4.1979 ],
[5.3048, 4.5708 ]])

       
fig, ax1 = plt.subplots()
fig.subplots_adjust(right=0.75)
ax2 = ax1.twinx()
ax3= ax1.twinx()
ax3.spines.right.set_position(("axes", 1.2))

p1=ax1.scatter(a[:,1],a[:,0],marker=markers())
p2=ax2.scatter(b[:,0],b[:,1],color='orange')
p3=ax3.scatter(c[:,0],c[:,1],color='green')

ax1.set_xlabel('pH')
ax1.set_ylabel('ionization')
ax2.set_ylabel('fluorescent Intensity')
ax3.set_ylabel('Rh')



ax1.yaxis.label.set_color('blue')
ax2.yaxis.label.set_color('orange')
ax3.yaxis.label.set_color('green')

ax1.tick_params(axis='y', colors='blue')
ax2.tick_params(axis='y', colors='orange')
ax3.tick_params(axis='y', colors='green')
ax1.tick_params(axis='x')
plt.title('Coil-globule PEAA')



#scaling the fluorescence and Rh between 0 and 1, assuming both edge cases represent it being in one state or the other
#will also scale the ionization data as I think there is apossibly 3 states or a continuum of pka's so will fit it up to pH 7.5

a_cut=a[a[:,1]<7.5]
a_scaled=np.transpose(np.vstack((a_cut[:,1],(((a_cut[:,0]-(np.min(a_cut[:,0])))/((np.max(a_cut[:,0]))-np.min(a_cut[:,0])))))))

b_scaled=np.transpose(np.vstack((b[:,0],(((b[:,1]-(np.min(b[:,1])))/((np.max(b[:,1]))-np.min(b[:,1])))))))

c_scaled=np.transpose(np.vstack((c[:,0],(((c[:,1]-(np.min(c[:,1])))/((np.max(c[:,1]))-np.min(c[:,1])))))))

c_scaled2=np.transpose(np.vstack((c[:,0],(((c[:,1]-(np.max(c[:,1])))/-((np.max(c[:,1]))-np.min(c[:,1])))))))



fig, ax1 = plt.subplots(num=2)
fig.subplots_adjust(right=0.75)
ax2 = ax1.twinx()
ax3= ax1.twinx()
ax3.spines.right.set_position(("axes", 1.2))

p1=ax1.scatter(a_scaled[:,0],a_scaled[:,1],color='blue')
p2=ax2.scatter(b_scaled[:,0],b_scaled[:,1],color='orange')
p3=ax3.scatter(c_scaled2[:,0],c_scaled2[:,1],color='green')

ax1.set_xlabel('pH')
ax1.set_ylabel('ionization')
ax2.set_ylabel('scaled fluorescent Intensity')
ax3.set_ylabel('scaled 1/Rh')

ax1.yaxis.label.set_color('blue')
ax2.yaxis.label.set_color('orange')
ax3.yaxis.label.set_color('green')


ax1.set_ylim(-0.05,1.05)
ax2.set_ylim(-0.05,1.05)
ax3.set_ylim(-0.05,1.05)
ax1.set_xlim(4.5,9)


ax1.tick_params(axis='y', colors='blue')
ax2.tick_params(axis='y', colors='orange')
ax3.tick_params(axis='y', colors='green')
ax1.tick_params(axis='x')
plt.title('Coil-globule PEAA')






def fh2(pH,Gh2,M):
    pka2=4.5
    r=2
    return fh(pH,Gh2,pka2,M,r)


fh3=np.vectorize(fh2)

pHrange=np.linspace(0,14,1401)

p0b=[4,10]

params_b, params_covariance_b = optimize.curve_fit(fh3, b_scaled[:,0],b_scaled[:,1],bounds=((3,1),(10,60)), p0=p0b)

params_c, params_covariance_c = optimize.curve_fit(fh3, c_scaled2[:,0],c_scaled2[:,1],bounds=((3,1),(10,60)), p0=p0b)


fig,ax=plt.subplots(num=3)

mcount=0
ax.scatter(b_scaled[:,0],b_scaled[:,1],marker=markers(),label='Fluorescence')
ax.scatter(c_scaled2[:,0],c_scaled2[:,1],marker=markers(),label='DLS')
ax.set_xlabel('pH')
ax.set_xlim((4,9))
ax.set_ylabel('${f_H}$')


ax.plot(pHrange, fh3(pHrange, params_b[0],params_b[1]))#,label='If:  M:'+str(round(params_b[1],2))+' E:'+str(round(params_b[0],2)))
ax2.legend()

ax.plot(pHrange, fh3(pHrange, params_c[0],params_c[1]))#,label='Rh:M:'+str(round(params_c[1],2))+' E:'+str(round(params_c[0],2)),)
ax.set_ylim((-0.05,1.05))
ax2.legend()

def theta(pH,Gh2,pka2,M,r):
    return (mp.power(10,(pH-pka2))/(1+mp.power(10,(pH-pka2))))*(1-fh(pH,Gh2,pka2,M,r))



theta2=np.vectorize(theta)

plt.figure(4)
plt.plot(pHrange, theta2(pHrange,params_b[0],4.5,params_b[1],2),label='If',color='orange')
ax2.legend()

plt.plot(pHrange, theta2(pHrange,params_c[0],4.5,params_c[1],2),label='Rh',color='green')
ax2.legend()

#data from https://doi.org/10.1021/ja00115a039

ta=np.array([[0.0335, 4.791   ],
[0.055, 5.1122   ],
[0.0724, 5.279   ],
[0.0898, 5.4458  ],
[0.1093, 5.5868  ],
[0.1338, 5.7532  ],
[0.1573, 5.8554  ],
[0.1849, 5.9445  ],
[0.2114, 6.0594  ],
[0.242, 6.1356   ],
[0.2706, 6.2504  ],
[0.29, 6.3012    ],
[0.3135, 6.4162  ],
[0.335, 6.5312   ],
[0.3595, 6.6979  ],
[0.385, 6.8128   ],
[0.4086, 6.9923  ],
[0.4341, 7.1587  ],
[0.4536, 7.3512  ],
[0.4705, 7.4601  ],
[0.4858, 7.5948  ],
[0.499, 7.7359   ],
[0.5165, 7.8705  ],
[0.5304, 8.0568  ],
[0.5472, 8.2301  ],
[0.5606, 8.4099  ],
[0.5764, 8.5381  ],
[0.5913, 8.7308  ],
[0.6082, 8.9104  ],
[0.6215, 9.0774  ],
[0.6374, 9.2185  ],
[0.6512, 9.3982  ],
[0.6644, 9.4944  ],
[0.6835, 9.687   ],
[0.7024, 9.8537  ],
[0.7264, 10.0202 ],
[0.7473, 10.1547 ],
[0.773, 10.3341  ],
[0.797, 10.4491  ],
[0.8159, 10.5899 ],
[0.8322, 10.7118 ],
[0.8496, 10.8013 ],
[0.8639, 10.8652 ],
[0.8833, 10.9353 ],
[0.8997, 11.0119 ]])


ax2=ax.twinx()
ax2.scatter(ta[:,1],ta[:,0],label="Titration",marker=markers())
ax2.set_ylabel(r'$\theta$')
fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)

fig,ax=plt.subplots(num=5)

mcount=0
ax.scatter(b_scaled[:,0],1-b_scaled[:,1],marker=markers(),label='Fluorescence')
ax.scatter(c_scaled2[:,0],1-c_scaled2[:,1],marker=markers(),label='DLS')
ax.set_xlabel('$pH$')
ax.set_xlim((4,9))
ax.set_ylabel('$1-f_H$')
ax.set_ylim((-0.05,1.05))


ax.plot(pHrange, 1-fh3(pHrange, params_b[0],params_b[1]))#,label='If:  M:'+str(round(params_b[1],2))+' E:'+str(round(params_b[0],2)))


ax.plot(pHrange, 1-fh3(pHrange, params_c[0],params_c[1]))#,label='Rh:M:'+str(round(params_c[1],2))+' E:'+str(round(params_c[0],2)),)


ax2=ax.twinx()
ax2.scatter(ta[:,1],ta[:,0],label="Titration",marker=markers())
ax2.set_ylabel(r'$\theta$')
ax2.set_ylim((-0.05,1.05))
fig.legend(loc="upper right", bbox_to_anchor=(1.01,0.4), bbox_transform=ax.transAxes,frameon=False, title='PEAA')
#plt.savefig('PEAA.pdf')

im = plt.imread(r'Figure_3_3_b_structure.png') # insert local path of the image.

newax = fig.add_axes([0.18,0.60,0.3,0.3], anchor='SW', zorder=1)
newax.imshow(im)
newax.axis('off')
plt.show()
#plt.savefig('PEAA.pdf')

#%%
plt.rcParams.update({'font.size':7})
plt.rcParams.update({'figure.autolayout': 1})#forcesa tight lauout
plt.rcParams['legend.fontsize'] =6
plt.rcParams["figure.figsize"] = (6/2.54,5.5/2.54)
plt.rcParams['legend.title_fontsize'] = 'small'
plt.rcParams['lines.markersize'] = 2
plt.rcParams['lines.linewidth'] = 1

fig9,(ax9) = plt.subplots(1,1,num=9)

mcount=0
ax9.scatter(b_scaled[:,0],1-b_scaled[:,1],marker=markers(),label='Fluorescence')
ax9.scatter(c_scaled2[:,0],1-c_scaled2[:,1],marker=markers(),label='DLS')
ax9.set_xlabel('$pH$')
ax9.set_xlim((4,9))
#ax9.set_ylabel('$1-f_H$')
ax9.set_ylim((-0.05,1.05))


ax9.plot(pHrange, 1-fh3(pHrange, params_b[0],params_b[1]))#,label='If:  M:'+str(round(params_b[1],2))+' E:'+str(round(params_b[0],2)))


ax9.plot(pHrange, 1-fh3(pHrange, params_c[0],params_c[1]))#,label='Rh:M:'+str(round(params_c[1],2))+' E:'+str(round(params_c[0],2)),)


#ax5=ax9.twinx()
ax9.scatter(ta[:,1],ta[:,0],label="Titration",marker=markers())
#ax9.set_ylabel(r'$\theta$')
ax9.set_ylim((-0.05,1.05))
h1, l1 = ax9.get_legend_handles_labels()
#h2, l2 = ax5.get_legend_handles_labels()
#h=h1+h2
#l=l1+l2
fig9.legend(loc="upper right", bbox_to_anchor=(1.01,0.4), bbox_transform=ax9.transAxes,frameon=False, title='PEAA',handles=h1,labels=l1)


im = plt.imread(r'PEAA.png') # insert local path of the image.

newax = fig9.add_axes([0.63,0.60,0.8,0.35], anchor='SW', zorder=1)
newax.imshow(im)
newax.axis('off')
plt.show()

# #%%for combined plotter
# plt.rcParams.update({'font.size':7})
# plt.rcParams.update({'figure.autolayout': 1})#forcesa tight lauout
# plt.rcParams['legend.fontsize'] =6
# plt.rcParams["figure.figsize"] = (12/2.54,5.5/2.54)
# plt.rcParams['legend.title_fontsize'] = 'small'
# plt.rcParams['lines.markersize'] = 2
# plt.rcParams['lines.linewidth'] = 1

# mcount=0
# ax9.scatter(b_scaled[:,0],1-b_scaled[:,1],marker=markers(),label='Fluorescence')
# ax9.scatter(c_scaled2[:,0],1-c_scaled2[:,1],marker=markers(),label='DLS')
# ax9.set_xlabel('$pH$')
# ax9.set_xlim((4,9))
# #ax9.set_ylabel('$1-f_H$')
# ax9.set_ylim((-0.05,1.05))


# ax9.plot(pHrange, 1-fh3(pHrange, params_b[0],params_b[1]))#,label='If:  M:'+str(round(params_b[1],2))+' E:'+str(round(params_b[0],2)))


# ax9.plot(pHrange, 1-fh3(pHrange, params_c[0],params_c[1]))#,label='Rh:M:'+str(round(params_c[1],2))+' E:'+str(round(params_c[0],2)),)


# #ax5=ax9.twinx()
# ax9.scatter(ta[:,1],ta[:,0],label="Titration",marker=markers())
# #ax9.set_ylabel(r'$\theta$')
# ax9.set_ylim((-0.05,1.05))
# h1, l1 = ax9.get_legend_handles_labels()
# #h2, l2 = ax5.get_legend_handles_labels()
# #h=h1+h2
# #l=l1+l2
# fig9.legend(loc="upper right", bbox_to_anchor=(1.01,0.4), bbox_transform=ax9.transAxes,frameon=False, title='PEAA',handles=h1,labels=l1)


# im = plt.imread(r'PEAA.png') # insert local path of the image.

# newax = fig9.add_axes([0.63,0.60,0.8,0.35], anchor='SW', zorder=1)
# newax.imshow(im)
# newax.axis('off')
# plt.show()
# #%%save image
# plt.savefig('combined_coil_globule.png', dpi=3000)
