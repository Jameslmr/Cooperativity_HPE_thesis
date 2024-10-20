# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 10:57:06 2022

@author: james
"""

#note this is all for the gao data, so for a basic system not an acidic one
#script to try and understand how theta works. I think the3 absolute values of the pka's do matter in this case, so it's really quite hard to fit stuff 

import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from cycler import cycler
import matplotlib.ticker as ticker
from scipy import optimize
from matplotlib.widgets import Slider, Button

plt.rcParams.update({'font.size':7})
plt.rcParams.update({'figure.autolayout': 0})#forcesa tight lauout
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


# also no disk state
#total/ionizable groups. so 2 is for a 1 to 1 monomer. set so the hydrophobic parameter is a group specific parameter and not dependent on composition. Easier to see if physical values come out.


#hydrophobic partition function

#f in this case (assuming we set pka3=pka means that f is the fraction of groups that can still ionize in the collapsed conformation)
def Ph(pH,pka1,M,f,r):
    return mp.power((1+(mp.power(10,(pH-pka1)))),(f*M))




#aqueous partition function

def Pa(pH,Gh,pka2,pka3,M,f2,r):
    return mp.power((mp.exp(-Gh*(r-1))),M)*mp.power((1+(mp.power(10,(pH-pka2)))),(f2*M))*mp.power(((1+(mp.power(10,(pH-pka3))))),((1-f2)*M))





#total partition function

def Pt(pH,pka1,Gh,pka2,pka3,M,f,f2,r):
    return Ph(pH,pka1,M,f,r)+Pa(pH,Gh,pka2,pka3,M,f2,r)


#defining the fraction of polymer in each state

def fh2(pH,pka1,Gh,pka2,pka3,M,f,f2,r):
    return np.array(Ph(pH,pka1,M,f,r)/Pt(pH,pka1,Gh,pka2,pka3,M,f,f2,r),dtype=float)

fh=np.vectorize(fh2)
    
def fa2(pH,pka1,Gh,pka2,pka3,M,f,f2,r):
    return np.array(Pa(pH,Gh,pka2,pka3,M,f2,r)/Pt(pH,pka1,Gh,pka2,pka3,M,f,f2,r),dtype=float)

fa=np.vectorize(fa2)

def theta2(pH,pka1,Gh,pka2,pka3,M,f,f2,r):
     return np.array(((mp.power(10,(pH))/(M*Pt(pH,pka1,Gh,pka2,pka3,M,f,f2,r)))*(((f*M)*mp.power(10,(-pka1))*mp.power((1+mp.power(10,(pH-pka1))),((f*M)-1)))+(mp.power(mp.exp(-Gh*(r-1)),M)*((mp.power((1+mp.power(10,(pH-pka3))),((1-f2)*M))*((f2*M)*mp.power(10,(-pka2))*mp.power((1+mp.power(10,(pH-pka2))),((f2*M)-1))))+(mp.power((1+mp.power(10,(pH-pka2))),((f2)*M))*(((1-f2)*M)*mp.power(10,(-pka3))*mp.power((1+mp.power(10,(pH-pka3))),(((1-f2)*M)-1)))))))),dtype=float)

thet=np.vectorize(theta2)

#we know need a graph where we can vary the quantities for pka's and Gh's and see the effect on the fractions. Slider would be ideal. Copied from the slider scripts





#defining our pH dimension
pH_range = np.linspace(2, 10, 401)

# Define initial parameters

init_Gh = 1.8
init_pka1 = 20
init_pka2 = 4.56
init_pka3 = 4.5
init_M = 20
init_f=0
init_f2=1
init_r=2





#naive fitting of what we have here

def Ph_n():
    return 1




#hydrophobic partition function

def Pa_n(pH,Gh2,pka2,M,r):
    return mp.power((mp.exp(-Gh2*(r-1))*(1+(mp.power(10,(pH-pka2))))),M)


def Pt_n(pH,Gh2,pka2,M,r):
    return Ph_n()+Pa_n(pH,Gh2,pka2,M,r)


#defining the fraction of polymer in each state

def fh_n(pH,Gh2,pka2,M,r):
    return Ph_n()/Pt_n(pH,Gh2,pka2,M,r)
    
def fa_n(pH,Gh2,pka2,M,r):
    return Pa_n(pH,Gh2,pka2,M,r)/Pt_n(pH,Gh2,pka2,M,r)


def fh_fit2(pH,Gh2,M):
    pka2=1.5
    r=2
    return np.array(fh_n(pH,Gh2,pka2,M,r),dtype=float)


fh_fit=np.vectorize(fh_fit2)

bounds=((0,0),(30,100))

p0=[5,40]

def fit_line(pH_range,pka1,Gh,pka2,pka3,M,f,f2,r):
    params_a, params_covariance_a = optimize.curve_fit(fh_fit, pH_range,fh(pH_range,pka1,Gh,pka2,pka3,M,f,f2,r),bounds=bounds, p0=p0)
    return params_a[0],params_a[1]








# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = plt.plot(pH_range, fh(pH_range, init_pka1, init_Gh,init_pka2,init_pka3,init_M,init_f,init_f2,init_r), lw=2,label='$f_{h}$')
line2, = plt.plot(pH_range, fa(pH_range, init_pka1, init_Gh,init_pka2,init_pka3,init_M,init_f,init_f2,init_r), lw=2, label='$f_{aq}$')
line3, = plt.plot(pH_range, thet(pH_range, init_pka1, init_Gh,init_pka2,init_pka3,init_M,init_f,init_f2,init_r), lw=2, label=r'$\theta$')

params_init=fit_line(pH_range, init_pka1, init_Gh,init_pka2,init_pka3,init_M,init_f,init_f2,init_r)

line4, = plt.plot(pH_range,fh_fit(pH_range,params_init[0],params_init[1]), lw=2, label="$\overline{f}_{h,fit}, M=$"+str(np.round(params_init[1],decimals=1)))

legend=ax.legend()



axcolor = 'lightgoldenrodyellow'
ax.margins(x=0) 

# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.3, bottom=0.3)

# Make a horizontal slider to control M.
axM = plt.axes([0.25, 0.17, 0.65, 0.03], facecolor=axcolor)
M_slider = Slider(
    ax=axM,
    label='M',
    valmin=1,
    valmax=100,
    valinit=init_M,
)

    
# Make a horizontal slider to control f.
axf = plt.axes([0.25, 0.13, 0.65, 0.03], facecolor=axcolor)
f_slider = Slider(
    ax=axf,
    label='f',
    valmin=0,
    valmax=1,
    valinit=init_f,
)

# Make a horizontal slider to control f.
axf2 = plt.axes([0.25, 0.09, 0.65, 0.03], facecolor=axcolor)
f2_slider = Slider(
    ax=axf2,
    label='f2',
    valmin=0,
    valmax=1,
    valinit=init_f2,
)

# Make a horizontal slider to control r.
axr = plt.axes([0.25, 0.04, 0.65, 0.03], facecolor=axcolor)
r_slider = Slider(
    ax=axr,
    label='r',
    valmin=1,
    valmax=5,
    valinit=init_r,
)


# Make a vertically oriented slider to control pka1
axpka1 = plt.axes([0.06, 0.25, 0.0225, 0.63], facecolor=axcolor)
pka1_slider = Slider(
    ax=axpka1,
    label="pKa1",
    valmin=0,
    valmax=20,
    valinit=init_pka1,
    orientation="vertical"
)

# Make a vertically oriented slider to control pka2
axpka2 = plt.axes([0.14, 0.25, 0.0225, 0.63], facecolor=axcolor)
pka2_slider = Slider(
    ax=axpka2,
    label="pKa2",
    valmin=1,
    valmax=20,
    valinit=init_pka2,
    orientation="vertical"
)

# Make a vertically oriented slider to control pka3
axpka3 = plt.axes([0.18, 0.25, 0.0225, 0.63], facecolor=axcolor)
pka3_slider = Slider(
    ax=axpka3,
    label="pKa3",
    valmin=1,
    valmax=20,
    valinit=init_pka3,
    orientation="vertical"
)





# Make a vertically oriented slider to control Gh2
axGh = plt.axes([0.1, 0.25, 0.0225, 0.63], facecolor=axcolor)
Gh_slider = Slider(
    ax=axGh,
    label="Gh",
    valmin=0,
    valmax=20,
    valinit=init_Gh,
    orientation="vertical"
)

















# The function to be called anytime a slider's value changes
def update(val):
    line.set_ydata(fh(pH_range, pka1_slider.val, Gh_slider.val,pka2_slider.val,pka3_slider.val,M_slider.val,f_slider.val,f2_slider.val,r_slider.val))
    line2.set_ydata(fa(pH_range, pka1_slider.val, Gh_slider.val,pka2_slider.val,pka3_slider.val,M_slider.val,f_slider.val,f2_slider.val,r_slider.val))
    line3.set_ydata(thet(pH_range, pka1_slider.val, Gh_slider.val,pka2_slider.val,pka3_slider.val,M_slider.val,f_slider.val,f2_slider.val,r_slider.val))
    
    params=fit_line(pH_range, pka1_slider.val, Gh_slider.val,pka2_slider.val,pka3_slider.val,M_slider.val,f_slider.val,f2_slider.val,r_slider.val)
    
    line4.set_ydata(fh_fit(pH_range,params[0],params[1]))
    line4.set_label('$f_{h,fit}, M=$'+str(np.round(params[1],decimals=1)))
    legend.remove()
    ax.legend()
    fig.canvas.draw_idle()


# register the update function with each slider

pka1_slider.on_changed(update)
Gh_slider.on_changed(update)
pka2_slider.on_changed(update)
pka3_slider.on_changed(update)
M_slider.on_changed(update)
f_slider.on_changed(update)
f2_slider.on_changed(update)
r_slider.on_changed(update)



# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.0, 0.0, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    Gh_slider.reset()
    pka1_slider.reset()
    pka2_slider.reset()
    pka3_slider.reset()
    M_slider.reset()
    f_slider.reset()
    f2_slider.reset()
    r_slider.reset()
    
button.on_clicked(reset)

plt.show()
    
#%%

#plotting everything in a separate graph so we can export it separately
plt.rcParams.update({'figure.autolayout': 1})#forcesa tight lauout

fig,(ax1,ax2)=plt.subplots(1,2,num=0)
ax1.text(-0.25, 0.95, '(a)', transform=ax1.transAxes, size=8)
ax2.text(-0.25, 0.95, '(b)', transform=ax2.transAxes, size=8)
fig.get_layout_engine().set(w_pad=0.5, h_pad=0.2, pad=0.2)

# Define parameters

Gh = 3
pka1 = 20
pka2 = 4.56
pka3 = 12
M = 40
f=0
f2=0.5
r=2



ax1.plot(pH_range, fh(pH_range, pka1, Gh,pka2,pka3,M,f,f2,r),label='$f_{H}$')
ax1.plot(pH_range, fa(pH_range, pka1, Gh,pka2,pka3,M,f,f2,r), label='$f_{aq}$')
ax1.plot(pH_range, thet(pH_range, pka1, Gh,pka2,pka3,M,f,f2,r), label=r'$\theta$')

params=fit_line(pH_range, pka1, Gh,pka2,pka3,M,f,f2,r)

ax1.plot(pH_range,fh_fit(pH_range,params[0],params[1]), label=r"$M_{eff}=$"+str(np.round(params[1],decimals=1)))

ax1.set_xlim(4,10)
ax1.set_xlabel('$pH$')
ax1.set_ylabel(r'$f$,$\theta$')
ax1.legend(frameon=0)



# Define parameters

Gh = 3
pka1 = 20
pka2 = 4.56
pka3 = 7.5
M = 40
f=0
f2=0.5
r=2



ax2.plot(pH_range, fh(pH_range, pka1, Gh,pka2,pka3,M,f,f2,r), label='$f_{H}$')
ax2.plot(pH_range, fa(pH_range, pka1, Gh,pka2,pka3,M,f,f2,r),  label='$f_{aq}$')
ax2.plot(pH_range, thet(pH_range, pka1, Gh,pka2,pka3,M,f,f2,r),  label=r'$\theta$')

params=fit_line(pH_range, pka1, Gh,pka2,pka3,M,f,f2,r)

ax2.plot(pH_range,fh_fit(pH_range,params[0],params[1]), label=r"$M_{eff}=$"+str(np.round(params[1],decimals=1)))

ax2.set_xlim(4,10)
ax2.set_xlabel('$pH$')
#ax2.set_ylabel(r'$f$,$\theta$')
ax2.legend(frameon=0)







