# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 23:01:14 2018

@author: anneb
"""


import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.colors as colors
import scipy as sp
import numpy as np
import scipy.interpolate
from numpy.ma import masked_array

num_x=200
num_y=200

deltaQ_min,step_deltaQ=0.001,2.1/num_y
KapaQ_min,step_KapaQ=0.0001,0.6/num_x
deltaQ_max=deltaQ_min+num_y*step_deltaQ
KapaQ_max=KapaQ_min+num_x*step_KapaQ

range_deltaQ=np.arange(deltaQ_min,deltaQ_max,step_deltaQ)
range_KapaQ=np.arange(KapaQ_min,KapaQ_max,step_KapaQ)

beta=0.0481589
Aex=20e-12 

epsilon1=np.zeros((num_y,num_x))
epsilon2=np.zeros((num_y,num_x))
Energy=np.zeros((num_y,num_x))
Energycor=np.zeros((num_y,num_x))
EnergyL=np.zeros((num_y,num_x))
r0=np.zeros((num_y,num_x))
r0L=np.zeros((num_y,num_x))
angle=np.zeros((num_y,num_x))
angleL=np.zeros((num_y,num_x))
limite=np.zeros((num_y,num_x))
vecteur1=np.zeros((num_y*num_x,2))
vecteur2=np.zeros((num_y*num_x,2))
number=np.zeros((num_y,num_x))

compteur=0  
for i in np.arange(0,num_y,1): 
    deltaQ=deltaQ_min+step_deltaQ*i
    for j in np.arange(0,num_x,1):
        KapaQ = KapaQ_min+step_KapaQ*j
#        if (((2*KapaQ+deltaQ)/np.sqrt(2))<1) :
        if ((1.081*KapaQ/deltaQ)<1) :
            epsilon2[i][j]=(13.581*KapaQ**2/deltaQ+3.876*deltaQ)
            number1=-beta*epsilon2[i][j]
            if (number1>(-1/np.exp(1))):
                r0[i][j]=1/(16*np.pi)*epsilon2[i][j]/(abs(np.log(beta*epsilon2[i][j]/(abs(np.log(beta*epsilon2[i][j]))))))
                r0L[i][j]=1/(16*np.pi)*epsilon2[i][j]/(abs(scipy.special.lambertw(number1,-1)))
                Energy[i][j]=1/(32*np.pi)*epsilon2[i][j]**2/abs(np.log(beta*epsilon2[i][j]))
                EnergyL[i][j]=(1/(32*np.pi)*epsilon2[i][j]**2/(scipy.special.lambertw(number1,-1))**2)*(-1*scipy.special.lambertw(number1,-1)-1/2)
                angle[i][j]=np.arccos(1.081*KapaQ/deltaQ)*180/np.pi 

        else :
            epsilon1[i][j]=(25.133*KapaQ-7.752*deltaQ)
            number2=-beta*epsilon1[i][j]
            if (number2>(-1/np.exp(1))):
                r0[i][j]=1/(16*np.pi)*epsilon1[i][j]/(abs(np.log(beta*epsilon1[i][j]/(abs(np.log(beta*epsilon1[i][j]))))))
                r0L[i][j]=1/(16*np.pi)*epsilon1[i][j]/(abs(scipy.special.lambertw(number2,-1)))
                Energy[i][j]=1/(32*np.pi)*epsilon1[i][j]**2/abs(np.log(beta*epsilon1[i][j]))
                EnergyL[i][j]=(1/(32*np.pi)*epsilon1[i][j]**2/(scipy.special.lambertw(number2,-1))**2)*(-1*scipy.special.lambertw(number2,-1)-1/2)
                angle[i][j]=0.001

    

line=np.sqrt(2)-2*range_KapaQ
    #%%
matplotlib.rc('xtick', labelsize=50) 
matplotlib.rc('ytick', labelsize=50) 


datar= masked_array(r0,r0==0)
datae= masked_array(Energy,Energy==0)
dataa= masked_array(angle,angle==0)

#



fig,ax = plt.subplots(figsize=(22, 18))
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)
plt.plot(range_KapaQ,line,"k--",linewidth=5)
plt.ylim((0,2.1))
#plt.xlim((0,0.75))
#plt.tight_layout()
CS=ax.contourf(range_KapaQ,range_deltaQ,datar,20,cmap=plt.cm.RdBu)
CS1=ax.contour(range_KapaQ,range_deltaQ,datar,20,colors=('k',),linewidths=0.75)
plt.setp(ax.yaxis.get_majorticklabels(), va="baseline" )
cbar = plt.colorbar(CS) 
cbar.add_lines(CS1)
#
#

fig,ax = plt.subplots(figsize=(22, 18))
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)
plt.plot(range_KapaQ,line,"k--",linewidth=5)
plt.ylim((0,2.1))
#plt.xlim((0,0.75))
CS=ax.contourf(range_KapaQ,range_deltaQ,datae,20,cmap=plt.cm.RdBu)
CS1=ax.contour(range_KapaQ,range_deltaQ,datae,20,colors=('k',),linewidths=0.75)
plt.setp(ax.yaxis.get_majorticklabels(), va="baseline" )
cbar = plt.colorbar(CS) 
cbar.add_lines(CS1)



fig,ax = plt.subplots(figsize=(22, 18))
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)
plt.plot(range_KapaQ,line,"k--",linewidth=5)
plt.ylim((0,2.1))
#plt.xlim((0,0.75))
#plt.setp( ax.xaxis.get_majorticklabels(), rotation=-45 ) 
CS=ax.contourf(range_KapaQ,range_deltaQ,dataa,20,cmap=plt.cm.RdBu)
CS1=ax.contour(range_KapaQ,range_deltaQ,dataa,20,colors=('k',),linewidths=0.75)
plt.setp(ax.yaxis.get_majorticklabels(), va="baseline" )
cbar = plt.colorbar(CS) 
cbar.add_lines(CS1)
