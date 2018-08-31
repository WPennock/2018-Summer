# -*- coding: utf-8 -*-
# Plotting Yingda

## Load Packages and Set Parameters
```python
from aide_design.play import *
from aide_design import floc_model as floc
from scipy.optimize import curve_fit    
import pdb

# Plotting settings
plt.rcParams['text.latex.preamble']=[r"\usepackage{txfonts}"]
params = {'text.usetex' : True,
          'font.size' : 14,
          'font.family' : 'serif',
          'text.latex.unicode': True,
          'axes.facecolor': 'white',
          'axes.labelcolor': 'black',
          'savefig.facecolor': 'white',
          'axes.edgecolor': 'black',
          'savefig.edgecolor': 'black',
          'xtick.color': 'black',
          'ytick.color': 'black'
          }
plt.rcParams.update(params)

u.define('NTU = 100/68*mg/L')

#k = 0.23 # had been 0.18
coag = np.array([0.53, 1.06, 1.59, 2.11, 2.56]) * u.mg/u.L
conc_humic_acid = np.array([0, 3, 6, 9, 12, 15]) * u.mg/u.L
floc.HumicAcid.Density = 1520
#dataset[0] is the 50NTU, dataset[1] is the 100NTU.
#within both subgroups, [0] is the pc.0, ranging evenly up to [5] which is the
# pc.15
dataset = np.array([[[0.634, 0.729, 0.891, 1.062, 1.205],
                     [0.563, 0.717, 0.903, 1.038, 1.193],
                     [0.136, 0.513, 0.793, 1.027, 1.095],
                     [0.109, 0.264, 0.749, 1.002, 1.089],
                     [0.084, 0.128, 0.647, 0.962, 1.057],
                     [0.061, 0.094, 0.308, 0.717, 0.928]
                     ],
                    [[0.746, 0.953, 1.191, 1.295, 1.414],
                     [0.563, 0.835, 1.085, 1.255, 1.403],
                     [0.185, 0.692, 0.971, 1.254, 1.390],
                     [0.105, 0.280, 0.956, 1.238, 1.361],
                     [0.097, 0.207, 0.740, 1.209, 1.316],
                     [0.084, 0.157, 0.566, 1.084, 1.314]
                     ]
                    ])
dataset[0][0]/dataset[1][0]
dataset_conc = 50*10**(-dataset[0])
dataset_conc
dataset2_conc = np.array([[12.44, 8.13, 5.86, 4.25, 3.24], # yingda's dupLicates
                          [14.65, 8.75, 5.97, 4.48, 3.25],
                          [37.41, 14.65, 7.69, 4.63, 3.93],
                          [39.63, 26.98, 8.51, 4.91, 4.04],
                          [40.83, 37.50, 11.32, 5.58, 4.29],
                          [42.95, 40.18, 25.18, 7.98, 5.04],
                           ])
dataset2_50 = -np.log10(dataset2_conc/50)
dataset2_50
dataset2 = np.array([[[ 0.60414962, 0.78887946, 0.93107239, 1.07058107, 1.18842499],
                      [ 0.53313238, 0.75696195, 0.92299567, 1.04769199, 1.18708664],
                      [ 0.1259823 , 0.53313238, 0.81304366, 1.03338901, 1.10457745],
                      [ 0.10094593, 0.26792806, 0.76904044, 1.00788851, 1.09258864],
                      [ 0.08799062, 0.12493874, 0.64512358, 0.95233581, 1.06651271],
                      [ 0.06600684, 0.09496007, 0.29791428, 0.79696711, 0.99653947]
                      ],
                     [[ 0.72153209,  0.97619048,  1.12111801,  1.28070175,  1.36749482],
                      [ 0.55383023,  0.85507246,  1.0942029 ,  1.25087719,  1.33126294],
                      [ 0.18012422,  0.70910973,  0.9989648 ,  1.25087719,  1.31055901],
                      [ 0.10766046,  0.28985507,  0.98861284,  1.22105263,  1.27846791],
                      [ 0.07246377,  0.20082816,  0.7494824 ,  1.1877193 ,  1.22049689],
                      [ 0.06935818,  0.14699793,  0.49585921,  1.1122807 ,  1.1873706 ]
                      ]
                      ])


coagGraph = np.arange(1 * 10**-5, 26.1 * 10**-4, 1 * 10**-4) * u.kg/u.m**3 # same as mathcad sheet
coag_graph = np.linspace(5*10**-5,50*10**-4,num=1000)*u.kg/u.m**3
coag_graph = coag_graph.to(u.mg/u.L)
enerDis = 4.514* u.mW/u.kg
# enerDis = 4.833 * u.mW/u.kg # same as mathcad sheet
temperature = 25 * u.degC # same as mathcad sheet
nu = pc.viscosity_kinematic(temperature)
resTime = 302 * u.s # same as mathcad sheet
tubeDiam = 3/8 * u.inch # same as mathcad sheet
velGrad = (np.sqrt(enerDis/nu)).to(u.s**(-1))
velGrad
# material properties also the same as in mathcad sheet


# trou.Leshooting discrepancy with mathcad
# problem was with odd calculation of phi w/in floc_model pc_viscous function.
@u.wraps(None, [u.W/u.kg, u.degK, u.s, u.m,
                u.kg/u.m**3, u.kg/u.m**3, u.kg/u.m**3, None,
                None, None, u.dimensionless, u.dimensionless], False)
def pc_viscous(energyDis, temp, time, diamTube,
               concClay, concAl, concNatOrgMat, natOrgMat,
               coag, material, fittingParam, ratioheightdiameter):
    return ((3/2)
            * np.log10((2/3) * np.pi * fittingParam * time
                       * (6.0/np.pi)**(2.0/3.0)
                       * np.sqrt(energyDis
                                 / (pc.viscosity_kinematic(temp).magnitude)
                                 )
                       * floc.alpha(diamTube, concClay, concAl, concNatOrgMat,
                               natOrgMat, coag, material, ratioheightdiameter)
                       * floc.frac_vol_floc_initial(concAl,concClay,coag,material)**(2/3)
                       + 1
                       )
            )
```
## fitting humic acid diameter and $k$
```python

# fitting k
def N_viscous (concClay,concNatOrgMat,concAl,material,natOrgMat,coag,diamTube,ratioheightdiameter,time,energyDis,temp):
    x = floc.alpha(diamTube, concClay, concAl, concNatOrgMat, natOrgMat, coag, material, ratioheightdiameter)*time*(np.sqrt(energyDis/pc.viscosity_kinematic(temp))).to(1/u.s)*floc.frac_vol_floc_initial(concAl, concClay, coag, material)**(2/3)
    return x.to(u.dimensionless)

def viscous_fit(N,k):
    return 1.5*np.log10(2.0/3.0*np.pi*k*N*(6.0/np.pi)**(2.0/3.0) + 1)
# # verify correctness
# pc_viscous(enerDis,temperature,resTime,tubeDiam,100*u.NTU,coag,0*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,kfit,floc.RATIO_HEIGHT_DIAM)
# n_test = N_viscous(100*u.NTU,coag,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
# viscous_fit(n_test,kfit)

N_fit50 = N_viscous(50*u.NTU,0*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
N_fit100 = N_viscous(100*u.NTU,0*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
N_graph = N_viscous(50*u.NTU,0*u.mg/u.L,coag_graph,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
# fit k for 50 NTU 0 mg/l ha data
kfit, kfitvar= curve_fit(viscous_fit,np.append(N_fit50,N_fit50),np.append(dataset[0][0],dataset2[0][0]))              
kfit

## fit k for combined 0 mg/l ha data
# kfit, kfitvar= curve_fit(viscous_fit,np.concatenate([N_fit50,N_fit100]),np.concatenate([dataset[0][0],dataset[1][0]]))              
kfit
## verify fit
N_graph50  = N_viscous(50*u.NTU,0*u.mg/u.L,coag_graph,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
N_graph100  = N_viscous(100*u.NTU,0*u.mg/u.L,coag_graph,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
N_plot = np.linspace(2.5,16,num=100)

plt.plot(np.append(N_fit50,N_fit50),np.append(dataset[0][0],dataset2[0][0]),'x')
plt.plot(N_graph50,pc_viscous(enerDis,temperature,resTime,tubeDiam,50*u.NTU,coag_graph,0*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,kfit,floc.RATIO_HEIGHT_DIAM),'k')
plt.show()
```
## fitting $d_{HA}$
### Brute Force Method
#### 50 NTU
```python
def SSE(a,b):
    '''function for determining sum-squared-error'''
    return np.sum((a-b)**2)

pC50NTU = np.concatenate((dataset[0],dataset2[0]))
pC100NTU = np.concatenate((dataset[1],dataset2[1]))
HArep = np.append(conc_humic_acid,conc_humic_acid)*u.mg/u.L
conv = 0.25 #criterion for pC*
dHA_range = np.arange(1,1001,1)*u.nm
dHA_50i = np.ones(len(HArep))
SSE50 = np.zeros([len(HArep),len(dHA_range)])
for i in range (0,len(HArep)):    
    for j in range(0,len(dHA_range)):
        floc.HumicAcid.Diameter = dHA_range[j].to(u.m).magnitude
        pCpred = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coag, HArep[i], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
        # pdb.set_trace()
        SSE50[i][j] = SSE(pCpred[np.where(pC50NTU[i]>conv)],pC50NTU[i][np.where(pC50NTU[i]>conv)])
    dHA_50i[i] = np.argmin(SSE50[i])    
dHA_50i

# Attempt to fit for all data at once
dHA_rangeA = np.arange(1,200,1)*u.nm
SSE50A = np.zeros(len(dHA_rangeA))
pC50NTUA = np.concatenate([np.concatenate(dataset[0][1:]),np.concatenate(dataset2[0][1:])])
np.shape(pC50NTUA)
coagA = np.concatenate([coag,coag,coag,coag,coag,coag,coag,coag,coag,coag])*u.mg/u.L
np.shape(coagA)
halfHA = np.zeros(25)
for i in range(0,len(halfHA)):
    halfHA[i] = conc_humic_acid[1+np.int(np.floor(i*5/len(halfHA)))].magnitude
HArepA = np.concatenate([halfHA,halfHA])*u.mg/u.L
np.shape(HArepA)
convA = 0.25
pCpredA = np.zeros([len(dHA_rangeA),len(coagA)])
np.shape(pCpredA)
for i in range(0,len(SSE50A)):
    # pdb.set_trace()
    floc.HumicAcid.Diameter = dHA_rangeA[i].to(u.m).magnitude
    for j in range(0,len(coagA)):
        pCpredA[i][j] = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagA[j], HArepA[j], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
    SSE50A[i] = SSE(pCpredA[i][np.where(pC50NTUA>convA)],pC50NTUA[np.where(pC50NTUA>convA)])

np.argmin(SSE50A)
np.shape(SSE50A)
```
#### 100 NTU
```python
# dHA_100i = np.ones(len(conc_humic_acid))
# SSE100 = np.zeros([len(conc_humic_acid),len(dHA_range)])
# for i in range (0,len(conc_humic_acid)):    
#     for j in range(0,len(dHA_range)):
#         floc.HumicAcid.Diameter = dHA_range[j].to(u.m).magnitude
#         pCpred = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coag, conc_humic_acid[i], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
#         # pdb.set_trace()
#         SSE100[i][j] = SSE(pCpred[np.where(dataset[1][i]>conv)],dataset[1][i][np.where(dataset[1][i]>conv)])
#         np.where(dataset[1][i]>conv)[0]
#     dHA_100i[i] = np.argmin(SSE100[i])    
#         # Return array of SSE
# dHA_100i
# dHA_fiti = np.mean([dHA_50i[2:],dHA_100i[2:]])
# dHA_fiti
# dHA_fit = dHA_range[int(np.round(dHA_fiti))]

# Using individual fits
# dHA_fit = dHA_range[int(np.round(np.mean(np.append(dHA_50i[2:6],dHA_50i[8:12]))))]
# floc.HumicAcid.Diameter = dHA_fit.to(u.m).magnitude

# Using collective fit
dHA_fit = dHA_rangeA[np.argmin(SSE50A)]
floc.HumicAcid.Diameter = dHA_fit.to(u.m).magnitude
floc.HumicAcid.Diameter
```
### Elegant method (curve_fit)
```python
# def pc_fit_dHA(dataset,dHA):
#     '''# For dataset, 0th row is coagulant, 1st is humic acid, 2nd is influent turbidity'''
#     G_Coag = floc.gamma_coag(dataset[2],dataset[0],floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM)
#     G_HA = np.minimum((((dataset[1].to(u.kg/u.m**3)).magnitude / (floc.conc_precipitate(dataset[0], floc.PACl).to(u.kg/u.m**3)).magnitude) * (floc.PACl.Density / floc.HumicAcid.Density) * (floc.PACl.Diameter / (4 * dHA)) ), np.ones(len(dataset[0])))
#     pred = np.zeros(len(dataset[0]))
#     for i in range(0,len(dataset[0])):
#         pred[i] = ((3/2) * np.log10((2/3) * np.pi * kfit * resTime.to(u.s).magnitude * (np.sqrt(enerDis.to(u.m**2/u.s**3) / (pc.viscosity_kinematic(temperature)) ).to(1/u.s)).magnitude * (2*(1-G_Coag[i])*(G_Coag[i]*(1-G_HA[i])) + (G_Coag[i]*(1-G_HA[i]))**2 + 2*(G_Coag[i]*(1-G_HA[i]))*(G_HA[i]*G_Coag[i]) ) * (np.pi/6)**(2/3) * (floc.Clay.Diameter / (floc.sep_dist_clay(dataset[2][i], floc.Clay).to(u.m)).magnitude ) ** 2 + 1 ) )
#     return pred    
#
# def pc_fit_dHA_alt(dataset,dHA):
#     '''# For dataset, 0th row is coagulant, 1st is humic acid'''
#     floc.HumicAcid.Diameter = dHA
#     return pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, dataset[0], dataset[1], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
#
#
# # 50 NTU Data            
# dHA_50 = np.zeros(6)
# dHAvar_50 = dHA_50
# for i in range(0,6):            
#     conv = 0.5
#     data = np.ones([3,len(np.where(dataset[0][i]>conv)[0])])
#     data[0] = coag[np.where(dataset[0][i]>conv)]
#     data[1] = data[1]*conc_humic_acid[i]
#     data[2] = (data[2]*50*u.NTU).to(u.mg/u.L)
#     data = data*u.mg/u.L
#     pC = dataset[0][i][np.where(dataset[0][i]>conv)]
#     dHA_50[i],dHAvar_50[i] = curve_fit(pc_fit_dHA,data,pC,bounds=(1e-9,1e-6),p0=1e-6,method='dogbox',maxfev=1000000)
# dHA_50    
# pc_fit_dHA(data,dHA_50[5])
# dataset[0][5]
#
# # Try fitting using already defined function.
# def pc_fit_dHA_alt(dataset,dHA):
#     '''# For dataset, 0th row is coagulant, 1st is humic acid'''
#     floc.HumicAcid.Diameter = dHA
#     return pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, dataset[0], dataset[1], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
#
# for i in range(0,6):
#         conv = 0.5
#         data = np.ones([2,len(np.where(dataset[0][i]>conv)[0])])
#         data[0] = coag[np.where(dataset[0][i]>conv)]
#         data[1] = data[1]*conc_humic_acid[i]
#         data = data*u.mg/u.L
#         pC = dataset[0][i][np.where(dataset[0][i]>conv)]
#         dHA_50[i],dHAvar_50[i] = curve_fit(pc_fit_dHA_alt,data,pC,bounds=(1e-9,1e-6),p0=1e-6,method='dogbox')
# # Plot to compare fitting functions
# plt.clf()
# plt.close('all')
# for i in range(0,6):
#     plt.plot(coag,dataset[0][i],'x')
#     dataplot= np.ones([3,len(coag_graph)])
#     dataplot[0] = coag_graph.to(u.mg/u.L)
#     dataplot[1] = (dataplot[1]*conc_humic_acid[i]).to(u.mg/u.L)
#     dataplot[2] = (dataplot[2]*50*u.NTU).to(u.mg/u.L)
#     dataplot = dataplot*u.mg/u.L
#     plt.plot(coag_graph.to(u.mg/u.L),pc_fit_dHA(dataplot,floc.HumicAcid.Diameter),'--')
#     plt.plot(coag_graph.to(u.mg/u.L),pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coag_graph, conc_humic_acid[i], floc.HumicAcid, floc.PACl, floc.Clay, k, floc.RATIO_HEIGHT_DIAM))    
# plt.show()    
# dHA_50
# # Trying with LMFIT
# from lmfit import Model
# pcmodel = Model(pc_fit_dHA)
# pcmodel.param_names
# pcmodel.independent_vars
# params = pcmodel.make_params(dHA = 40e-9)
# result = pcmodel.fit(pC,data,dHA=40e-9)
# # def pc_fit_k_dHA(data,k,dHA):
# #     '''# For data, 0th row is coagulant 1st is humic acid, 2nd is influent turbidity'''
# #     G_Coag = floc.gamma_coag(data[2]*u.mg/u.L,data[0]*u.mg/u.L,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM)
# #     G_HA = np.minimum(((data[1]*u.mg/u.L / floc.conc_precipitate(data[0]*u.mg/u.L, floc.PACl))
# #                 * (floc.PACl.Density / floc.HumicAcid.Density)
# #                 * (floc.PACl.Diameter*u.m / (4 * dHA*u.nm))
# #                 ).to(u.dimensionless),
# #                np.ones(len(data[0])))
# #     pred = np.zeros(len(data[0]))
# #     for i in range(0,len(data[0])):
# #         pred[i] = ((3/2) * np.log10((2/3) * np.pi * k * resTime * np.sqrt(enerDis.to(u.m**2/u.s**3) / (pc.viscosity_kinematic(temperature)) ) * ( 2*(1-G_Coag[i])*(G_Coag[i]*(1-G_HA[i])) + (G_Coag[i]*(1-G_HA[i]))**2 + 2*(G_Coag[i]*(1-G_HA[i]))*(G_HA[i]*G_Coag[i]) ) * (np.pi/6)**(2/3) * (floc.Clay.Diameter*u.m / floc.sep_dist_clay(data[2][i], floc.Clay) ) ** 2 + 1 ) )
# #     return pred    
# #
# # fit_50 = np.zeros([6,1])
# # fitvar_50 = fit_50
# #
# # for i in range(0,6):            
# #    data = np.ones((3,len(np.where(dataset[0][i]>0.5)[0])))
# #    data[0] = coag[np.where(dataset[0][i]>0.5)]
# #    data[1] = data[1]*conc_humic_acid[i]
# #    data[2] = data[2]*50*u.NTU
# #    pC = dataset[0][i][np.where(dataset[0][i]>0.5)]
# #    fit_50[i],fitvar_50[i] = curve_fit(pc_fit_k_dHA,data,pC)
# # fit, fitvar = curve_fit(pc_fit_k_dHA,data,pC)
# # fit
# # pc_fit_k_dHA(data,0.0691)
#
# # 100 NTU Data            
# dHA_100 = np.zeros(6)
# dHAvar_100 = dHA_100
# for i in range(0,6):            
#    data = np.ones((3,len(np.where(dataset[1][i]>0.5)[0])))
#    data[0] = coag[np.where(dataset[1][i]>0.5)]
#    data[1] = data[1]*conc_humic_acid[i]
#    data[2] = data[2]*100*u.NTU
#    pC = dataset[0][i][np.where(dataset[0][i]>0.5)]
#    dHA_50[i],dHAvar_50[i] = curve_fit(pc_fit_dHA,data,pC,p0=100)
# dHA_100
```
## Evaluate Goodness of Fit (By humic acid concentration)
### 50 NTU
```python
def RMSE(A,B):
    return np.sqrt(np.mean((A-B)**2))

RMSE50 = np.zeros(len(conc_humic_acid))
for i in range (0,len(conc_humic_acid)):    
    RMSE50[i] = RMSE(pC50NTU[i],pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coag, HArep[i], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM))
RMSE50
RMSE50mean = np.mean(RMSE50)
RMSE50mean
```
### 100 NTU
```python
RMSE100 = np.zeros(len(conc_humic_acid))
for i in range (0,len(conc_humic_acid)):    
    RMSE100[i] = RMSE(pC100NTU[i],pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coag, conc_humic_acid[i], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM))
RMSE100
RMSE100mean = np.mean(RMSE100)
RMSE100mean
```
## Evaluate Goodness of Fit (Whole data set)
### 50 NTU
```python
coagAll = np.concatenate([coag,coag,coag,coag,coag,coag,coag,coag,coag,coag,coag,coag])*u.mg/u.L
np.shape(coagAll)
halfHAAll = np.zeros(30)
for i in range(0,len(halfHAAll)):
    halfHAAll[i] = conc_humic_acid[np.int(np.floor(i*6/len(halfHAAll)))].magnitude
HArepAll = np.concatenate([halfHAAll,halfHAAll])*u.mg/u.L
np.shape(HArepAll)
pC50NTUAll = np.concatenate([np.concatenate(dataset[0]),np.concatenate(dataset2[0])])
np.shape(pC50NTUAll)
pC50predAll = np.zeros(len(coagAll))
for i in range(0,len(coagAll)):
    pC50predAll[i] = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagAll[i], HArepAll[i], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
RMSE50All = RMSE(pC50NTUAll,pC50predAll)
RMSE50All

```
### 100 NTU
```python
pC100NTUAll = np.concatenate([np.concatenate(dataset[1]),np.concatenate(dataset2[1])])
pC100predAll = np.zeros(len(coagAll))
for i in range (0,len(coagAll)):    
    pC100predAll[i] = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagAll[i], HArepAll[i], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
RMSE100All = RMSE(pC100NTUAll,pC100predAll)
RMSE100All
```

#Begin graphing the 50NTU datasets
```python
plt.clf()
plt.close('all')

plt.subplot(121)
plt.title('50 NTU Graph')
plt.ylabel('pC*')
plt.xlabel('coagulant dosage (mg/L)')

plt.plot(coag, dataset[0][0], 'r.', coag, dataset[0][1], 'b.', coag, dataset[0][2], 'g.',
         coag, dataset[0][3], 'm.', coag, dataset[0][4], 'c.', coag, dataset[0][5], 'y.')

#I wish there was a cleaner way to assign these but I can't think
# of what it would be.
# Until floc_model is fixed, use locally defined pc_viscous
line0mg50 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagGraph, 0 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line3mg50 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagGraph, 3 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line6mg50 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagGraph, 6 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line9mg50 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagGraph, 9 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line12mg50 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagGraph, 12 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line15mg50 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagGraph, 15 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)

x = coagGraph.to(u.mg/u.L)
plt.plot(x, line0mg50, 'r', x, line3mg50, 'b', x, line6mg50, 'g',
         x, line9mg50, 'm', x, line12mg50, 'c', x, line15mg50, 'y')

#Begin graphing the 100NTU datasets
plt.subplot(122)
plt.title('100 NTU Graph')
plt.ylabel('pC*')
plt.xlabel('coagulant dosage (mg/L)')

plt.plot(coag, dataset[1][0], 'r.', coag, dataset[1][1], 'b.', coag, dataset[1][2], 'g.',
         coag, dataset[1][3], 'm.', coag, dataset[1][4], 'c.', coag, dataset[1][5], 'y.')

line0mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 0 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line3mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 3 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line6mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 6 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line9mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 9 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line12mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 12 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
line15mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 15 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)

x = coagGraph.to(u.mg/u.L)
plt.plot(x, line0mg100, 'r', x, line3mg100, 'b', x, line6mg100, 'g',
         x, line9mg100, 'm', x, line12mg100, 'c', x, line15mg100, 'y')

#And now we display our graph!

#nua = 15*u.mm**2/u.s
#EDR1 = 1.5*u.m**2/u.s**3
#EDR2 = 6*u.m**2/u.s**3
#eta1 = (nua**3/EDR1)**(1/4)
#eta2 = (nua**3/EDR2)**(1/4)
#eta1.to(u.mm)
#eta2.to(u.mm)
plt.savefig('Yingda.png',format='png')
plt.show()
```
# Okay, but let's make some publication quality graphs!


## 50 NTU Data
```python
# Create function for old model (no HA)
@u.wraps(None, [u.W/u.kg, u.degK, u.s, u.m,
               u.kg/u.m**3, u.kg/u.m**3,
               None, None, u.dimensionless, u.dimensionless], False)
def pc_viscous_old(EnergyDis, Temp, Time, DiamTube,
              ConcClay, ConcAl, coag, material,
              FittingParam, RatioHeightDiameter):
   return ((3/2)
           * np.log10((2/3) * np.pi * FittingParam * Time
                      * np.sqrt(EnergyDis
                                / (pc.viscosity_kinematic(Temp).magnitude)
                                )
                      * (2*floc.gamma_coag(ConcClay,ConcAl,coag,material,DiamTube,RatioHeightDiameter)-floc.gamma_coag(ConcClay,ConcAl,coag,material,DiamTube,RatioHeightDiameter)**2)
                      * floc.frac_vol_floc_initial(ConcAl,ConcClay,coag,material)**(2/3)
                      + 1
                      )
           )

def N_old(ConcClay,ConcAl,coag,material,DiamTube,RatioHeightDiameter,Time,EnergyDis,Temp):
	x = (2*floc.gamma_coag(ConcClay, ConcAl, coag, material, DiamTube, RatioHeightDiameter)-floc.gamma_coag(ConcClay, ConcAl, coag, material, DiamTube, RatioHeightDiameter)**2)*Time*(np.sqrt(EnergyDis/pc.viscosity_kinematic(Temp))).to(1/u.s)*floc.frac_vol_floc_initial(ConcAl, ConcClay, coag, material)**(2/3)
	return x.to(u.dimensionless)

N_fit_old = N_old(50*u.NTU,coag,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
N_graph_old = N_old(50*u.NTU,coag_graph,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)

def viscous_fit(N,k):
    return 1.5*np.log10(2.0/3.0*np.pi*k*N*(6.0/np.pi)**(2.0/3.0) + 1)

k_old, kvar_old = curve_fit(viscous_fit,np.append(N_fit_old,N_fit_old),np.append(dataset[0][0],dataset2[0][0]))              

plt.clf()
plt.close('all')
plt.figure(0)
plt.plot(coag,dataset[0][0],'--rx', label=r'0 mg/L HA')
plt.plot(coag,dataset2[0][0],'--rx')
plt.plot(coag,dataset[0][1],'--b+', label=r'3 mg/L HA')
plt.plot(coag,dataset2[0][1],'--b+')
plt.plot(coag,dataset[0][2],'--gs', markerfacecolor='none', label=r'6 mg/L HA')
plt.plot(coag,dataset2[0][2],'--gs', markerfacecolor='none')
plt.plot(coag,dataset[0][3],'--mD', markerfacecolor='none', label=r'9 mg/L HA')
plt.plot(coag,dataset2[0][3],'--mD', markerfacecolor='none')
plt.plot(coag,dataset[0][4],'--co', markerfacecolor='none', label=r'12 mg/L HA')
plt.plot(coag,dataset2[0][4],'--co', markerfacecolor='none')
plt.plot(coag,dataset[0][5],'--^',c='xkcd:brown', markerfacecolor='none', label=r'12 mg/L HA')
plt.plot(coag,dataset2[0][5],'--^',c='xkcd:brown', markerfacecolor='none')
plt.plot(coag_graph.to(u.mg/u.L)[90:500],viscous_fit(N_graph_old,k_old)[90:500],'k',label=r'Pennock et al. (2018)')
#plt.plot(x[4:],pc_viscous_old(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagGraph[4:], floc.PACl, floc.Clay, k50, floc.RATIO_HEIGHT_DIAM),'k',label=r'Pennock et al. (2018)')
plt.xlabel(r'Coagulant Dose (mg/L)')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0.5, 2.6, 0, 1.6])
plt.legend(loc='upper left',ncol=2,borderpad=0.1,handletextpad=0.1,labelspacing=0,columnspacing=0.1,edgecolor='inherit')
plt.tight_layout()
plt.savefig('50NTU.png',format='png')
plt.savefig('50NTU.eps',format='eps')
plt.show()
```

## 50 NTU Fit
```python
plt.clf()
plt.close('all')
plt.figure()
plt.plot(coag,dataset[0][0],'rx', label=r'0 mg/L HA Data')
plt.plot(coag,dataset2[0][0],'rx')
plt.plot(coag,dataset[0][1],'b+', label=r'3 mg/L HA Data')
plt.plot(coag,dataset2[0][1],'b+')
plt.plot(coag,dataset[0][2],'gs', markerfacecolor='none', label=r'6 mg/L HA Data')
plt.plot(coag,dataset2[0][2],'gs', markerfacecolor='none')
plt.plot(coag,dataset[0][3],'mD', markerfacecolor='none', label=r'9 mg/L HA Data')
plt.plot(coag,dataset2[0][3],'mD', markerfacecolor='none')
plt.plot(coag,dataset[0][4],'co', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.plot(coag,dataset2[0][4],'co', markerfacecolor='none')
plt.plot(coag,dataset[0][5],'^',c='xkcd:brown', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.plot(coag,dataset2[0][5],'^',c='xkcd:brown', markerfacecolor='none')
plt.plot(x, line0mg50, 'r', label=r'0 mg/L HA Model')
plt.plot(x, line3mg50, 'b', label=r'3 mg/L HA Model')
plt.plot(x, line6mg50, 'g', label=r'6 mg/L HA Model')
plt.plot(x, line9mg50, 'm', label=r'9 mg/L HA Model')
plt.plot(x, line12mg50, 'c', label=r'12 mg/L HA Model')
plt.plot(x, line15mg50, 'xkcd:brown', label=r'15 mg/L HA Model')
plt.xlabel(r'Coagulant Dose (mg/L)')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0,2.6,0,1.25])
plt.legend(loc=2,bbox_to_anchor=(0,-0.6,1,0.4),ncol=2,borderpad=0.1,handletextpad=0.2,labelspacing=0,columnspacing=0.2,edgecolor='white')
#plt.tight_layout()

plt.savefig('50NTUfit.png',format='png', bbox_inches = "tight")
plt.savefig('50NTUfit.eps',format='eps', bbox_inches = "tight")
plt.show()
```

## 100 NTU Fit
```python
plt.clf()
plt.close('all')
plt.figure()
plt.plot(coag,dataset[1][0],'rx', label=r'0 mg/L HA Data')
plt.plot(coag,dataset2[1][0],'rx')
plt.plot(coag,dataset[1][1],'b+', label=r'3 mg/L HA Data')
plt.plot(coag,dataset2[1][1],'b+')
plt.plot(coag,dataset[1][2],'gs', markerfacecolor='none', label=r'6 mg/L HA Data')
plt.plot(coag,dataset2[1][2],'gs', markerfacecolor='none')
plt.plot(coag,dataset[1][3],'mD', markerfacecolor='none', label=r'9 mg/L HA Data')
plt.plot(coag,dataset2[1][3],'mD', markerfacecolor='none')
plt.plot(coag,dataset[1][4],'co', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.plot(coag,dataset2[1][4],'co', markerfacecolor='none')
plt.plot(coag,dataset[1][5],'^',c='xkcd:brown', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.plot(coag,dataset2[1][5],'^',c='xkcd:brown', markerfacecolor='none')
plt.plot(x, line0mg100, 'r', label=r'0 mg/L HA Model')
plt.plot(x, line3mg100, 'b', label=r'3 mg/L HA Model')
plt.plot(x, line6mg100, 'g', label=r'6 mg/L HA Model')
plt.plot(x, line9mg100, 'm', label=r'9 mg/L HA Model')
plt.plot(x, line12mg100, 'c', label=r'12 mg/L HA Model')
plt.plot(x, line15mg100, 'xkcd:brown', label=r'15 mg/L HA Model')
plt.xlabel(r'Coagulant Dose (mg/L)')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0,2.6,0,1.5])
plt.legend(loc=2,bbox_to_anchor=(0,-0.6,1,0.37),ncol=2,borderpad=0.1,handletextpad=0.2,labelspacing=0,columnspacing=0.2,edgecolor='white')
#plt.tight_layout()

plt.savefig('100NTUfit.png',format='png', bbox_inches = "tight")
plt.savefig('100NTUfit.eps',format='eps', bbox_inches = "tight")
plt.show()
```
## Collision Potential
```python
N50 = np.zeros([len(conc_humic_acid),len(coag)])
for i in range(0,len(N50)):
    N50[i] = N_viscous(50*u.NTU,conc_humic_acid[i],coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)

plt.clf()
plt.close('all')
plt.figure()
plt.semilogx(N50[0],dataset[0][0],'rx', label=r'0 mg/L HA Data')
plt.semilogx(N50[0],dataset2[0][0],'rx')
plt.semilogx(N50[1],dataset[0][1],'b+', label=r'3 mg/L HA Data')
plt.semilogx(N50[1],dataset2[0][1],'b+')
plt.semilogx(N50[2],dataset[0][2],'gs', markerfacecolor='none', label=r'6 mg/L HA Data')
plt.semilogx(N50[2],dataset2[0][2],'gs', markerfacecolor='none')
plt.semilogx(N50[3],dataset[0][3],'mD', markerfacecolor='none', label=r'9 mg/L HA Data')
plt.semilogx(N50[3],dataset2[0][3],'mD', markerfacecolor='none')
plt.semilogx(N50[4],dataset[0][4],'co', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.semilogx(N50[4],dataset2[0][4],'co', markerfacecolor='none')
plt.semilogx(N50[5],dataset[0][5],'^',c='xkcd:brown', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.semilogx(N50[5],dataset2[0][5],'^',c='xkcd:brown', markerfacecolor='none')
plt.semilogx(N_graph,viscous_fit(N_graph,kfit),'xkcd:navy',label=r'Model')
plt.xlabel(r'$\alpha\overline{G}\theta\phi^{2/3}$')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0.2, 20, 0, 1.7])
plt.legend(loc='upper left',ncol=2,borderpad=0.1,handletextpad=0.1,labelspacing=0,columnspacing=0.1,edgecolor='inherit')
plt.tight_layout()
plt.savefig('Collision.png',format='png')
plt.savefig('Collision.eps',format='eps')
plt.show()
```
### Collision Potential with 100 NTU
```python
N100 = np.zeros([len(conc_humic_acid),len(coag)])
for i in range(0,len(N100)):
    N100[i] = N_viscous(100*u.NTU,conc_humic_acid[i],coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)

plt.clf()
plt.close('all')
plt.figure()
# 50 NTU data
plt.semilogx(N50[0],dataset[0][0],'rx', label=r'50 NTU, 0 mg/L HA Data')
plt.semilogx(N50[0],dataset2[0][0],'rx')
plt.semilogx(N50[1],dataset[0][1],'b+', label=r'50 NTU, 3 mg/L HA Data')
plt.semilogx(N50[1],dataset2[0][1],'b+')
plt.semilogx(N50[2],dataset[0][2],'gs', markerfacecolor='none', label=r'50 NTU, 6 mg/L HA Data')
plt.semilogx(N50[2],dataset2[0][2],'gs', markerfacecolor='none')
plt.semilogx(N50[3],dataset[0][3],'mD', markerfacecolor='none', label=r'50 NTU, 9 mg/L HA Data')
plt.semilogx(N50[3],dataset2[0][3],'mD', markerfacecolor='none')
plt.semilogx(N50[4],dataset[0][4],'co', markerfacecolor='none', label=r'50 NTU, 12 mg/L HA Data')
plt.semilogx(N50[4],dataset2[0][4],'co', markerfacecolor='none')
plt.semilogx(N50[5],dataset[0][5],'^',c='xkcd:brown', markerfacecolor='none', label=r'50 NTU, 15 mg/L HA Data')
plt.semilogx(N50[5],dataset2[0][5],'^',c='xkcd:brown', markerfacecolor='none')
# Model
plt.semilogx(N_graph,viscous_fit(N_graph,kfit),'xkcd:navy',label=r'Model')
# 100 NTU data
plt.semilogx(N100[0],dataset[1][0],'kx', label=r'100 NTU, 0 mg/L HA Data' )
plt.semilogx(N100[0],dataset2[1][0],'kx')
plt.semilogx(N100[1],dataset[1][1], 'k+', label=r'100 NTU, 3 mg/L HA Data')
plt.semilogx(N100[1],dataset2[1][1],'k+')
plt.semilogx(N100[2],dataset[1][2],'ks', markerfacecolor='none', label=r'100 NTU, 6 mg/L HA Data')
plt.semilogx(N100[2],dataset2[1][2],'ks', markerfacecolor='none')
plt.semilogx(N100[3],dataset[1][3],'kD', markerfacecolor='none', label=r'100 NTU, 9 mg/L HA Data')
plt.semilogx(N100[3],dataset2[1][3],'kD', markerfacecolor='none')
plt.semilogx(N100[4],dataset[1][4],'ko', markerfacecolor='none', label=r'100 NTU, 12 mg/L HA Data')
plt.semilogx(N100[4],dataset2[1][4],'ko', markerfacecolor='none')
plt.semilogx(N100[5],dataset[1][5],'k^', markerfacecolor='none', label=r'100 NTU, 15 mg/L HA Data')
plt.semilogx(N100[5],dataset2[1][5],'k^', markerfacecolor='none')

plt.xlabel(r'$\overline{\alpha}\overline{G}\theta\phi^{2/3}$')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0.2, 20, 0, 1.7])
plt.legend(loc=2,bbox_to_anchor=(-0.125,-0.6,1,0.37),ncol=2,borderpad=0.1,handletextpad=0.2,labelspacing=0,columnspacing=0.2,edgecolor='white')
# plt.tight_layout()
plt.savefig('Collision100.png',format='png', bbox_inches='tight')
plt.savefig('Collision100.eps',format='eps', bbox_inches='tight')
plt.show()
```
## $\overline{\Gamma}_\mathrm{HA}$
```python
plt.clf()
plt.close('all')
plt.figure()
plt.plot(coag_graph.to(u.mg/u.L),floc.gamma_humic_acid_to_coag(coag_graph,conc_humic_acid[0],floc.HumicAcid,floc.PACl),'r',label=r'0 mg/L HA')
plt.plot(coag_graph.to(u.mg/u.L),floc.gamma_humic_acid_to_coag(coag_graph,conc_humic_acid[1],floc.HumicAcid,floc.PACl),'b',label=r'3 mg/L HA')
plt.plot(coag_graph.to(u.mg/u.L),floc.gamma_humic_acid_to_coag(coag_graph,conc_humic_acid[2],floc.HumicAcid,floc.PACl),'g',label=r'6 mg/L HA')
plt.plot(coag_graph.to(u.mg/u.L),floc.gamma_humic_acid_to_coag(coag_graph,conc_humic_acid[3],floc.HumicAcid,floc.PACl),'m',label=r'9 mg/L HA')
plt.plot(coag_graph.to(u.mg/u.L),floc.gamma_humic_acid_to_coag(coag_graph,conc_humic_acid[4],floc.HumicAcid,floc.PACl),'c',label=r'12 mg/L HA')
plt.plot(coag_graph.to(u.mg/u.L),floc.gamma_humic_acid_to_coag(coag_graph,conc_humic_acid[5],floc.HumicAcid,floc.PACl),'xkcd:brown',label=r'15 mg/L HA')
plt.xlabel(r'Coagulant Dose (mg/L)')
plt.ylabel(r'$\overline{\Gamma}_\mathrm{HA-PACl}$')
plt.axis([0, 3, 0, 1])
plt.legend(loc='upper right',borderpad=0.1,handletextpad=0.1,labelspacing=0,columnspacing=0.1,edgecolor='inherit')
plt.tight_layout()
plt.savefig('Gamma_HA.png',format='png')
plt.savefig('Gamma_HA.eps',format='eps')
plt.show()
```
## $\alpha$ vs. Dose
```python
plt.clf()
plt.close('all')
plt.figure()
# Red is 0 mg/L
plt.plot(coag_graph,floc.alpha_pacl_clay(tubeDiam,50*u.NTU,coag_graph,0*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'r:',label=r'$\alpha_\mathrm{PACl-Clay}$ 0 mg/L HA')
plt.plot(coag_graph,floc.alpha_pacl_pacl(tubeDiam,50*u.NTU,coag_graph,0*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'r-.',label=r'$\alpha_\mathrm{PACl-PACl}$ 0 mg/L HA')
plt.plot(coag_graph,floc.alpha_pacl_nat_org_mat(tubeDiam,50*u.NTU,coag_graph,0*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'r--',label=r'$\alpha_\mathrm{HA-PACl}$ 0 mg/L HA')
plt.plot(coag_graph,floc.alpha(tubeDiam,50*u.NTU,coag_graph,0*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'r-',label=r'$\alpha_\mathrm{Total}$ 0 mg/L HA')
# Magenta is 9 mg/L
plt.plot(coag_graph,floc.alpha_pacl_clay(tubeDiam,50*u.NTU,coag_graph,9*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'m:',label=r'$\alpha_\mathrm{PACl-Clay}$ 9 mg/L HA')
plt.plot(coag_graph,floc.alpha_pacl_pacl(tubeDiam,50*u.NTU,coag_graph,9*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'m-.',label=r'$\alpha_\mathrm{PACl-PACl}$ 9 mg/L HA')
plt.plot(coag_graph,floc.alpha_pacl_nat_org_mat(tubeDiam,50*u.NTU,coag_graph,9*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'m--',label=r'$\alpha_\mathrm{HA-PACl}$ 9 mg/L HA')
plt.plot(coag_graph,floc.alpha(tubeDiam,50*u.NTU,coag_graph,9*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'m-',label=r'$\alpha_\mathrm{Total}$ 9 mg/L HA')
# Brown in 15 mg/L
plt.plot(coag_graph,floc.alpha_pacl_clay(tubeDiam,50*u.NTU,coag_graph,15*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),':',color='xkcd:brown',label=r'$\alpha_\mathrm{PACl-Clay}$ 15 mg/L HA')
plt.plot(coag_graph,floc.alpha_pacl_pacl(tubeDiam,50*u.NTU,coag_graph,15*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'-.',color='xkcd:brown',label=r'$\alpha_\mathrm{PACl-PACl}$ 15 mg/L HA')
plt.plot(coag_graph,floc.alpha_pacl_nat_org_mat(tubeDiam,50*u.NTU,coag_graph,15*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'--',color='xkcd:brown',label=r'$\alpha_\mathrm{HA-PACl}$ 15 mg/L HA')
plt.plot(coag_graph,floc.alpha(tubeDiam,50*u.NTU,coag_graph,15*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'-',color='xkcd:brown',label=r'$\alpha_\mathrm{Total}$ 15 mg/L HA')
plt.xlabel(r'Coagulant Dose (mg/L)')
plt.ylabel(r'$\alpha$')
plt.axis([0, 2.7, 0, 0.5])
plt.legend(loc=2,bbox_to_anchor=(-0.24,-0.6,1.5,0.4),ncol=3,borderpad=0.1,handletextpad=0.2,labelspacing=0,columnspacing=0.2,edgecolor='white')
# plt.tight_layout()
plt.savefig('alphas.png',format='png',bbox_inches='tight')
plt.savefig('alphas.eps',format='eps',bbox_inches='tight')
plt.show()
```
## $\alpha$ vs. Dose (15 mg/L only)
```python
plt.clf()
plt.close('all')
plt.figure()

# Brown in 15 mg/L
plt.plot(coag_graph,floc.alpha_pacl_clay(tubeDiam,50*u.NTU,coag_graph,15*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),':',color='r',label=r'$\alpha_\mathrm{PACl-Clay}$ 15 mg/L HA')
plt.plot(coag_graph,floc.alpha_pacl_pacl(tubeDiam,50*u.NTU,coag_graph,15*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'-.',color='r',label=r'$\alpha_\mathrm{PACl-PACl}$ 15 mg/L HA')
plt.plot(coag_graph,floc.alpha_pacl_nat_org_mat(tubeDiam,50*u.NTU,coag_graph,15*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'--',color='r',label=r'$\alpha_\mathrm{HA-PACl}$ 15 mg/L HA')
plt.plot(coag_graph,floc.alpha(tubeDiam,50*u.NTU,coag_graph,15*u.mg/u.L,floc.HumicAcid,floc.PACl,floc.Clay,floc.RATIO_HEIGHT_DIAM),'-',color='r',label=r'$\alpha_\mathrm{Total}$ 15 mg/L HA')
plt.xlabel(r'Coagulant Dose (mg/L)')
plt.ylabel(r'$\alpha$')
plt.axis([0, 2.7, 0, 0.3])
plt.legend()
# plt.legend(loc=2,bbox_to_anchor=(-0.24,-0.6,1.5,0.4),ncol=3,borderpad=0.1,handletextpad=0.2,labelspacing=0,columnspacing=0.2,edgecolor='white')
# plt.tight_layout()
plt.savefig('alpha15.png',format='png',bbox_inches='tight')
plt.savefig('alpha15.eps',format='eps',bbox_inches='tight')
plt.show()
```
## Collision Potential vs. Effluent concentration
```python
def CPvy(L,L0):
  CP = (floc.sep_dist_clay(L, floc.Clay)**(2)-floc.sep_dist_clay(L0, floc.Clay)**(2))/(floc.Clay.Diameter*u.m)**2
  return CP

L = np.arange(0.1,100,0.1)*u.NTU

#  Creating y-axis values for data
def CPvx(ConcClay,ConcNatOrgMat,ConcAl,material,NatOrgMat,coag,DiamTube,RatioHeightDiameter,Time,EnergyDis,Temp,k):
    return 2*np.pi/3*k*floc.alpha(DiamTube, ConcClay, ConcAl, ConcNatOrgMat, NatOrgMat, coag, material, RatioHeightDiameter)*Time*(np.sqrt(EnergyDis/pc.viscosity_kinematic(Temp))).to(1/u.s)

def Eff(pC,Inf):
    return Inf*10**(-pC)    

Inf = np.array([50,100])*u.NTU

# fit k for C/C0
def C_C0_fit(N,k):
    return (2*np.pi/3*(6/np.pi)**(2/3)*k*N + 1)**(-3/2)

kfit_C_C0, kfitvar_C_C0 = curve_fit(C_C0_fit,np.append(N_fit50,N_fit50),np.append(10**(-dataset[0][0]),10**(-dataset2[0][0])))              
kfit_C_C0
plt.plot(np.append(N_fit50,N_fit50),np.append(10**(-dataset[0][0]),10**(-dataset2[0][0])),'x')
plt.plot(N_graph50,C_C0_fit(N_graph50,kfit_C_C0),'k')
plt.show()
# 50 NTU values
N_50_0 = CPvx(50*u.NTU,0*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_50_3 = CPvx(50*u.NTU,3*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_50_6 = CPvx(50*u.NTU,6*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_50_9 = CPvx(50*u.NTU,9*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_50_12 = CPvx(50*u.NTU,12*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_50_15 = CPvx(50*u.NTU,15*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)

# 100 NTU values
N_100_0 = CPvx(100*u.NTU,0*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_100_3 = CPvx(100*u.NTU,3*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_100_6 = CPvx(100*u.NTU,6*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_100_9 = CPvx(100*u.NTU,9*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_100_12 = CPvx(100*u.NTU,12*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)
N_100_15 = CPvx(100*u.NTU,15*u.mg/u.L,coag,floc.Clay,floc.HumicAcid,floc.PACl,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature,kfit_C_C0)

# Make plot
plt.clf()
plt.close('all')
plt.figure()
#50 NTU Model
plt.loglog(L,CPvy(L,50*u.NTU),label=r'$C_0=$ 50 NTU Model')
#50 NTU Data
plt.loglog(Eff(dataset[0][0],Inf[0]),N_50_0,'rx', label=r'50 NTU, 0 mg/L HA Data')
plt.loglog(Eff(dataset2[0][0],Inf[0]),N_50_0,'rx')
plt.loglog(Eff(dataset[0][1],Inf[0]),N_50_3,'b+', label=r'50 NTU, 3 mg/L HA Data')
plt.loglog(Eff(dataset2[0][1],Inf[0]),N_50_3,'b+')
plt.loglog(Eff(dataset[0][2],Inf[0]),N_50_6,'gs', markerfacecolor='none', label=r'50 NTU, 6 mg/L HA Data')
plt.loglog(Eff(dataset2[0][2],Inf[0]),N_50_6,'gs', markerfacecolor='none')
plt.loglog(Eff(dataset[0][3],Inf[0]),N_50_9,'mD', markerfacecolor='none', label=r'50 NTU, 9 mg/L HA Data')
plt.loglog(Eff(dataset2[0][3],Inf[0]),N_50_9,'mD', markerfacecolor='none')
plt.loglog(Eff(dataset[0][4],Inf[0]),N_50_12,'co', markerfacecolor='none', label=r'50 NTU, 12 mg/L HA Data')
plt.loglog(Eff(dataset2[0][4],Inf[0]),N_50_12,'co', markerfacecolor='none')
plt.loglog(Eff(dataset[0][5],Inf[0]),N_50_15,'^',c='xkcd:brown', markerfacecolor='none', label=r'50 NTU, 15 mg/L HA Data')
plt.loglog(Eff(dataset2[0][5],Inf[0]),N_50_15,'^',c='xkcd:brown', markerfacecolor='none')
# 100 NTU Model
plt.loglog(L,CPvy(L,100*u.NTU),label=r'$C_0=$ 100 NTU Model')
# 100 NTU Data
plt.loglog(Eff(dataset[1][0],Inf[1]),N_100_0,'kx', label=r'100 NTU, 0 mg/L HA Data' )
plt.loglog(Eff(dataset2[1][0],Inf[1]),N_100_0,'kx')
plt.loglog(Eff(dataset[1][1],Inf[1]),N_100_3, 'k+', label=r'100 NTU, 3 mg/L HA Data')
plt.loglog(Eff(dataset2[1][1],Inf[1]),N_100_3,'k+')
plt.loglog(Eff(dataset[1][2],Inf[1]),N_100_6,'ks', markerfacecolor='none', label=r'100 NTU, 6 mg/L HA Data')
plt.loglog(Eff(dataset2[1][2],Inf[1]),N_100_6,'ks', markerfacecolor='none')
plt.loglog(Eff(dataset[1][3],Inf[1]),N_100_9,'kD', markerfacecolor='none', label=r'100 NTU, 9 mg/L HA Data')
plt.loglog(Eff(dataset2[1][3],Inf[1]),N_100_9,'kD', markerfacecolor='none')
plt.loglog(Eff(dataset[1][4],Inf[1]),N_100_12,'ko', markerfacecolor='none', label=r'100 NTU, 12 mg/L HA Data')
plt.loglog(Eff(dataset2[1][4],Inf[1]),N_100_12,'ko', markerfacecolor='none')
plt.loglog(Eff(dataset[1][5],Inf[1]),N_100_15,'k^', markerfacecolor='none', label=r'100 NTU, 15 mg/L HA Data')
plt.loglog(Eff(dataset2[1][5],Inf[1]),N_100_15,'k^', markerfacecolor='none')
# Settings
plt.axis([3E0, 1E2, 1E2, 5E3])
plt.xlabel(r'Final Concentration (NTU)')
plt.ylabel(r'$\frac{2}{3}k\pi \overline{\alpha}\overline{G}\theta$')
plt.legend(loc=2,bbox_to_anchor=(-0.125,-0.6,1,0.37),ncol=2,borderpad=0.1,handletextpad=0.2,labelspacing=0,columnspacing=0.2,edgecolor='white')
# plt.tight_layout()
plt.savefig('performance.png',format='png',bbox_inches='tight')
plt.savefig('performance.eps',format='eps',bbox_inches='tight')
plt.show()
```
## Comparing number of HA macromolecules to number of PACl precipitates
```python
N_PACl = floc.num_nanoclusters(coag,floc.PACl)
N_PACl.to(1/u.cm**3)
N_HA = conc_humic_acid/(floc.HumicAcid.Density*u.kg/u.m**3*np.pi/6*(75*u.nm)**3)
N_HA.to(1/u.cm**3)
RatioHP = np.zeros([len(conc_humic_acid),len(coag)])
for i in range(0,len(conc_humic_acid)):
    for j in range(0,len(coag)):
        RatioHP[i][j] = N_HA[i]/N_PACl[j]
RatioHP

# Calculating max RatioHP
RatioHP_max = np.pi*(90*u.nm)**2/(np.pi/4*(75*u.nm)**2)
RatioHP_max
```
