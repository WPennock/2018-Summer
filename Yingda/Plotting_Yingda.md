# -*- coding: utf-8 -*-
# Plotting Yingda

## Load Packages and Set Parameters
```python
from aide_design.play import *
from aide_design import floc_model as floc
from scipy.optimize import curve_fit    
import pdb

u.NTU = 100/68*u.mg/u.L

#k = 0.23 # had been 0.18
coag = np.array([0.53, 1.06, 1.59, 2.11, 2.56]) * u.mg/u.L
conc_humic_acid = np.array([0, 3, 6, 9, 12, 15]) * u.mg/u.L
floc.HumicAcid.Density = 1520
#dataset[0] is the 50NTU, dataset[1] is the 100NTU.
#Within both subgroups, [0] is the pC.0, ranging evenly up to [5] which is the
# pC.15
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
dataset_pC = 50*10**(-dataset[0])
coagGraph = np.arange(1 * 10**-4, 26.1 * 10**-4, 1 * 10**-4) * u.kg/u.m**3 # Same as Mathcad sheet
coag_graph = np.linspace(5*10**-5,26*10**-4,num=100)*u.kg/u.m**3
enerDis = 4.833 * u.mW/u.kg # Same as Mathcad sheet
temperature = 25 * u.degC # Same as Mathcad sheet
resTime = 302 * u.s # Same as Mathcad sheet
tubeDiam = 3/8 * u.inch # Same as Mathcad sheet
# Material properties also the same as in Mathcad sheet

# Troubleshooting discrepancy with Mathcad
# Problem was with odd calculation of phi w/in floc_model pc_viscous function.
@u.wraps(None, [u.W/u.kg, u.degK, u.s, u.m,
                u.kg/u.m**3, u.kg/u.m**3, u.kg/u.m**3, None,
                None, None, u.dimensionless, u.dimensionless], False)
def pc_viscous(EnergyDis, Temp, Time, DiamTube,
               ConcClay, ConcAl, ConcNatOrgMat, NatOrgMat,
               coag, material, FittingParam, RatioHeightDiameter):
    return ((3/2)
            * np.log10((2/3) * np.pi * FittingParam * Time
                       * np.sqrt(EnergyDis
                                 / (pc.viscosity_kinematic(Temp).magnitude)
                                 )
                       * floc.alpha(DiamTube, ConcClay, ConcAl, ConcNatOrgMat,
                               NatOrgMat, coag, material, RatioHeightDiameter)
                       * floc.frac_vol_floc_initial(ConcAl,ConcClay,coag,material)**(2/3)
                       + 1
                       )
            )

```
## Fitting Humic Acid Diameter and $k$
```python

# Fitting k
def N_viscous (ConcClay,ConcAl,coag,material,DiamTube,RatioHeightDiameter,Time,EnergyDis,Temp):
    x = floc.alpha(DiamTube, ConcClay, ConcAl, 0*u.mg/u.L, floc.HumicAcid, coag, material, RatioHeightDiameter)*Time*(np.sqrt(EnergyDis/pc.viscosity_kinematic(Temp))).to(1/u.s)*floc.frac_vol_floc_initial(ConcAl, ConcClay, coag, material)**(2/3)
    return x.to(u.dimensionless)

def viscous_fit(N,k):
    return 1.5*np.log10(2.0/3.0*np.pi*k*N*(6.0/np.pi)**(2.0/3.0) + 1)

N_fit50 = N_viscous(50*u.NTU,coag,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
N_fit100 = N_viscous(100*u.NTU,coag,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
N_graph = N_viscous(50*u.NTU,coag_graph,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
## Fit k for combinded 0 mg/L HA data
kfit, kfitvar= curve_fit(viscous_fit,np.concatenate([N_fit50,N_fit100]),np.concatenate([dataset[0][0],dataset[1][0]]))              
kfit
## verify fit
N_graph50  = N_viscous(50*u.NTU,coag_graph,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
N_graph100  = N_viscous(100*u.NTU,coag_graph,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
plt.plot(np.concatenate([N_fit50,N_fit100]),np.concatenate([dataset[0][0],dataset[1][0]]),'x')
N_plot = np.linspace(2.5,16,num=100)
plt.plot(N_plot,viscous_fit(N_plot,kfit),'k')
plt.show()
```
## Fitting $d_HA$
### Brute Force method
#### 50 NTU
```python
def SSE(A,B):
    '''Function for determining sum-squared-error'''
    return np.sum((A-B)**2)

conv = 0.0 #criterion for pC*
dHA_range = np.arange(1,501,1)*u.nm
dHA_50i = np.ones(len(conc_humic_acid))
SSE50 = np.zeros([len(conc_humic_acid),len(dHA_range)])
for i in range (0,len(conc_humic_acid)):    
    for j in range(0,len(dHA_range)):
        floc.HumicAcid.Diameter = dHA_range[j].to(u.m).magnitude
        pCpred = pc_viscous(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coag, conc_humic_acid[i], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
        # pdb.set_trace()
        SSE50[i][j] = SSE(pCpred[np.where(dataset[0][i]>conv)],dataset[0][i][np.where(dataset[0][i]>conv)])
        np.where(dataset[0][i]>conv)[0]
    dHA_50i[i] = np.argmin(SSE50[i])    
        # SSE of pC*
        # Return array of SSE
dHA_50i
SSE50[2]
# Return index of min of SSE
```
#### 100 NTU
```python
dHA_100i = np.ones(len(conc_humic_acid))
SSE100 = np.zeros([len(conc_humic_acid),len(dHA_range)])
for i in range (0,len(conc_humic_acid)):    
    for j in range(0,len(dHA_range)):
        floc.HumicAcid.Diameter = dHA_range[j].to(u.m).magnitude
        pCpred = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coag, conc_humic_acid[i], floc.HumicAcid, floc.PACl, floc.Clay, kfit, floc.RATIO_HEIGHT_DIAM)
        # pdb.set_trace()
        SSE100[i][j] = SSE(pCpred[np.where(dataset[1][i]>conv)],dataset[1][i][np.where(dataset[0][i]>conv)])
        np.where(dataset[1][i]>conv)[0]
    dHA_100i[i] = np.argmin(SSE100[i])    
        # Return array of SSE
dHA_100i
dHA_fiti = np.mean([dHA_50i[2:],dHA_100i[2:]])
dHA_fiti
dHA_fit = dHA_range[int(np.round(dHA_fiti))]
floc.HumicAcid.Diameter = dHA_fit.to(u.m).magnitude
floc.HumicAcid.Diameter
kfit
```
### Elegant method (curve_fit)
```python
# def pc_fit_dHA(dataset,dHA):
#     '''# For dataset, 0th row is coagulant, 1st is humic acid, 2nd is influent turbidity'''
#     G_Coag = floc.gamma_coag(dataset[2],dataset[0],floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM)
#     G_HA = np.minimum((((dataset[1].to(u.kfitg/u.m**3)).magnitude / (floc.conc_precipitate(dataset[0], floc.PACl).to(u.kfitg/u.m**3)).magnitude) * (floc.PACl.Density / floc.HumicAcid.Density) * (floc.PACl.Diameter / (4 * dHA)) ), np.ones(len(dataset[0])))
#     pred = np.zeros(len(dataset[0]))
#     for i in range(0,len(dataset[0])):
#         pred[i] = ((3/2) * np.log10((2/3) * np.pi * kfit * resTime.to(u.s).magnitude * (np.sqrt(enerDis.to(u.m**2/u.s**3) / (pc.viscosity_kfitinematic(temperature)) ).to(1/u.s)).magnitude * (2*(1-G_Coag[i])*(G_Coag[i]*(1-G_HA[i])) + (G_Coag[i]*(1-G_HA[i]))**2 + 2*(G_Coag[i]*(1-G_HA[i]))*(G_HA[i]*G_Coag[i]) ) * (np.pi/6)**(2/3) * (floc.Clay.Diameter / (floc.sep_dist_clay(dataset[2][i], floc.Clay).to(u.m)).magnitude ) ** 2 + 1 ) )
#     return pred    
# help(floc.gamma_humic_acid_to_coag)
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
```

#Begin graphing the 100NTU datasets
```python
plt.subplot(122)
plt.title('100 NTU Graph')
plt.ylabel('pC*')
plt.xlabel('coagulant dosage (mg/L)')

plt.plot(coag, dataset[1][0], 'r.', coag, dataset[1][1], 'b.', coag, dataset[1][2], 'g.',
         coag, dataset[1][3], 'm.', coag, dataset[1][4], 'c.', coag, dataset[1][5], 'y.')

line0mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 0 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, k100, floc.RATIO_HEIGHT_DIAM)
line3mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 3 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, k100, floc.RATIO_HEIGHT_DIAM)
line6mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 6 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, k100, floc.RATIO_HEIGHT_DIAM)
line9mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 9 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, k100, floc.RATIO_HEIGHT_DIAM)
line12mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 12 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, k100, floc.RATIO_HEIGHT_DIAM)
line15mg100 = pc_viscous(enerDis, temperature, resTime, tubeDiam, 100 * u.NTU, coagGraph, 15 * u.mg/u.L, floc.HumicAcid, floc.PACl, floc.Clay, k100, floc.RATIO_HEIGHT_DIAM)

x = coagGraph.to(u.mg/u.L)
plt.plot(x, line0mg100, 'r', x, line3mg100, 'b', x, line6mg100, 'g',
         x, line9mg100, 'm', x, line12mg100, 'c', x, line15mg100, 'y')
```


#And now we display our graph!
```python
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
```python
plt.rcParams['text.latex.preamble']=[r"\usepackage{txfonts}"]
params = {'text.usetex' : True,
          'font.size' : 14,
          'font.family' : 'serif',
          'text.latex.unicode': True,
          'axes.facecolor': 'white',
          'savefig.facecolor': 'white',
          'axes.edgecolor': 'black',
          'savefig.edgecolor': 'black'
          }
plt.rcParams.update(params)
```

## 50 NTU Data
```python
# Create function for old model (no HA)
#@u.wraps(None, [u.W/u.kg, u.degK, u.s, u.m,
#                u.kg/u.m**3, u.kg/u.m**3,
#                None, None, u.dimensionless, u.dimensionless], False)
#def pc_viscous_old(EnergyDis, Temp, Time, DiamTube,
#               ConcClay, ConcAl, coag, material,
#               FittingParam, RatioHeightDiameter):
#    return ((3/2)
#            * np.log10((2/3) * np.pi * FittingParam * Time
#                       * np.sqrt(EnergyDis
#                                 / (pc.viscosity_kinematic(Temp).magnitude)
#                                 )
#                       * (2*floc.gamma_coag(ConcClay,ConcAl,coag,material,DiamTube,RatioHeightDiameter)-floc.gamma_coag(ConcClay,ConcAl,coag,material,DiamTube,RatioHeightDiameter)**2)
#                       * floc.frac_vol_floc_initial(ConcAl,ConcClay,coag,material)**(2/3)
#                       + 1
#                       )
#            )

def N_old(ConcClay,ConcAl,coag,material,DiamTube,RatioHeightDiameter,Time,EnergyDis,Temp):
	x = (2*floc.gamma_coag(ConcClay, ConcAl, coag, material, DiamTube, RatioHeightDiameter)-floc.gamma_coag(ConcClay, ConcAl, coag, material, DiamTube, RatioHeightDiameter)**2)*Time*(np.sqrt(EnergyDis/pc.viscosity_kinematic(Temp))).to(1/u.s)*floc.frac_vol_floc_initial(ConcAl, ConcClay, coag, material)**(2/3)
	return x.to(u.dimensionless)

N_fit_old = N_old(50*u.NTU,coag,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)
N_graph_old = N_old(50*u.NTU,coag_graph,floc.PACl,floc.Clay,tubeDiam,floc.RATIO_HEIGHT_DIAM,resTime,enerDis,temperature)

def viscous_fit(N,k):
    return 1.5*np.log10(2.0/3.0*np.pi*k*N*(6.0/np.pi)**(2.0/3.0) + 1)

k_old, kvar_old = curve_fit(viscous_fit,N_fit_old,dataset[0][0])              

plt.clf()
plt.close('all')
plt.figure(0)
plt.plot(coag,dataset[0][0],'--rx', label=r'0 mg/L HA')
plt.plot(coag,dataset[0][1],'--b+', label=r'3 mg/L HA')
plt.plot(coag,dataset[0][2],'--gs', markerfacecolor='none', label=r'6 mg/L HA')
plt.plot(coag,dataset[0][3],'--mD', markerfacecolor='none', label=r'9 mg/L HA')
plt.plot(coag,dataset[0][4],'--co', markerfacecolor='none', label=r'12 mg/L HA')
plt.plot(coag,dataset[0][5],'--^',c='xkcd:brown', markerfacecolor='none', label=r'12 mg/L HA')
plt.plot(coag_graph.to(u.mg/u.L)[19:],viscous_fit(N_graph_old,k_old)[19:],'k',label=r'Pennock et al. (2018)')
#plt.plot(x[4:],pc_viscous_old(enerDis, temperature, resTime, tubeDiam, 50 * u.NTU, coagGraph[4:], floc.PACl, floc.Clay, k50, floc.RATIO_HEIGHT_DIAM),'k',label=r'Pennock et al. (2018)')
plt.xlabel(r'Coagulant Dose (mg/L)')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0, 3, 0, 1.7])
plt.legend(loc='upper left',ncol=2,borderpad=0.1,handletextpad=0.1,labelspacing=0,columnspacing=0.1,edgecolor='inherit')
plt.tight_layout()
dataset[0][0]
plt.savefig('50NTU.png',format='png')
plt.savefig('50NTU.eps',format='eps')
plt.show()
```

## Figure 5
```python
plt.clf()
plt.close('all')
plt.figure(5)
plt.plot(coag,dataset[0][0],'rx', label=r'0 mg/L HA Data')
plt.plot(coag,dataset[0][1],'b+', label=r'3 mg/L HA Data')
plt.plot(coag,dataset[0][2],'gs', markerfacecolor='none', label=r'6 mg/L HA Data')
plt.plot(coag,dataset[0][3],'mD', markerfacecolor='none', label=r'9 mg/L HA Data')
plt.plot(coag,dataset[0][4],'co', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.plot(coag,dataset[0][5],'^',c='xkcd:brown', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.plot(x, line0mg50, 'r', label=r'0 mg/L HA Model')
plt.plot(x, line3mg50, 'b', label=r'3 mg/L HA Model')
plt.plot(x, line6mg50, 'g', label=r'6 mg/L HA Model')
plt.plot(x, line9mg50, 'm', label=r'9 mg/L HA Model')
plt.plot(x, line12mg50, 'c', label=r'12 mg/L HA Model')
plt.plot(x, line15mg50, 'xkcd:brown', label=r'15 mg/L HA Model')
plt.xlabel(r'Coagulant Dose (mg/L)')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0,3,0,1.5])
plt.legend(loc=2,bbox_to_anchor=(0,-0.6,1,0.4),ncol=2,borderpad=0.1,handletextpad=0.2,labelspacing=0,columnspacing=0.2,edgecolor='white')
#plt.tight_layout()

plt.savefig('DuFig5.png',format='png', bbox_inches = "tight")
plt.savefig('DuFig5.eps',format='eps', bbox_inches = "tight")
```

## Figure 9
```python
plt.clf()
plt.close('all')
plt.figure(9)
plt.plot(coag,dataset[1][0],'rx', label=r'0 mg/L HA Data')
plt.plot(coag,dataset[1][1],'b+', label=r'3 mg/L HA Data')
plt.plot(coag,dataset[1][2],'gs', markerfacecolor='none', label=r'6 mg/L HA Data')
plt.plot(coag,dataset[1][3],'mD', markerfacecolor='none', label=r'9 mg/L HA Data')
plt.plot(coag,dataset[1][4],'co', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.plot(coag,dataset[1][5],'^',c='xkcd:brown', markerfacecolor='none', label=r'12 mg/L HA Data')
plt.plot(x, line0mg100, 'r', label=r'0 mg/L HA Model')
plt.plot(x, line3mg100, 'b', label=r'3 mg/L HA Model')
plt.plot(x, line6mg100, 'g', label=r'6 mg/L HA Model')
plt.plot(x, line9mg100, 'm', label=r'9 mg/L HA Model')
plt.plot(x, line12mg100, 'c', label=r'12 mg/L HA Model')
plt.plot(x, line15mg100, 'xkcd:brown', label=r'15 mg/L HA Model')
plt.xlabel(r'Coagulant Dose (mg/L)')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0,3,0,1.5])
plt.legend(loc=2,bbox_to_anchor=(0,-0.6,1,0.37),ncol=2,borderpad=0.1,handletextpad=0.2,labelspacing=0,columnspacing=0.2,edgecolor='white')
#plt.tight_layout()

plt.savefig('DuFig9.png',format='png', bbox_inches = "tight")
plt.savefig('DuFig9.eps',format='eps', bbox_inches = "tight")
```
