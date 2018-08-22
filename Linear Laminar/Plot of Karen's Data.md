# Plot of Karen's data

## Import Packages
```python
from aide_design.play import *
from aide_design import floc_model as floc
from scipy.optimize import curve_fit

u.define('NTU = 100/68*mg/L')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes', facecolor='w',edgecolor='k')  
```
## Import Data
```python
# Karen's Data
Lam_5_1 = pd.read_csv("Karen's Data/Karen_5_1.csv")
Lam_15_1 = pd.read_csv("Karen's Data/Karen_15_1.csv")
Lam_50_1 = pd.read_csv("Karen's Data/Karen_50_1.csv")
Lam_150_1 = pd.read_csv("Karen's Data/Karen_150_1.csv")
Lam_500_1 = pd.read_csv("Karen's Data/Karen_500_1.csv")
Lam_5_2 = pd.read_csv("Karen's Data/Karen_5_2.csv")
Lam_15_2 = pd.read_csv("Karen's Data/Karen_15_2.csv")
Lam_50_2 = pd.read_csv("Karen's Data/Karen_50_2.csv")
Lam_150_2 = pd.read_csv("Karen's Data/Karen_150_2.csv")
Lam_500_2 = pd.read_csv("Karen's Data/Karen_500_2.csv")
Aggregated = pd.read_csv("Karen's Data/Aggregated.csv")
Lam = [Lam_5_1,Lam_15_1,Lam_50_1,Lam_150_1,Lam_500_1,Lam_5_2,Lam_15_2,Lam_50_2,Lam_150_2,Lam_500_2]

# Set up for her original assumptions.
KPACl = floc.Chemical('KPACl', (180 * u.nm).to(u.m).magnitude, 1138, 1.039,
                      'KPACl', AluminumMPM=13)
KClay = floc.Material('KClay', 2 * 10**-6, 2650, None)
u.define('KNTU = 1.4*mg/L')

# Process aggregated and individual sets with original values
Aggregated["KGamma"] = floc.gamma_coag(Aggregated["Influent (NTU)"].values*u.KNTU, Aggregated["Dose (mg/L)"].values*u.mg/u.L, KPACl, KClay, 0.375*u.inch, 0.1)
Aggregated["Old Composite Parameter"] = Aggregated["KGamma"].values*Aggregated["G (Hz)"].values*Aggregated["theta (s)"].values*Aggregated["phi"].values**(2/3)
Aggregated["Old alpha Composite Parameter"] = (2*Aggregated["KGamma"].values-Aggregated["KGamma"].values**2)*Aggregated["G (Hz)"].values*Aggregated["theta (s)"].values*Aggregated["phi"].values**(2/3)
for i in range(0,len(Lam)):
    Lam[i]["KGamma"] = floc.gamma_coag(Lam[i]["Influent (NTU)"].values*u.KNTU, Lam[i]["Dose (mg/L)"].values*u.mg/u.L, KPACl, KClay, 0.375*u.inch, 0.1)
for i in range(0,len(Lam)):
    Lam[i]['Old Composite Parameter'] = Lam[i]["KGamma"].values*Lam[i]["G (Hz)"].values*Lam[i]["theta (s)"].values*Lam[i]["phi"].values**(2/3)
for i in range(0,len(Lam)):
    Lam[i]['Old alpha Composite Parameter'] = (2*Lam[i]["KGamma"].values-Lam[i]["KGamma"].values**2)*Lam[i]["G (Hz)"].values*Lam[i]["theta (s)"].values*Lam[i]["phi"].values**(2/3)
# Process aggregated and individual sets with updated values for clay dimensions
Aggregated['MGamma'] = floc.gamma_coag(Aggregated["Influent (NTU)"].values*u.NTU, Aggregated["Dose (mg/L)"].values*u.mg/u.L, KPACl, floc.Clay, 0.375*u.inch, 0.1)
Aggregated["Mid Composite Parameter"] = Aggregated["MGamma"].values*Aggregated["G (Hz)"].values*Aggregated["theta (s)"].values*Aggregated["phi"].values**(2/3)
Aggregated["Mid alpha Composite Parameter"] = (2*Aggregated["MGamma"].values-Aggregated["MGamma"].values**2)*Aggregated["G (Hz)"].values*Aggregated["theta (s)"].values*Aggregated["phi"].values**(2/3)
for i in range(0,len(Lam)):
    Lam[i]['MGamma'] = floc.gamma_coag(Lam[i]["Influent (NTU)"].values*u.NTU, Lam[i]["Dose (mg/L)"].values*u.mg/u.L, KPACl, floc.Clay, 0.375*u.inch, 0.1)
for i in range(0,len(Lam)):
    Lam[i]["Mid Composite Parameter"] = Lam[i]["MGamma"].values*Lam[i]["G (Hz)"].values*Lam[i]["theta (s)"].values*Lam[i]["phi"].values**(2/3)
for i in range(0,len(Lam)):
    Lam[i]["Mid alpha Composite Parameter"] = (2*Lam[i]["MGamma"].values-Lam[i]["MGamma"].values**2)*Lam[i]["G (Hz)"].values*Lam[i]["theta (s)"].values*Lam[i]["phi"].values**(2/3)
```
## Plot Data
Set up intermediate values for plotting.
```python
# Set up markers
markers = ['kx','ks','kD','ko','k^','k+','ks','kD','ko','k^']
fills = ['full','none','none','none','none','full','full','full','full','full']
labels = [r'5 NTU, $\theta = $ 800 s',r'15 NTU, $\theta = $ 800 s',r'50 NTU, $\theta = $ 800 s',
         r'150 NTU, $\theta = $ 800 s',r'500 NTU, $\theta = $ 800 s',r'5 NTU, $\theta = $ 1200 s',
         r'15 NTU, $\theta = $ 1200 s',r'50 NTU, $\theta = $ 1200 s',r'150 NTU, $\theta = $ 1200 s',
         r'500 NTU, $\theta = $ 1200 s']

# Set up model plot         
N = np.arange(0.1,100,0.1)
def ModelV(N,k):
    return 3/2*np.log10(2/3*(6/np.pi)**(2/3)*k*np.pi*N+1)

# Fit k to Data
k_vO,k_vOvar  = curve_fit(ModelV,Aggregated["Old Composite Parameter"],Aggregated["pC*"])
k_vM,k_vMvar  = curve_fit(ModelV,Aggregated["Mid Composite Parameter"],Aggregated["pC*"])
k_vOa,k_vOavar  = curve_fit(ModelV,Aggregated["Old alpha Composite Parameter"],Aggregated["pC*"])
k_vMa,k_vMavar  = curve_fit(ModelV,Aggregated["Mid alpha Composite Parameter"],Aggregated["pC*"])
```
### Gamma
Look at Karen's data with all original values in terms of Gamma.
```python
# Plot Data according to original values with Gamma

plt.clf()
plt.figure(0,figsize=(7,7))
plt.xlabel(r'$\overline{\Gamma}\overline{G}\theta\phi_0^{2/3}$')
plt.ylabel(r'$\mathrm{p}C^*$')
plt.tight_layout()
plt.axis([0.1,100,0,2])
for i in range(0,len(Lam)):
    plt.semilogx(Lam[i]['Old Composite Parameter'],Lam[i]['pC*'],markers[i],fillstyle=fills[i],label=labels[i])
plt.semilogx(N,ModelV(N,k_vO),'r',label='Viscous Model')
plt.legend()
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(6, 6)
fig.savefig('Karen_Old.png',format='png')
plt.show()
plt.close()
```
Look at Karen's data with modified clay dimensions in terms of Gamma.
```python
# Plot Data according to modified values for clay in terms of Gamma
plt.clf()
plt.figure(1,figsize=(7,7))
plt.xlabel(r'$\overline{\Gamma}\overline{G}\theta\phi_0^{2/3}$')
plt.ylabel(r'$\mathrm{p}C^*$')
plt.tight_layout()
plt.axis([0.1,100,0,2])
for i in range(0,len(Lam)):
    plt.semilogx(Lam[i]['Mid Composite Parameter'],Lam[i]['pC*'],markers[i],fillstyle=fills[i],label=labels[i])
plt.semilogx(N,ModelV(N,k_vM),'r',label='Viscous Model')
plt.legend()
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(6, 6)
fig.savefig('Karen_Mid.png',format='png')
plt.show()
plt.close()
```
### Alpha
Look at Karen's data with all original values in terms of alpha.
```python
# Plot Data according to original values with alpha

plt.clf()
plt.figure(2,figsize=(7,7))
plt.xlabel(r'$\overline{\alpha}\overline{G}\theta\phi_0^{2/3}$')
plt.ylabel(r'$\mathrm{p}C^*$')
plt.tight_layout()
plt.axis([0.1,100,0,2])
for i in range(0,len(Lam)):
    plt.semilogx(Lam[i]['Old alpha Composite Parameter'],Lam[i]['pC*'],markers[i],fillstyle=fills[i],label=labels[i])
plt.semilogx(N,ModelV(N,k_vOa),'r',label='Viscous Model')
plt.legend()
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(6, 6)
fig.savefig('Karen_Old_alpha.png',format='png')
fig.savefig('Karen_Old_alpha.eps',format='eps')
plt.show()
plt.close()
```
Look at Karen's data with modified clay dimensions in terms of alpha.
```python
# Plot Data according to modified values for clay in terms of alpha
plt.clf()
plt.figure(3,figsize=(7,7))
plt.xlabel(r'$\overline{\alpha}\overline{G}\theta\phi_0^{2/3}$')
plt.ylabel(r'$\mathrm{p}C^*$')
plt.tight_layout()
plt.axis([0.1,100,0,2])
for i in range(0,len(Lam)):
    plt.semilogx(Lam[i]['Mid alpha Composite Parameter'],Lam[i]['pC*'],markers[i],fillstyle=fills[i],label=labels[i])
plt.semilogx(N,ModelV(N,k_vMa),'r',label='Viscous Model')
plt.legend()
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(6, 6)
fig.savefig('Karen_Mid_alpha.png',format='png')
plt.show()
plt.close()
```
## Fit Statistics
```python
def RMSE(A,B):
	return np.sqrt(np.mean((A-B)**2))

# Viscous Fit
RMSEv = RMSE(Aggregated["pC*"],ModelV(Aggregated["Composite Parameter"],k_v))
RMSEv
```
