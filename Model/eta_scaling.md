# Scaling relative velocity by $\eta$
If the length scale controlling relative velocity is $\eta$, the Kolmogorov microscale, then both the viscous and inertial equations will integrate to the same form:

$$a$$

## Load Packages
``` python
from aide_design.play import *
from aide_design import floc_model as floc
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes', facecolor='w',edgecolor='k')
```

## Import Data
```python
Turb = pd.read_csv("Turbulent.csv")
# Adjust Gamma
u.define('NTU = 100/68*mg/L')
# Adjust Gamma for new coagulant diameter.
Turb["Gamma"] = floc.gamma_coag(Turb["Influent Turbidity (NTU)"].values*u.NTU, Turb["PACl Dose (mg/L)"].values*u.mg/u.L, floc.PACl, floc.Clay, 1.25*u.inch, 0.1)         
Turb["eta Parameter"] = (2*Turb['Gamma'].values-Turb['Gamma'].values**2)*6/(floc.Clay.Diameter*u.m)*Turb["phi"].values*(Turb["Kinematic Viscosity (mm^2/s)"].values*u.mm**2/u.s)**(1/4)*(Turb["Energy Dissipation Rate (mW/kg)"].values*u.mW/u.kg)**(1/4)*Turb["Residence Time (s)"].values*u.s
Turb["eta Parameter"].values
```
## Plotting
```python
def eta_model(param,k):
    return -np.log10(-k*param+1)
eta_N = np.arange(0,5e5,1e3)
plt.clf()
plt.figure(0,figsize=(12, 8))
plt.plot(Turb["eta Parameter"],Turb["pC*"],'kx',label='Turbulent Data')
plt.plot(eta_N,eta_model(eta_N,0.5),'r',label=r'$\eta$ Model')
plt.legend(loc=2)
plt.text(-1.5,1,'k = '+str(ut.sig(k_i[0],3)))
plt.xlabel(r'$\log_{10}\left[\frac{8}{9}k\pi \varepsilon^{1/3}\theta\alpha\right]$')
plt.ylabel(r'$\log_{10}\left[\frac{\Lambda^{3n}-\Lambda_0^{3n}}{D^{3n-2/3}}\right]$')
plt.savefig('Linear_Plot_Inertial.png',format='png')
plt.show()
```
$$\rlap{\kern.08em--}V$$
