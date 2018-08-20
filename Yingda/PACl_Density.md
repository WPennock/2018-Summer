# Finding PACl Density
This is an attempt to find the right value of PACl density to use for calculations.

## Load Packages
``` python
from aide_design.play import *
from aide_design import floc_model as floc
from periodictable import formulas as form
```
## Find molecular weights
```python


Al = form.formula("Al")
MW_Al = (Al.molecular_mass*u.g*u.N_A).to(u.g/u.mol)

# Based on Casey et al., 2005
W30 = 1*u.nm # Width of Al30
H30 = 2*u.nm # Length of Al30
V30 = W30**2*H30 # Volume of Al30
Al30 = form.formula("Al2O8Al28(OH)56(H2O)26")
rhoAl30 = Al30.molecular_mass*u.g/V30
rhoAl30.to(u.kg/u.m**3)
MW_Al30 = (Al30.molecular_mass*u.g*u.N_A).to(u.g/u.mol)

Al13 = form.formula("Al13O4(OH)24") # Has 12 bound waters
MW_Al13 = (Al13.molecular_mass*u.g*u.N_A).to(u.g/u.mol)
MW_Al13

PACl = form.formula("(Al13(OH)35Cl4)")
MW_PACl = (PACl.molecular_mass*u.g*u.N_A).to(u.g/u.mol)
```
## Find Density of PACl
```python
# From PACl Density.xmcd in Spring 2017.

C_Stock = np.array([635, 3177, 9037, 11014, 71, 70600, 7090])*u.mg/u.L # Concentration of stock in mg/L as Al
rho_Stock = np.array([0.998, 1.01, 1.02, 1.056, 0.998, 1.271, 1.267])*u.g/u.cm**3 # Measured density of stock solutions
n = 13 # Number of Al per Al13

C_asAl13 = C_Stock/(n*MW_Al)*MW_Al13
# C_asAl13 = C_Stock*MW_Al/MW_Al13
C_asAl13
```
$$\rho_\mathrm{Al_{13}}=\frac{\rho_\mathrm{H_2O}}{1-\frac{\rho_\mathrm{Stock}-\rho_\mathrm{H_2O}}{C_\mathrm{PACl\:as\:Al_{13}}}}$$

```python
Temp = (22*u.degC).to(u.K) # Assume from measurements in turbulent experiments
rho_H2O = pc.density_water(Temp).to(u.g/u.cm**3)
rho_Al13 = rho_H2O/(1-((rho_Stock-rho_H2O)/(C_asAl13)))
rho_Al13.to(u.kg/u.m**3)
((rho_Stock-rho_H2O)/C_asAl13).to(u.dimensionless)
```
