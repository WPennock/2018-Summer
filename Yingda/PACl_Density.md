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
rhoAl30 = Al30.molecular_mass*u.g/V30 # Density given assumed volume of Al30
rhoAl30.to(u.kg/u.m**3)
MW_Al30 = (Al30.molecular_mass*u.g*u.N_A).to(u.g/u.mol)

Al13 = form.formula("Al13O4(OH)24") # Has 12 bound waters
MW_Al13 = (Al13.molecular_mass*u.g*u.N_A).to(u.g/u.mol)
MW_Al13

PACl = form.formula("(Al2(OH)3Cl3)")
MW_PACl = (PACl.molecular_mass*u.g*u.N_A).to(u.g/u.mol)

AlOH3 = form.formula("Al(OH)3")
MW_AlOH3 = (AlOH3.molecular_mass*u.g*u.N_A).to(u.g/u.mol)

Al2O3 = form.formula("Al2O3")
MW_Al2O3 = (Al2O3.molecular_mass*u.g*u.N_A).to(u.g/u.mol)
```
## Convert concentration as Al to as precipitate molecule
```python
Temp = (20*u.degC).to(u.K) # Assume from measurements in turbulent experiments
rho_H2O = pc.density_water(Temp).to(u.g/u.cm**3)

# From PACl Density.xmcd in Spring 2017.

C_Stock = np.array([635, 3177, 9037, 11014, 71, 70600, 70900])*u.mg/u.L # Concentration of stock in mg/L as Al
rho_Stock = np.array([0.998, 1.01, 1.02, 1.056, 0.998, 1.271, 1.267])*u.g/u.cm**3 # Measured density of stock solutions

# Using Holland data sheet from 2/13/18
rho_Holland = 1.275*rho_H2O
C_Holland = 0.106*rho_Holland*((2*MW_Al)/MW_Al2O3)
C_Holland.to(u.g/u.L)

# Numbers of Al per molecule
n_Al13 = 13 # Number of Al per Al13
n_Al30 = 30 # Number of Al per Al30
n_PACl = 2 # Number of Al per PACl
n_AlOH3 = 1 # Number of Al per AlOH3

# Transform Concentrations
def C_Al_as(C,n_Al,MW):
    return C/(n_Al*MW_Al)*MW

C_asAl13 = C_Al_as(C_Stock,n_Al13,MW_Al13)
np.mean(C_asAl13/C_Stock)
C_asAl30 = C_Al_as(C_Stock,n_Al30,MW_Al30)# Concentration as Al30
np.mean(C_asAl30/C_Stock)
C_asAlOH3 = C_Al_as(C_Stock,n_AlOH3,MW_AlOH3) # Concentration as AlOH3
np.mean(C_asAlOH3/C_Stock)
C_asPACl = C_Al_as(C_Stock,n_PACl,MW_PACl) # Concentration as PACl
np.mean(C_asPACl/C_Stock)
```
## Solve for density of precipitate molecule
$$\rho_\mathrm{Al_{13}}=\frac{\rho_\mathrm{H_2O}}{1-\frac{\rho_\mathrm{Stock}-\rho_\mathrm{H_2O}}{C_\mathrm{PACl\:as\:Al_{13}}}}$$

```python
def rho_part(rho_H2O,rho_bulk,C):
    return (rho_H2O/(1-((rho_bulk-rho_H2O)/(C)))).to(u.kg/u.m**3)

rho_Al13 = rho_part(rho_H2O,rho_Stock,C_asAl13)    
np.mean(rho_Al13.to(u.kg/u.m**3))

rho_Al30 = rho_part(rho_H2O,rho_Stock,C_asAl30)    
np.mean(rho_Al30.to(u.kg/u.m**3))

rho_AlOH3 = rho_part(rho_H2O,rho_Stock,C_asAlOH3)    
np.mean(rho_AlOH3.to(u.kg/u.m**3))

rho_PACl = rho_part(rho_H2O,rho_Stock,C_asPACl)    
np.mean(rho_PACl.to(u.kg/u.m**3))
```
## Compare with original values
```python
# Deconstructing original values
rho_PACl_xmcd = np.array([1.11E3,1.009E3,1.142E3,1.085E3,1.207E3,1.143E3,1.14E3])*u.kg/u.m**3    

C_xmcd = ((rho_Stock-rho_H2O)/(1-(rho_H2O/rho_PACl_xmcd))).to(u.mg/u.L)
C_xmcd
(C_xmcd/C_Stock).to(u.dimensionless)    
```
## Adjust density of $\mathrm{H_2O}$ for dissolved ions
```python
Base = 0.7 # Basicity given by Holland
n_OH = 3*n_PACl*Base # From basicity eqn. for PACl
R_OH = 3*Base # Obtain stoichiometric ratio of OH to Al
n_Cl = 3*n_PACl*(1 - Base) # From basicity eqn. for PACl
R_Cl = n_Cl/n_PACl # Stoichiometric ratio of Cl to Al
# Calculate molarity of Aluminum
M_Al = (0.106*rho_Holland/MW_Al2O3).to(u.mol/u.L)
# Calculate molarity of Cl
M_Cl = (R_Cl*M_Al).to(u.mol/u.L)
# Based on molarity, Cl is insignificant compared to either Na or Ca in determining density.
# Calculate molarity of hydroxide
M_OH = (R_OH*M_Al).to(u.mol/u.L)

M_Na = M_OH # Molarity of Na if NaOH used to titrate
M_Na
rho_if_Na = 1.11*u.g/u.cm**3 # Approximated from https://www.engineeringtoolbox.com/density-aqueous-solution-inorganic-sodium-salt-concentration-d_1957.html

rho_PACl_if_Na = rho_part(rho_if_Na,rho_Holland,C_Al_as(C_Holland,n_PACl,MW_PACl))
rho_PACl_if_Na.to(u.g/u.cm**3)

M_Ca = M_OH/8 # Molarity of Ca if calcium aluminate used to titrate (Based on Li et al., 2010)
M_Ca
rho_if_Ca = 1.02*u.g/u.cm**3 # Approximated from https://www.engineeringtoolbox.com/density-aqueous-solution-inorganic-chlorides-salt-concentration-d_1955.html

rho_Al13_if_Ca = rho_part(rho_if_Ca,rho_Holland,C_Al_as(C_Holland,n_Al13,MW_Al13))
rho_Al13_if_Ca

rho_Al13_if_H2O = rho_part(rho_H2O,rho_Holland,C_Al_as(C_Holland,n_Al13,MW_Al13))
rho_Al13_if_H2O
```
## Fractal $\mathrm{Al}_{13}$ approach
Assuming nanoglob reaches average size of 90 nm with fractal dimension of 2.3 from $\mathrm{Al}_{13}$ molecules of diameter 1.2 nm and density 2.42 g/cm^3.

```python
d_0 = 1.2*u.nm # Diameter of Al13
d = 90*u.nm # Diameter of nanoglob
Df = 2.3 # Fractal dimension
i = (d/d_0)**Df # Number of Al13 in nanoglob (from Weber-Shirk & Lion, 2010)

m_Al13 = Al13.molecular_mass*u.g # Mass of single Al13

m_nano = i*m_Al13 # Mass of nanoglob
V_nano = np.pi*d**3/6 # Volume of nanoglob
rho_nano = m_nano/V_nano
rho_nano.to(u.g/u.cm**3)

# Incorporating water
V_water = V_nano - i*np.pi*d_0**3/6
V_water/V_nano
M_nano = m_nano + V_water*rho_H2O
Rho_nano = M_nano/V_nano
Rho_nano.to(u.g/u.cm**3)

```
