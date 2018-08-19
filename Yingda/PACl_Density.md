# Finding PACl Density
This is an attempt to find the right value of PACl density to use for calculations.

## Load Packages
``` python
from aide_design.play import *
from aide_design import floc_model as floc
```
## Find density of $\mathrm{Al}_{30}$
Based on Casey et al., 2005
```python
W30 = 1*u.nm # Width of Al30
H30 = 2*u.nm # Length of Al30
V30 = W30**2*H30 # Volume of Al30

from periodictable import formulas as form
Al30 = form.formula("Al2O8Al28(OH)56(H2O)26")
rhoAl30 = Al30.molecular_mass*u.g/V30
rhoAl30.to(u.kg/u.m**3)
MWAl30 = (Al30.molecular_mass*u.g*u.N_A)
MWAl30.to(u.g/u.mol)

Al13 = form.formula("Al13O4(OH)24")
(Al13.molecular_mass*u.g*u.N_A).to(u.g/u.mol)

Keggin = form.formula("Al13O4(OH)24")
(Keggin.molecular_mass*u.g*u.N_A).to(u.g/u.mol)

PACl = form.formula("(Al13(OH)35Cl4)")
(PACl.molecular_mass*u.g*u.N_A).to(u.g/u.mol)
```
