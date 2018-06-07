# Supporting Calculations for Yingda's Paper
```python
from aide_design.play import *
from aide_design import floc_model as floc
```
## Calculation of $G$
```python
Q = 6*u.mL/u.s
ID = 9.52*u.mm
A = np.pi*ID**2/4
R_c = 15*u.cm
Temperature = 25*u.degC
nu = pc.viscosity_kinematic(Temperature)
L = 25.45*u.m
T = (L*A/Q).to(u.s)
# Using what is
G_Straight = floc.g_straight(Q,ID)
G_Coil = floc.g_coil(Q,ID,R_c,Temperature)
Ratio_Coil = G_Coil/G_Straight
G_s_Yingda = G_Coil/Ratio_Coil
# Attempting from scratch
Re_Yingda = (Q/(np.pi*ID**2/4)*ID/(nu)).to(u.dimensionless)
f_Yingda = pc.fric(Q,ID,nu,0.001)
f_Coil = Ratio_Coil*f_Yingda
hf_Coil = (L*f_Coil*8/np.pi**2/u.g_0*Q**2/ID**5).to(u.m)
EDR_Coil = (hf_Coil*u.g_0/T).to(u.mW/u.kg)
G_Calc = np.sqrt(EDR_Coil/nu).to(1/u.s)
G_Calc
G_Coil
```
