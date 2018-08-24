# Supporting Calculations for Yingda's Paper
```python
from aide_design.play import *
from aide_design import floc_model as floc
import pdb
```
## Calculation of $G$
### Parameters
```python
Q = 6*u.mL/u.s
ID = 9.52*u.mm
A = np.pi*ID**2/4
V = (Q/A).to(u.m/u.s)
V.to(u.mm/u.s)
R_c = 15*u.cm
Temperature = 25*u.degC
nu = pc.viscosity_kinematic(Temperature)
L = 25.45*u.m
T = (L/V).to(u.s)
Re = (V*ID/(nu)).to(u.dimensionless)
De = (Re*(ID/R_c)**(1/2)).to(u.dimensionless)
```
### Navier-Stokes
```python
G_Straight = floc.g_straight(Q,ID) # Navier-Stokes
G_Straight
G_Coil = floc.g_coil(Q,ID,R_c,Temperature)
G_Coil
Ratio_Coil = G_Coil/G_Straight
Ratio_Coil
G_s_Yingda = G_Coil/Ratio_Coil
G_s_Yingda
```
### Darcy-Weisbach
```python
f_Yingda = pc.fric(Q,ID,nu,0.001)
f_Yingda
f_Coil = Ratio_Coil*f_Yingda
f_Coil
hf_Coil = (L*f_Coil_Calc*8/np.pi**2/u.g_0*Q**2/ID**5).to(u.m)
hf_Coil
hf_Coil_Calc = Ratio_Coil_C*(32*nu*L*V/(u.g_0*ID**2))
hf_Coil_Calc.to(u.m)
EDR_Coil = (hf_Coil*u.g_0/T).to(u.mW/u.kg)
EDR_Coil
G_Calc = np.sqrt(EDR_Coil/nu).to(1/u.s)
G_Calc
G_Coil
```
### Monroe's Equation
```python
G_M = 4*np.sqrt(2)*V/ID*(1+0.033*(np.log10(De))**4)**(1/2)
G_M.to(u.s**(-1))

# Trying to match Monroe's result with Darcy-Weisbach
Ratio_Coil_C= 1+0.033*(np.log10(De))**4
f_Coil_C= f_Coil*Ratio_Coil_C
hf_Coil_C= Ratio_Coil_C*(32*nu*L*V/(u.g_0*ID**2))
EDR_Coil_C = (hf_Coil_C*u.g_0/T).to(u.mW/u.kg)
EDR_Coil_C
G_Calc_C = np.sqrt(EDR_Coil_C/nu)
G_Calc_C.to(1/u.s)
```
### Calculation in Paper
```python
hf_Straight = (32*nu*L*V/(u.g_0*ID**2))
hf_Straight.to(u.m)
EDR_Straight = u.g_0*hf_Straight/T
EDR_Straight.to(u.mW/u.kg)
G_Straight = np.sqrt(EDR_Straight/nu)
G_Straight.to(1/u.s)
G_Coil_Paper = G_Straight*np.sqrt(1+0.033*np.log10(De)**4)
G_Coil_Paper.to(1/u.s)
```

## Analyze Goodness of Extracted 100 NTU Points
```python
dataset = np.array([[0.746, 0.953, 1.191, 1.295, 1.414],
 [0.563, 0.835, 1.085, 1.255, 1.403],
 [0.185, 0.692, 0.971, 1.254, 1.390],
 [0.105, 0.280, 0.956, 1.238, 1.361],
 [0.097, 0.207, 0.740, 1.209, 1.316],
 [0.084, 0.157, 0.566, 1.084, 1.314]
 ]
)
extract = pd.read_csv('Modified Axes.csv',header=-1)
extract
extract3 = pd.read_csv('Modified_12.csv',header=-1)

points = np.zeros([5,12])
points[0] = extract[1][extract[0].between(0.4,0.6)]
points[1] = extract[1][extract[0].between(1.0,1.1)]
points[2] = extract[1][extract[0].between(1.5,1.6)]
points[3] = extract3[1]
points[4] = extract[1][extract[0].between(2.6,2.7)]

data = np.transpose(dataset)

filtered = np.zeros(np.shape(data))

def find_nearest(guesses, data):
    for i in range(0,len(data)):
        for j in range(0,len(data[0])):
            filtered[i][j] = (np.abs(guesses[i]- data[i][j])).argmin()
    return filtered

exclude = find_nearest(points,data)    
exclude
data[4]
points[4]
exclude[4] = [0, 2, 4, 7, 8, 11]
exclude

difference = np.zeros(np.shape(data))
for i in range(0,len(data)):
    for j in range(0,len(data[0])):
        # pdb.set_trace()
        difference[i][j] = (points[i][int(exclude[i][j])]-data[i][j])/data[i][j]
difference*100

replicates = np.zeros(np.shape(data))
for i in range(0,len(replicates)):
    # pdb.set_trace()
    replicates[i] = np.sort(np.delete(points[i],exclude[i]))[::-1]
replicates = np.transpose(replicates)
replicates
```
## Calculation of Settler Flow
```python
D_S = 2.66*u.cm
L_S = 1.37*u.m
alpha_S = 60*u.deg
V_c = 0.12*u.mm/u.s
Q_S = np.pi/4*D_S**2*V_c*(L_S/D_S*np.cos(alpha_S)+np.sin(alpha_S))
Q_S.to(u.mL/u.s)
L_S.to(u.inch)
(4.5*u.ft).to(u.inch)
D_S.to(u.inch)

# What was actual capture velocity if Q_S = 1.5 mL/s?
Q_Actual = 1.5*u.mL/u.s
V_c_Actual = Q_Actual*4/np.pi/D_S**2/(L_S/D_S*np.cos(alpha_S)+np.sin(alpha_S))
V_c_Actual.to(u.mm/u.s)
```
