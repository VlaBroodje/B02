import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import pandas as pd

# custom functions and constants
from constants import *
from functions import *

# XFLR data from postive y onwards at aoa 0
df = pd.read_csv('./Alpha0.txt', on_bad_lines='skip')
Alpha0_arr = df.to_numpy()

y_span0 = Alpha0_arr[19:,0]
chord0 = Alpha0_arr[19:,1]
aInduced0 = Alpha0_arr[19:,2]
Cl0 = Alpha0_arr[19:,3]
Cdi0 = Alpha0_arr[19:,5]
Cm0 = Alpha0_arr[19:,7]

# XFLR data from postive y onwards at aoa 10
df1 = pd.read_csv('./Alpha10.txt', on_bad_lines='skip')
Alpha10_arr = df1.to_numpy()

y_span10 = Alpha10_arr[19:,0]
chord10 = Alpha10_arr[19:,1]
aInduced10 = Alpha10_arr[19:,2]
Cl10 = Alpha10_arr[19:,3]
Cdi10 = Alpha10_arr[19:,5]
Cm10 = Alpha10_arr[19:,7]

#Cl equation determination
Cl_intercept = CL0
Cl_slope = (Cl10-Cl0)/(CL10-CL0)

#Cm equation determination
Cm_slope = (Cm10-Cm0)/(CM10-CM0)
CM_a = (CM10-CM0)/10
CL_a = (CL10-CL0)/10


def plot_deflection():
    v = deflection(y_linspace)
    print("Deflection at tip:", v[-1], "m")
    plt.plot(y_linspace,v)
    plt.title("Deflection along span")
    plt.xlabel("Span (m)")
    plt.ylabel("Deflection (m)")
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


#Number of stringers on top and bottom
top_stringers = 10    
bottom_stringers = 10
N_stringers = top_stringers + bottom_stringers
I_xx = I_zz = I_xz = .0

#geometry inputs: areas of top and bottom stringers
A_top = .00024
A_bottom = .0002

"""
x,z and area off all centroids of parts. One part has the same index number on the arrays.
I.e. the first stringer has x position x_arr[0], z position z_arr[0] and area A_arr[0]
The coordinates will be defined relative to the LE first and x,z will be relative to the chord (x=0 is LE and x=1 is TE)
"""

x_arr = np.zeros(N_stringers + 4,dtype=np.float64) # +4 accounts for the big spars
z_arr = np.zeros_like(x_arr) #zero_like generates an array of zeros with the same shape as x_arr
A_arr = np.zeros_like(x_arr)


# Define coordinates of the 4 big spars
spars = np.empty(4, dtype=object) #Array containing the 4 coordinates of the saprs
spars[0] = np.array([x_front,zmin_front])  #Bottom left
spars[1] = np.array([x_front,zmax_front])   #Top left
spars[2] = np.array([x_rear,zmax_rear])    #Top right
spars[3] = np.array([x_rear,zmin_rear])   #Bottom right

spars_MoI_results = spars_MoI(spars, t)
I_xx += spars_MoI_results[0]
I_zz += spars_MoI_results[1]
I_xz += spars_MoI_results[2]

# Assign big spars to last indeces in arrays
x_arr[N_stringers:] = spars_MoI_results[3]
z_arr[N_stringers:] = spars_MoI_results[4]
A_arr[N_stringers:] = spars_MoI_results[5]

# Assign areas to area array corresponding to stringers on top and bottom spar
A_arr[0:top_stringers] = A_top                  # Top stringers                
A_arr[top_stringers:N_stringers] = A_bottom     # Bottom stringers

# Generate equidistant points to function as stringers
x_arr[0:top_stringers] = np.linspace(x_front,x_rear,top_stringers)
x_arr[top_stringers:N_stringers] = np.linspace(x_front,x_rear,top_stringers)

z_arr[0:top_stringers] = np.linspace(zmax_front,zmax_rear,top_stringers)
z_arr[top_stringers:N_stringers] = np.linspace(zmin_front,zmin_rear,bottom_stringers)

# Calulate centroid 
Cx,Cz = determine_centroid(x_arr,z_arr,A_arr)
print(f"The centroid is {Cx},{Cz}")

# Assume centroid is also Shear Centre for now
SCx,SCz = Cx,Cz

# Determine new coordinate system relative to centroid (xp = x' and zp = z')
xp_arr = x_arr - Cx 
zp_arr = z_arr - Cz

# Calculate final moments of inertia using parallel axis theorem
I_xx += parallel_axis(A_arr,zp_arr,zp_arr)
I_zz += parallel_axis(A_arr,xp_arr,xp_arr)
I_xz += parallel_axis(A_arr,xp_arr,zp_arr)

print(f"For the entire wingbox:\nI_xx = {I_xx}\nI_zz = {I_zz}\nI_xz={I_xz}")


# Cl(y), CL, Cm(y), C(y) at alpha = 0
funcCl = sp.interpolate.CubicSpline(y_span0,Cl0)
funcCm = sp.interpolate.CubicSpline(y_span0,Cm0)

def funcChord(y):
    return Cr + (Ct - Cr)/(b/2) * y

funcCd = sp.interpolate.interp1d(y_span0, Cdi0,kind='cubic',fill_value="extrapolate")
# Extrapolated functions of dCl/dCL and dCm/dCM at every point.
funcdCldCL = sp.interpolate.CubicSpline(y_span0,Cl_slope)
funcdCmdCM = sp.interpolate.CubicSpline(y_span0,Cm_slope)
    

y_linspace = np.linspace(0,b/2,500)

# Distributed aerodynamic coefficient functions
def distributed_Cm(y,CM=CM0):
    return funcCm(y) + funcdCmdCM(y) * (CM - CM0)

def distributed_Cl(y,CL=CL0):
    return funcCl(y) + funcdCldCL(y) * (CL - CL0)

# Distributed second moments of area functions, scaled by local chord^4
# Because Ixx and Iyy are defined for chord = 1m

def distributed_Ixx(y):
    return I_xx * (funcChord(y))**4

def distributed_Izz(y):
    return I_zz * (funcChord(y))**4

def distributed_Ixz(y):
    return I_xz * (funcChord(y))**4


""" 
These function below should be changed to account for the changing local Cl as a function of big CL
(See above functions: distributed_Cm() and distributed_Cl())
"""

def lift(y,CL=CL0,q=qCruise):
    return distributed_Cl(y,CL) * funcChord(y) * q

Cd0 = 0.0488

def drag(y):
    return((funcCd(y)+Cd0/len(y_linspace))*qCruise*funcChord(y))


def moment(y,CL=CL0,CM=CM0, q=qCruise):
    print("ЕБАСИ",CM)
    c = funcChord(y)
    # Moment due to Cm
    Cm_m = distributed_Cm(y,CM) * (c**2) * q
    print(distributed_Cm(y,CM))
    # Moment due to Cl, assuming AC is at c/4
    Cl_m = distributed_Cl(y,CL) * c * q * (SCx - c/4)
    # print(Cm_m)
    return Cm_m + Cl_m


# Set up force arrays (V = shear, M = moment, D = dra)
V_arr = np.empty_like(y_linspace)
M_arr = np.empty_like(y_linspace)
D_arr = np.empty_like(y_linspace)
c_arr = funcChord(y_linspace)

#(integrated section volume/entire volue) * Mfuel
# volume equals area integrated over span
# area varies with the square of the chord
# chord varies linearly over the span

def Chord(y):
    return (Ct + (Cr-Ct)/(b/2)*(b/2-y))

def Area(y):
    return (Chord(y)**2) * sum(A_arr)

def IdilFuel(y):
   return((-funcChord(y)**2)*(-0.005007189*y+0.184802636)*775)

print("Fuel Weight",IdilFuel(b/2)-IdilFuel(0))
print(Mfuel/2)

#Critical Load Case
nFactor = -1 * 1.5 
vCrit = 231.4

def force_distribution(y_linspace,CL=CL0,q=qCruise,W_F_bool=False, nFactor =nFactor):
    L_arr = np.empty_like(y_linspace)
    Gear_W_array = np.empty_like(y_linspace)
    Fuel_W_array = np.empty_like(y_linspace)
    struct_W_array = np.empty_like(y_linspace)

    for i,y in enumerate(y_linspace):
        L_arr[i] = sp.integrate.quad(lift,y,b/2,args=(CL, q))[0]
        Fuel_W_array[i] = - (IdilFuel(b/2)-IdilFuel(y)) * 9.81*nFactor
        struct_W_array[i] = -sp.integrate.quad(Area, y, b/2)[0] * rho * 9.81*nFactor
        Gear_W_array[i] = - Mgear * 9.81*nFactor if (y<GearDist) else 0*nFactor

    V_arr = L_arr + Gear_W_array + Fuel_W_array*W_F_bool+ struct_W_array

    V_arr_func = sp.interpolate.CubicSpline(y_linspace,V_arr)

    return V_arr_func,V_arr

LiftTot = sp.integrate.quad(lift, 0, b/2) #integration of lift for one wing

def BendingMoment(x):
    return funcCl(x)*funcChord(x)*qCruise*x

TotalBendingMoment = sp.integrate.quad(BendingMoment,0,b/2)

# Interpolation verification plots
""" 
y_arr = np.linspace(0,b/2)
interp_Cm10 = np.zeros_like(y_arr)
interp_Cm5 = np.zeros_like(y_arr)
interp_Cm0 = np.zeros_like(y_arr)

interp_Cl10 = np.zeros_like(y_arr)
interp_Cl5 = np.zeros_like(y_arr)
interp_Cl0 = np.zeros_like(y_arr)

for i,y in enumerate(y_arr):
    interp_Cm10[i] = distributed_Cm(y,CM10)
    interp_Cm5[i] = distributed_Cm(y,(CM10 + CM0)/2)
    interp_Cm0[i] = distributed_Cm(y,CM0)
    interp_Cl10[i] = distributed_Cl(y,CL10)
    interp_Cl5[i] = distributed_Cl(y,(CL10 + CL0)/2)
    interp_Cl0[i] = distributed_Cl(y,CL0)

plt.subplot(211)
plt.plot(y_arr,interp_Cm10,label='interpolated cm10')
plt.plot(y_arr,interp_Cm5,label='interpolated cm5')
plt.plot(y_arr, interp_Cm0,label='interpolated cm0')
plt.plot(y_span0,Cm10,label='Cm10')
plt.plot(y_span0,Cm0,label='Cm0')
plt.title("Cm distribution verification")
plt.legend()


plt.subplot(212)
plt.plot(y_arr,interp_Cl10,label='interpolated cl10')
plt.plot(y_arr,interp_Cl5,label='interpolated cl5')
plt.plot(y_arr, interp_Cl0,label='interpolated cl0')
plt.plot(y_span0,Cl10,label='Cl10')
plt.plot(y_span0,Cl0,label='Cl0')
plt.title("Cl distribution verification")
plt.legend()

plt.show() 
"""

Load_CL = determine_CL(vCrit,nFactor, W = (MTOM-Mfuel)*9.81)
print(f"Critical Load Case CL: {Load_CL} at V = {vCrit} m/s and n = {nFactor}")
print("ALPHA",(Load_CL-0.180034)/((1.042209-0.180034)/10))

V_func,V_arr = force_distribution(y_linspace,CL=Load_CL,q=.5*vCrit**2*rhoISA)

for i,y in enumerate(y_linspace):
    M0 = sp.integrate.quad(V_func,0,b/2)
    M = sp.integrate.quad(V_func,y,b/2)
    M_arr[i]=-M[0]

# Create interpolated function from M_arr for better precision
M_func = sp.interpolate.CubicSpline(y_linspace, M_arr)

def deflection(y_linspace):
    # Correct integration: EI * d²v/dy² = M(y)
    # Integrate twice from root (y=0) with boundary conditions: v(0)=0, dv/dy(0)=0
    # First integration: dv/dy = ∫[0 to y] M(s)/(EI(s)) ds, with dv/dy(0) = 0
    ddvdy = M_func(y_linspace) / (Emod * distributed_Ixx(y_linspace))
    dvdy = sp.integrate.cumulative_trapezoid(ddvdy, y_linspace, initial=0)
    # Second integration: v = ∫[0 to y] slope(s) ds, with v(0) = 0
    v = sp.integrate.cumulative_trapezoid(dvdy, y_linspace, initial=0)
    return v


A_enclosed = 0
for i in range(len(spars)-1):
    x0, z0 = spars[i]
    x1, z1 = spars[i+1]
    A_enclosed += (x0 * z1) - (x1 * z0)

A_enclosed = abs(A_enclosed) / 2

# print("Enclosed area:", A_enclosed)
# print(f"Enclosed area at cr = {Cr} m is {A_enclosed * Cr**2} m^2")

""" 
# Print polynomial coefficients for shear, moment and deflection
print("DEFLECTION: ",np.polyfit(y_linspace, v, 5))
print("SHEAR: ",np.polyfit(y_linspace, V_arr, 8))
print("MOMENT: ",np.polyfit(y_linspace, M_arr, 8))
 """

plt.plot(y_linspace, V_arr, label='Shear Force')
plt.plot(y_linspace, M_arr, label='Bending Moment')
plt.ylabel('Force (N) / Moment (Nm)')
plt.xlabel('Spanwise Location (m)')
plt.title('Shear Force and Bending Moment Distribution along Span')
plt.legend()
plt.show()
