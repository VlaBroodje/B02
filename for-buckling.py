import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from package4 import A_enclosed, y_linspace, t, b, zmax_front, zmax_rear, zmin_front, zmin_rear, Emod, x_arr, z_arr, A_arr, funcChord, I_xx, distributed_Ixx, M_func, V_func
from bendingtorsion import funcT
from constants import *

#variables
k_s = 10
k_v = 2
k_c = 4
K = 3
pi = np.pi

# geometric spar parameters
a = .04    # height and width of L stringer section
t_s = .006   # thickness of L stringer section

def L_area(a,t):
    return 2*a*t - t*t

def L_centroid(a,t):
    C = t * (a*t + a*a - t*t)/ (2*L_area(a,t))
    return C

C = L_centroid(a,t)

def L_Imin(a,t):
    return t*a**3 / 12 + 2*(t*a)*(2*(.5*a+C)**2) 

def tau_cr_spar(thickness, b_plate,k_s=10):
    return ((pi**2)*k_s*Emod)/(12*(1-poisson**2)) * ((thickness/b_plate)**2)

def column_buckling_stress(K,I_min,A,L):
    return (K*(pi**2)*Emod*I_min)/(L*L*A)

#input:
nribs = 5 #placeholder
ribspacing = (b/2)/nribs

#forces
# V_arr_func, V_arr = force_distribution(y_linspace)
V_arr = V_func(y_linspace)
Torque_arr = np.empty_like(y_linspace)

# Torque integrated to get internal torque distribution
for i,y in enumerate(y_linspace):
    T = sp.integrate.quad(funcT,y,b/2)[0]
    Torque_arr[i] = T

#second moment of inertia for I shape
A_spar = L_area(a,t_s)
I_min = L_Imin(a,t_s)

#web buckling front spar
thickness = t * funcChord(y_linspace)
b_plate_front = min(ribspacing, min((zmax_front-zmin_front)*funcChord(y_linspace)))
tau_critical_front = tau_cr_spar(thickness, b_plate_front, k_s)
print(f'Critical stress for web of the front spar {tau_critical_front}')

#web buckling rear spar
b_plate_rear = min(ribspacing, min((zmax_rear-zmin_rear)*funcChord(y_linspace)))
tau_critical_rear = tau_cr_spar(thickness, b_plate_rear, k_s)
print(f'Critical stress for web of the rear spar {tau_critical_rear}')

#average web shear  
tau_average = V_arr/((zmax_front-zmin_front)*funcChord(y_linspace)*thickness + (zmax_rear-zmin_rear)*funcChord(y_linspace)*thickness)
tau_max_average = k_v * tau_average

# shear due to torsion
q = Torque_arr/(2*A_enclosed)

tau_max_ave_front = tau_max_average + q/thickness
tau_max_ave_rear = tau_max_average - q/thickness
#print(f'Maximum applied stress for spars {round(max(tau_max_ave_front,2)), tau_max_ave_rear}')

#skin buckling
sigma_crit_top = ((pi**2)*k_c*Emod)/(12*(1-poisson**2))*((t/b_plate_front)**2)
print(f'Critical stress for skin (top) {sigma_crit_top}')
sigma_crit_bottom = ((pi**2)*k_c*Emod)/(12*(1-poisson**2))*((t/b_plate_rear)**2)
print(f'Critical stress for skin (bottom) {sigma_crit_bottom}')

nult = 3.11*1.5

sigma_max_compressive = nult*max(z_arr)*funcChord(y_linspace)*M_func(y_linspace) / distributed_Ixx(y_linspace)
print(f'Applied stress compressive {sigma_max_compressive}')
sigma_max_tensile = min(z_arr)*funcChord(y_linspace)*M_func(y_linspace) / distributed_Ixx(y_linspace)
print(f'Applied stress tensile {sigma_max_tensile}')
# https://calcresource.com/moment-of-inertia-angle.html



#column buckling
sigma_crit_buckling = column_buckling_stress(K,I_min,A_spar,ribspacing)
print(f'Critical column buckling stress: {round(sigma_crit_buckling/1e6,2)} MPa')

#margins of safety 
margin_of_safety_web_front = tau_critical_front / tau_max_average
margin_of_safety_web_rear = tau_critical_rear / tau_max_average

margin_of_safety_skin_front_compressive = sigma_crit_top/ sigma_max_compressive 
margin_of_safety_skin_front_tensile = sigma_crit_top / sigma_max_tensile

margin_of_safety_skin_rear_compressive = sigma_crit_bottom / sigma_max_compressive
margin_of_safety_skin_rear_tensile = sigma_crit_bottom / sigma_max_tensile

margin_of_safety_stringer_compressive = sigma_crit_buckling / sigma_max_compressive
margin_of_safety_stringer_tensile = sigma_crit_buckling / sigma_max_tensile

#diagramms
plt.plot(y_linspace, margin_of_safety_web_front, label='Margin of safety (web, front)')
plt.plot(y_linspace, margin_of_safety_web_rear, label='Margin of safety (web, rear)')
plt.ylim(-5,5)
plt.ylabel('Safety Margin')
plt.xlabel('Spanwise Location (m)')
plt.title('Safety Margin in Spar Web Distribution along Span')
plt.legend()
plt.show()
