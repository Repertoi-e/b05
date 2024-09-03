import numpy as np
import math
import matplotlib.pyplot as plt

R_AIR = 287      # J/kg/K
GRAV = 9.80665   # m/s^2
GAMMA = 1.4      # the thermodynamic degree of freedom thing

PI = math.pi
def SQRT(x): return np.sqrt(x)

T_sealevel_isa = 288.15   # K
T_tropopause_isa = 217    # K 

p_sealevel_isa = 101325   # Pa
p_tropopause_isa = 22600  # Pa

# Use T_offset=15 for e.g. a hot-day 
def isa(h_in_meters, T_offset=0):
    assert(h_in_meters >= 0)

    if h_in_meters < 11_000:
        l = -0.0065  # temperature lapse in K/m
        T = T_sealevel_isa + T_offset + l * h_in_meters    
        p = p_sealevel_isa * (T/T_sealevel_isa) ** (-GRAV/(l * R_AIR))
    if h_in_meters >= 11_000 and h_in_meters < 20_000:
        T = T_tropopause_isa + T_offset
        p = p_tropopause_isa * math.exp(-GRAV * (h_in_meters - 11_000) / (R_AIR * T))
    
    rho = p / (R_AIR * T)
    
    return T, p, rho

def speed_of_sound(T):
    return math.sqrt(GAMMA * R_AIR * T)

# TODO: Test
print("ISA 1600", isa(1600))

#
# CONSTANTS FOR THE PLANE
#

bypass_ratio = 12  # TODO placeholder value
num_engines = 2    # TODO placeholder value

M_cruise = 0.82
altitude_cr = 11887.2 # m, 39000 ft

T_cruise, _, _ = isa(altitude_cr)
V_cr = M_cruise * speed_of_sound(T_cruise)

specific_energy = 44 # kerosene, MJ/kg  TODO placeholder value

jet_efficiency = V_cr / (22 * bypass_ratio**(-0.19)*10**(-6)) / (specific_energy*10**6)

aspect_ratio = 9.4   # TODO placeholder value

c_f_equivalent = 0.00264 # estimated from Figure 6.3 in reader, TODO replace with analysis?
S_wet_ratio    = 6       # estimated from Figure 6.3 in reader, TODO replace with analysis?

c_d0_cruise = c_f_equivalent * S_wet_ratio # 6.15

oswald_phi = 0.97    # Figure 6.4 
oswald_psi = 0.0075  # Figure 6.4 
oswald_eff_cruise = 1/(PI * aspect_ratio * oswald_psi + 1/oswald_phi) # 6.17

flap_deflection_takeoff = 15    # in degrees  TODO placeholder value
flap_deflection_landing = 35    # in degrees  TODO placeholder value

def c_d0(gear_extended, flap_deflection):
    res = c_d0_cruise + 0.0013 * flap_deflection # 7.57
    if gear_extended: res += 0.02 # 7.58 Reader says to est. a value between 0.0100 and 0.0250 
    return res

def oswald_eff(flap_deflection):
    return oswald_eff_cruise + 0.0026 * flap_deflection # 7.56 wing-mounted engines 

max_LD = 1/2 * SQRT(PI * aspect_ratio * oswald_eff_cruise / c_d0_cruise)

cl_max_cruise  = 1.3  # TODO: placeholder
cl_max_takeoff = 1.7  # TODO: placeholder
cl_max_landing = 1.9  # TODO: placeholder
l_lf = 1856  # m 
l_to = 2790  # m 

mass_fr_cr = 0.95           # TODO: placeholder
mass_fr_climb = mass_fr_cr
mass_fr_landing = 0.79      # TODO: placeholder
mass_fr_app = mass_fr_landing

def cruise_speed(wing_load):
    T_cr, p_cr, rho_cr  = isa(altitude_cr)
    thr_lapse_cr = thrust_lapse(T_cr, 0.82, p_cr)

    T_W = mass_fr_cr/thr_lapse_cr * (c_d0_cruise * 0.5 * rho_cr * V_cr ** 2)/(mass_fr_cr * wing_load)
    + (mass_fr_cr * wing_load)/(PI * aspect_ratio * 0.5 * rho_cr * V_cr ** 2) # Equation 7.22
    return T_W

def min_speed(approach_speed):
    rho_app = p_sealevel_isa / (T_sealevel_isa * R_AIR)
    m_speed = (1/mass_fr_app) * (rho_app/2) * (approach_speed/1.23) ** 2 * cl_max_landing # Equation 7.6
    return m_speed

def climb_rate(wing_load, altitude_climb, rho, climb_lift_coeff): # climb rate req. is a variable right now since we don't know the altitude requirement
    # calculate climb_lift_coeff
    T_climb, p_climb, rho_climb = isa(altitude_climb)
    sound_speed_climb = speed_of_sound(T_climb)
    climb_speed = SQRT(wing_load * (2 / rho) * (1 / climb_lift_coeff))
    M_climb_rt = (climb_speed, T_climb)    
    thr_lapse_climb = thr_lapse(T, M_climb_rt, p) # add Mach number
    climb_rt_req = ... # TBD or calculate?
    climb_rate_TW = (mass_fr_climb/thr_lapse) * (SQRT((climb_rt_req ** 2/(mass_fr_climb * wing_load))
    * (rho / 2) * SQRT(c_d0_cruise * PI * aspect_ratio * oswald_eff_cruise)) + 2
    * SQRT(c_d0_cruise/(PI * aspect_ratio * oswald_eff_cruise))) # Equation 7.43

    return climb_rate_TW

def climb_gradient(thr_lapse, mass_fr, zero_lift_drag_coeff, oswald_factor, climb_grad):
    # calculate thr lapse 
    
    climb_grad_TW = (mass_fr/thr_lapse) * (climb_grad/100 + 2 * SQRT(zero_lift_drag_coeff/(PI * oswald_factor * aspect_ratio)))
    return climb_grad_TW

        
def landing_field(rho_lf, cl_lf):
    W_S = (1/mass_fr_landing) * (l_lf/cl_lf) * ((rho_lf * cl_max_landing)/2) # Equation 7.15
    return W_S       

def takeoff_field(wing_load, altitude, h2):
    T, p, rho = isa(altitude)

    cl =(1/1.13)**2 * cl_max_takeoff
    M = SQRT(wing_load * 2/rho * 1/cl)/ speed_of_sound(T)

    takeoff_field_TW = 1.15 * thrust_lapse(T, M, p)[0] * SQRT((wing_load)/(l_to* 0.85 * rho *GRAV* math.pi * aspect_ratio * oswald_eff(flap_deflection_takeoff)))+((4*h2)/(l_to)) # Equation 7.67 
    return takeoff_field_TW   


def mach(speed, T): 
    return speed / speed_of_sound(T)

def total_temperature(M, T_static): # 7.24
    return T_static * (1 + ((GAMMA - 1) / 2) * M**2)

def total_pressure(M, p_static): # 7.25
    return p_static * (1 + ((GAMMA - 1) / 2) * M**2)**(GAMMA/(GAMMA-1))

def thrust_lapse(T, M, p): # 7.29 - 7.32 
    delta_t = total_pressure(M, p) / p_sealevel_isa
    theta_t = total_temperature(M, T) / T_sealevel_isa
    
    if np.isscalar(theta_t):
        theta_t = np.array([theta_t])  # Convert scalar to 1-element array for processing

    theta_t_break = 1.07   # According to Ref. [4], modern engines have values of θt break ranging between 1.06 and 1.08.

    assert(bypass_ratio > 0 and bypass_ratio < 15)

    result = np.zeros_like(theta_t)
    if bypass_ratio < 5:
        mask = theta_t <= theta_t_break
        if mask.any():
            result[mask] = delta_t
        if (~mask).any():
            result[~mask] = delta_t * (1 - 2.1 * (theta_t[~mask] - theta_t_break) / theta_t[~mask])
    elif bypass_ratio >= 5 and bypass_ratio < 15:
        mask = theta_t <= theta_t_break
        if mask.any():
            result[mask] = delta_t * (1 - (0.43 + 0.014 * bypass_ratio) * np.sqrt(M))
        if (~mask).any():
            result[~mask] = delta_t * (1 - (0.43 + 0.014 * bypass_ratio) * np.sqrt(M) - (3 * (theta_t[~mask] - theta_t_break)) / (M + 1.5))
    return result

W_S = np.arange(0, 9000, step=1000)
W_S_min_speed = min_speed(approach_speed=78) # TODO approach speed is placeholder, should be TLAR

_, _, landing_field_rho = isa(h_in_meters=0) # TODO assuming airport is at 0 altitude
W_S_landing_field_length = landing_field(landing_field_rho, cl_lf=0.45) # Example 7.4 shows using cl_lf=0.45 for jet aircraft

W_S_cruise_speed = cruise_speed(W_S)
# W_S_climb_rate = climb_rate(W_S, altitude_climb=7400, thr_lapse=1) # TODO placeholder altitude climb


W_S_takeoff_field = takeoff_field(W_S, altitude=500, h2=11) # TODO placeholder altitude and h2

plt.vlines([W_S_min_speed, W_S_landing_field_length])
plt.plot(W_S, W_S_cruise_speed)
plt.plot(W_S, W_S_takeoff_field)
plt.title("Matching diagram")
plt.xlabel("W/S [N/m^2]")
plt.ylabel("T/W [N/N]")
plt.show()

#
# Climb rate requirements (from ADSE)
#                            25,119	25.121(a)	25.121(b)	25.121(c)	25.121(d)
# Altitude		               	  0	        0	        0	        0	        0
# No. of operating engines		  2	        1	        1	        1	        1
# Mass fraction			          1	        1	        1	        1	     0,79
# Temperature Difference		  0	        0	        0	        0	        0
# Temperature			     288,15	   288,15	    288,15	    288,15	   288,15

speed_of_sound_climb = speed_of_sound(288.15)
# W_S_climb_gradient = climb_gradient(...)
