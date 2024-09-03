import numpy as np
import math

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
    assert(h_in_meters > 0)

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

jet_efficiency = V_cr / (22 * bypass_ratio**(-0,19)*10**(-6)) / (specific_energy*10**6)

aspect_ratio = 9   # TODO placeholder value

c_f_equivalent = 0.00264 # estimated from Figure 6.3 in reader, TODO replace with analysis?
S_wet_ratio    = 6       # estimated from Figure 6.3 in reader, TODO replace with analysis?

c_d0_cruise = c_f_equivalent * S_wet_ratio # 6.15

oswald_phi = 0.97    # Figure 6.4 
oswald_psi = 0.0075  # Figure 6.4 
oswald_eff_cruise = 1/(PI * aspect_ratio * oswald_psi + 1/oswald_phi) # 6.17

flap_deflection_takeoff = 15    # in degrees  TODO placeholder value
flap_deflection_landing = 35    # in degrees  TODO placeholder value

def c_d0(gear_exteded, flap_deflection):
    res = c_d0_cruise + 0.0013 * flap_deflection # 7.57
    if gear_exteded: res += 0.02 # 7.58 Reader says to est. a value between 0.0100 and 0.0250 
    return res

def oswald_eff(flap_deflection):
    return oswald_eff_cruise + 0.0026 * flap_deflection # 7.56 wing-mounted engines 

max_LD = 1/2 * SQRT(PI * aspect_ratio * oswald_eff_cruise / c_d0_cruise)

cl_max_cruise  = 1.3  # TODO: placeholder
cl_max_takeoff = 1.7  # TODO: placeholder
cl_max_landing = 1.9  # TODO: placeholder
l_lf = 1856  # m 
l_to = 2790  # m 

mass_fr_cr = ...
mass_fr_climb = ...
mass_fr_landing = ...
mass_fr_app = mass_fr_landing

def cruise_speed(wing_load):
    T_cr, p_cr, rho_cr  = isa(altitude_cr)
    thr_lapse_cr = thrust_lapse(T_cr, 0.82, p_cr)
    zero_lift_drag = c_d0_cruise

    T_W = mass_fr_cr/thr_lapse_cr * (zero_lift_drag * 0.5 * rho_cr * V_cr ** 2)/(mass_fr_cr*wing_load)
    + (mass_fr_cr * wing_load)/(PI * aspect_ratio * 0.5 * rho_cr * V_cr ** 2) # Equation 7.22
    return T_W

def min_speed(approach_speed):
    rho_app = p_sealevel_isa / (T_sealevel_isa * R_AIR)
    m_speed = (1/mass_fr_app) * (rho_app/2) * (approach_speed/1.23) ** 2 * cl_max_landing # Equation 7.6
    return m_speed

def climb_rate(wing_load, thr_lapse, climb_rt_req, rho, zero_lift_drag):

    climb_rate_TW = (mass_fr_climb/thr_lapse) * (SQRT((climb_rt_req ** 2/(mass_fr_climb * wing_load))
    * (rho / 2) * SQRT(zero_lift_drag * PI * aspect_ratio * oswald_eff_cruise)) + 2
    * SQRT(zero_lift_drag/(PI * aspect_ratio * oswald_eff_cruise))) # Equation 7.43

    return climb_rate_TW

def climb_gradient(thr_lapse, mass_fr, zero_lift_drag_coeff, oswald_factor, climb_grad):
    climb_grad_TW = (mass_fr/thr_lapse) * (climb_grad/100 + 2 * SQRT(zero_lift_drag_coeff/(PI * oswald_factor * aspect_ratio)))
    return climb_grad_TW

        
def landing_field(rho_lf, cl_lf):
    W_S = (1/mass_fr_landing) * (l_lf/cl_lf) * ((rho_lf * cl_max_landing)/2) # Equation 7.15
    return W_S       

def takeoff_field(rho_to, W_S, h2, T_to, p_to, M_to):
    takeoff_field_TW = 1,15 * thrust_lapse(T_to, M_to, p_to) * SQRT((W_S)/(l_to* 0,85 * rho_to *GRAV* math.pi * aspect_ratio * oswald_eff(flap_deflection_takeoff)))+((4*h2)/(l_to)) # Equation 7.67 
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

    theta_t_break = 1.07   # According to Ref. [4], modern engines have values of θt break ranging between 1.06 and 1.08.

    assert(bypass_ratio > 0 and bypass_ratio < 15)
    if bypass_ratio < 5:
        if theta_t <= theta_t_break:
            return delta_t
        else:
            return delta_t * (1 - 2.1 * (theta_t - theta_t_break)/theta_t)
    if bypass_ratio >= 5 and bypass_ratio < 15:
        if theta_t <= theta_t_break:
            return delta_t * (1 - (0.43 + 0.014 * bypass_ratio) * SQRT(M))
        else:
            return delta_t * (1 - (0.43 + 0.014 * bypass_ratio) * SQRT(M) - (3*(theta_t - theta_t_break))/(1.5 + M))


W_S = np.arrange(0, 9000, step=1000)
