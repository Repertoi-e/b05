import math

R_AIR = 287      # J/kg/K
GRAV = 9.80665   # m/s^2
GAMMA = 1.4      # the thermodynamic degree of freedom thing

T_sealevel_isa = 288.15   # K
T_tropopause_isa = 217    # K 

p_sealevel_isa = 101325   # Pa
p_tropopause_isa = 22600  # Pa


# constants for the plane
cl_max_cruise = ...
cl_max_takeoff = ...
cl_max_landing = ...
l_lf = 1856 # m 
l_to = 2790 # m 
bypass_ratio = ...
num_engines = ...

def SQRT(x): return math.sqrt(x)

def cruise_speed(V_cr, rho_cr, thr_lapse_cr, wing_load, mass_fr_cr, zero_lift_drag, aspect_ratio):
    T_W = mass_fr_cr/thr_lapse_cr * (zero_lift_drag * 0.5 * rho_cr * V_cr ** 2)/(mass_fr_cr*wing_load)
    + (mass_fr_cr * wing_load)/(math.pi * aspect_ratio * 0.5 * rho_cr * V_cr ** 2)
    return T_W

def min_speed(cl_max, mass_fr, approach_speed, rho):
    m_speed = (1/mass_fr) * (rho/2) * (approach_speed/1.23) ** 2 * cl_max
    return m_speed

def climb_rate(wing_load, mass_fr, thr_lapse, climb_rt_req, rho, zero_lift_drag, aspect_ratio, oswald_eff_fr):
    climb_rate_TW = (mass_fr/thr_lapse) * (SQRT((climb_rt_req ** 2/(mass_fr * wing_load))
    * (rho / 2) * SQRT(zero_lift_drag * math.pi * aspect_ratio * oswald_eff_fr)) + 2
    * SQRT(zero_lift_drag/(math.pi * aspect_ratio * oswald_eff_fr)))
    return climb_rate_TW

def climb_gradient(thr_lapse, mass_fr, zero_lift_drag_coeff, oswald_factor, climb_grad, aspect_ratio):
    climb_grad_TW = (mass_fr/thr_lapse) * (climb_grad/100 + 2 * SQRT(zero_lift_drag_coeff/(math.pi * oswald_factor * aspect_ratio)))
    return climb_grad_TW

        
def landing_field(mass_fr_lf, rho_lf):
    W_S = (1/mass_fr_lf) * (l_lf/cl_lf) * ((rho_lf * cl_max_landing)/2)
    return W_S       

def takeoff_field(w_to, rho_to, W_S, T_to, p_to):
    cl_2 = ((1/1.13) ** 2) * cl_max_takeoff
    v2 = math.sqrt(W_S*(2/rho_to)*(1/cl_2))
    M = (v2)/(math.sqrt(GAMMA*R_AIR*T_to))
    T_t = T_to*(1+0,2*(M ** 2))
    theta_t = T_t/T_sealevel_isa
    p_t = p_to * (1+0,2*(M ** 2)) ** 3,5
    delta_t = p_t/p_sealevel_isa
    
    
    
    
    

# Use T_offset=15 for e.g. a hot-day 
def isa(h_in_meters, T_offset=0):
    assert(h_in_meters > 0)

    l = -0.0065  # temperature lapse in K/m
    

    if h_in_meters < 11_000:
        T = T_sealevel_isa + T_offset + l * h_in_meters    
        p = p_sealevel_isa * (T/T_sealevel_isa) ** (-GRAV/(l * R_AIR))
    if h_in_meters >= 11_000 and h_in_meters < 20_000:
        T = T_tropopause_isa + T_offset
        p = p_tropopause_isa * math.exp(-GRAV * (h_in_meters - 11_000) / (R_AIR * T))
    
    rho = p / (R_AIR * T)
    
    return T, p, rho

# TODO: Test
print("ISA 1600", isa(1600))

def speed_of_sound(T):
    return math.sqrt(GAMMA * R_AIR * T)

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
            return delta_t * (1 - (0.43 + 0.014 * bypass_ratio) * math.sqrt(M))
        else:
            return delta_t * (1 - (0.43 + 0.014 * bypass_ratio) * SQRT(M) - (3*(theta_t - theta_t_break))/(1.5 + M))


