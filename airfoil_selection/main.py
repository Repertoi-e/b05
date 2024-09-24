import numpy as np


def isa(h_in_meters, T_offset=0):
    T_sealevel_isa = 288.15   # K
    T_tropopause_isa = 216.6  # K 
    p_sealevel_isa = 101325   # Pa
    p_tropopause_isa = 22600  # Pa
    assert(h_in_meters >= 0)
    if h_in_meters < 11_000:
        l = -0.0065  # temperature lapse in K/m
        T = T_sealevel_isa + T_offset + l * h_in_meters
        p = p_sealevel_isa * np.pow((T_sealevel_isa + l * h_in_meters)/T_sealevel_isa, -g/(l * R))
    if h_in_meters >= 11_000 and h_in_meters < 20_000:
        T = T_tropopause_isa + T_offset
        p = p_tropopause_isa * np.exp(-g * (h_in_meters - 11_000) / (R * T_tropopause_isa))
    
    rho = p / (R * T)
    
    return T, p, rho

#Constants
R = 287
g = 9.80665
M = 0.82 
altitude =  11887.2  # m
T_cr, p_cr, rho_cr = isa(altitude)

sweep_angle = (15*np.pi)/180  # deg
S = 375  # m^2
W_start = 206000 * g  # N  # MTOW
W_fuel = 73775 * g  # N  
W_end = W_start - W_fuel  # N 
a = np.sqrt(T_cr * R * 1.4) #speed of sound [m/s]
V_cr = M * a #cruise speed [m/s]



def lift_coeff(): # lec 1 slide 69, 70
    q = 1/2 * rho_cr * V_cr ** 2
    V_eff = V_cr * np.cos(sweep_angle)
    q_eff = 1/2 * rho_cr * V_eff ** 2
    cL_design = 1.1 * (1/q) * (0.5) * (W_start/S + W_end/S)
    cl_airfoil = (cL_design * q)/q_eff
    return cL_design, cl_airfoil
    

