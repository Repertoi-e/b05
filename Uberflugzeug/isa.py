import math

R_AIR = 287      # J/kg/K
GRAV = 9.80665   # m/s^2

T_sealevel_isa = 288.15   # K
T_tropopause_isa = 216.6  # K 

p_sealevel_isa = 101325   # Pa
p_tropopause_isa = 22600  # Pa

def isa(h_in_meters, T_offset=0):
    'Use T_offset=15 for a hot-day example'

    assert(h_in_meters >= 0)

    if h_in_meters < 11_000:
        l = -0.0065  # temperature lapse in K/m
        T = T_sealevel_isa + T_offset + l * h_in_meters
        p = p_sealevel_isa * math.pow((T_sealevel_isa + l * h_in_meters)/T_sealevel_isa, -GRAV/(l * R_AIR))
    if h_in_meters >= 11_000 and h_in_meters < 20_000:
        T = T_tropopause_isa + T_offset
        p = p_tropopause_isa * math.exp(-GRAV * (h_in_meters - 11_000) / (R_AIR * T_tropopause_isa))
    
    rho = p / (R_AIR * T)
    
    return T, p, rho