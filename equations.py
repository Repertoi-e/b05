import numpy as np
import math
import matplotlib.pyplot as plt

R_AIR = 287      # J/kg/K
GRAV = 9.80665   # m/s^2
GAMMA = 1.4      # Ratio of specific heats for air

PI = math.pi
def SQRT(x): return np.sqrt(x)

T_sealevel_isa = 288.15   # K
T_tropopause_isa = 216.6  # K 

p_sealevel_isa = 101325   # Pa
p_tropopause_isa = 22600  # Pa

# Use T_offset=15 for e.g. a hot-day 
def isa(h_in_meters, T_offset=0):
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

def speed_of_sound(T):
    return math.sqrt(GAMMA * R_AIR * T)

#
# CONSTANTS FOR THE PLANE
#

bypass_ratio = 12  # From Class I weight estimation
num_engines = 2    # 

climb_rt_req = 12.7        # From REQ-OPER-0
altitude_climb_rt_req = 0  # From REQ-OPER-0

M_cruise = 0.82       #              From REQ-CRUISE-01
altitude_cr = 11887.2 # m, 39000 ft, From REQ-CRUISE-01

approach_speed = 80 # m/s   From REQ-LAND-05

T_cruise, _, _ = isa(altitude_cr)
V_cr = M_cruise * speed_of_sound(T_cruise)

specific_energy = 43 # kerosene, MJ/kg   From Class I weight estimation

jet_efficiency = V_cr / (22 * bypass_ratio**(-0.19)*10**(-6)) / (specific_energy*10**6)

aspect_ratio = 9.4   # From Class I weight estimation

c_f_equivalent = 0.0026    # From Class I weight estimation, estimated from Figure 6.3 in reader
S_wet_ratio    = 6.0       # From Class I weight estimation, estimated from Figure 6.3 in reader

c_d0_cruise = c_f_equivalent * S_wet_ratio # 6.15

oswald_phi = 0.97    # Figure 6.4 
oswald_psi = 0.0075  # Figure 6.4 
oswald_eff_cruise = 1/(PI * aspect_ratio * oswald_psi + 1/oswald_phi) # 6.17

flap_deflection_takeoff = 15    # in degrees  
flap_deflection_landing = 35    # in degrees  

def c_d0(gear_extended, flap_deflection):
    res = c_d0_cruise + 0.0013 * flap_deflection # 7.57
    if gear_extended: res += 0.02 # 7.58 Reader says to est. a value between 0.0100 and 0.0250 
    return res

def oswald_eff(flap_deflection):
    return oswald_eff_cruise + 0.0026 * flap_deflection # 7.56 wing-mounted engines 

max_LD = 1/2 * SQRT(PI * aspect_ratio * oswald_eff_cruise / c_d0_cruise)

cl_max_cruise  = 1.2  
cl_max_takeoff = 1.6  
cl_max_landing = 1.8  

l_to = 2790  # m, From REQ-APA-02
l_lf = 1856  # m, From REQ-APA-03

weight_takeoff = 200900.0 * GRAV  # From Class I weight estimation

mass_fr_cr = 0.95           
mass_fr_climb = mass_fr_cr
mass_fr_landing = 0.79      
mass_fr_app = mass_fr_landing

def min_speed():
    rho_app = p_sealevel_isa / (T_sealevel_isa * R_AIR)
    m_speed = (1/mass_fr_app) * (rho_app/2) * (approach_speed/1.23) ** 2 * cl_max_landing # Equation 7.6
    return m_speed

def landing_field(h_landing):
    cl_lf=0.45 # For jet aircraft

    _, _, rho_l  = isa(h_landing, T_offset=15)

    W_S = (1/mass_fr_landing) * (l_lf/cl_lf) * ((rho_l * cl_max_landing)/2) # Equation 7.15
    return W_S       

def cruise_speed(wing_load):
    T_cr, p_cr, rho_cr  = isa(altitude_cr)
    thr_lapse_cr = thrust_lapse(T_cr, M_cruise, p_cr)

    T_W = mass_fr_cr/thr_lapse_cr * (c_d0_cruise * 0.5 * rho_cr * V_cr ** 2)/(mass_fr_cr * wing_load)
    + (mass_fr_cr * wing_load)/(PI * aspect_ratio * oswald_eff_cruise * 0.5 * rho_cr * V_cr ** 2) # Equation 7.22
    return T_W

def climb_rate(wing_load):
    climb_lift_coeff = SQRT(c_d0_cruise * PI * oswald_eff_cruise * aspect_ratio)
    
    T_climb, p_climb, rho_climb = isa(altitude_climb_rt_req, T_offset=15)

    climb_speed = SQRT(wing_load * (2 / rho_climb) * (1 / climb_lift_coeff))
    M_climb_rt = climb_speed / speed_of_sound(T_climb)    

    climb_rate_TW = (mass_fr_climb/np.vectorize(thrust_lapse)(T_climb, M_climb_rt, p_climb)) * (SQRT((climb_rt_req ** 2/(mass_fr_climb * wing_load))
    * (rho_climb / 2) * SQRT(c_d0_cruise * PI * aspect_ratio * oswald_eff_cruise)) + 2
    * SQRT(c_d0_cruise/(PI * aspect_ratio * oswald_eff_cruise))) # Equation 7.43

    return climb_rate_TW

def climb_gradient(wing_load, altitude_climb_rate, zero_lift_drag_coeff, oswald_factor, climb_grad_percent, oei_condition):
    climb_lift_coeff = SQRT(zero_lift_drag_coeff * PI * oswald_factor * aspect_ratio)
    
    T_climb, p_climb, rho_climb = isa(altitude_climb_rate)
    
    climb_speed = SQRT(wing_load * (2 / rho_climb) * (1 / climb_lift_coeff))
    M_climb_rt = climb_speed / speed_of_sound(T_climb)    
    
    if oei_condition:
        oei_factor = num_engines / (num_engines - 1)
    else:
        oei_factor = 1

    climb_grad_TW = oei_factor * (mass_fr_climb/np.vectorize(thrust_lapse)(T_climb, M_climb_rt, p_climb)) * (climb_grad_percent/100 + 2 * SQRT(zero_lift_drag_coeff/(PI * oswald_factor * aspect_ratio)))
    return climb_grad_TW

def takeoff_field(wing_load, h_takeoff, oei_condition):
    h2 = 11 # m, this is for CS-25 planes

    T, p, rho = isa(h_takeoff, T_offset=15)
    cl_2 =(1/1.13)**2 * cl_max_takeoff

    v2 = SQRT(wing_load * (2 / rho) * (1 / cl_2))
    M = v2 / speed_of_sound(T)

    if oei_condition:
        oei_factor = num_engines / (num_engines - 1)
    else:
        oei_factor = 1

    k_t = 0.85 # Assumed for jet airplanes

    takeoff_field_TW = 1.15 * np.vectorize(thrust_lapse)(T, M, p) * SQRT((oei_factor * wing_load)/(l_to * k_t * rho *GRAV* math.pi * aspect_ratio * oswald_eff(flap_deflection_takeoff))) + oei_factor * ((4*h2)/(l_to)) # Equation 7.67 
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
        if  theta_t <= theta_t_break:
            return delta_t
        else:
            return delta_t * (1 - 2.1 * (theta_t - theta_t_break) / theta_t)
    elif bypass_ratio >= 5 and bypass_ratio < 15:

        if theta_t <= theta_t_break:
            return delta_t * (1 - (0.43 + 0.014 * bypass_ratio) * np.sqrt(M))
        else:
            return delta_t * (1 - (0.43 + 0.014 * bypass_ratio) * np.sqrt(M) - (3 * (theta_t - theta_t_break)) / (M + 1.5))

def calculate_takeoff_field_length_from_design_point(h_takeoff, oei_condition, wing_area, thrust_takeoff):
    h2 = 11 # m, this is for CS-25 planes

    T, p, rho = isa(h_takeoff, T_offset=15)

    cl_2 =(1/1.13)**2 * cl_max_takeoff

    print("cl_2", cl_2)
    Tv2 = thrust_takeoff * thrust_lapse(T, approach_speed / speed_of_sound(T), p)
    print("Tv2", Tv2)

    if oei_condition:
        oei_factor = (num_engines - 1) / num_engines
    else:
        oei_factor = 1

    print("oei_factor", oei_factor)

    c_d0_to = c_d0(gear_extended=True, flap_deflection=flap_deflection_takeoff)
    e_to = oswald_eff(flap_deflection_takeoff)

    print("oei_factor", oei_factor)

    k_t = 0.85 # Assumed for jet airplanes
    return weight_takeoff**2 / (rho * GRAV * wing_area * cl_2 * k_t * Tv2) + 2 * h2 * (oei_factor * Tv2 / weight_takeoff - (c_d0_to / cl_2 + cl_2 / (PI * aspect_ratio * e_to)))



#Wing Planform formulas

def wing_area(WingSpan, ChordLengthAtTip, ChordLengthAtRoot):
    S = WingSpan(ChordLengthAtTip + ChordLengthAtRoot)/2
    return S

def calculate_aspect_ratio(WingSpan, WingArea):
    AR = WingSpan^2/WingArea
    return AR

def calculate_taper_ratio(ChordLengthAtTip, ChordLengthAtRoot):
    TR = ChordLengthAtTip / ChordLengthAtRoot
    return TR

def angle_of_quarter_chord_line(MachAtCruise):
    Angle = math.acos(1.16/(MachAtCruise+0.5))
    return Angle

def angle_of_leading_edge(angle_of_quarter_chord_line, ChordLengthAtRoot, WingSpan, TaperRatio):
    Angle = math.tan(angle_of_quarter_chord_line) - ChordLengthAtRoot*(TaperRatio-1)/(2*WingSpan)

def wing_span(aspect_ratio, S):
    b = math.sqrt(aspect_ratio*360)
    return b 

#must be less than 80m


def taper_ratio(angle_of_quarter_chord_line):
    taper_ratio = -0.0083 * angle_of_quarter_chord_line + 0.4597
    return taper_ratio

def root_chord(taper_ratio, S, b):
    root_chord = (2*360)/((1+taper_ratio)*wing_span)
    return root_chord

def tip_chord(taper_ratio, root_chord):
    tip_chord = taper_ratio * root_chord
    return tip_chord

def mean_geometric_chord(root_chord, taper_ratio):
    MGC = (2/3)*root_chord*((1 + taper_ratio + (taper_ratio)^2)/( 1 + taper_ratio))
    return MGC

def YLEMGC(b, taper_ratio):
    YLEMGC = (wing_span/6) * ( 1 + 2*taper_ratio)/(1 + taper_ratio)
    return YLEMGC

def XLEMGC(YLEMGC, AngleOfLeadingEdge):
    XLEMGC = YLEMGC * math.tan(AngleOfLeadingEdge)
    return XLEMGC


def chordlength(LengthFromRoot, ChordLengthAtRoot, ChordLengthAtTip, WingSpan): 

    Percentage = LengthFromRoot/WingSpan
    
    Length = ChordLengthAtRoot -  Percentage*(ChordLengthAtRoot - ChordLengthAtTip)
    return Length

def mean_aerodynamic_chord(ChordLength, WingArea):
    pass