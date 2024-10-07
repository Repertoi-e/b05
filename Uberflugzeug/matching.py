from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QFileDialog, QTabWidget,
                             QVBoxLayout, QHBoxLayout, QWidget, QPlainTextEdit, QLineEdit, QLabel, QMessageBox, QToolBar, QFormLayout, QSizePolicy)

from isa import *
from utils import *

import pandas as pd

import matplotlib.ticker as plticker

import tab

class MatchingDiagram:
    def __init__(self, airplane):
        self.airplane = airplane
        for key, value in vars(airplane).items():
            setattr(self, key, value)

    def min_speed(self, cl_max_landing, approach_speed):
        rho_app = p_sealevel_isa / (T_sealevel_isa * R_AIR)
        m_speed = (1/self.mass_fr_app) * (rho_app/2) * (approach_speed/1.23) ** 2 * cl_max_landing # Equation 7.6
        return m_speed

    def landing_field(self, h_landing, cl_max_landing):
        cl_lf=0.45 # For jet aircraft

        _, _, rho_l  = isa(h_landing, T_offset=15)

        W_S = (1/self.mass_fr_land) * (self.l_lf/cl_lf) * ((rho_l * cl_max_landing)/2) # Equation 7.15
        return W_S       

    def cruise_speed(self, wing_load):
        T_cr, p_cr, rho_cr  = isa(self.altitude_cr)
        thr_lapse_cr = self.thrust_lapse(T_cr, self.M_cruise, p_cr)

        T_W = self.mass_fr_cr/thr_lapse_cr * ((self.c_d0_cruise * 0.5 * rho_cr * self.V_cr ** 2)/(self.mass_fr_cr * wing_load)
        + (self.mass_fr_cr * wing_load)/(PI * self.aspect_ratio * self.oswald_eff_cruise * 0.5 * rho_cr * self.V_cr ** 2)) # Equation 7.22
        return T_W

    def climb_rate(self, wing_load):
        climb_lift_coeff = SQRT(self.c_d0_cruise * PI * self.oswald_eff_cruise * self.aspect_ratio)
        
        T_climb, p_climb, rho_climb = isa(self.altitude_climb_rt_req, T_offset=15)

        climb_speed = SQRT(wing_load * (2 / rho_climb) * (1 / climb_lift_coeff))
        M_climb_rt = climb_speed / speed_of_sound(T_climb)    

        climb_rate_TW = (self.mass_fr_climb/np.vectorize(self.thrust_lapse)(T_climb, M_climb_rt, p_climb)) * (SQRT((self.climb_rt_req ** 2/(self.mass_fr_climb * wing_load))
        * (rho_climb / 2) * SQRT(self.c_d0_cruise * PI * self.aspect_ratio * self.oswald_eff_cruise)) + 2
        * SQRT(self.c_d0_cruise/(PI * self.aspect_ratio * self.oswald_eff_cruise))) # Equation 7.43

        return climb_rate_TW

    def climb_gradient(self, wing_load, altitude_climb_rate, zero_lift_drag_coeff, oswald_factor, climb_grad_percent, oei_condition):
        climb_lift_coeff = SQRT(zero_lift_drag_coeff * PI * oswald_factor * self.aspect_ratio)
        
        T_climb, p_climb, rho_climb = isa(altitude_climb_rate)
        
        climb_speed = SQRT(wing_load * (2 / rho_climb) * (1 / climb_lift_coeff))
        M_climb_rt = climb_speed / speed_of_sound(T_climb)    
        
        if oei_condition:
            oei_factor = self.num_engines / (self.num_engines - 1)
        else:
            oei_factor = 1

        climb_grad_TW = oei_factor * (self.mass_fr_climb/np.vectorize(self.thrust_lapse)(T_climb, M_climb_rt, p_climb)) * (climb_grad_percent/100 + 2 * SQRT(zero_lift_drag_coeff/(PI * oswald_factor * self.aspect_ratio)))
        return climb_grad_TW

    def takeoff_field(self, wing_load, h_takeoff, oei_condition, cl_max):
        h2 = 11 # m, this is for CS-25 planes

        T, p, rho = isa(h_takeoff, T_offset=15)
        cl_2 =(1/1.13)**2 * cl_max
        
        v2 = SQRT(wing_load * (2 / rho) * (1 / cl_2))
        M = v2 / speed_of_sound(T)

        if oei_condition:
            oei_factor = self.num_engines / (self.num_engines - 1)
        else:
            oei_factor = 1

        k_t = 0.85 # Assumed for jet airplanes

        takeoff_field_TW = 1.15 * np.vectorize(self.thrust_lapse)(T, M, p) * SQRT((oei_factor * wing_load)/(self.l_to * k_t * rho *GRAV* math.pi * self.aspect_ratio * self.airplane.oswald_eff(self.flap_deflection_takeoff))) + oei_factor * ((4*h2)/(self.l_to)) # Equation 7.67 
        return takeoff_field_TW   


    def thrust_lapse(self, T, M, p): # 7.29 - 7.32 
        delta_t = total_pressure(M, p) / p_sealevel_isa
        theta_t = total_temperature(M, T) / T_sealevel_isa
        
        theta_t_break = 1.07   # According to Ref. [4], modern engines have values of θt break ranging between 1.06 and 1.08.

        if not (self.bypass_ratio > 0 and self.bypass_ratio < 15):
            raise AirplaneValueError("Bypass ratio must be between 0 and 15")

        if self.bypass_ratio < 5:
            if theta_t <= theta_t_break:
                return delta_t
            else:
                return delta_t * (1 - 2.1 * (theta_t - theta_t_break) / theta_t)
        elif self.bypass_ratio >= 5 and self.bypass_ratio < 15:
            if theta_t <= theta_t_break:
                return delta_t * (1 - (0.43 + 0.014 * self.bypass_ratio) * np.sqrt(M))
            else:
                return delta_t * (1 - (0.43 + 0.014 * self.bypass_ratio) * np.sqrt(M) - (3 * (theta_t - theta_t_break)) / (M + 1.5))
            
    def calculate_takeoff_field_length_from_design_point(self, h_takeoff, oei_condition, wing_area, thrust_takeoff, approach_speed, cl_max_takeoff):
        h2 = 11 # m, this is for CS-25 planes

        T, p, rho = isa(h_takeoff, T_offset=15)

        cl_2 =(1/1.13)**2 * cl_max_takeoff

        Tv2 = thrust_takeoff * self.thrust_lapse(T, mach(approach_speed, T), p)

        if oei_condition:
            oei_factor = (self.num_engines - 1) / self.num_engines
        else:
            oei_factor = 1

        c_d0_to = self.airplane.c_d0(gear_extended=True, flap_deflection=self.flap_deflection_takeoff)
        e_to = self.airplane.oswald_eff(self.flap_deflection_takeoff)

        k_t = 0.85 # Assumed for jet airplanes
        return self.MTOW**2 / (rho * GRAV * wing_area * cl_2 * k_t * Tv2) + 2 * h2 * (oei_factor * Tv2 / self.MTOW - (c_d0_to / cl_2 + cl_2 / (PI * self.aspect_ratio * e_to)))


class Tab(tab.Tab):
    def __init__(self, uber):
        super().__init__(uber, "Matching Diagram")
        self.second_pass = None

    def init_ui_for_airplane(self):
        super().init_ui_for_airplane()

        self.ax2 = self.ax.twinx()

        self.create_input_box("Approach Speed", "approach_speed")
        self.create_input_box("Number of engines", "num_engines")

        self.create_input_box("Mass fraction cruise", "mass_fr_cr")
        self.create_input_box("Mass fraction climb", "mass_fr_climb")
        self.create_input_box("Mass fraction approach", "mass_fr_app")
        self.create_input_box("Mass fraction landing", "mass_fr_land")

        self.separator(title="TLAR")
        
        self.create_input_box("Climb rate [m/s]", "climb_rt_req")
        self.create_input_box("Altitude for climb rate", "altitude_climb_rt_req")

        self.create_input_box("Take-off field length [m]", "l_to")
        self.create_input_box("Landing field length [m]", "l_lf")


def recalculate(tab):
    a = tab.uber.airplane

    d = MatchingDiagram(a)

    v_stall = a.approach_speed / 1.23
    rho = isa(0, T_offset=15)[2]

    a.q_stall = 1/2 * rho * v_stall**2

    ax, ax2 = tab.ax, tab.ax2
    ax.clear()
    ax2.clear()

    a.WS = 5900 # Dummy for now, calculated later!
    
    W_S = np.arange(0, 10000, step=100)

    a.CL_max_takeoff = 1.7 # This will be recalced once the design point WS is chosen later!

    W_S_takeoff_field = d.takeoff_field(W_S[4:], h_takeoff=1500, oei_condition=True, cl_max=a.CL_max_takeoff) # This is not critical! 

    W_S_cruise_speed = d.cruise_speed(W_S[4:])  # First few elements are divide by 0 so skip them 
    W_S_climb_rate = d.climb_rate(W_S[4:])

    W_S_climb_gradient_25_119  = d.climb_gradient(W_S, altitude_climb_rate=0, zero_lift_drag_coeff=a.c_d0(gear_extended=True, flap_deflection=a.flap_deflection_landing), oswald_factor=a.oswald_eff(a.flap_deflection_landing), climb_grad_percent=3.2, oei_condition=False) 
    W_S_climb_gradient_25_121a = d.climb_gradient(W_S, altitude_climb_rate=0, zero_lift_drag_coeff=a.c_d0(gear_extended=True, flap_deflection=a.flap_deflection_takeoff), oswald_factor=a.oswald_eff(a.flap_deflection_takeoff), climb_grad_percent=0, oei_condition=True) 
    W_S_climb_gradient_25_121b = d.climb_gradient(W_S, altitude_climb_rate=0, zero_lift_drag_coeff=a.c_d0(gear_extended=False, flap_deflection=a.flap_deflection_takeoff), oswald_factor=a.oswald_eff(a.flap_deflection_takeoff), climb_grad_percent=2.4, oei_condition=True)
    W_S_climb_gradient_25_121c = d.climb_gradient(W_S, altitude_climb_rate=0, zero_lift_drag_coeff=a.c_d0_cruise, oswald_factor=a.oswald_eff_cruise, climb_grad_percent=1.2, oei_condition=True)
    W_S_climb_gradient_25_121d = d.climb_gradient(W_S, altitude_climb_rate=0, zero_lift_drag_coeff=a.c_d0(gear_extended=False, flap_deflection=a.flap_deflection_landing), oswald_factor=a.oswald_eff(a.flap_deflection_landing), climb_grad_percent=2.1, oei_condition=True)

    # Get our reference aircraft database to plot them
    sheets_url = "https://docs.google.com/spreadsheets/d/1HbuwOGN-pKfY95QqNByAydbRJVicJ3oztc3d0fufQT8/export?format=csv&gid=0#gid=0"

    df = pd.read_csv(sheets_url, header=[0,1], nrows=8)

    # This will be recalced once the design point WS is chosen later!
    cl_max_landing = 1/a.q_stall * a.WS

    W_S_min_speed = d.min_speed(cl_max_landing=cl_max_landing, approach_speed=a.approach_speed)

    W_S_landing_field_length = d.landing_field(h_landing=0, cl_max_landing=cl_max_landing)

    ax.vlines([W_S_landing_field_length], ymin=0, ymax=1, color="orange", label="REQ-LAND-01 Landing Field Length")
    ax.vlines([W_S_min_speed], ymin=0, ymax=1, color="b", label="REQ-LAND-05 Minimum Speed")

    cl_max = 1/a.q_stall * W_S
    ax2.plot(W_S, cl_max, label=f"$C_{{L_{{max}}}}$, $V_{{app}} = {a.approach_speed:.2f}\ m/s$", linestyle="dotted")

    ax.plot(W_S[4:], W_S_takeoff_field, label="REQ-OFF-01 Take-off Field Length", color="green")

    ax.plot(W_S, W_S_climb_gradient_25_121a, label="REQ-OFF-02 Climb Gradient Take-off", color="lightblue")
    ax.plot(W_S, W_S_climb_gradient_25_121b, label="REQ-OFF-03 Climb Gradient Take-off", color="orange")
    ax.plot(W_S, W_S_climb_gradient_25_121c, label="REQ-OFF-04 Climb Gradient Take-off", color="lightgreen")

    ax.plot(W_S[4:], W_S_climb_rate, color="yellow", label="REQ-OPER-01 Climb Rate")

    ax.plot(W_S[4:], W_S_cruise_speed, color="purple", label="REQ-CRU-01 Cruise Speed")

    ax.plot(W_S, W_S_climb_gradient_25_119, color="red", label="REQ-LAND-02 Climb Gradient Landing")

    ax.plot(W_S, W_S_climb_gradient_25_121d, color="cyan", label="REQ-LAND-03 Climb Gradient Approach")

    ax.scatter(df["W/S"].values, df["T_TO/W_TO"].values, c="blue", s=20, label="Reference Aircraft")

    # Choose design point as min W_S_cruise_speed, since that's the limiting case
    index_approximate_min = np.argmin(W_S_cruise_speed)
    min_ws = W_S[index_approximate_min - 10]
    max_ws = W_S[index_approximate_min + 10]

    min_TW, min_WS = None, None
    for ws_search in np.arange(min_ws, max_ws, 10):
        if min_TW is None or min_TW > d.cruise_speed(ws_search):
            min_TW = d.cruise_speed(ws_search)
            min_WS = ws_search

    if min_WS > W_S_landing_field_length:
        min_WS = W_S_landing_field_length

    if min_WS > W_S_min_speed:
        min_WS = W_S_min_speed

    min_TW = d.cruise_speed(min_WS)

    a.WS = min_WS
    a.TW = min_TW
    a.CL_max = 1/a.q_stall * a.WS

    # Recalc with proper WS
    if tab.second_pass is None:
        tab.second_pass = True
        recalculate(tab)
        return

    tab.second_pass = None

    a.S = a.MTOW / a.WS
    a.T = a.TW * a.MTOW 

    print(f"Design thrust to weight: {a.TW:.4f} N/N")
    print(f"Thrust required: {a.T/1000:.0f} kN")
    print()
    print(f"Design wing loading: {a.WS:.2f} N/m^2")
    print(f"Design C_L_max: {a.CL_max:.3f}")

    to_length = None
    a.CL_max_takeoff = 0.2
    while to_length is None or to_length > a.l_to:
        to_length = d.calculate_takeoff_field_length_from_design_point(h_takeoff=1500, oei_condition=True, wing_area=a.S, thrust_takeoff=a.T, approach_speed=a.approach_speed, cl_max_takeoff=a.CL_max_takeoff)
        a.CL_max_takeoff += 0.1

    print(f"Take-off C_L_max required: {a.CL_max_takeoff:.2f}")

    ax.scatter([a.WS], [a.TW], c="red", s=60, label="Design Point")

    ax.set_ylim((0, 0.6))
    ax.set_ylabel("$T/W\ [N/N]$")
    ax2.set_ylabel("$C_{L_{max}}$")
    ax.set_xlabel("$W/S\ [N/m^2]$")
    
    plt.xlim(1200, 10000)
    ax2.set_ylim((0, 4.85))

    ax.xaxis.set_major_locator(plticker.MultipleLocator(1000)) 
    ax.yaxis.set_major_locator(plticker.MultipleLocator(0.1)) 

    ax.grid()
    ax2.grid()

    # ax.legend() Too busy!
    ax2.legend()

    # ax.savefig('matching_diagram.svg', format='svg')

    tab.canvas.draw()

