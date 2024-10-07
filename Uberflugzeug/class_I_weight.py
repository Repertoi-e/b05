from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QFileDialog, QTabWidget,
                             QVBoxLayout, QHBoxLayout, QWidget, QPlainTextEdit, QLineEdit, QLabel, QMessageBox, QToolBar, QFormLayout, QSizePolicy)

import numpy as np
from utils import *

import tab

class Tab(tab.Tab):
    def __init__(self, uber):
        super().__init__(uber, "Class I Weight Estimation")

    def init_ui_for_airplane(self):
        super().init_ui_for_airplane()

        self.create_input_box("Bypass Ratio", "bypass_ratio")
        self.create_input_box("Fuel specific energy [MJ/kg]", "fuel_specific_energy")

        self.create_input_box("Fuel contingency for range", "f_con")

        self.create_input_box("Divert range [km]", "R_div")
        self.create_input_box("Loiter time [s]", "t_loiter")

        self.create_input_box("Mass fraction OE to MTO", "m_oe_to_mto")

        self.separator(title="TLAR")

        self.create_input_box("Mach cruise", "M_cruise")
        self.create_input_box("Cruise altitude [m]", "altitude_cr")

        self.create_input_box("Max payload req. [kg]", "m_max_pl")

        self.create_input_box("Design payload req. [kg]", "m_pl_des")
        self.create_input_box("Design range req. [km]", "R_des")

        self.create_input_box("Ferry range req. [km]", "R_ferry_req")

        self.separator(title="Additional TLAR @ MTOW req")

        self.create_input_box("Payload [kg]", "m_pl_max_mtow_and_full_fuel")
        self.create_input_box("Range [km]", "R_max_mtow_and_full_fuel")

        self.separator(title="Constants for drag polar estimation")

        self.create_input_box("Wet area ratio", "S_wet_ratio")
        self.create_input_box("c_f_equivalent", "c_f_equivalent")
        self.create_input_box("oswald_phi", "oswald_phi")
        self.create_input_box("oswald_psi", "oswald_psi")

def recalculate(tab):
    a = tab.uber.airplane

    a.added_fuel = 0
    a.max_pl_for_fuel_exchange = a.m_max_pl

    # True air speed at cruise
    a.T_cruise, _, a.rho_cr = isa(a.altitude_cr)
    a.V_cr = a.M_cruise * speed_of_sound(a.T_cruise)

    a.TSFC = 22 * a.bypass_ratio**(-0.19)*10**(-6)

    a.jet_efficiency = a.V_cr / (a.TSFC) / (a.fuel_specific_energy*10**6)

    # Rough cd0 estimation
    a.c_d0_cruise = a.c_f_equivalent * a.S_wet_ratio # 6.15
    a.oswald_eff_cruise = 1/(PI * a.aspect_ratio * a.oswald_psi + 1/a.oswald_phi) # 6.17

    a.max_LD = 1/2 * SQRT(PI * a.aspect_ratio * a.oswald_eff_cruise / a.c_d0_cruise)

    while True:
        r_lost = 1/0.7 * a.max_LD * (a.altitude_cr + a.V_cr**2/(2 * GRAV)) / 1000
        
        def r_eq(range):
            return a.V_cr * a.t_loiter/1000 + 1.2 * a.R_div + (1+ a.f_con) * (range + r_lost)

        r_eq_des = r_eq(a.R_des)

        def fuel_mass_fr(r_eq):
            return 1 - np.exp(-r_eq * 1000 / (a.jet_efficiency * (a.fuel_specific_energy * 1e6 / GRAV) * a.max_LD))

        beta_f_des = fuel_mass_fr(r_eq_des)
        m_mto_des = a.m_pl_des / (1 - a.m_oe_to_mto - beta_f_des) + a.added_fuel
        m_oe_des = a.m_oe_to_mto * (m_mto_des - a.added_fuel)
        m_f_des = beta_f_des * m_mto_des 

        r_aux = r_eq_des - a.R_des
        
        def range(m_oe, m_pl, m_f):
            return 1/1000 * a.jet_efficiency * a.max_LD * (a.fuel_specific_energy * 1e6 / GRAV) * np.log((m_oe + m_pl + m_f)/(m_oe + m_pl)) - r_aux

        m_f_max_pl = m_mto_des - m_oe_des - a.m_max_pl
        r_max_pl = range(m_oe_des, a.m_max_pl, m_f_max_pl)

        m_f_ferry = m_mto_des - m_oe_des - (a.m_max_pl - a.max_pl_for_fuel_exchange)
        r_ferry = range(m_oe_des, 0, m_f_ferry)

        r_des = round(range(m_oe_des, a.m_pl_des, m_f_des), 0)

        if r_ferry > a.R_ferry_req:
            # print("Ferry range too big. Removing 10 kg from payload fuel exchange and trying again...")
            a.max_pl_for_fuel_exchange -= 10
            continue

        m_f_req_req = m_f_des - (a.m_pl_des - a.m_pl_max_mtow_and_full_fuel)
        req_req_R = range(m_oe_des, a.m_pl_max_mtow_and_full_fuel, m_f_req_req)
        if req_req_R < a.R_max_mtow_and_full_fuel:
            # print("Range not enough for max MTOW and full fuel. Adding 10 kg of fuel and trying again...", req_req_R, a.R_max_mtow_and_full_fuel)
            a.added_fuel += 10
            continue
        break

    a.beta_f_des = beta_f_des
    a.m_fuel_design = m_f_des
    a.MTOW = m_mto_des * GRAV
    a.MOE = m_oe_des * GRAV

    a.m_f_ferry = m_f_ferry
    a.R_des = r_des
    a.R_ferry = r_ferry
    a.R_max_pl = r_max_pl
        
    print(f".... Total added fuel to reach critical TLAR {a.added_fuel:.0f} kg")
    print(f"Design fuel mass fraction: {a.beta_f_des:.4f}")
    print()
    print(f"MTOW: {a.MTOW / GRAV:.0f} kg")
    print(f"MOE: {a.MOE / GRAV:.0f} kg")
    print()
    print(f"Max payload: {a.m_max_pl:.0f} kg")
    print(f"Max range @ max payload: {a.R_max_pl:.0f} km")
    print()
    print(f"Design payload: {a.m_pl_des:.0f} kg")
    print(f"Range @ design payload: {a.R_des:.0f} km")
    print(f"m_fuel, design: {a.m_fuel_design:.0f} kg")
    print()
    print(f"Ferry range: {a.R_ferry:.0f} km")
    print(f"m_fuel, max: {a.m_fuel_design + a.max_pl_for_fuel_exchange:.0f} kg")

    ax = tab.ax
    ax.clear()

    ax.plot([0, r_max_pl], [a.m_max_pl, a.m_max_pl], color="red", label="TLAR Max Payload")
    ax.hlines([a.m_pl_des], [0], [r_des], linestyles="dotted", label="TLAR Design Payload", color="red")
    ax.scatter([a.R_max_mtow_and_full_fuel, a.R_ferry_req], [a.m_pl_max_mtow_and_full_fuel, 0], color="red", label="TLAR Missions")

    ax.plot([r_max_pl, r_des, r_ferry], [a.m_max_pl, a.m_pl_des, 0], label="Trade-off Payload vs. Fuel")
    
    ax.grid()
    ax.set_xlabel("Range [km]")
    ax.legend()
    ax.set_ylabel("Payload [kg]")

    tab.canvas.draw()

