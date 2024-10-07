from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QFileDialog, QTabWidget, QComboBox,
                             QVBoxLayout, QHBoxLayout, QWidget, QPlainTextEdit, QLineEdit, QLabel, QMessageBox, QToolBar, QFormLayout, QSizePolicy)

import math
from utils import *

from scipy import integrate

import tab

def HLD_surface_available(le_angle, te_angle, root_chord, start_x, end_x, wing_surface):
    S_start = ((math.tan(te_angle*math.pi/180)-math.tan(le_angle*math.pi/180))/2)*(start_x**2) + root_chord*start_x
    S_end = ((math.tan(te_angle*math.pi/180)-math.tan(le_angle*math.pi/180))/2)*(end_x**2) + root_chord*end_x
    S_max = 2*(S_end - S_start)  # FULL WING !
    Swf_over_S = S_max/wing_surface
    return S_max, Swf_over_S

def spanwise_positioning_HLD(le_angle, te_angle, root_chord, start_x, wing_surface, s_wf_over_s):
    a = ((math.tan(te_angle*math.pi/180)-math.tan(le_angle*math.pi/180))/2)
    b = root_chord
    c = -(((math.tan(te_angle*math.pi/180)-math.tan(le_angle*math.pi/180))/2)*(start_x**2) + root_chord*start_x) - 0.5*wing_surface*s_wf_over_s
    Delta = b**2 - 4*a*c
    x_1 = (-b+math.sqrt(Delta))/(2*a)
    x_2 = (-b-math.sqrt(Delta))/(2*a)
    return x_1, x_2

def Delta_CL_Max_obtained(Swf_S_te, Swf_S_le, deltaCl_te, deltaCl_le, hinge_sweep_angle_TE, hinge_sweep_angle_LE):
    'Function to calculate the delta Cl obtained by the combination of LE and TE HLDs'
    delta_CL_Max_TE = 0.9 * deltaCl_te * Swf_S_te*math.cos(hinge_sweep_angle_TE * (math.pi / 180))
    delta_CL_Max_LE = 0.9 * deltaCl_le * Swf_S_le*math.cos(hinge_sweep_angle_LE * (math.pi / 180))
    delta_CL_Max = delta_CL_Max_LE + delta_CL_Max_TE
    return delta_CL_Max, delta_CL_Max_TE, delta_CL_Max_LE

def check_validity_engine(ending_x_engine, maximum_x):
    if ending_x_engine <= maximum_x:
        print("→ all checks out")
    else:
        print("!!! HLDs do not fit with engine, fuselage and ailerons !!! ")

class Tab(tab.Tab):
    def __init__(self, uber):
        super().__init__(uber, "Ailerons && HLD")

    def init_ui_for_airplane(self):
        super().init_ui_for_airplane()

        self.separator(title="Ailerons")

        self.create_input_box("Aileron max deflection [deg]", "aileron_max_deflection")

        self.create_input_box("Aileron start fraction b/2", "aileron_start")
        self.create_input_box("Aileron end fraction b/2", "aileron_end")
        self.create_input_box("cf/c aileron", "aileron_cf_ratio")

        self.create_input_box("Aileron effectiveness", "aileron_effectiveness")
        self.create_input_box("Roll required in 1.5s [deg]", "bank_angle_required")

        self.separator(title="HLD")

        self.create_input_box("Flap deflection take-off [deg]", "flap_deflection_takeoff")
        self.create_input_box("Flap deflection landing [deg]", "flap_deflection_landing")

        self.create_input_box("Take-off angle of attack [deg]", "take_off_angle_of_attack")
        self.create_input_box("Landing angle of attack [deg]", "landing_angle_of_attack")

        self.create_input_box("cf/c TE", "cf_c_TE")

        self.create_input_box("Swf/S TE", "SwfS_TE")
        self.create_input_box("Swf/S LE", "SwfS_LE")

        self.separator(title="Subsystems and margins")

        self.create_input_box("Y fuselage [m]", "Y_fuselage")
        self.create_input_box("Margin fuselage [m]", "margin_fuselage")

        self.create_input_box("Margin ailerons [m]", "margin_ailerons")

        self.create_input_box("Margin wing tip [m]", "margin_wing_tip")

        self.create_input_box("Y engine [m]", "Y_engine")
        self.create_input_box("dy nacelle [m]", "dy_nacelle")

def recalculate(tab):
    a = tab.uber.airplane

    # WING DATA
    LE_sweep = math.degrees(a.sweep_angle_le) # [deg]
    TE_sweep = math.degrees(a.sweep_angle_at(1))  # [deg]

    r_c = a.c_r  # [m]

    print("        ### AILERONS ###")

    aileron_b_0 = a.aileron_start * a.b / 2
    aileron_b_1 = a.aileron_end * a.b / 2

    C_LA = ((2 * a.clalpha_airfoil * a.aileron_effectiveness) / (a.S * a.b)) * integrate.quad(lambda y: y * a.chord(y), aileron_b_0, aileron_b_1)[0]
    C_LP = -(4 * (a.clalpha_airfoil + a.cd0_airfoil_approach)) / (a.S * a.b**2) * integrate.quad(lambda y: y**2 * a.chord(y), 0, a.b/2)[0]
    P_val = -(C_LA / C_LP) * math.radians(a.aileron_max_deflection) * ((2 * a.approach_speed) / a.b)
    roll_t = math.radians(a.bank_angle_required) / P_val

    print(f"Aileron C_LA: {C_LA:4f}")
    print(f"Aileron C_LP: {C_LP:4f}")
    print(f"Aileron P: {P_val:.2f}")
    print(f"Aileron time to roll: {roll_t:.2f} s")

    if roll_t > 1.5:
        print("!!! Roll requirement not met in time !!!")
    print()

    # SUBSYSTEMS AND MARGINS
    Y_start = round(a.Y_fuselage + a.margin_fuselage, 3)  # [m]
    Y_end_max_TE = round(aileron_b_0 - a.margin_ailerons, 3)  # [m]
    Y_end_max_LE = round(0.5*a.b - a.margin_wing_tip, 3)  # [m]

    # SURFACE CALCULATIONS
    # TRAILING EDGE
    max_HLD_Surface_TE, max_S_wf_over_S_TE = HLD_surface_available(LE_sweep, TE_sweep, r_c, Y_start, Y_end_max_TE, a.S)
    # LEADING EDGE
    max_HLD_Surface_LE, max_S_wf_over_S_LE = HLD_surface_available(LE_sweep, TE_sweep, r_c, Y_start, Y_end_max_LE, a.S)
    
    Fowler_TE_HLD_spanwise_pos = spanwise_positioning_HLD(LE_sweep, TE_sweep, r_c, Y_start, a.S, a.SwfS_TE)
    Fowler_LE_HLD_spanwise_pos = spanwise_positioning_HLD(LE_sweep, TE_sweep, r_c, Y_start, a.S, a.SwfS_LE)
    Y_LE_end = Fowler_LE_HLD_spanwise_pos[0]  # end y span-wise position half wing
    Y_TE_end = Fowler_TE_HLD_spanwise_pos[0]  # end y span-wise position half wing

    # Surface covered to Nacelle
    surface_to_nacelle = HLD_surface_available(LE_sweep, TE_sweep, r_c, Y_start, a.Y_engine-0.5*a.dy_nacelle, a.S)
    surface_fuselage_to_nacelle_ratio = surface_to_nacelle[1]
    # Left to cover (Nacelle to y_end)
    # LEADING EDGE
    Fowler_LE_HLD_spanwise_pos_w_engine = spanwise_positioning_HLD(LE_sweep, TE_sweep, r_c, a.Y_engine+0.5*a.dy_nacelle, a.S, a.SwfS_LE-surface_fuselage_to_nacelle_ratio)
    Y_LE_end_w_engine = Fowler_LE_HLD_spanwise_pos_w_engine[0]
    # TRAILING EDGE
    Fowler_TE_HLD_spanwise_pos_w_engine = spanwise_positioning_HLD(LE_sweep, TE_sweep, r_c, a.Y_engine+0.5*a.dy_nacelle, a.S, a.SwfS_TE-surface_fuselage_to_nacelle_ratio)
    Y_TE_end_w_engine = Fowler_TE_HLD_spanwise_pos_w_engine[0]

    print()
    print("        ### INITIAL HLD CALCULATIONS ###")
    print("S available for TE HLD's is:", round(max_HLD_Surface_TE, 2), "m^2")
    print("Max TE HLD S_wf/S is:       ", round(max_S_wf_over_S_TE, 3))
    print("S available for LE HLD's is:", round(max_HLD_Surface_LE, 2), "m^2")
    print("Max LE HLD S_wf/S is:       ", round(max_S_wf_over_S_LE, 3))
    print()
    print("LE HLD for the Fowler flap from y =", Y_start, "m to", round(Fowler_LE_HLD_spanwise_pos[0], 3), "m")
    # print("For the LEADING edge we DISREGARD the solution y=", round(Fowler_LE_HLD_spanwise_pos[1], 2), "m")
    print("TE HLD for the Fowler flap from y =", Y_start, "m to", round(Fowler_TE_HLD_spanwise_pos[0], 3), "m")
    # print("For the TRAILING edge we DISREGARD the solution y=", round(Fowler_TE_HLD_spanwise_pos[1], 2), "m")
    print()
    print()

    if (max_S_wf_over_S_TE - a.SwfS_TE) > 0.1:
        print("!!! Swf/S TE Fowler can be reduced !!! ")

    if (max_S_wf_over_S_LE - a.SwfS_LE) > 0.1:
        print("!!! Swf/S LE Fowler can be reduced !!! ")

    TE_y0 = Y_start
    TE_y1 = a.Y_engine-0.5*a.dy_nacelle
    TE_y2 = a.Y_engine+0.5*a.dy_nacelle
    TE_y3 = round(Y_TE_end_w_engine, 2)

    LE_y0 = Y_start
    LE_y1 = a.Y_engine-0.5*a.dy_nacelle
    LE_y2 = a.Y_engine+0.5*a.dy_nacelle
    LE_y3 = round(Y_LE_end_w_engine, 2)

    print("        ### HLD WITH ENGINE ###")
    print("TE HLD with engine y = [", TE_y0, ",", TE_y1, "] ∪ [", TE_y2, ",", TE_y3, "] m")
    # print("For the TRAILING edge we DISREGARD the solution y=", Fowler_TE_HLD_spanwise_pos_w_engine[1], "m")
    check_validity_engine(Y_TE_end_w_engine, Y_end_max_TE)
    print("LE HLD with engine y = [", LE_y0, ",", LE_y1, "] ∪ [", LE_y2, ",", LE_y3, "] m")
    check_validity_engine(Y_LE_end_w_engine, Y_end_max_LE)
    # print("For the TRAILING edge we DISREGARD the solution y=", Fowler_TE_HLD_spanwise_pos_w_engine[1], "m")
    print()

    CL_takeoff_clean = a.CLalpha_wing * math.radians(a.take_off_angle_of_attack)
    CL_landing_clean = a.CLalpha_wing * math.radians(a.landing_angle_of_attack)
    
    delta_CL_take_off_required = a.CL_max_takeoff - CL_takeoff_clean
    delta_CL_approach_required = a.CL_max - CL_landing_clean

    print(f'Required delta CL at take-off: {delta_CL_take_off_required:.3f}')
    print(f'Required delta CL at landing: {delta_CL_approach_required:.3f}')

    to_ratio = a.flap_deflection_takeoff/a.flap_deflection_landing

    # Values from ADSEE formula sheet
    def interpolate(x, x1=10, y1=0.45, x2=40, y2=0.6):
        return y1 + (x - x1) / (x2 - x1) * (y2 - y1)

    delta_c_cf_to = interpolate(a.flap_deflection_takeoff)
    delta_c_cf_land = interpolate(a.flap_deflection_landing)

    c_prime_c_to = (1/a.cf_c_TE + delta_c_cf_to) * a.cf_c_TE
    c_prime_c_land = (1/a.cf_c_TE + delta_c_cf_land) * a.cf_c_TE

    # Parameters to calculate total increment in Cl by the HLDs of design option 1 (TE: fowler flaps + LE: kruger flap)
    deltaCl_airfoil_TE_design_option1_landing = 1.3 * (c_prime_c_to) # values taken from type of TE HLDs
    deltaCl_airfoil_TE_design_option1_take_off = 1.3 * (c_prime_c_land) 
    deltaCl_airfoil_LE_design_option1 = 0.3 # value taken from type of LE HLDs

    # Parameters to calculate total increment in Cl by the HLDs of design option 2 (TE: single-slotted flaps + LE: kruger flap)
    deltaCl_airfoil_TE_design_option2_landing = 1.3 # value taken from type of TE HLDs
    deltaCl_airfoil_TE_design_option2_take_off = deltaCl_airfoil_TE_design_option2_landing * to_ratio # only a fr. of landing values, as TE HLD's are not fully deployed
    deltaCl_airfoil_LE_design_option2 = 0.3 # value taken from type of LE HLDs

    # Parameters to calculate size of TE HLDs
    # ratio_c_option1 = deltaCl_airfoil_TE_design_option1_landing / 1.3
    # c_prime_option1 = ratio_c_option1 * a.mac #[m]
    
    # cf_TE = 0.38 * a.mac # [m]
    # ratio_deltac_cf_option2 = 0.31
    # delta_c_option2 = ratio_deltac_cf_option2 * cf_TE
    # c_prime_option2 = delta_c_option2 + a.mac

    # Calculate area of TE HLDs
    # ratio_S_prime_S_option1 = 1 + max_S_wf_over_S_TE * (c_prime_option1 / a.mac - 1)
    # ratio_S_prime_S_option2 = 1 + max_S_wf_over_S_TE * (c_prime_option2 / a.mac - 1)
    # area_TE_HLD_option1 = a.S * (ratio_S_prime_S_option1 - 1) / 2 #Area of each Fowler flap
    # area_TE_HLD_option2 = a.S * (ratio_S_prime_S_option2 - 1) / 2 #Area of each Single-Sloted flap

    delta_CL_max_obtained_option1_landing, _, _ = Delta_CL_Max_obtained(Swf_S_te=a.SwfS_TE, Swf_S_le=a.SwfS_LE, deltaCl_te=deltaCl_airfoil_TE_design_option1_landing, deltaCl_le=deltaCl_airfoil_LE_design_option1, hinge_sweep_angle_TE=TE_sweep, hinge_sweep_angle_LE=LE_sweep)
    delta_CL_max_obtained_option1_take_off, _, _ = Delta_CL_Max_obtained(Swf_S_te=a.SwfS_TE, Swf_S_le=a.SwfS_LE, deltaCl_te=deltaCl_airfoil_TE_design_option1_take_off, deltaCl_le=deltaCl_airfoil_LE_design_option1, hinge_sweep_angle_TE=TE_sweep, hinge_sweep_angle_LE=LE_sweep)

    delta_CL_max_obtained_option2_landing, _, _ = Delta_CL_Max_obtained(Swf_S_te=a.SwfS_TE, Swf_S_le=a.SwfS_LE, deltaCl_te=deltaCl_airfoil_TE_design_option2_landing, deltaCl_le=deltaCl_airfoil_LE_design_option2, hinge_sweep_angle_TE=TE_sweep, hinge_sweep_angle_LE=LE_sweep)
    delta_CL_max_obtained_option2_take_off, _, _ = Delta_CL_Max_obtained(Swf_S_te=a.SwfS_TE, Swf_S_le=a.SwfS_LE, deltaCl_te=deltaCl_airfoil_TE_design_option2_take_off, deltaCl_le=deltaCl_airfoil_LE_design_option2, hinge_sweep_angle_TE=TE_sweep, hinge_sweep_angle_LE=LE_sweep)

    print(f'delta CL, combination 1 of HLDs, take-off: {delta_CL_max_obtained_option1_take_off:.3f}')
    print(f'delta CL, combination 1 of HLDs, landing: {delta_CL_max_obtained_option1_landing:.3f}')
    print(f'delta CL, combination 2 of HLDs, take-off: {delta_CL_max_obtained_option2_take_off:.3f}')
    print(f'delta CL, combination 2 of HLDs, landing: {delta_CL_max_obtained_option2_landing:.3f}')

    if delta_CL_max_obtained_option1_landing > delta_CL_approach_required:
        print('Combination 1 of HLDs works for landing')
    if delta_CL_max_obtained_option1_take_off > delta_CL_take_off_required:
        print('Combination 1 of HLDs works for take_off')
    if delta_CL_max_obtained_option1_landing < delta_CL_approach_required:
        print('Combination 1 of HLDs does not work for landing')
    if delta_CL_max_obtained_option1_take_off < delta_CL_take_off_required:
        print('Combination 1 of HLDs does not work for take_off')

    if delta_CL_max_obtained_option2_landing > delta_CL_approach_required:
        print('Combination 2 of HLDs works for landing')
    if delta_CL_max_obtained_option2_take_off > delta_CL_take_off_required:
        print('Combination 2 of HLDs works for take_off')
    if delta_CL_max_obtained_option2_landing < delta_CL_approach_required:
        print('Combination 2 of HLDs does not work for landing')
    if delta_CL_max_obtained_option2_take_off < delta_CL_take_off_required:
        print('Combination 2 of HLDs does not work for take_off')

    # print(f'The area of each TE HLD corresponding to option 1 is: {round(area_TE_HLD_option1,2)}')
    # print(f'The area of each TE HLD corresponding to option 2 is: {round(area_TE_HLD_option2,2)}')

    x = [0, a.b/2, a.b/2, 0,0]
    y = [0, math.tan(a.sweep_angle_le)*a.b/2,  math.tan(a.sweep_angle_le)*a.b/2 + a.c_t, a.c_r, 0]

    xmac = [a.ymac, a.ymac]
    ymac = [math.tan(a.sweep_angle_le) * a.ymac, math.tan(a.sweep_angle_le) * a.ymac + a.mac]

    ax = tab.ax
    ax.clear()
    
    ax.grid()
    ax.plot(x, y, label = "Wing Outline")
    ax.plot(xmac, ymac, label = 'MAC')

    ax.plot([a.Y_fuselage, a.Y_fuselage], [0, 20], label='Fuselage')

    for label, y0, y1 in [("TE HLD", TE_y0, TE_y1), (None, TE_y2, TE_y3)]:
        xf_0 = math.tan(math.radians(TE_sweep)) * y0 + a.c_r
        xf_1 = math.tan(math.radians(TE_sweep)) * y1 + a.c_r
        ax.plot([y0, y1, y1, y0, y0], [xf_0, xf_1, xf_1 - a.chord(y1) * a.cf_c_TE, xf_0 - a.chord(y0) * a.cf_c_TE, xf_0], label=label, color="red")


    for label, y0, y1 in [("LE HLD", LE_y0, LE_y1), (None, LE_y2, LE_y3)]:
        ratio_clf = 0.15
        
        xlf_0 = math.tan(math.radians(TE_sweep)) * y0 + a.c_r
        xlf_1 = math.tan(math.radians(TE_sweep)) * y1 + a.c_r
        ax.plot([y0, y1, y1, y0, y0], [xlf_0 - a.chord(y0) + a.chord(y0) * ratio_clf, xlf_1 - a.chord(y1) + a.chord(y1) * ratio_clf,
                xlf_1 - a.chord(y1), xlf_0 - a.chord(y0), xlf_0 - a.chord(y0) + a.chord(y0) * ratio_clf], label=label, color="purple")
    
    y0 = aileron_b_0
    y1 = aileron_b_1
    al_0 = math.tan(math.radians(TE_sweep)) * y0 + a.c_r
    al_1 = math.tan(math.radians(TE_sweep)) * y1 + a.c_r
    ax.plot([y0, y1, y1, y0, y0], [al_0, al_1, al_1 - a.chord(y1) * a.aileron_cf_ratio, al_0 - a.chord(y0) * a.aileron_cf_ratio, al_0], label="Aileron", color="cyan")

    ax.set_xlabel("[m]")
    ax.set_ylabel("[m]")
    ax.invert_yaxis()

    ax.legend()

    tab.canvas.draw()
