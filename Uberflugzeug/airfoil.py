from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QFileDialog, QTabWidget, QComboBox,
                             QVBoxLayout, QHBoxLayout, QWidget, QPlainTextEdit, QLineEdit, QLabel, QMessageBox, QToolBar, QFormLayout, QSizePolicy)

import math
import os
import glob

import pandas as pd

from xfoil import XFoil
from xfoil.model import Airfoil

from scipy.optimize import fsolve

from utils import *

import tab

class Tab(tab.Tab):
    def __init__(self, uber):
        super().__init__(uber, "Airfoil")

    def init_ui_for_airplane(self):
        super().init_ui_for_airplane(plot_nrows=1, plot_ncols=2)

        dat_files = glob.glob(os.path.join(os.path.dirname(__file__), 'airfoils/*.dat'))

        self.airfoils = []

        combobox = QComboBox()

        for file_path in dat_files:
            with open(file_path, 'r') as file:
                name = os.path.basename(file_path)
                combobox.addItem(name)
                self.airfoils.append((name, pd.read_table(file_path, sep='\s+', skiprows=[0], names=['x','y', 'tc'], index_col=False)))

        def select_airfoil(i):
            airfoil_name, airfoil = self.airfoils[i]
            self.airfoil_name = airfoil_name
            self.airfoil = airfoil
            self.uber.recalculate()

        combobox.setCurrentIndex(3)
        combobox.currentIndexChanged.connect(select_airfoil)
        self.airfoil_name, self.airfoil = self.airfoils[3]

        self.separator("Select airfoil")
        self.addRow(combobox)

        self.separator("DATCOM")

        self.create_input_box("CL/cl ratio", "CL_cl_ratio")
        self.create_input_box("Technology f-r for DATCOM lift slope", "technology_factor_for_slope")


def beta(M_inf):
    return math.sqrt(1 - M_inf ** 2)

def laitone(M_inf, x):
    correction_term = M_inf ** 2 * (1 + ((GAMMA - 1) / 2) * M_inf ** 2) / (2 * beta(M_inf))
    return x / (beta(M_inf) + correction_term * x)

def karmantsien(M_inf, x):
    correction_term = (M_inf ** 2) / (1 + beta(M_inf))
    return x / (beta(M_inf) + correction_term * (x/2))

def prandtlglauert(M_inf, x):
    return x / (beta(M_inf))

def recalculate(tab):
    a = tab.uber.airplane

    q_cr = 0.5 * a.rho_cr * a.V_cr**2

    _CLdes = 1.1 / (q_cr) * (a.WS * a.mass_fr_cr + a.WS * a.mass_fr_land) / 2

    cosLE = math.cos(a.sweep_angle_le)

    V_eff = cosLE * a.V_cr
    Cldes = _CLdes * a.V_cr**2 / V_eff**2

    a.Re_cr = Re(a.altitude_cr, V_eff, a.mac)
    a.Re_app = Re(0, cosLE * a.approach_speed, a.mac)

    a.CLdes_airplane = _CLdes
    a.Cldes_airfoil = Cldes

    print("        @ CRUISE")
   
    print(f"C_L_des (airplane): {_CLdes:.4f}")
    print(f"C_l_des (airfoil):  {Cldes:.4f}")
    print()
    print("Re_cr:", a.Re_cr)
    print()

    xf_cr = XFoil()
    xf_cr.airfoil=Airfoil(x=tab.airfoil.x.values, y=tab.airfoil.y.values)
    xf_cr.Re = a.Re_cr
    xf_cr.M = 0.0
    xf_cr.max_iter = 40

    aoa, cd, cm, min_cp = xf_cr.cl(Cldes)
    print(f"{tab.airfoil_name} @ Cldes:\n    aoa={aoa:.1f}, cd={cd:.6f}, min cp={min_cp:.2f}, cm={cm:.4f}")

    a.design_aoa = aoa
    a.airfoil_cd = laitone(a.M_cruise, cd)

    t_c = tab.airfoil.iloc[0].tc / 100

    if 0:
        def find_m_cr(M_cr, corr_func):
            lhs = corr_func(M_cr, min_cp)
            rhs = (2 / (GAMMA * M_cr**2)) * (((1 + ((GAMMA - 1) / 2) * M_cr**2) / (1 + (GAMMA - 1) / 2))**(GAMMA / (GAMMA - 1)) - 1)
            return lhs - rhs
    
        # These are not reliable for non-zero AOA cause pressure coefficient is really wrong from XFOIL 

        M_cr_laitone = fsolve(lambda M_cr: find_m_cr(M_cr, laitone), 0.7)[0] / cosLE
        M_cr_karmantsien = fsolve(lambda M_cr: find_m_cr(M_cr, karmantsien), 0.7)[0] / cosLE
        M_cr_prandtl = fsolve(lambda M_cr: find_m_cr(M_cr, prandtlglauert), 0.7)[0] / cosLE

        M_dd_laitone = M_cr_laitone + (0.1/80)**(1/3)
        M_dd_karmantsien = M_cr_karmantsien + (0.1/80)**(1/3)
        M_dd_prandtl = M_cr_prandtl + (0.1/80)**(1/3)

        print(f"    M_cr_laitone={M_cr_laitone:.2f}, M_cr_karmantsien={M_cr_karmantsien:.2f}, M_cr_prandtl={M_cr_prandtl:.2f}")
        print(f"    M_dd_laitone={M_dd_laitone:.2f}, M_dd_karmantsien={M_dd_karmantsien:.2f}, M_dd_prandtl={M_dd_prandtl:.2f}")

    k_a = 0.95

    M_dd = (k_a/cosLE) - (t_c/(cosLE**2)) - (_CLdes/(10 * (cosLE**3)))
    M_cr = M_dd - (0.1/80)**(1/3)
    print(f"    M_cr_est={M_cr:.2f}, M_dd_est={M_dd:.2f}")

    print()
    print(f"Design AOA: {aoa:2f}")
    print(f"Airfoil Cd, laitone corrected: {a.airfoil_cd:4f}")

    a.M_dd_airfoil = M_dd
    a.M_cr_airfoil = M_cr

    ax_curves = tab.ax[0]
    ax_curves.clear()
     
    df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'airfoils/Polar_Graph_1.csv'), header=0)
    for i, (airfoil_name, airfoil) in enumerate(tab.airfoils):
        df_a = df.iloc[:, [3 * i, 3 * i + 1]].apply(pd.to_numeric, errors='coerce').dropna()
        aoa = df_a.iloc[:, 0]
        cl = df_a.iloc[:, 1]

        ax_curves.plot(aoa, cl, label=' '.join(cl.name.split()[:3]))

    ax_curves.set_xlabel("Alpha")
    ax_curves.set_ylabel("Cl")
    ax_curves.set_title(f"Lift curves (approach)")
    ax_curves.set_ylim((-0.5, 2.5))
    ax_curves.set_xlim((-10, 25))
    ax_curves.grid()
    # ax_curves.legend()

    #
    # Chosen airfoil @ approach conditions
    #
    xf_app = XFoil()
    xf_app.airfoil=Airfoil(x=tab.airfoil.x.values, y=tab.airfoil.y.values)
    xf_app.Re = a.Re_app
    xf_app.M = mach(a.approach_speed, 288.15)
    xf_app.max_iter = 40

    print()
    print()
    print("        @ APPROACH")
    print("Re_approach:", a.Re_app)
    print()
    # Find alpha L=0
    a.alpha_l0, _, _, _ = xf_app.cl(0)
    print(f"Airfoil alpha,L=0: {a.alpha_l0:2f} deg")    

    # Find airfoil Cl_curve 
    approach_aoa, approach_cl, approach_cd, _, _ = xf_app.aseq(0, 5, 0.5)
    a.clalpha_airfoil, a.cl_at_aoa0 = np.polyfit(np.radians(approach_aoa), approach_cl, 1)
    print(f"Airfoil Cl_alpha: {a.clalpha_airfoil:2f} [1/rad]")    
    print(f"Airfoil Cl_aoa=0: {a.cl_at_aoa0:2f}")    

    a.cd0_airfoil_approach = approach_cd[0]

    # Find max Cl
    _, approach_cl, _, _, _ = xf_app.aseq(15, 20, 0.5)
    approach_cl = approach_cl[~np.isnan(approach_cl)]
    a.Clmax_airfoil = max(approach_cl)

    print(f"Airfoil Cl_max: {a.Clmax_airfoil:2f}")    
    print()

    a.CLmax_airplane = a.Clmax_airfoil * a.CL_cl_ratio
    print(f"DATCOM wing CL_max: {a.CLmax_airplane:2f}")    

    beta = np.sqrt(1 - a.M_cruise**2)
    a.CLalpha_wing = 2 * np.pi * a.aspect_ratio / (2 + np.sqrt(4 + (a.aspect_ratio * beta / a.technology_factor_for_slope)**2 * (1 + math.tan(a.sweep_angle_at(0.5))**2 / (beta**2))))
    print(f"DATCOM wing CL_alpha: {a.CLalpha_wing:2f} [1/rad]")    

    a.stall_aoa = math.degrees(a.CLmax_airplane / a.CLalpha_wing) + a.alpha_l0 + 3
    print(f"DATCOM wing stall AOA: {a.stall_aoa:2f} [deg]")    
    print()

    a.e = 2 / (2 - a.aspect_ratio + math.sqrt(4 + a.aspect_ratio**2 * (1 + math.tan(a.sweep_angle_at(0.5))**2)))
    print(f"Estimated wing Oswald efficiency: {a.e:2f}")

    ax_wing = tab.ax[1]
    ax_wing.clear()

    ax_wing.set_xlabel("Alpha")
    ax_wing.set_ylabel("Cl")

    ax_wing.plot([-10, a.alpha_l0, 10], [np.radians(-10 - a.alpha_l0) * a.CLalpha_wing, 0, np.radians(10 - a.alpha_l0) * a.CLalpha_wing], label="Wing [DATCOM]")
    
    ax_wing.plot([-10, a.alpha_l0, 10], [np.radians(-10 - a.alpha_l0) * a.clalpha_airfoil, 0, np.radians(10 - a.alpha_l0) * a.clalpha_airfoil], label="Airfoil [XFOIL]")

    ax_wing.scatter([a.stall_aoa], [a.CLmax_airplane], label="Stall point [DATCOM]", color="red")
    ax_wing.set_ylim((-0.5, 2.5))
    ax_wing.set_xlim((-10, 25))

    ax_wing.hlines([0], -10, 25, color="black")
    ax_wing.vlines([0], -10, 25, color="black")

    ax_wing.set_title(f"Lift curve {' '.join(tab.airfoil_name.split()[0:2])} (approach)")

    ax_wing.legend()

    ax_wing.grid()

    tab.canvas.draw()