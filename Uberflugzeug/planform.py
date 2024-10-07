import math
import matplotlib.pyplot as plt

from utils import *

from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QFileDialog, QTabWidget,
                             QVBoxLayout, QHBoxLayout, QWidget, QPlainTextEdit, QLineEdit, QLabel, QMessageBox, QToolBar, QFormLayout, QSizePolicy)

import tab

class Tab(tab.Tab):
    def __init__(self, uber):
        super().__init__(uber, "Planform")

    def init_ui_for_airplane(self):
        super().init_ui_for_airplane()

        self.create_input_box("Aspect Ratio", "aspect_ratio")

def recalculate(tab):
    a = tab.uber.airplane

    a.sweep_angle_c4 = math.acos(1.16/(a.M_cruise + 0.5))
    a.taper_ratio = 0.2 * (2 - a.sweep_angle_c4)

    a.b = math.sqrt(a.aspect_ratio * a.S)
    if a.b > 80:
        raise AirplaneValueError("Wing span must be < 80")

    a.c_r = (2 * a.S) / ((1+a.taper_ratio) * a.b)
    a.c_t = a.taper_ratio * a.c_r

    a.dihedral_angle_deg = 3 - 0.1 * math.degrees(a.sweep_angle_c4) + 2

    a.sweep_angle_le = math.tan(a.sweep_angle_c4) - a.c_r*(a.taper_ratio-1)/(2*a.b)

    a.mac = (2/3) * a.c_r * ((1 + a.taper_ratio + (a.taper_ratio)**2)/(1 + a.taper_ratio))

    a.ymac = (a.b / 6) * (1 + 2*a.taper_ratio)/(1 + a.taper_ratio)
    a.xlemac = a.ymac * math.tan(a.sweep_angle_le)

    AR_max = math.floor(7.72*(2-a.taper_ratio)*math.exp(-0.043  * a.sweep_angle_c4))
    if a.aspect_ratio > AR_max:
        raise AirplaneValueError(f"Aspect ratio must be < {AR_max}")

    print(f"Wing area: {a.S:.1f} m^2")
    print(f"Wing span: {a.b:.2f} m")
    print(f"Aspect ratio: {a.aspect_ratio:.2f}")
    print(f"Taper ratio: {a.taper_ratio:.2f}")
    print(f"Root chord: {a.c_r:.1f} m")
    print(f"Tip chord: {a.c_t:.1f} m")
    print()
    print(f"Sweep angle at c/4: {math.degrees(a.sweep_angle_c4):.2f} deg")
    print(f"Sweep angle at LE: {math.degrees(a.sweep_angle_le):.2f} deg")
    print(f"Dihedral angle: {a.dihedral_angle_deg:.1f} deg")
    print(f"MAC: {a.mac:.1f} m, X_LEMAC: {a.xlemac:.1f} m, YMAC: {a.ymac:.1f} m")

    x = [0, a.b/2, a.b/2, 0,0]
    y = [0, math.tan(a.sweep_angle_le)*a.b/2,  math.tan(a.sweep_angle_le)*a.b/2 + a.c_t, a.c_r, 0]

    xmac = [a.ymac, a.ymac]
    ymac = [math.tan(a.sweep_angle_le) * a.ymac, math.tan(a.sweep_angle_le) * a.ymac + a.mac]

    ax = tab.ax
    ax.clear()
    
    ax.grid()
    ax.plot(x, y, label = "Wing Outline")
    ax.plot(xmac, ymac, label = 'MAC')
    ax.set_xlabel("[m]")
    ax.set_ylabel("[m]")
    ax.invert_yaxis()

    tab.canvas.draw()
