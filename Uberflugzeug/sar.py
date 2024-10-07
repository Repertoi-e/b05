from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QFileDialog, QTabWidget, QComboBox,
                             QVBoxLayout, QHBoxLayout, QWidget, QPlainTextEdit, QLineEdit, QLabel, QMessageBox, QToolBar, QFormLayout, QSizePolicy)

import math
from utils import *

import tab

class Tab(tab.Tab):
    def __init__(self, uber):
        super().__init__(uber, "SAR")

    def init_ui_for_airplane(self):
        super().init_ui_for_airplane()

def recalculate(tab):
    a = tab.uber.airplane

    cdi = a.CLdes_airplane**2 / (math.pi * a.aspect_ratio * a.e) 
    print("cdi:", cdi)

    D = (a.airfoil_cd + cdi) * 0.5 * a.rho_cr * a.V_cr**2 * a.S
    SAR_drag_wing_only = a.V_cr / (D * a.TSFC)

    print(f"SAR (drag from wing only): {SAR_drag_wing_only:.2f} m/kg")