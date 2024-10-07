import numpy as np
import math
import matplotlib.pyplot as plt
import json
from isa import *

GAMMA = 1.4      # Ratio of specific heats for air

PI = math.pi
def SQRT(x): 
    return np.sqrt(x)

def speed_of_sound(T):
    return math.sqrt(GAMMA * R_AIR * T)

def mach(speed, T): 
    return speed / speed_of_sound(T)

def total_temperature(M, T_static): # 7.24
    return T_static * (1 + ((GAMMA - 1) / 2) * M**2)

def total_pressure(M, p_static): # 7.25
    return p_static * (1 + ((GAMMA - 1) / 2) * M**2)**(GAMMA/(GAMMA-1))

def Re(altitude, V, L):
    T, _, rho = isa(altitude)

    mu_0 = 1.716e-5
    C = 110.4
    T0 = 273.15

    mu = mu_0 * (T/T0)**(3/2) * (T0 + C) / (T + C)
    return rho * V * L / mu


class AirplaneValueError(BaseException):
    def __init__(self, msg):
        self.msg = msg