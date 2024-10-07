from isa import *
from utils import *

class Airplane:
    def __init__(self, data):
        self.loaded_from_data = data
        for key, value in data.items():
            if not key.startswith("__comment_"):
                setattr(self, key, value)

    def __repr__(self):
        return f"<Airplane: {self.__dict__}>"

    def c_d0(self, gear_extended, flap_deflection):
        "Roughly esimates cd0 + contribution from flaps and landing gear"
        res = self.c_d0_cruise + 0.0013 * flap_deflection # 7.57
        if gear_extended: res += 0.02 # 7.58 Reader says to est. a value between 0.0100 and 0.0250 
        return res

    def oswald_eff(self, flap_deflection):
        "Roughly esimates Oswald efficiency factor + contribution from flaps"
        return self.oswald_eff_cruise + 0.0026 * flap_deflection # 7.56 wing-mounted engines 
    
    def sweep_angle_at(self, x_c):
        """
        x_c should be from 0 to 1
        0 at LE and 1 on TE, 0.25 at c/4, etc.
        """
        return math.atan((math.tan(self.sweep_angle_le) - x_c * (2*self.c_r / self.b) * (1 - self.taper_ratio)))

    def chord(self, y):
        """
        Returns the length of the chord located at 'y' spanwise (from 0 to b/2)
        """
        return self.c_r*(((self.taper_ratio-1)/(self.b/2))*y+1)
