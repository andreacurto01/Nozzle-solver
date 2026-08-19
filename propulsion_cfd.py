import numpy as np
from scipy.optimize import fsolve

class RocketNozzleCFD:
    def __init__(self, Pc, Tc, gamma, R, D_throat, D_exit, L_convergent, L_divergent):
        """
        Initialize engine and nozzle parameters.
        Pc: Chamber pressure (Pa)
        Tc: Chamber temperature (K)
        gamma: Specific heat ratio (Cp/Cv)
        R: Specific gas constant (J/kg*K)
        """
        self.Pc = Pc
        self.Tc = Tc
        self.gamma = gamma
        self.R = R
        
        # Nozzle Geometry
        self.At = np.pi * (D_throat / 2)**2
        self.Ae = np.pi * (D_exit / 2)**2
        self.L_conv = L_convergent
        self.L_div = L_divergent
        
        # Arrays to be populated
        self.x = None
        self.A = None
        self.Mach = None
        self.P = None
        self.T = None
        
    def generate_geometry(self, num_points=100):
        """Generates the area profile along the x-axis of the nozzle"""
        # Create x-axis: from -L_conv to +L_div (throat is at x = 0)
        x_conv = np.linspace(-self.L_conv, 0, num_points // 2)
        x_div = np.linspace(0, self.L_div, num_points // 2 + 1)[1:] # Avoid duplicating x=0
        self.x = np.concatenate([x_conv, x_div])
        
        # Linear approximation for diameters (simplified geometric profile)
        D_throat = 2 * np.sqrt(self.At / np.pi)
        D_exit = 2 * np.sqrt(self.Ae / np.pi)
        D_inlet = D_throat * 2.5 # Assumption for chamber inlet size
        
        D_x = np.zeros_like(self.x)
        D_x[self.x <= 0] = D_inlet + (D_throat - D_inlet) * (self.x[self.x <= 0] + self.L_conv) / self.L_conv
        D_x[self.x > 0] = D_throat + (D_exit - D_throat) * self.x[self.x > 0] / self.L_div
        
        self.A = np.pi * (D_x / 2)**2
        return self.x, self.A

    def _area_mach_equation(self, M, A_target, region):
        """Isentropic area-Mach relation residual function to be solved for M"""
        g = self.gamma
        term = (2 / (g + 1)) * (1 + 0.5 * (g - 1) * M**2)
        exponent = (g + 1) / (2 * (g - 1))
        
        # Physical constraints to guide fsolve and prevent unphysical solutions
        if region == 'subsonic' and M > 1.0:
            return 1e6 * (M - 1.0) # Numerical penalty for supersonic values in convergent
        if region == 'supersonic' and M < 1.0:
            return 1e6 * (1.0 - M) # Numerical penalty for subsonic values in divergent
            
        return (1/M) * (term**exponent) - (A_target / self.At)

    def solve_flow(self):
        """Solves the quasi-1D flow calculating Mach, static Pressure, and static Temperature"""
        self.Mach = np.zeros_like(self.x)
        
        for i, x_val in enumerate(self.x):
            current_area = self.A[i]
            
            if x_val <= 0:
                initial_guess = 0.3
                region = 'subsonic'
            else:
                initial_guess = 2.0
                region = 'supersonic'
                
            # Solves the non-linear area-Mach equation
            M_sol = fsolve(self._area_mach_equation, x0=initial_guess, args=(current_area, region))
            self.Mach[i] = M_sol[0]
            
        # Calculate static T and P using isentropic stagnation relations
        self.T = self.Tc / (1 + 0.5 * (self.gamma - 1) * self.Mach**2)
        self.P = self.Pc / (1 + 0.5 * (self.gamma - 1) * self.Mach**2)**(self.gamma / (self.gamma - 1))
        
        return self.Mach, self.P, self.T

# --- MODULE TEST ---
if __name__ == "__main__":
    # Input data for a typical liquid rocket engine (e.g., LOX/Kerosene)
    Pc_input = 7e6          # 70 bar in Pascal
    Tc_input = 3300.0       # 3300 Kelvin
    gamma_input = 1.22
    R_input = 360.0         # Specific gas constant (J/kg*K)
    
    # Geometry data
    D_throat_input = 0.08   # 8 cm throat diameter
    D_exit_input = 0.24     # 24 cm exit diameter
    L_conv_input = 0.15     # 15 cm convergent length
    L_div_input = 0.40      # 40 cm divergent length

    # Instantiate and run the solver
    nozzle = RocketNozzleCFD(Pc_input, Tc_input, gamma_input, R_input, D_throat_input, D_exit_input, L_conv_input, L_div_input)
    x, A = nozzle.generate_geometry()
    M, P, T = nozzle.solve_flow()
    
    print("=== CFD CORE RESULTS ===")
    print(f"Inlet  -> Mach: {M[0]:.2f}, Pressure: {P[0]/1e5:.1f} bar, Temp: {T[0]:.0f} K")
    print(f"Throat -> Mach: {M[len(M)//2]:.2f}, Pressure: {P[len(M)//2]/1e5:.1f} bar, Temp: {T[len(M)//2]:.0f} K")
    print(f"Exit   -> Mach: {M[-1]:.2f}, Pressure: {P[-1]/1e5:.2f} bar, Temp: {T[-1]:.0f} K")
