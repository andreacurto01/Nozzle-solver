import numpy as np

class RocketNozzleThermal:
    """
    Thermal Module for the AeroRocket-Nozzle-Solver.
    Calculates the gas-side heat transfer coefficient using Bartz's Equation
    and solves the 1D radial thermal resistance network iteratively with variable wall thickness.
    """

    def __init__(self, cfd_solver, fluid_properties, thermal_properties):
        # Extract mesh and physical results from CFD
        self.x = cfd_solver.x
        self.A = cfd_solver.A
        self.M = cfd_solver.Mach
        self.T_static = cfd_solver.T
        
        # Extract chamber and flow parameters from CFD
        self.gamma = cfd_solver.gamma
        self.P_c = cfd_solver.Pc      # Chamber pressure [Pa]
        self.T_0 = cfd_solver.Tc      # Stagnation (Chamber) temperature [K]
        self.At = cfd_solver.At       # Throat area [m^2]
        
        # Fluid / Combustion transport properties 
        self.Pr = fluid_properties['Pr']          
        self.mu_0 = fluid_properties['mu_0']      
        self.Cp = fluid_properties['Cp']          
        self.c_star = fluid_properties['c_star']  
        
        # Structural Wall and Cooling Jacket parameters
        self.k_w = thermal_properties['k_w']      
        self.h_c = thermal_properties['h_c']      
        self.T_cool = thermal_properties['T_cool']  
        
        # GENERA IL PROFILO DI SPESSORE VARIABILE (Dinamico lungo l'asse x)
        self.t_w = np.zeros_like(self.x)
        for i, x_val in enumerate(self.x):
            if x_val < -0.05:
                self.t_w[i] = 0.0040  # 4.0 mm in Camera/Inlet per reggere la pressione
            elif x_val >= -0.05 and x_val <= 0.05:
                dist_from_throat = abs(x_val)
                # Restringimento lineare simmetrico fino a 1.2 mm esatti al collo (x=0)
                self.t_w[i] = 0.0012 + (0.0040 - 0.0012) * (dist_from_throat / 0.05)
            else:
                self.t_w[i] = 0.0020  # 2.0 mm nel Divergente di scarico
        
        # Reconstruct local diameters from the area profile
        self.D = 2.0 * np.sqrt(self.A / np.pi)
        self.D_t = 2.0 * np.sqrt(self.At / np.pi)
        
        # Allocate empty arrays to store the calculated thermal profiles
        self.T_aw = np.zeros_like(self.x)
        self.h_g = np.zeros_like(self.x)
        self.T_wi = np.zeros_like(self.x)
        self.T_wo = np.zeros_like(self.x)
        self.q_flux = np.zeros_like(self.x)

    def calculate_adiabatic_wall_temperature(self):
        r = self.Pr**(1.0 / 3.0)
        self.T_aw = self.T_static * (1.0 + r * ((self.gamma - 1.0) / 2.0) * self.M**2)

    def solve_thermal_node(self, idx, max_iter=100, tol=1e-4):
        T_wi_guess = 0.5 * (self.T_cool + self.T_aw[idx])
        
        bartz_base = (0.026 / (self.D_t**0.2)) * \
                     ((self.mu_0**0.2 * self.Cp) / (self.Pr**0.6)) * \
                     ((self.P_c / self.c_star)**0.8) * \
                     ((self.D_t / self.D[idx])**1.8)
        
        # Recuperiamo lo spessore locale per questo specifico nodo
        t_w_local = self.t_w[idx]
        
        for iteration in range(max_iter):
            g_factor = (self.gamma - 1.0) / 2.0
            sub_bracket_1 = (T_wi_guess / (2.0 * self.T_0)) * (1.0 + g_factor * self.M[idx]**2) + 0.5
            sub_bracket_2 = 1.0 + g_factor * self.M[idx]**2
            
            sigma = (sub_bracket_1**(-0.68)) * (sub_bracket_2**(-0.12))
            h_g_current = bartz_base * sigma
            
            # Calcolo della resistenza usando lo SPESSORE LOCALE PUNTUALE
            R_total = (1.0 / h_g_current) + (t_w_local / self.k_w) + (1.0 / self.h_c)
            
            q_current = (self.T_aw[idx] - self.T_cool) / R_total
            
            T_wi_new = self.T_aw[idx] - (q_current / h_g_current)
            T_wo_new = self.T_cool + (q_current / self.h_c)
            
            if abs(T_wi_new - T_wi_guess) < tol:
                return h_g_current, T_wi_new, T_wo_new, q_current
            
            T_wi_guess = T_wi_new
            
        return h_g_current, T_wi_new, T_wo_new, q_current

    def run_solver(self):
        self.calculate_adiabatic_wall_temperature()
        for i in range(len(self.x)):
            h_g_i, T_wi_i, T_wo_i, q_i = self.solve_thermal_node(i)
            self.h_g[i] = h_g_i
            self.T_wi[i] = T_wi_i
            self.T_wo[i] = T_wo_i
            self.q_flux[i] = q_i
            
        return {
            'x': self.x,
            'h_g': self.h_g,
            'T_wi': self.T_wi,
            'T_wo': self.T_wo,
            'q_flux': self.q_flux,
            'T_aw': self.T_aw
        } 
