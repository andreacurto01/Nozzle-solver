import numpy as np

class RocketNozzleStructural:
    """
    Structural Module for the AeroRocket-Nozzle-Solver.
    Calculates mechanical and thermal stresses along the nozzle profile
    using local thick shell theory and variable wall thickness.
    """
    
    def __init__(self, cfd_solver, thermal_solver, material_properties):
        # Mesh and Geometry from previous modules
        self.x = cfd_solver.x
        self.A = cfd_solver.A
        self.P_gas = cfd_solver.P          # Local static pressure from CFD [Pa]
        self.D_ext = thermal_solver.D      # Outer diameter of the inner wall [m]
        
        # Ereditiamo l'esatto profilo di spessore variabile calcolato dal modulo termico
        self.t_w = thermal_solver.t_w      
        
        # Calcolo dei raggi locali (Interno ed Esterno) punto per punto
        self.R_ext = self.D_ext / 2.0
        self.R_int = self.R_ext - self.t_w
        
        # Thermal profiles from Step 2
        self.T_wi = thermal_solver.T_wi    # Inner wall temperature [K]
        self.T_wo = thermal_solver.T_wo    # Outer wall temperature [K]
        
        # Material Mechanical Properties
        self.E = material_properties['E']          
        self.alpha = material_properties['alpha']  
        self.nu = material_properties['nu']        
        self.sigma_y = material_properties['sigma_y'] 
        
        self.P_amb = 101325.0              # Standard atmospheric pressure [Pa]
        
        # Pre-allocate stress arrays
        self.sigma_hoop_p = np.zeros_like(self.x)
        self.sigma_hoop_t_i = np.zeros_like(self.x)
        self.sigma_hoop_t_o = np.zeros_like(self.x)
        
        self.sigma_axial_p = np.zeros_like(self.x)
        self.sigma_axial_t_i = np.zeros_like(self.x)
        self.sigma_axial_t_o = np.zeros_like(self.x)
        
        self.sigma_vm_i = np.zeros_like(self.x)
        self.sigma_vm_o = np.zeros_like(self.x)
        self.margin_of_safety = np.zeros_like(self.x)

    def calculate_pressure_stresses(self):
        delta_P = self.P_gas - self.P_amb
        for i in range(len(self.x)):
            ri = self.R_int[i]
            ro = self.R_ext[i]
            
            # Hoop Stress from pressure (Lamé - Maximum at inner surface)
            self.sigma_hoop_p[i] = delta_P[i] * (ro**2 + ri**2) / (ro**2 - ri**2)
            # Axial Stress from pressure (Closed-end cylinder assumption)
            self.sigma_axial_p[i] = delta_P[i] * (ri**2) / (ro**2 - ri**2)

    def calculate_thermal_stresses(self):
        delta_T = self.T_wi - self.T_wo
        # Real-time calculation of thermal stress based on localized temperature gradient
        factor = (self.E * self.alpha * delta_T) / (2.0 * (1.0 - self.nu))
        
        self.sigma_hoop_t_i = -factor
        self.sigma_hoop_t_o = factor
        self.sigma_axial_t_i = -factor
        self.sigma_axial_t_o = factor

    def evaluate_von_mises_and_margins(self):
        self.calculate_pressure_stresses()
        self.calculate_thermal_stresses()
        
        for i in range(len(self.x)):
            # 1. Total Stresses at INNER surface
            s_hoop_i = self.sigma_hoop_p[i] + self.sigma_hoop_t_i[i]
            s_axial_i = self.sigma_axial_p[i] + self.sigma_axial_t_i[i]
            s_radial_i = -self.P_gas[i]
            
            vm_i_sq = 0.5 * ((s_hoop_i - s_axial_i)**2 + 
                             (s_axial_i - s_radial_i)**2 + 
                             (s_radial_i - s_hoop_i)**2)
            self.sigma_vm_i[i] = np.sqrt(vm_i_sq)
            
            # 2. Total Stresses at OUTER surface
            s_hoop_p_o = 2.0 * (self.P_gas[i] - self.P_amb) * (self.R_int[i]**2) / (self.R_ext[i]**2 - self.R_int[i]**2)
            s_hoop_o = s_hoop_p_o + self.sigma_hoop_t_o[i]
            s_axial_o = self.sigma_axial_p[i] + self.sigma_axial_t_o[i]
            s_radial_o = -self.P_amb
            
            vm_o_sq = 0.5 * ((s_hoop_o - s_axial_o)**2 + 
                             (s_axial_o - s_radial_o)**2 + 
                             (s_radial_o - s_hoop_o)**2)
            self.sigma_vm_o[i] = np.sqrt(vm_o_sq)
            
            # 3. Margin of Safety (MoS) based on the worst-case stress
            max_vm = max(self.sigma_vm_i[i], self.sigma_vm_o[i])
            self.margin_of_safety[i] = (self.sigma_y / max_vm) - 1.0
            
        return {
            'x': self.x,
            'sigma_vm_inner': self.sigma_vm_i,
            'sigma_vm_outer': self.sigma_vm_o,
            'margin_of_safety': self.margin_of_safety
        }