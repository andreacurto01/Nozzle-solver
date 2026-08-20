# =====================================================================
# AEROSOCKET-NOZZLE-SOLVER: MAIN MULTIPHYSICS INTEGRATION TEST
# =====================================================================

import numpy as np

from propulsion_cfd import RocketNozzleCFD
from thermal_bartz import RocketNozzleThermal
from RocketNozzleStructural import RocketNozzleStructural

if __name__ == "__main__":
    print("====================================================")
    print("  LAUNCHING AEROROCKET MULTIPHYSICS SOLVER (2026)   ")
    print("====================================================\n")

    # -----------------------------------------------------------------
    # STEP 1: INPUT DATA SETUP
    # -----------------------------------------------------------------
    Pc_input = 70e5          # Chamber Pressure: 70 bar [Pa]
    Tc_input = 3300.0        # Chamber Temperature [K]
    gamma_input = 1.22       # Specific heat ratio
    R_input = 360.0          # Specific gas constant [J/kg*K]
    
    D_throat_input = 0.08    # Throat diameter: 8 cm [m]
    D_exit_input = 0.24      # Exit diameter: 24 cm [m]
    L_conv_input = 0.15      # Convergent length: 15 cm [m]
    L_div_input = 0.40       # Divergent length: 40 cm [m]

    fluid_props = {
        'Pr': 0.8,           
        'mu_0': 8.5e-5,      
        'Cp': 2000.0,        
        'c_star': 1650.0     
    }

    thermal_props = {       
        'k_w': 320.0,        
        'h_c': 25000.0,      
        'T_cool': 350.0      
    }

    material_props = {
        'E': 110e9,          
        'alpha': 17e-6,      
        'nu': 0.33,          
        'sigma_y': 420e6     
    }

    # -----------------------------------------------------------------
    # STEP 2: EXECUTE MULTIPHYSICS PIPELINE
    # -----------------------------------------------------------------
    
    # 1. Run CFD Module
    print("--> Running Quasi-1D Isentropic CFD Solver...")
    cfd = RocketNozzleCFD(
        Pc_input, Tc_input, gamma_input, R_input, 
        D_throat_input, D_exit_input, L_conv_input, L_div_input
    )
    cfd.generate_geometry(num_points=100)
    cfd.solve_flow()
    
    # 2. Run Thermal Module
    print("--> Running 1D Radial Thermal Resistance Network (Bartz)...")
    thermal = RocketNozzleThermal(cfd, fluid_props, thermal_props)
    thermal.run_solver()
    
    # 3. Run Structural Module
    print("--> Running Structural Stress & Margin of Safety Evaluation...")
    structural = RocketNozzleStructural(cfd, thermal, material_props)
    struct_results = structural.evaluate_von_mises_and_margins()
    
    print("\n====================================================")
    print("               SIMULATION RESULTS SUMMARY           ")
    print("====================================================")
    
    idx_inlet = 0
    idx_throat = np.argmin(np.abs(cfd.x))
    idx_exit = -1
    
    sections = [("INLET", idx_inlet), ("THROAT", idx_throat), ("EXIT", idx_exit)]
    
    for name, idx in sections:
        print(f"\n[{name} SECTION (x = {cfd.x[idx]:.3f} m)]")
        print(f"  Gas Pressure     : {cfd.P[idx]/1e5:.1f} bar")
        print(f"  Wall Temp (In/Out): {thermal.T_wi[idx]:.1f} K / {thermal.T_wo[idx]:.1f} K")
        print(f"  Von Mises Inner  : {structural.sigma_vm_i[idx]/1e6:.1f} MPa")
        print(f"  Von Mises Outer  : {structural.sigma_vm_o[idx]/1e6:.1f} MPa")
        
        mos = structural.margin_of_safety[idx]
        status = "PASSED" if mos >= 0 else "FAILED (Plastic Yielding!)"
        print(f"  Margin of Safety : {mos:.3f} --> {status}")

    print("\n====================================================")
    print("Multiphysics Analysis Complete. Modular File Integration Successful!")
    print("====================================================")
