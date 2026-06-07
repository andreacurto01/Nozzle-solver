import numpy as np
import matplotlib.pyplot as plt
from propulsion_cfd import RocketNozzleCFD
from thermal_bartz import RocketNozzleThermal
from RocketNozzleStructural import RocketNozzleStructural

def generate_multiphysics_plots():
    Pc_input = 70e5          
    Tc_input = 3300.0        
    gamma_input = 1.22       
    R_input = 360.0          
    D_throat_input = 0.08    
    D_exit_input = 0.24      
    L_conv_input = 0.15      
    L_div_input = 0.40       

    fluid_props = {'Pr': 0.8, 'mu_0': 8.5e-5, 'Cp': 2000.0, 'c_star': 1650.0}
    thermal_props = {'t_w': 0.004, 'k_w': 320.0, 'h_c': 25000.0, 'T_cool': 350.0}
    material_props = {'E': 110e9, 'alpha': 17e-6, 'nu': 0.33, 'sigma_y': 420e6}

    cfd = RocketNozzleCFD(Pc_input, Tc_input, gamma_input, R_input, D_throat_input, D_exit_input, L_conv_input, L_div_input)
    cfd.generate_geometry(num_points=100)
    cfd.solve_flow()
    
    thermal = RocketNozzleThermal(cfd, fluid_props, thermal_props)
    thermal.run_solver()
    
    structural = RocketNozzleStructural(cfd, thermal, material_props)
    structural.evaluate_von_mises_and_margins()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("AeroRocket Nozzle Solver - Multiphysics Structural Assessment", fontsize=14, fontweight='bold')

    sigma_vm_i_MPa = structural.sigma_vm_i / 1e6
    sigma_vm_o_MPa = structural.sigma_vm_o / 1e6
    sigma_y_MPa = material_props['sigma_y'] / 1e6

    ax1.plot(cfd.x, sigma_vm_i_MPa, label="Von Mises - Inner Wall (Gas Side)", color='tab:red', linewidth=2)
    ax1.plot(cfd.x, sigma_vm_o_MPa, label="Von Mises - Outer Wall (Coolant Side)", color='tab:orange', linewidth=2, linestyle='--')
    ax1.axhline(y=sigma_y_MPa, color='black', linestyle=':', label="Material Yield Strength ($\sigma_y$)", linewidth=1.5)
    
    ax1.axvline(x=0, color='gray', linestyle='-', alpha=0.5)
    ax1.text(0.005, sigma_y_MPa * 1.1, "Throat (x=0)", color='gray', fontweight='bold')

    ax1.set_ylabel("Equivalent Stress [MPa]", fontsize=11)
    ax1.set_title("Von Mises Stress Distribution vs. Yield Limit", fontsize=12, loc='left')
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper right')

    mos = structural.margin_of_safety
    
    ax2.axhline(y=0.0, color='black', linestyle='-', linewidth=1.2)
    
    ax2.plot(cfd.x, mos, color='tab:blue', linewidth=2, label="Calculated MoS")
    ax2.fill_between(cfd.x, mos, 0, where=(mos >= 0), facecolor='tab:green', alpha=0.3, label="Safe Region (MoS > 0)")
    ax2.fill_between(cfd.x, mos, 0, where=(mos < 0), facecolor='tab:red', alpha=0.3, label="Yielding Region (MoS < 0)")
    
    ax2.axvline(x=0, color='gray', linestyle='-', alpha=0.5)
    
    ax2.set_xlabel("Nozzle Axis Position (x) [m]", fontsize=11)
    ax2.set_ylabel("Margin of Safety (MoS) [-]", fontsize=11)
    ax2.set_title("Structural Margin of Safety Along Nozzle Profile", fontsize=12, loc='left')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='lower right')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    generate_multiphysics_plots()