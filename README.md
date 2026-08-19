NOZZLE SOLVER

A coupled quasi-1D multiphysics solver designed to evaluate the thermo-structural performance of liquid rocket engine nozzles. The framework integrates aerodynamic expansion, high-temperature heat transfer (Bartz relation), and structural mechanics (Von Mises yield criterion) to validate the nozzle integrity under realistic operational loads.

Project Overview

This project simulates a high-performance regeneratively cooled copper-alloy (Cu-Cr-Zr / GRCop-84) rocket nozzle operating at a chamber pressure of 70 bar and a combustion temperature of 3300 K.

The solver couples three core engineering disciplines sequentially:

Gas Dynamics (CFD): Solves quasi-1D isentropic flow equations to extract pressure, temperature, and Mach number profiles along the nozzle axis.
Thermal Analysis: Computes the gas-side heat transfer coefficient using the Bartz Relation, solving the radial thermal gradient across the nozzle wall.
Structural Mechanics: Evaluates equivalent Von Mises stresses on both the inner (gas-side) and outer (coolant-side) walls, accounting for combined pressure loads and thermal expansion strains.

Repository Structure

The architecture is modularized into 5 interconnected Python scripts:

propulsion_cfd.py: Defines the nozzle geometry and handles the aerodynamic flow solver.
thermal_bartz.py: Computes the convective heat transfer coefficients and wall temperature profiles.
RocketNozzleStructural.py: Performs the stress analysis and calculates the engineering Margin of Safety (MoS).
nozzle_solver.py: The main execution script that couples the physics modules together.
plot_multiphysics.py: Post-processing script that generates the visual data plots.

Material Properties & Input Data

The simulation utilizes properties baseline to advanced aerospace-grade copper liners:

**Young's Modulus (E): 110 GPa
**Coefficient of Thermal Expansion (alpha): 17*10^-6 K^(-1)
**Poisson's Ratio (nu): 0.33
**Yield Strength (sigma_y): 420 MPa
**Thermal Conductivity (k_w): 320 W/mK

Results

![Von Mises stress and Margin of Safety distribution along the nozzle axis](images/nozzle_stress_mos.png)
*Fig. 1 – Von Mises stress distribution (top) and structural Margin of Safety (bottom) along the nozzle axis. The vertical line marks the throat (x=0). Note the local stress minimum and MoS minimum at the throat, driven by the minimum-thickness constraint, and the sharp discontinuity past x≈0.05 m where the wall thickness switches back to its default value in the divergent section — an artifact of the thin-wall quasi-1D assumption, which does not capture the local stress concentration a real step-change in thickness would produce.*

Thermal vs. Mechanical Stress Coupling In RocketNozzleStructural.py, the solver evaluates the total stress by superimposing mechanical hoop stress (tensile) and thermal stress (clamped expansion). The simulation plots clearly show that the outer wall stress (sigma_vm_o) is higher than the inner wall stress (sigma_vm_i).
Mechanical pressure from cfd.P generates a uniform tensile stress across the wall. However, the high thermal gradient computed in thermal_bartz.py induces a massive expansion delta. The hot inner wall is constrained by the colder outer structure, inducing a compressive thermal stress. On the inner wall, these loads partially cancel out (sigma_{total} = sigma_{pressure} - |sigma_{thermal}|). On the outer wall, both pressure and thermal expansion act as tensile loads, reinforcing each other (sigma_{total} = sigma_{pressure} + sigma_{thermal}$), which pushes the equivalent Von Mises stress closer to the sigma_y limit of 420 MPa.

Throat Thickness Optimization At the throat (x = 0), the code in thermal_bartz.py implements a conditional piecewise function that aggressively drops the wall thickness t_w to a minimum of 1.2 mm.
The Bartz relation predicts a massive spike in the gas-side convective heat transfer coefficient (h_g) at the throat due to high mass flux. If the wall were thick, the conductive thermal resistance (R_{th} = t_w / k_w) would cause a catastrophic temperature delta across the material. By forcing t_w = 0.0012 at the throat, the code minimizes R_{th}. Even though a thinner wall provides less structural material to fight the 70 bar chamber pressure, reducing the thermal stress contribution outweighs the mechanical loss. This keeps the combined Von Mises stress well below the material yield strength, as shown by the stable valley in the stress plot at x = 0.

Divergent Section Over-design (CFD-Structural Mismatch) In the divergent section, the structural script reports a rapidly increasing Margin of Safety , eventually peaking past 7.0 near the nozzle exit.
In propulsion_cfd.py, the supersonic expansion logic causes the static pressure cfd.P to drop exponentially towards ambient atmospheric levels. Concurrently, the geometry loop defaults the wall thickness back to a static t_w = 20 mm in the else block. Because the driving mechanical load (gas pressure) crolls down to almost zero while the structural thickness remains a robust 2.0 mm, the nozzle becomes massively over-designed in the divergent section. A Margin of Safety 7 implies that the material is working at less than 15% of its structural capacity. This insight exposes a clear engineering optimization path: a variable, tapering wall thickness function could be implemented in the divergent section to shave off unnecessary copper mass, directly improving the rocket's payload capacity.
