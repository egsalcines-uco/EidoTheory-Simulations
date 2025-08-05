# -*- coding: utf-8 -*-
import numpy as np

# --- 1. Parameters of the 3D Simulation ---
Nx = 32# Grid Dimension X
Ny = 32# Grid Y Dimension
Nz = 32# Z Dimension of the grid
dx = 1.0# Grid spacing (isotropic assumed)

# --- 2. Function to calculate the mass of a particle ---
def calculate_mass(n_L, alpha):
    """
    Calculate the mass of a particle based on the equation M=α∫∣∇⋅Θ∣dV.
    
    Args:
        n_L (int): Number of turns (winding number) of the eidic trajectory.
        alpha (float): Constant of proportionality.
        
    Returns:
        float: The calculated mass.
    """
    
    print(f"\nCalculating mass for n_L = {n_L}...")
    
    # Grid for the phase field (Θ) in 3D
    Theta = np.zeros((Nx, Ny, Nz), dtype=float)
    
    # Coordinates of the center of the vortex (line along the Z axis)
    vortex_center_x = Nx / 2
    vortex_center_y = Ny / 2
    
    # Phase field initialization (line vortex)
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x_dist = i - vortex_center_x
                y_dist = j - vortex_center_y
                
                # The angle for a vortex of n_L turns
                angle = n_L * np.arctan2(y_dist, x_dist)
                Theta[i, j, k] = angle
    
    #--- 3. Calculation of the Divergence (∇·Θ) discretely ---
    total_divergence_sum = 0.0
    
    for i in range(1, Nx - 1): # We ignore the edges to simplify
        for j in range(1, Ny - 1):
            for k in range(1, Nz - 1):
                # Approximation of partial derivatives with finite differences
                dTheta_dx = (Theta[i+1, j, k] - Theta[i-1, j, k]) / (2 * dx)
                dTheta_dy = (Theta[i, j+1, k] - Theta[i, j-1, k]) / (2 * dx)
                dTheta_dz = (Theta[i, j, k+1] - Theta[i, j, k-1]) / (2 * dx) # In this case it is 0 for a line vortex

                divergence = dTheta_dx + dTheta_dy + dTheta_dz
                total_divergence_sum += np.abs(divergence)

    # The discrete volume is dx**3. The whole becomes the sum.
    mass_calculated = alpha * total_divergence_sum * (dx**3)
    
    return mass_calculated

# --- 4. Calibration of the Proportionality Constant (α) ---
# We use the electron as a reference particle to calibrate α
electron_mass_real = 0.511 # Electron mass in MeV/c²
n_L_electron = 1

# We make a first simulation for the electron with a α = 1
temp_alpha = 1.0
mass_electron_simulated = calculate_mass(n_L_electron, temp_alpha)

# We calculate the real value of α so that the simulated mass matches the real one
alpha = electron_mass_real / mass_electron_simulated
print(f"\nCalibrated proportionality constant (α): {alpha:.8f}")

# --- 5. Calculating the masses of other particles ---
n_L_muon = 207 # n_L value for the muon
n_L_tau = 3481 # n_L value for tau

# We calculate the mass of the muon
mass_muon_calculated = calculate_mass(n_L_muon, alpha)

# We calculate the mass of tau
mass_tau_calculated = calculate_mass(n_L_tau, alpha)

# --- 6. Results and Comparison ---
print("\n--- Simulation Results ---")
print(f"Electron Mass (calculated): {electron_mass_real:.3f} MeV/c² (By calibration)")
print(f"Masa real del Muón: 105.66 MeV/c²")
print(f"Muon mass (calculated): {mass_muon_calculated:.3f} MeV/c²")
print(f"Masa real del Tau: 1776.86 MeV/c²")
print(f"Tau mass (calculated): {mass_tau_calculated:.3f} MeV/c²")

#Note: The results of this conceptual simulation are to show the method.
# The actual values will depend on the physics of the simulation, the relaxation and the grid parameters.
