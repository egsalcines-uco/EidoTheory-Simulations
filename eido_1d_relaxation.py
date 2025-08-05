# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Fundamental Matrices ---
I = np.identity(2)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

# --- 2. Auxiliary Functions ---
def su2_exp(angle, axis):
    """
    Creates an SU(2) matrix from a rotation angle and an axis vector.
    U = exp(-i * (angle/2) * (axis . sigma))
    """
    if np.linalg.norm(axis) == 0:
        return I
    axis = axis / np.linalg.norm(axis)
    axis_sigma = axis[0]*sigma_x + axis[1]*sigma_y + axis[2]*sigma_z
    return np.cos(angle/2) * I - 1j * np.sin(angle/2) * axis_sigma

def project_to_su2(matrix):
    """
    Projects a 2x2 matrix back onto the SU(2) group.
    Approximation by normalization.
    """
    det = np.linalg.det(matrix)
    if det == 0:
        return I
    return matrix / np.sqrt(det)

def get_n_charge(U_field):
    """
    Calculates the winding number (n) of a 1D SU(2) matrix field robustly.
    Sums the phase change of the (0,0) element of each neighboring matrix.
    """
    charge_sum = 0
    for i in range(len(U_field)):
        prev_i = (i - 1 + len(U_field)) % len(U_field)
        
        # Get the phase of the (0,0) elements of the neighboring matrices
        phase_prev = np.angle(U_field[prev_i, 0, 0])
        phase_current = np.angle(U_field[i, 0, 0])
        
        # Calculate the phase difference, ensuring it's in [-pi, pi]
        d_phase = phase_current - phase_prev
        if d_phase > np.pi:
            d_phase -= 2 * np.pi
        elif d_phase < -np.pi:
            d_phase += 2 * np.pi
            
        charge_sum += d_phase
        
    # The winding number is the total sum of phase changes divided by 2pi
    return charge_sum / (2 * np.pi)

# --- 3. Simulation Parameters ---
Nx = 100       
num_steps = 5000  
dt = 0.1       

# --- 4. Field Initialization and n=1 Vortex ---
U = np.zeros((Nx, 2, 2), dtype=complex)
x = np.linspace(-10, 10, Nx)
n_history = []
angle_history = []

print("Initializing field with n=1 vortex...")
vortex_axis = np.array([0, 0, 1])
for i in range(Nx):
    # CORRECTION: The angle must go from 0 to 4pi for the SU(2) exponent
    # to produce a 2pi phase winding (an n=1 vortex).
    theta_i = 4 * np.pi * i / Nx
    U[i] = su2_exp(theta_i, vortex_axis)

initial_n_charge = get_n_charge(U)
print(f"Initial winding number: {initial_n_charge:.4f}")

initial_angles = np.array([np.angle(U[i, 0, 0]) for i in range(Nx)])

# --- 5. Relaxation Loop (Energy Dynamics) ---
print("Starting relaxation simulation...")
for step in range(num_steps):
    U_new = np.zeros_like(U)
    for i in range(Nx):
        prev_i = (i - 1 + Nx) % Nx
        next_i = (i + 1) % Nx
        
        force_term = U[prev_i] + U[next_i]
        
        U_new[i] = U[i] + dt * force_term
        
        U_new[i] = project_to_su2(U_new[i])

    U = U_new.copy()
    
    n_history.append(get_n_charge(U))

print("Simulation completed.")

final_n_charge = n_history[-1]
print(f"Final winding number: {final_n_charge:.4f}")

final_angles = np.array([np.angle(U[i, 0, 0]) for i in range(Nx)])

# --- 6. Results Visualization ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot the reference lines first
ax1.axhline(y=1, color='r', linestyle='--', label='n=1 (Stable)')
ax1.axhline(y=-1, color='r', linestyle='--', label='n=-1 (Stable)')
# Then plot the blue line with a higher zorder to be on top
ax1.plot(np.arange(num_steps) * dt, n_history, label='n(t)', zorder=10)

ax1.set_title('Evolution of Winding Number (n)')
ax1.set_xlabel('Relaxation Time')
ax1.set_ylabel('n')
ax1.grid(True)
ax1.legend()

# Plot the final configuration first (solid line)
ax2.plot(x, final_angles, label='Final Configuration', color='b', zorder=5) 
# Plot the initial configuration after, with a different color (red dashed)
# This ensures it's visible even if it overlaps
ax2.plot(x, initial_angles, label='Initial Configuration', color='r', linestyle='--', zorder=10) 

ax2.set_title('Angular Field Configuration (Phase of element U[0,0])')
ax2.set_xlabel('Position x')
ax2.set_ylabel('Phase')
ax2.legend()
ax2.grid(True)

# Adjust space between subplots
plt.subplots_adjust(hspace=0.4)

plt.tight_layout()
plt.show()