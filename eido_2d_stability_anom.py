# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os

# --- 1. Fundamental Matrices ---
I = np.identity(2)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

# --- 2. Auxiliary Functions ---
def project_to_su2(matrix):
    det = np.linalg.det(matrix)
    if det == 0:
        return I
    return matrix / np.sqrt(det)

def calculate_spin_from_position_angle(cx, cy, r, Nx, Ny, n_L):
    """
    Calculates the spin based on the change in position angle.
    """
    cx = int(round(cx)); cy = int(round(cy))
    pts = []
    # Corrected contour loop
    for j in range(cy - r, cy + r + 1): pts.append(((cx - r) % Nx, j % Ny))
    for i in range(cx - r + 1, cx + r + 1): pts.append((i % Nx, (cy + r) % Ny))
    for j in range(cy + r - 1, cy - r - 1, -1): pts.append(((cx + r) % Nx, j % Ny))
    for i in range(cx + r - 1, cx - r - 1, -1): pts.append((i % Nx, (cy - r) % Ny))
    
    total_phase_change = 0.0
    for k in range(len(pts)):
        i0, j0 = pts[k]
        i1, j1 = pts[(k + 1) % len(pts)]
        
        phi0 = np.arctan2(j0 - cy, i0 - cx)
        phi1 = np.arctan2(j1 - cy, i1 - cx)
        
        d_phi = phi1 - phi0
        if d_phi > np.pi: d_phi -= 2 * np.pi
        elif d_phi < -np.pi: d_phi += 2 * np.pi
        
        total_phase_change += d_phi
    
    return n_L * total_phase_change / (2 * np.pi)

def compute_energy_xy_vectorized(U):
    U_ip1 = np.roll(U, -1, axis=0)
    U_jp1 = np.roll(U, -1, axis=1)
    
    term_x = np.trace(U_ip1 @ np.conj(U).transpose(0,1,3,2), axis1=2, axis2=3)
    term_y = np.trace(U_jp1 @ np.conj(U).transpose(0,1,3,2), axis1=2, axis2=3)
    
    return float((2*np.prod(U.shape[:2])) - np.real((term_x + term_y).sum()))

def plot_and_save_phase(U, title, filename):
    plt.figure(figsize=(8, 8))
    phase = np.angle(U[:, :, 0, 0])
    plt.imshow(phase, cmap='hsv', origin='lower')
    plt.colorbar(label='Phase')
    plt.title(title)
    plt.savefig(filename)
    plt.close()

# --- 3. Main Simulation Function ---
def run_simulation(name, n_L, damping_factor, num_steps, Nx, Ny, noise_amplitude):
    print(f"\n--- Starting simulation for {name} ---")

    # --- Field Initialization with n_L Vortex ---
    U = np.zeros((Nx, Ny, 2, 2), dtype=complex)
    cx, cy = Nx / 2 - 0.5, Ny / 2 - 0.5

    print(f"Initializing 2D field with vortex charge={n_L}...")
    i_idx, j_idx = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')
    x_dist, y_dist = i_idx - cx, j_idx - cy

    for i in range(Nx):
        for j in range(Ny):
            phi = np.arctan2(y_dist[i,j], x_dist[i,j])
            psi_real = np.cos(n_L * phi)
            psi_imag = np.sin(n_L * phi)
            U[i,j] = np.array([[psi_real - 1j*psi_imag, 0], [0, psi_real + 1j*psi_imag]])

    r2 = x_dist**2 + y_dist**2
    U[r2 < 1.0] = np.zeros((2, 2))
    
    # Saves the initial phase image
    plot_and_save_phase(U, f"Initial Phase ({name})", f"{name}_initial_phase.png")

    # --- Relaxation Loop ---
    print("Starting relaxation simulation...")
    energy_history = []
    
    initial_spin = calculate_spin_from_position_angle(cx, cy, r=min(Nx, Ny)//4, Nx=Nx, Ny=Ny, n_L=n_L)
    initial_energy = compute_energy_xy_vectorized(U)
    print(f"Step 0/{num_steps} | Energy: {initial_energy:.6f} | Spin: {initial_spin:.4f}")

    for step in range(1, num_steps + 1):
        U_left = np.roll(U, 1, axis=0)
        U_right = np.roll(U, -1, axis=0)
        U_up = np.roll(U, 1, axis=1)
        U_down = np.roll(U, -1, axis=1)
        
        neighbor_avg = (U_left + U_right + U_up + U_down)
        U_new = (1.0 - damping_factor) * U + damping_factor * (neighbor_avg/4.0)

        for i in range(Nx):
            for j in range(Ny):
                U_new[i, j] = project_to_su2(U_new[i, j])

        noise = np.exp(1j * np.random.uniform(-noise_amplitude, noise_amplitude, size=(Nx, Ny)))
        U = U_new * noise.reshape(Nx, Ny, 1, 1)

        if step > num_steps // 2:
            energy_history.append(compute_energy_xy_vectorized(U))
        
        if step % (num_steps / 10) == 0:
            energy = compute_energy_xy_vectorized(U)
            spin = calculate_spin_from_position_angle(cx, cy, r=min(Nx, Ny)//4, Nx=Nx, Ny=Ny, n_L=n_L)
            print(f"Step {step}/{num_steps} | Energy: {energy:.6f} | Spin: {spin:.4f}")

    # --- Final Results ---
    final_energy = compute_energy_xy_vectorized(U)
    final_spin = calculate_spin_from_position_angle(cx, cy, r=min(Nx, Ny)//4, Nx=Nx, Ny=Ny, n_L=n_L)
    residual_variance = np.var(energy_history)
    
    # Saves the final phase image
    plot_and_save_phase(U, f"Final Phase ({name})", f"{name}_final_phase.png")

    print("------------------------------")
    print(f"Simulation of {name} completed.")
    print(f"Final energy: {final_energy:.6f}")
    print(f"Spin: {final_spin:.4f}")
    print(f"Residual Energy Variance: {residual_variance:.6f}")

# --- 4. Main Execution Block ---
if __name__ == "__main__":
    leptons = [
        {"name": "Electron", "n_L": 1, "damping_factor": 0.05},
        {"name": "Muon", "n_L": 207, "damping_factor": 0.8},
        {"name": "Tau", "n_L": 3477, "damping_factor": 1.0},
    ]

    for lepton in leptons:
        run_simulation(
            name=lepton["name"],
            n_L=lepton["n_L"],
            damping_factor=lepton["damping_factor"],
            num_steps=50000,
            Nx=100,
            Ny=100,
            noise_amplitude=1e-4
        )
