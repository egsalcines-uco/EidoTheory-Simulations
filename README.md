# EidoTheory-Simulations

This repository contains the Python simulation code accompanying the paper "From Plausibility to Predictability: Angular Field Dynamics and the Lepton Mass Hierarchy in Eido Theory" by Enrique García Salcines.

This work explores a new physical ontology based on a fundamental field of discrete angular orientations (the Eido field). The simulations here serve as a proof-of-concept for the theory's key hypotheses:
1.  **Topological Stability:** That elementary particles, viewed as topological configurations (eidos), are robust and stable.
2.  **Quantitative Predictability:** That the mass of leptons can be predicted from the topological properties of the Eido field.

## Project Structure

* `eido_1d_relaxation.py`: A simulation that demonstrates the stability of a topological vortex (an eido) in a simplified one-dimensional angular field. It shows that the winding number `n` is a conserved quantity.
* `eido_3d_mass_prediction.py`: A conceptual simulation that calculates the mass hierarchy of leptons (muon and tau) by approximating a volume integral over a 3D angular field. It uses the electron's mass to calibrate a coupling constant and successfully predicts the other masses with high accuracy.

## Key Results from the Simulations

* The 1D simulation confirms that non-trivial configurations of the Eido field are stable under a relaxation dynamic.
* The 3D simulation, using the electron mass as a single calibration point, predicts the masses of the muon and tau with an accuracy greater than 99.8%.
* The results suggest a fundamental relationship: $M \approx n_L \cdot M_e$, where $M$ is the particle mass, $M_e$ is the electron mass, and $n_L$ is a topological winding number.

## Getting Started

### Prerequisites

You need a working Python environment with the following libraries installed:
* `numpy`
* `matplotlib`

You can install them using pip:
```bash
pip install numpy matplotlib
```
How to Run
  1. Clone this repository to your local machine:
     ```bash
     git clone https://github.com/[Tu_Usuario]/EidoTheory-Simulations.git
     ```
  2. Navigate to the repository directory.
  3. To run the 1D simulation:
     ```bash
     python eido_1d_relaxation.py
     ```
  4.To run the 3D mass prediction simulation:
    ```bash
    python eido_3d_mass_prediction.py
    ```
### License 
This project is licensed under the MIT License - see the LICENSE file for details.

### Contact
For any questions or feedback regarding the code or the theory, please contact Enrique García Salcines at egsalcines@uco.es.
