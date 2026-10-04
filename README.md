# Neural Network Approximation of eMPC on ESP32

Bachelor's thesis project: a neural network approximating the explicit MPC
control law of an inverted pendulum on a cart with Robust Tube MPC, deployed on an ESP32 for
real-time execution.

## Results

<p align="center">
  <img width="700" alt="sim_nn_plant_dist" src="https://github.com/user-attachments/assets/19a19334-0b79-4535-ab6c-b430433272a0" />
  <br/>
  <em>Closed-loop simulation of the inverted pendulum on a cart: plant state under the NN-approximated tube MPC, with disturbance.</em>
</p>

<p align="center">
  <img width="700" alt="NN on ESP32" src="https://github.com/user-attachments/assets/09ce9747-5e93-401c-9323-a70eff5cd6ff" />
  <br/>
  <em>Neural network controller deployed on the ESP32, running in the loop.</em>
</p>

<p align="center">
  <img width="700" alt="Benchmark timing" src="https://github.com/user-attachments/assets/75da46c0-143e-408c-89de-e5a4ec740aed" />
  <br/>
  <em>Benchmark of per-step computation time on the ESP32.</em>
</p>

Demo video: [https://youtu.be/v6TLr3-HuwE]

## Structure
- `matlab/` – controller design, NN training, plots
- `firmware/pendulum_nn/` – ESP32 sketch and exported NN weights
- `python/` – plant simulation, perturbation tests, data plotting
- `data/` – logged microcontroller data
- `docs/` – thesis and presentation

## Stack
MATLAB, C/C++ (ESP32), Python
