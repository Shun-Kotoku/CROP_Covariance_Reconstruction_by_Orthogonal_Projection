# Reconstructing the information processing capacity of physical systems from noisy observations 
This repository contains the code for CROP (Covariance Reconstruction by Orthogonal Projection), a method for estimating noise-free information processing capacity (IPC) from noisy reservoir states. 

## Requirements 

- C++ compiler
- GNU make
- [Armadillo] (https://arma.sourceforge.net/)

### Armadillo 
Armadillo is required for the linear algebra operations used in this code. 
Please install Armadillo according to your operating system: https://arma.sourceforge.net/download.html

## Project Structure
This repository consists of four main files:

- **`MC.cpp`** — Entry point. Runs the memory capacity calculation for a single configuration (parameters are set in `Initialization.h`).
- **`Data.h`** — Defines the `Data` struct, which holds all state (matrices, file paths, parameters) passed between functions.
- **`Initialization.h`** — Configures run parameters and file paths (`setOutputFile`), and loads the input data files (`loadData`).
- **`Memory.h`** — Core algorithm: computes the IPC via Legendre polynomial reconstruction, and estimates the noise-free IPC using CROP.

## Input Data Format
Two binary files are required:

1. **Input (target) signal** — a column vector of length >= `inputLength`, used to compute the Legendre-polynomial reconstruction targets. This must be a one-dimensional, i.i.d. signal drawn from Uniform[-1, 1]. If you use a different input distribution (e.g. i.i.d. Gaussian), the orthogonal polynomial family in `Memory.h` must be changed accordingly (e.g. Hermite instead of Legendre).
2. **State matrix** — shape `(inputLength - startTargetRow) x (nodeNum + 1)`. The first column must be a constant bias term (all ones); the remaining `nodeNum` columns are the reservoir node states.

By default, the code expects:
- Input signal at: `<folderName>/input_1.01e6.bin`
- State matrix at: `<folderName>/stateData_<filenameBase>.bin`

## Adapting to Your Own Data
Edit `setOutputFile()` in `Initialization.h` to match your data:

- `data.folderName` — directory containing the input signal and state matrix file
- `data.inputfileName` — name of the input signal file
- `data.inputLength` — length of the input signal / number of rows used from the state matrix
- `data.startTargetRow` — number of initial rows to discard (transient). This must be larger than the maximum value in `stepsForMemory`, since earlier target values are indexed relative to this offset.
- `data.nodeNum` — number of reservoir nodes (excluding the bias column)
- `data.filenameBase` — naming pattern used to locate your state matrix file
- `data.stepsForMemory` — max delay steps to evaluate, for each polynomial degree
