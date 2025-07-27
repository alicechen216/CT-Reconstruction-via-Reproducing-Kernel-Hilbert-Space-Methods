# CT Reconstruction via Reproducing Kernel Hilbert Space Methods

This repository reimplements the numerical experiments from the paper:

> **"Regularization in a functional reproducing kernel Hilbert space"**  
> by Ruoming Wang and Yuesheng Xu, *Journal of Complexity*, Vol. 66, 2021, 101567.  
> [DOI: 10.1016/j.jco.2021.101567](https://doi.org/10.1016/j.jco.2021.101567)
> 
![Result](comparison.png)

We focus on applying kernel-based regularization methods to reconstruct images from limited-angle Radon transform measurements using:

- **FR** (Functional Reproducing kernel method)
- **RR** (Reproducing kernel method with Gaussian basis)
- **LR-G** (L2-regularized Gaussian basis method)

## 📁 Project Structure

```
.
├── main_experiment.m                 # Main script to run the experiment
├── calculate_metrics.m              # MSE / PSNR evaluation
├── generate_phantoms.m              # Phantom image generator (Shepp-Logan, bullseye, etc.)
├── reconstruct_FR.m                 # FR method reconstruction
├── reconstruct_RR.m                 # RR method reconstruction
├── reconstruct_LR_G.m               # LR-G method reconstruction
├── radon_transform_FR_kernel.m      # Analytic Radon transform of FR kernel
├── radon_transform_Gaussian_basis.m# Analytic Radon transform of Gaussian basis
├── evaluate_FR_kernel_on_grid.m     # Evaluation of FR kernel on pixel grid
└── README.md                        # This file
```

## 🧠 Methods Overview

Each method solves a regularized linear inverse problem using closed-form or analytic expressions for Radon transforms of kernel functions. The reconstruction is performed over a 2D domain discretized into a pixel grid.

- **FR** solves:  
  \[(F_R + \lambda I)\mathbf{c} = \mathbf{y}\]
  where \(F_R\) is computed analytically.

- **RR** and **LR-G** use Gaussian kernels with differing inner product structures.

All methods generate reconstructions over a grid using the learned coefficients and precomputed kernel evaluations.

## ▶️ Running the Experiment

1. Open MATLAB and navigate to the project folder.
2. Run:
   ```matlab
   main_experiment
   ```

This will:
- Generate a phantom image (e.g., Shepp-Logan)
- Generate synthetic Radon data with random projection directions
- Add optional Gaussian noise
- Reconstruct the image using all three methods
- Display and compute metrics (MSE and PSNR)

## 🖼 Example Output

The following figure is generated showing original and reconstructed images:

```
+----------------+--------------------+
| Original       | FR (PSNR: XX dB)   |
|                |                    |
+----------------+--------------------+
| RR (PSNR: XX)  | LR-G (PSNR: XX)    |
+----------------+--------------------+
```

> ⚠️ Note: Due to changes in resolution and noise levels, your results may differ from the paper unless hyperparameters are retuned.

## 🧪 Parameters and Tuning

- Resolution: `img_size = 256` (set to `512` for paper-quality results)
- Number of Radon samples: `m = 4096` (set to `16384` for paper-quality results)
- Regularization parameters (in `main_experiment.m`):
  ```matlab
  params_FR.lambda = 3.2258e-5;
  params_FR.delta = 3.5484;
  params_FR.gamma = 0.0007;

  params_RR.lambda = 1.0e-6;
  params_RR.gamma = 0.0020;

  params_LRG.lambda = 9.6865;
  params_LRG.gamma = 0.0010;
  ```

> You are encouraged to experiment with these values if using smaller images or fewer projections.

## 🛠 Requirements

- MATLAB R2018b or later
- No toolboxes beyond base MATLAB
