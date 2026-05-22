# AD-BOMA: Amplitude-Dependent Bayesian Operational Modal Analysis

**Amplitude-Dependent Bayesian Operational Modal Analysis for Tall Buildings under Typhoon Excitation: Field Validation and Identifiability Limits**

Jun Hu and Heung-fai Lam

*Submitted to Mechanical Systems and Signal Processing*

---

## Overview

This repository provides the code and processed data accompanying the paper on **AD-BOMA**, a Bayesian OMA framework that embeds parametric amplitude–frequency and amplitude–damping relationships directly into the spectral-matrix likelihood.

Key components:

| Module | Description |
|--------|-------------|
| `src/ad_boma.py` | AD-BOMA likelihood, power-law and log-linear amplitude models |
| `src/two_step.py` | Two-step estimation pipeline (per-window peak ID + curve fitting) |
| `src/sat_mcmc.py` | Sequential Adaptive Transport MCMC with staleness monitoring |
| `src/multi_lorentzian.py` | Trace PSD multi-Lorentzian fitting for closely-spaced modes |
| `src/damage_detection.py` | Baseline-only damage detection (Methods A/B/B'/C) |
| `src/plotting.py` | Figure generation utilities |
| `data/synthetic/` | Synthetic benchmark generator (6-mode, 15-channel) |
| `data/field/` | Anonymised processed modal parameters from the Typhoon Mangkhut 2018 field study |

## Key Results

- **Baseline-only fitting** achieves 100% detection at 5% FPR (200 MC repetitions), matching full AD-BOMA — the fitting strategy, not the inference method, is decisive.
- **αγ non-identifiability**: the power-law curvature exponent γ is practically unidentifiable from single-event data (~1.9 amplitude decades).
- **Closely-spaced modes** (35–40 mHz separation) recovered to within 0.3% via multi-Lorentzian trace PSD fitting.
- Full 141-window, 6-mode pipeline runs in ~3 min (two-step) or ~8.3 h (full Bayesian) on an RTX 4090.

## Data Availability

The **synthetic data generator** is fully open. The **field monitoring data** cannot be shared publicly due to confidentiality agreements with the building owner. Processed modal parameter estimates (per-window frequencies, damping ratios, and RMS amplitudes) are included in `data/field/` in anonymised form using window indices rather than site-identifying metadata.

## Installation

```bash
git clone https://github.com/junhu22/SAT-MCMC.git
cd SAT-MCMC
conda create -n sat-mcmc python=3.10
conda activate sat-mcmc
pip install -r requirements.txt
```

### Dependencies

- Python ≥ 3.10
- PyTorch ≥ 2.0 (CUDA recommended)
- NumPy, SciPy, Matplotlib
- nflows (normalizing flow library)

## Quick Start

```python
# Generate synthetic data
python src/generate_synthetic.py --modes 6 --windows 40 --amplitude-range 8 633

# Run two-step estimation
python src/two_step.py --data data/synthetic/benchmark_6mode.npz

# Run full SAT-MCMC (GPU recommended)
python src/sat_mcmc.py --data data/field/mangkhut_processed.npz --prior physical --device cuda
```

## Reproducing Paper Figures

```bash
# All synthetic benchmark figures (Figs 1–5)
python notebooks/run_synthetic_benchmarks.py

# Field study figures (Figs 6–13)
python notebooks/run_field_study.py

# Damage detection comparison (Fig 3, Table 5)
python notebooks/run_damage_detection.py --mc-repetitions 200
```

## Citation

```bibtex
@article{hu2026adboma,
  title   = {Amplitude-Dependent {Bayesian} Operational Modal Analysis for
             Tall Buildings under Typhoon Excitation: Field Validation and
             Identifiability Limits},
  author  = {Hu, Jun and Lam, Heung-fai},
  journal = {Mechanical Systems and Signal Processing},
  year    = {2026},
  note    = {Submitted}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Contact

Jun Hu — junhu22@whut.edu.cn — ORCID: [0000-0003-2104-2886](https://orcid.org/0000-0003-2104-2886)
