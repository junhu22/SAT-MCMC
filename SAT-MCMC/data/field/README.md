# Field Data: Typhoon Mangkhut 2018

## Data Description

This directory contains **anonymised processed modal parameter estimates** from the
Typhoon Mangkhut 2018 field study. The raw accelerometer data cannot be shared
publicly due to confidentiality agreements with the building owner.

## Anonymisation

- All data are indexed by **window number** (0–140), not by timestamp or location.
- No site-identifying metadata (building name, address, coordinates) is included.
- Sensor IDs (TM38, etc.) are internal labels that do not identify the building.
- Mode frequencies, damping ratios, and RMS amplitudes are physical measurements
  that do not uniquely identify the structure.

## File Format

`mangkhut_processed.npz` contains:

| Key | Shape | Description |
|-----|-------|-------------|
| `window_index` | (141,) | Clean window indices (4 outliers removed from 145 total) |
| `rms_amplitude` | (141,) | 15-channel mean RMS acceleration (µg) |
| `freq_tx1` | (141,) | Per-window TX1 frequency (Hz) |
| `freq_ty1` | (141,) | Per-window TY1 frequency (Hz) |
| `freq_r1` | (141,) | Per-window R1 frequency (Hz) |
| `freq_tx2a` | (141,) | Per-window TX2a frequency (Hz) |
| `freq_tx2b` | (141,) | Per-window TX2b frequency (Hz) |
| `freq_ty2` | (141,) | Per-window TY2 frequency (Hz) |
| `damping_*` | (141,) | Corresponding damping ratios |
| `fit_params` | dict | Two-step power-law fit parameters per mode |
