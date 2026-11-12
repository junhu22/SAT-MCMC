# AD-BOMA: Amplitude-Dependent Bayesian Operational Modal Analysis

**Amplitude-Dependent Bayesian Operational Modal Analysis for Tall Buildings under Typhoon Excitation: Field Evidence and Identifiability Limits**

Jun Hu and Heung-fai Lam  
*Mechanical Systems and Signal Processing* (2026), article 115163  
DOI: https://doi.org/10.1016/j.ymssp.2026.115163

## Scope

This repository provides the supplementary comparison, calibration and
processed-data materials associated with the article. It preserves the
supporting files supplied with the accepted manuscript rather than an
earlier prototype implementation.

The materials include processed modal estimates, calibration and
comparison scripts, a comparison protocol and two sensitivity studies
for Window 144. The comparison protocol is not presented as a completed
raw-signal head-to-head benchmark against all established methods.

## Contents

- supplementary/comparison/README.md: package instructions and scope
- supplementary/comparison/VERSIONS.md: processed-data version definitions
- supplementary/comparison/COMPARISON_PROTOCOL.md: comparison design and metrics
- supplementary/comparison/results/: two versions of anonymised field estimates
- supplementary/comparison/baselines/: spectral identification utilities
- supplementary/comparison/.audit/: calibration/validity scripts and JSON outputs
- supplementary/window144_impact.json: Window 144 sensitivity
- supplementary/window144_routeA.json: Window 144 sensitivity for Route A

The comparison package retains the original folder structure and
file relationships submitted with the manuscript. Start with its
README and protocol. Scripts depend on a scientific Python environment;
NumPy, SciPy and pandas are used by common utilities.

## Scope of the scientific results

Operational field estimates in the article come from two-step spectral
fitting and amplitude-law regression. Baseline-only amplitude fitting
accounts for the synthetic damage-detection improvement. The plug-in
transport-MCMC runs did not converge and do not supply field posterior
estimates. A separately converged exact-joint field block was used
solely to assess amplitude-only model adequacy.

Raw accelerometer signals and building-identifying data are withheld
under the field-data confidentiality agreement. Only the processed,
anonymised modal estimates are released.

## Citation

Hu, J., and Lam, H.-f. (2026). Amplitude-Dependent Bayesian Operational
Modal Analysis for Tall Buildings under Typhoon Excitation: Field
Evidence and Identifiability Limits. Mechanical Systems and Signal
Processing, article 115163.
https://doi.org/10.1016/j.ymssp.2026.115163

## License

Repository code and documentation are provided under the MIT License.
Third-party dependencies retain their own licenses.
