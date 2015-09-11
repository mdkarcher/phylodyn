phylodyn
========

The purpose of `phylodyn` is to facilitate phylodynamic inference and analysis in an approachable R package.

## Installation

1. Install (if necessary) package dependencies and helpers `ape`, `INLA`, `spam` and `devtools` using `install.packages`.

2. Load `devtools` using `library(devtools)`.

3. Install `phylodyn` using 

    a. `install_github("phylodyn", username="mdkarcher")`, or

    b. `install_github("phylodyn", username="mdkarcher", build_vignettes = TRUE)` if you want some illustrative vignettes (note: using `build_vignettes = TRUE` will make the install take longer).

## Vignettes

1. **SimpleBNPR**: A short example showing how to use BNPR and BNPR-PS on simulated data.

2. **SimplePhyloinfer**: A short example comparing BNPR with a split HMC MCMC sampler approach.

3. **LongPhyloinfer**: A longer example comparing BNPR with multiple MCMC samplers, including split HMC as in SimplePhyloinfer.

## References

1. J. A. Palacios and V. N. Minin.
Integrated nested Laplace approximation for Bayesian nonparametric phylodynamics.
In *Proceedings of the Twenty-Eighth International Conference on Uncertainty in Artificial Intelligence*, pages 726–735, 2012.

2. M. S. Gill, P. Lemey, N. R. Faria, A. Rambaut, B. Shapiro, and M. A. Suchard.
Improving Bayesian population dynamics inference: a coalescent-based model for multiple loci.
*Molecular biology and evolution*, 30(3):713–724, 2013.

3. Lan, S., Palacios, J. A., Karcher, M., Minin, V. N., & Shahbaba, B. (2014).
An Efficient Bayesian Inference Framework for Coalescent-Based Nonparametric Phylodynamics.
*arXiv preprint arXiv*:1412.0158.