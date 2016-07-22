\name{phylodyn-package}
\alias{phylodyn-package}
\alias{phylodyn}
\docType{package}
\title{Tools for phylodynamics and population genetics}
\description{
Tools for analyzing genealogy and genome sequence data using INLA and MCMC.
}
\details{
\tabular{ll}{
Package: \tab phylodyn\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-12\cr
License: \tab GPL-3\cr
}
An overview of how to use the package ***
}
\author{
Michael Karcher
Julia Palacios
Shiwei Lan
Vladimir Minin

Maintainer: Michael Karcher <mkarcher@uw.edu>
}
\references{
1. A. Rambaut, O. G. Pybus, M. I. Nelson, C. Viboud, J. K. Taubenberger, E. C. Holmes
[The genomic and epidemiological dynamics of human influenza A
virus](http://www.nature.com/doifinder/10.1038/nature06945).
*Nature*, 453(7195): 615-619, 2008.

2. J. A. Palacios and V. N. Minin.
[Integrated nested Laplace approximation for Bayesian nonparametric phylodynamics](http://www.auai.org/uai2012/papers/310.pdf).
In *Proceedings of the Twenty-Eighth International Conference in Uncertainty in Artificial Intelligence*, pages 726-735, 2012.

3. D. Zinder, T. Bedford, E. B. Baskerville, R. J. Woods, M. Roy, M. Pascual.
[Seasonality in the migration and establishment of H3N2 Influenza lineages with epidemic growth and decline](http://bmcevolbiol.biomedcentral.com/articles/10.1186/s12862-014-0272-2).
*BMC Evolutionary Biology*, 14(1): 272, 2014.

4. S. Lan, J. A. Palacios, M. Karcher, V. N. Minin, and B. Shahbaba
[An Efficient Bayesian Inference Framework for Coalescent-Based Nonparametric Phylodynamics](http://bioinformatics.oxfordjournals.org/content/31/20/3282),
*Bioinformatics*, 31(20): 3282-3289, 2015.

5. M. D. Karcher, J. A. Palacios, T. Bedford, M. A. Suchard, and V. N. Minin.
[Quantifying and mitigating the effect of preferential sampling on phylodynamic inference](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004789).
*PLOS Computational Biology*, 12:e1004789, 2016.

6. J.A Palacios, J. Wakeley,  and S. Ramachandran. [Bayesian nonparametric inference of population size changes from sequential genealogies.](http://www.genetics.org/content/early/2015/07/28/genetics.115.177980) *Genetics* Vol. 201:281-304, 2015.
}

\keyword{ package }
\seealso{
\code{\link[INLA:INLA-package]{INLA}}
}
\examples{
library(phylodyn)

traj = exp_traj
gene = coalsim(samp_times = 0, n_sampled = 100, traj = traj)
res_BNPR    = BNPR(data = gene)
res_BNPR_PS = BNPR_PS(data = gene)

par(mfrow=c(1,2))
plot_BNPR(res_BNPR, traj = traj, main="BNPR")
plot_BNPR(res_BNPR_PS, traj = traj, main="BNPR-PS")
}
