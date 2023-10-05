
# Multivariate MAPIT Documentation <img src="man/figures/logo.png" align="right" alt="" width="120"/>

[![R CMD check](https://github.com/lcrawlab/mvMAPIT/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/lcrawlab/mvMAPIT/actions/workflows/check-standard.yaml)
[![Docker Image CI](https://github.com/lcrawlab/mvMAPIT/actions/workflows/docker-image.yml/badge.svg)](https://github.com/lcrawlab/mvMAPIT/actions/workflows/docker-image.yml)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/mvMAPIT)](https://cranlogs.r-pkg.org/badges/grand-total/mvMAPIT)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mvMAPIT)](https://cran.r-project.org/package=mvMAPIT)

Find the full package documentation including examples and articles here: [Multivariate MAPIT Documentation](https://lcrawlab.github.io/mvMAPIT/).


## The multivariate MArginal ePIstasis Test (mvMAPIT)
This R package is a generalization of the [MAPIT
implementation](https://github.com/lorinanthony/MAPIT) by Crawford et
al. (2017)[^1] for any number of traits as described by Stamp et al. (2023)[^2].
The univariate MAPIT test for marginal epistasis is implemented as the special
case of running multivariate MAPIT with a single trait.

mvMAPIT is implemented as a set of R and C++ routines, which can be
carried out within an R environment.

### Introduction
Epistasis, commonly defined as the interaction between genetic loci, is known to
play an important role in the phenotypic variation of complex traits. As a
result, many statistical methods have been developed to identify genetic variants
that are involved in epistasis, and nearly all of these approaches carry out
this task by focusing on analyzing one trait at a time. However, because of the
large combinatorial search space of interactions, most epistasis mapping
methods face enormous computational challenges and often suffer from low
statistical power.

Previous studies have shown that jointly modeling multiple phenotypes can often
dramatically increase statistical power for association mapping. Therefore, here
we present the **multivariate MArginal ePIstasis Test (mvMAPIT)** – a
multi-outcome generalization of a recently proposed epistatic detection method
which seeks to detect *marginal epistasis* or the combined pairwise interaction
effects between a given variant and all other variants. By searching for marginal
epistatic effects, one can identify genetic variants that are involved in
epistasis without the need to identify the exact partners with which the variants
interact – thus, potentially alleviating much of the statistical and computational
burden associated with conventional explicit search based methods. Our proposed
mvMAPIT builds upon this strategy by leveraging correlation structures between
traits to improve the identification of variants involved in epistasis. We
formulate mvMAPIT as a multivariate linear mixed model and develop a multi-trait
variance component estimation algorithm for efficient parameter inference and
*P*-value computation. Together with reasonable model approximations, our proposed
approach is scalable to moderately sized GWA studies.


### The Method
The **multivariate MArginal ePIstasis Test** is a multi-outcome extension of the
statistical framework MAPIT which aims to identify variants that are involved in
epistatic interactions by leveraging the correlation structure of non-additive
genetic variation that is shared between multiple traits. The key idea behind the
concept of marginal epistasis is to identify variants that are involved in
epistasis while avoiding the need to explicitly conduct an exhaustive search over
all possible pairwise interactions. As an overview of mvMAPIT and its
corresponding software implementation, we will assume that we have access to a
GWA study on `N` individuals denoted as `D = {X,Y}` where `X` is an `N x J` matrix
of genotypes with `J` denoting the number of SNPs (each of which is encoded as
`{0,1,2}` copies of a reference allele at each locus `j`) and `Y` denoting a `N x D`
matrix holding `D` different traits that are measured for each of the `N`
individuals.

The goal of mvMAPIT is to identify variants that have non-zero interaction effects
with any other variant in the data. To accomplish this, we examine each SNP in
turn and assess the null hypothesis that its corresponding variance component is zero. In
practice, we use a computationally efficient method of moments algorithm called MQS from Zhou (2017)[^3]
to estimate model parameters and to carry out calibrated statistical tests within
mvMAPIT.

## Installation

The package needs compilation but the released version can be installed from
CRAN.

```R
install.packages("mvMAPIT")
```

### The R Environment
R is a widely used, free, and open source software environment for
statistical computing and graphics. The most recent version of R can be
downloaded from the [Comprehensive R Archive Network
(CRAN)](https://cran.r-project.org/). CRAN provides precompiled binary
versions of R for Windows, macOS, and select Linux distributions that
are likely sufficient for many users' needs. Users can also install R
from source code; however, this may require a significant amount of
effort. For specific details on how to compile, install, and manage R
and R-packages, refer to the manual [R Installation and
Administration](https://cran.r-project.org/doc/manuals/r-release/R-admin.html).

### R Packages Required for mvMAPIT
mvMAPIT requires the installation of the following R libraries:

- [checkmate](https://cran.r-project.org/package=checkmate)
- [CompQuadForm](https://cran.r-project.org/package=CompQuadForm)
- [dplyr](https://cran.r-project.org/package=dplyr)
- [foreach](https://cran.r-project.org/package=foreach)
- [harmonicmeanp](https://cran.r-project.org/package=harmonicmeanp)
- [logging](https://cran.r-project.org/package=logging)
- [mvtnorm](https://cran.r-project.org/package=mvtnorm)
- [Rcpp](https://cran.r-project.org/package=Rcpp)
- [RcppAlgos](https://cran.r-project.org/package=RcppAlgos)
- [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo)
- [RcppParallel](https://cran.r-project.org/package=RcppParallel)
- [RcppSpdlog](https://cran.r-project.org/package=RcppSpdlog)
- [tidyr](https://cran.r-project.org/package=tidyr)

The easiest method to install these packages is with the following
example command entered in an R shell:

``` {.R}
install.packages(c( 'checkmate', 
                    'CompQuadForm', 
                    'dplyr', 
                    'foreach', 
                    'harmonicmeanp', 
                    'logging', 
                    'mvtnorm', 
                    'Rcpp', 
                    'RcppAlgos', 
                    'RcppArmadillo', 
                    'RcppParallel', 
                    'RcppSpdlog', 
                    'testthat', 
                    'tidyr'), 
                    dependencies = TRUE);
```

Alternatively, one can also [install R packages from the
command-line](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

### Installing mvMAPIT from Sources
The easiest way to install the package from sources is to change into
the directory of mvMAPIT and run `R CMD INSTALL . --preclean`. The
`--preclean` flag makes sure that the latest state is run.

### C++ Functions Required for MAPIT

The code in this repository assumes that basic Fortran and C++ libraries and compilers are already set up on the running personal computer or
cluster. If not, the mvMAPIT functions and necessary Rcpp packages will
not work properly. A simple option is to use
[gcc](https://gcc.gnu.org/). macOS users may use this collection by
installing the [Homebrew package manager](https://brew.sh/index.html) and
then typing the following into the terminal:

``` {.bash}
brew install gcc
```
### OpenMP
Note that mvMAPIT takes advantage of [OpenMP](https://www.openmp.org/), an
API for multi-platform shared-memory parallel programming in C/C++. This
is to speed up the computational time of the modeling algorithm.
Unfortunately, macOS does not currently support OpenMP under the default
compiler. A work around to use OpenMP in R on macOS can be found
[here](https://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/).
mvMAPIT can be compiled without OpenMP, but we recommend using it if
applicable for scalability.

### Known Issues
- When your compiler changes, some R package dependencies might need to be recompiled. This is likely the case if the compilation error explicitly names an R package in the local library.

- On macOS, you might need to run `brew reinstall z3` to fix `'libz3.4.11.dylib' (no such file)` related errors ([clang issues](https://github.com/Homebrew/discussions/discussions/3920)).

- For extra tips on how to run C++ on macOS, please visit
<https://seananderson.ca/2013/11/18/rcpp-mavericks.html>.

- For tips on how to avoid errors dealing with `-lgfortran` or `-lquadmath`, please visit
<https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/>.

------------------------------------------------------------------------

## Questions and Feedback
For questions or concerns with the MAPIT functions, please contact
[Lorin Crawford](mailto:lcrawford@microsoft.com) or
[Julian Stamp](mailto:julian_stamp@brown.edu).

We appreciate any feedback you may have with our repository and instructions.

## References
[^1]: L. Crawford, P. Zeng, S. Mukherjee, X. Zhou (2017). Detecting
    epistasis with the marginal epistasis test in genetic mapping
    studies of quantitative traits. *PLoS Genet*. **13**(7): e1006869.
    <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006869>

[^2]: J. Stamp, A. DenAdel, D. Weinreich, L. Crawford (2023). Leveraging the
    Genetic Correlation between Traits Improves the Detection of Epistasis in
    Genome-wide Association Studies. *G3 Genes|Genomes|Genetics*, **13**(8), jkad118. doi:
    <https://doi.org/10.1093/g3journal/jkad118>

[^3]: X. Zhou (2017). A unified framework for variance component estimation with summary statistics
  in genome-wide association studies. *Ann Appl Stat*. **11**(4): 2027-2051.
  <https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/issue-4/A-unified-framework-for-variance-component-estimation-with-summary-statistics/10.1214/17-AOAS1052.full>
