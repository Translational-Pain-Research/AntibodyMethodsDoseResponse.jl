# AntibodyMethodsDoseResponse

## About

This package implements methods to analyze rate-constant distributions for (antibody-binding) dose-response curves, keeping the dependencies as minimal as possible. This assures the most general compatibility and allows to define the analysis from scratch.

!!! info "Convenience package"
	For most users, this minimal approach and the additional flexibility may not be very useful. Instead, the [`AntibodyMethodsDoseResponseConvenience`](@ref api_convenience) package should be preferred, as it pre-defines the recommended analyses and implements the necessary functions.

In essence, the dose-response data is modelled by the following accumulation model (see [https://arxiv.org/abs/2407.06052](https://arxiv.org/abs/2407.06052) for further details):
```math
r(a) = \int_0^\infty g(K_\tau) \left(1- e^{-\frac{a}{K_\tau}} \right) \ dK_\tau
```
 where ``K_\tau \sim \frac{1}{k_{\text{on}}}`` is proportionally inverse to the binding rate ``k_{\text{on}}`` and ``g(K_\tau)`` is the ``K_\tau``-density (amount) of epitopes that exhibit the rate constant ``K_\tau``. Alternatively, the Langmuir isotherm (equilibrium model) can be used. In either case, the analysis consists of a model fit, estimating the density ``g(K_\tau)`` from measurement data:

```@example index
using AntibodyMethodsDoseResponseConvenience, Measures # hide
results, data, replicates =  load_results("examples/default_concentrations") # hide
dr_plot, density_plot = bin_analysis_plot(results,data,replicates) # hide
plot(dr_plot, density_plot, layout = (1,2), size = (800,300), margins = 4mm) # hide
```

In this documentation, both the terms **density plot** and **histogram** will be used interchangeably to refer to the plot of the estimated density ``g(K_\tau)``.

## About the tutorials

* The **quick start guide** 
    + A minimal tutorial, covering only the necessary steps to analyze dose-response data without any in-depth explanations. 
    + From here, the [AnitbodyMethodsDoseResponseConvenience API](@ref api_convenience) is highly recommended to explore the additional options of the convenience functions that are not covered in the minimal tutorial.
    + The recommendation for most users.

* The **Detailed explanations**
    + Covers the background and the internals of this package. 
    + Assumes (at some points) to be read in order from top to bottom.
    + The convenience functions (from the quick start guide) use the same data types and methods internally.
    + Intended for developers that need access to the internal methods and that want to take full control over the analysis.

## Installation

The `AntibodyMethods` packages are not part of the official registry. Thus, the installation differs slightly from the usual package installation process.

### Simple installation


Open the Julia repl and press `]` to enter the package manger.

```raw

               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version ****
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |


pkg>
```

For a fresh install of Julia, make sure to add the official registry `General` first:

```raw
pkg> registry add General
```

Next, add the `AntibodyPackagesRegistry`:

```raw
pkg> registry add https://github.com/AntibodyPackages/AntibodyPackagesRegistry
```

Finally, add `AntibodyMethodsDoseResponseConvenience`, which install the necessary dependencies:

```raw
pkg> add AntibodyMethodsDoseResponseConvenience
```

### Manual installation

It is also possible to install the dependencies for `AntibodyMethodsDoseResponse` manually, with the following commands (exactly in the order below):

```julia
using Pkg
Pkg.add(url="https://github.com/AntibodyPackages/FittingObjectiveFunctions")
Pkg.add(url="https://github.com/AntibodyPackages/AdaptiveDensityApproximation")
Pkg.add(url="https://github.com/AntibodyPackages/AntibodyMethodsDoseResponse")
```

To install the `AntibodyMethodsDoseResponseConvenience` package manually, two additional packages need to be installed:

```julia
Pkg.add(url="https://github.com/AntibodyPackages/AntibodyMethodsDoseResponseRecipes")
Pkg.add(url="https://github.com/AntibodyPackages/AntibodyMethodsDoseResponseConvenience")
```

!!! info
	`AntibodyMethodsDoseResponseRecipes` contains plotting instructions for `Plots.jl`, allowing to plot the data structures defined/used in `AntibodyMethodsDoseResponse`.


## How to cite the package

If you would like to cite this package for scientific purposes, you might also want to cite the corresponding paper [https://arxiv.org/abs/2407.06052](https://arxiv.org/abs/2407.06052).
