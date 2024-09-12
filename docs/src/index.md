# AntibodyMethodsDoseResponse

## About

[`AntibodyMethodsDoseResponse.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponse.jl), [`AntibodyMethodsDoseResponseConvenience.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponseConvenience.jl) and [`AntibodyMethodsDoseResponseRecipes.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponseRecipes.jl) are [Julia](https://julialang.org/) packages for the analysis of (antibody-binding) dose-response curves.

In essence, the dose-response data is modelled by the following accumulation model (see [https://arxiv.org/abs/2407.06052](https://arxiv.org/abs/2407.06052) for further details and applications):
```math
r(a) = \int_0^\infty g(K_\tau) \left(1- e^{-\frac{a}{K_\tau}} \right) \ dK_\tau
```
where ``K_\tau \sim \frac{1}{k_{\text{on}}}`` is proportionally inverse to the binding rate ``k_{\text{on}}`` and ``g(K_\tau)`` is the ``K_\tau``-density of epitopes. The packages also implement the Langmuir isotherm (equilibrium model). In either case, the analysis consists of a model fit, estimating the density ``g(K_\tau)`` from measurement data.

```@example index
using AntibodyMethodsDoseResponseConvenience, Measures # hide
results, data, replicates =  load_results("examples/default_concentrations") # hide
dr_plot, density_plot = bin_analysis_plot(results,data,replicates) # hide
plot(dr_plot, density_plot, layout = (1,2), size = (800,300), margins = 4mm) # hide
```

In this documentation, both the terms **density plot** and **histogram** will be used interchangeably to refer to the plot of the estimated density ``g(K_\tau)``.

!!! info "The different packages"
    | Package | Short Description |
    | :------ | :-------- |
    | [`AntibodyMethodsDoseResponseConvenience.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponseConvenience.jl) | Convenience package for the analysis of dose-response curves. Contains predefined analysis and plotting methods. Recommended for the analysis of dose-response data. |
    | [`AntibodyMethodsDoseResponse.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponse.jl) | Minimal package defining the underlying models and methods for the analysis of dose-response curves. Intended for the development of analysis methods from scratch. Requires fewer dependencies. |
    | [`AntibodyMethodsDoseResponseRecipes.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponseRecipes.jl) |  [`Plots.jl`](https://docs.juliaplots.org/stable/) recipes for `AntibodyMethodsDoseResponse` objects. |



## About the tutorials

* The **quick start guide** 
    + A minimal tutorial, covering only the necessary steps to analyze dose-response data without any in-depth explanations. 
    + From here, the [AnitbodyMethodsDoseResponseConvenience API](@ref api_convenience) is highly recommended to explore the additional options of the convenience functions that are not covered in the minimal tutorial.
    + The recommendation for most users.

* The **Detailed explanations**
    + Covers the background and the internals of the packages.
    + Assumes (at some points) to be read in order from top to bottom.
    + The convenience functions (from the quick start guide) use the same data types and methods internally.
    + Intended for developers that need access to the internal methods and that want to take full control over the analysis.

## Installation

First, add the registry `Translational-Pain-Julia-Registry`:

```julia
using Pkg
Pkg.Registry.add()
Pkg.Registry.add(RegistrySpec(url = "https://github.com/Translational-Pain-Research/Translational-Pain-Julia-Registry"))
```

Then, the packages can be installed as usual. E.g. [`AntibodyMethodsDoseResponseConvenience.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponseConvenience.jl), which installs all necessary dependencies for the analysis of dose-response curves:

```julia
using Pkg
Pkg.add("AntibodyMethodsDoseResponseConvenience")
```

## How to cite the package

If you would like to cite this package for scientific purposes, you might also want to cite the corresponding paper [https://arxiv.org/abs/2407.06052](https://arxiv.org/abs/2407.06052).
