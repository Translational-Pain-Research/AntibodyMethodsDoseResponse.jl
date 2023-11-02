# FittingObjectiveFunctions

This is a lightweight package without dependencies to create objective functions for model-fitting.
To avoid dependencies, this package does not include optimizers / samplers.


## Installation

This package is not in the general registry and needs to be installed from the GitHub repository:

```@julia
using Pkg
Pkg.add(url="https://github.com/AntibodyPackages/FittingObjectiveFunctions")
```


## Simple example

After the installation, the package can be used like any other package:
```@example 1
using FittingObjectiveFunctions
```

Consider the following example data-set:

```@example 1
using Plots

X = collect(1:10)
Y = [1.0, 1.78, 3.64, 3.72, 5.33, 2.73, 7.52, 9.19, 6.73, 8.95]
ΔY = [0.38, 0.86, 0.29, 0.45, 0.66, 2.46, 0.39, 0.45, 1.62, 1.54]

scatter(X,Y, yerror = ΔY, legend=:none, xlabel = "X", ylabel="Y")
```

Before objective functions can be created, the data needs to be summarized in a `FittingData` object:

```@example 1
data = FittingData(X,Y,ΔY)
nothing #hide
```

Information about the model needs to be summarized in a `ModelFunctions` object. Here we choose a simple linear model ``m(x,\lambda) = \lambda x ``:


```@example 1
model = ModelFunctions((x,λ) -> λ*x) 
nothing #hide
```

A weighted least squares objective can be be constructed as follows:

```@example 1
lsq = lsq_objective(data,model)

plot(λ-> lsq(λ), legend = :none, xlabel = "λ", ylabel = "lsq")
```