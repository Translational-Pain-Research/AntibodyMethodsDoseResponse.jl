# Models

## [Binding models](@id binding_models)

Two binding model kernels are available:

```math
\text{Accumulation:}\quad  x(a) = g(1-e^{-\frac{a}{K_\tau}})\qquad \text{and}\qquad \text{Langmuir isotherm:}\quad x(a) = \frac{g}{1+\frac{K_d}{a}}
```

* The Langmuir isotherm describes the equilibrium of binding (with rate ``k_{\text{on}}``) and unbinding (with rate ``k_{\text{off}}``) of antibodies, which is characterized by the **dissociation constant** ``K_d = \frac{k_{\text{off}}}{k_{\text{on}}}``.
* The accumulation model describes the accumulation of antibodies (with rate ``k_{\text{on}}``) during the incubation time ``\tau``, which is characterized by the **accessibility constant** ``K_\tau = \frac{1}{k_{\text{on}}\cdot \tau}``.

!!! tip "Jovanovic isotherm"
	The accumulation model resembles the *Jovanovic isotherm* structurally. However, the Jovanovic isotherm describes an equilibrium of binding and unbinding processes, whereas the accumulation model describes the accumulation over time, which is stopped at a finite time ``\tau`` before reaching the saturation. Nevertheless, because of the structural similarity, the accumulation model can be used as drop-in replacement to analyze equilibrium data with the Jovanovic isotherm.



Both model kernels are used for the rate constant distribution approach [(Svitel et al. 2003)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1302986/):

```math
\begin{aligned}
\text{Accumulation model:} \quad x(a) &= \int_0^\infty g(k)(1-e^{-\frac{a}{k}})\ dk\\
\text{Langmuir model:} \quad  x(a) &= \int_0^\infty \frac{g(k)}{1+\frac{k}{a}} \ dk\ .
\end{aligned}
```
In both cases, the density ``g(k)`` describes the number of epitopes that exhibit the respective constant ``k`` (accessibility / affinity depending on the model). The density approach models the superposition of different classes of epitopes being present in the system at the same time (e.g. in complex cellular systems).



!!! info "Default model"
	The accumulation model is set as default model for all methods in this package. Accordingly this documentation describes most terms, results, interpretations from the perspective of the accumulation model.


## Integral approximation

The analysis of dose-response data with the models above essentially means the estimation of the density ``g(k)`` from the dose-response data. To simplify the inference problem, ``g(k)`` is approximated by a sum of constant functions over a disjunct set of intervals ``\{I_j\}_{j=1}^m``:

```math
g(k) \approx \sum_{j=1}^m g_j \cdot \chi_{I_j}(x) \qquad \text{with} \qquad \chi_{I_j}(x) = \left\{ \begin{array}{ll} 1 & \ , x \in I_j\\ 0 &\ , \text{else} \end{array} \right.
```
With this approximation, the accumulation model becomes:

```math
x(a) \approx \sum_{j=1}^m g_j \int_{I_j}(1-e^{-\frac{a}{k}})\ dk\ .
```
Thus, the inference problem is the estimation of the parameters ``g_1,\ldots, g_m``.

The intervals ``I_j`` and the weights ``g_j`` are implemented as one-dimensional grid with [`AdaptiveDensityApproximation.jl`](https://translational-pain-research.github.io/AdaptiveDensityApproximation-documentation/).

```@example Models
using AntibodyMethodsDoseResponseConvenience #hide
using AdaptiveDensityApproximation
grid = create_grid([1,2,3,5])
```

The example above created a grid corresponding to the intervals ``\{[1,2],[2,3],[3,5]\}`` with weights ``g_j = 1``.  To view the grid properties, [`AdaptiveDensityApproximationRecipes.jl`](https://github.com/Translational-Pain-Research/AdaptiveDensityApproximationRecipes.jl) can be used:
```@example Models
using AdaptiveDensityApproximationRecipes, Plots
plot(grid)
```

New weights can be set with [`import_weights!`](https://translational-pain-research.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.import_weights!):
```@example Models
import_weights!(grid, [1,2,0.5])
plot(grid)
```

When the intervals have different lengths, the overall impact of the weights ``g_j`` becomes skewed, since the term ``\int_{I_j}(1-e^{-\frac{a}{k}})\ dk`` depends on the interval length. For the estimation of parameters it can be beneficial to use length-normalized parameters:

```math
g_j = \frac{\lambda_j}{\text{length}(I_j)}
```

In other words, while ``g_j`` is the density value, ``\lambda_j`` is the number of epitopes with ``K_\tau\in I_j``.

Finally, the analytical solution of ``\int_{I_j}(1-e^{-\frac{a}{k}})\ dk`` requires the *exponential integral* function, not implemented in Julia Base. To avoid additional dependencies, this integral is approximated by
```math
\int_{I_j}(1-e^{-\frac{a}{k}})\ dk \approx \text{length}(I_j) \cdot (1-e^{-\frac{a}{\text{center}(I_j)}})
```
```math
\Rightarrow \quad x(a) \approx \sum_{j=1}^m g_j \int_{I_j}(1-e^{-\frac{a}{k}})\ dk\ \approx \sum_{j=1}^m g_j \cdot \text{length}(I_j) \cdot (1-e^{-\frac{a}{\text{center}(I_j)}}) = \sum_{j=1}^m \lambda_j \cdot (1-e^{-\frac{a}{\text{center}(I_j)}})
```

!!! info "Weights of the grid"
	The weights of `OneDimGrid` objects are always treated as ``\lambda_j`` for the analysis of dose-response data. Thus, the weights describe the number of epitopes with ``K_\tau`` in the given interval, not the ``K_\tau``-density value!

!!! remark "Inverse constant"
	To solve the integrals analytically, the inverse constant ``\widetilde{K} = \frac{1}{K_\tau}`` can be used: 
	```math
	\int (1-e^{-a\cdot \widetilde{k}}) \ d \widetilde{k} = \frac{e^{-a\cdot \widetilde{k}}}{a} + \widetilde{k} + \text{constant}\ .
	```
	While the use of the inverse constant does not change the discrete superposition (countable sum), the integral approximation uses a different  density ``\widetilde{g}`` if the inverse constant is used.

## Obtain model functions

Having specified the intervals with a `grid`, the model function can be obtained with the following model generators: [`accumulation_model`](@ref), [`langmuir_model`](@ref), [`accumulation_inv_const_model`](@ref) and [`langmuir_inv_const_model`](@ref).

```@example Models
model, init_params, centers, volumes = accumulation_model(grid, offset = 10) 
nothing #hide
```

`model` is a [`ModelFunctions`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.ModelFunctions) object from [`FittingObjectiveFunctions.jl`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/), containing both the model function and the partial derivatives w.r.t. to the parameters.

`init_params` contains the weights of the grid, and as last element the `offset` if `offset != nothing`:

```@example Models
println(init_params)
```

`centers` contains the center points of the intervals of the `grid`:

```@example Models
println(centers)
```

and `volumes` contains the lengths of the intervals of the grid:

```@example Models
println(volumes)
```


## Tips for girds

Choosing unequal interval sizes in the example above was not just for demonstration purposes. In fact, a rule of thumb for the ``K_\tau`` domain is to use the concentration/dilution range of the dose-response curve, which often spans multiple orders of magnitude. Consider for example the domain ``[10^{-8},10^{-2}]``:

```@example Models
plot(create_grid(LinRange(1e-8,1e-2,50)))
```

In a linear scale, this interval discretization seems to resolve the domain well enough. But plotting the same grid in a logarithmic scale leads to:

```@example Models
plot(create_grid(LinRange(1e-8,1e-2,50)), xaxis = :log, xticks = [10.0^-i for i in 2:8])
```

Equally sized intervals poorly resolve the smaller orders of magnitude. A single interval covers the domain ``[10^{-8},10^{-4}]``, while almost all intervals subdivide the domain ``[10^{-3},10^{-2}]``. A solution for this problem is to use logarithmically sized intervals:

```@example Models
plot(create_grid(LogRange(1e-8,1e-2,50)), xaxis = :log, xticks = [10.0^-i for i in 2:8])
```

!!! tip "LogRange"
	Julia provides the `LinRange` function to create equally spaced points in a given range. [`AntibodyMethodsDoseResponse.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponse.jl) adds [`LogRange`](@ref) to create logarithmically spaced points in a given range, following the same general syntax: `LogRange(start, stop, n_points, [base = 10])`.
