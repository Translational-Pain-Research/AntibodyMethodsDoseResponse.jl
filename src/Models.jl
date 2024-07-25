####################################################################################################
# Models for dose-response-curve fitting.
####################################################################################################



# Construction of model functions
####################################################################################################
export accumulation_model, accumulation_inv_const_model, langmuir_model, langmuir_inv_const_model  


# Internal function to test if a grid satisfies the properties of a K_τ-density.
# Only positive, real numbers (not ±Inf and not NaN).
function density_grid_tests(grid::AdaptiveDensityApproximation.OneDimGrid)
	for block in AdaptiveDensityApproximation.values(grid)
		if isnan(block.interval.left) || isnan(block.interval.right)
			throw(ArgumentError("K_τ values must not be NaN. Grid domain contains NaN."))
		elseif isinf(block.interval.left) || isinf(block.interval.right)
			throw(ArgumentError("K_τ values must not be ±Inf. Grid domain contains ±Inf."))
		elseif block.interval.left < 0
			throw(DomainError(block.interval.left,"K_τ values must not be negative. Grid domain contains negative values."))
		elseif !(typeof(block.weight) <: Real)
			throw(ArgumentError("Weight type of the grid must be <: Real, but is $(typeof(block.weight))."))
		end
	end
end



# Common steps to create return-values for the models.
function generate_model(model,∂model,centers, volumes, parameters, offset)
	# Offset requires to redefine the methods (additional element for the parameters array).
	# Different names to avoid infinite recursion at the redefinition for the offset case.
	if isnothing(offset)
		f,λ, ∂f = model, parameters, ∂model
	else
		# Parameters array gets an additional element, the offset.
		# The offset does not belong to the grid, s.t. the model needs to be redefined to only consider parameters[1:end-1] as K_τ-grid values.
		λ = vcat(parameters,offset)
		f = @inline function(a,λ) return model(a,λ[1:end-1]) + λ[end] end
		∂f = vcat([@inline function(a,λ) return ∂(a,λ[1:end-1]) end for ∂ in ∂model], @inline function(a,λ) return 1 end)
	end

	return ModelFunctions(f,partials = ∂f), λ, centers, volumes
end






@doc raw"""
	accumulation_model(grid::OneDimGrid; offset = nothing)
Create a multi-epitope accumulation model. Returns `(model,λ,centers,volumes)` where 

* `model` is a [`ModelFunctions`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.ModelFunctions) object.
* `λ` is an initial parameter array (the weights of the grid and the offset if `offset != nothing`). If `offset != nothing`, the last element is the offset parameter `λ[end] = offset`.
* `centers` and `volumes` are the remaining properties of the grid, see [`export_all`](https://antibodypackages.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.export_all)

**Model function**

The following model function and partial derivatives are used:

```math 
\text{model}(a,\lambda) = \lambda_e + \sum_i \lambda_i \left(1-e^{-\frac{a}{c_i}}\right)  \approx  \lambda_e + \sum_i \int_{l_i}^{u_i} \frac{\lambda_i}{u_i-l_i}\left(1-e^{-\frac{a}{k}}\right) \ dk 
```
```math 
\partial_{\lambda_j} \text{model}(a,\lambda) = 1-e^{-\frac{a}{c_i}} \quad, \qquad \partial_{\lambda_e} \text{model}(a,\lambda) = 1
```
where ``a`` is the antibody concentration, ``c_i`` are the centers of the `grid` intervals ``[u_i,l_i]`` and ``\lambda_e`` is the offset (if `offset != nothing`).
"""
function accumulation_model(grid::AdaptiveDensityApproximation.OneDimGrid ; offset::Union{Nothing,R} = nothing) where R <: Real
	density_grid_tests(grid)

	centers, volumes, parameters = export_all(grid)

	model = let centers = centers
		@inline function(a,λ)
			return sum(λ[i]*(1-exp(-a/centers[i])) for i in 1:length(volumes))
		end
	end

	∂model = Function[]
	for i in eachindex(centers)
		∂ = let center = centers[i]
			@inline function(a,λ)
				return (1-exp(-a/center))
			end
		end
		push!(∂model,∂)
	end

	return generate_model(model,∂model,centers,volumes, parameters,offset)
end






@doc raw"""
	accumulation_inv_const_model(grid::OneDimGrid; offset = nothing)
Create a multi-epitope accumulation model with `1/K_τ = k_a * τ` as constant domain. Returns `(model,λ,centers,volumes)` where 

* `model` is a [`ModelFunctions`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.ModelFunctions) object.
* `λ` is an initial parameter array (the weights of the grid and the offset if `offset != nothing`). If `offset != nothing`, the last element is the offset parameter `λ[end] = offset`.
* `centers` and `volumes` are the remaining properties of the grid, see [`export_all`](https://antibodypackages.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.export_all)

**Model function**

The following model function and partial derivatives are used:

```math 
\text{model}(a,\lambda) = \lambda_e + \sum_i \lambda_i \left(1+\frac{1}{a\cdot(u_i-l_i)}\left(e^{-a u_i}-e^{-a l_i} \right) \right) =   \lambda_e + \sum_i \int_{l_i}^{u_i} \frac{\lambda_i}{u_i-l_i}\left(1-e^{-a k}\right) \ dk 
```
```math 
\partial_{\lambda_j} \text{model}(a,\lambda) = 1+\frac{1}{a\cdot(u_i-l_i)}\left(e^{-a u_i}-e^{-a l_i} \right) \quad, \qquad \partial_{\lambda_e} \text{model}(a,\lambda) = 1
```
where ``a`` is the antibody concentration, ``c_i`` are the centers of the `grid` intervals ``[u_i,l_i]`` and ``\lambda_e`` is the offset (if `offset != nothing`).
"""
function accumulation_inv_const_model(grid::AdaptiveDensityApproximation.OneDimGrid ; offset::Union{Nothing,R} = nothing) where R <: Real
	density_grid_tests(grid)

	centers, volumes, parameters = export_all(grid)
	left_edges = centers .- volumes ./ 2
	right_edges = centers .+ volumes ./ 2

	model = let volumes = volumes, left_edges = left_edges, right_edges = right_edges
		@inline function(a,λ)
			return sum(λ[i]*(1+1/(a*volumes[i])*(exp(-a*right_edges[i])-exp(-a*left_edges[i]))) for i in 1:length(volumes))
		end
	end

	∂model = Function[]
	for i in eachindex(centers)
		∂ = let volumes = volumes, left_edges = left_edges, right_edges = right_edges
			@inline function(a,λ)
				return 1+1/(a*volumes[i])*(exp(-a*right_edges[i])-exp(-a*left_edges[i]))
			end
		end
		push!(∂model,∂)
	end

	return generate_model(model,∂model,centers,volumes, parameters,offset)
end

















@doc raw"""
	langmuir_model(grid::OneDimGrid; offset = nothing)
Create a multi-epitope Langmuir model. Returns `(model,λ,centers,volumes)` where 

* `model` is a [`ModelFunctions`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.ModelFunctions) object.
* `λ` is an initial parameter array (the weights of the grid and the offset if `offset != nothing`). If `offset != nothing`, the last element is the offset parameter `λ[end] = offset`.
* `centers` and `volumes` are the remaining properties of the grid, see [`export_all`](https://antibodypackages.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.export_all).


**Model function**

The following model function and partial derivatives are used:

```math 
\text{model}(a,\lambda) = \lambda_e + \sum_i \frac{\lambda_i\cdot a}{ (u_i - l_i)} \ln\left(\frac{a+u_i}{a+l_i}\right)  =  \lambda_e + \sum_i \int_{l_i}^{u_i} \frac{\frac{\lambda_i}{u_i - l_i}}{1+\frac{k}{a}} \ dk 
```
```math 
\partial_{\lambda_j} \text{model}(a,\lambda) = \frac{a}{u_j-l_j} \ln\left(\frac{a+u_j}{a+l_j}\right)\quad, \qquad \partial_{\lambda_e} \text{model}(a,\lambda) = 1
```
where ``a`` is the antibody concentration, ``[l_i,u_i]`` are the intervals of the `grid` and ``\lambda_e`` is the offset (if `offset != nothing`).
"""
function langmuir_model(grid::AdaptiveDensityApproximation.OneDimGrid ; offset::Union{Nothing,R} = nothing) where R <: Real
	density_grid_tests(grid)

	centers, volumes, parameters = export_all(grid)
	left_edges = centers .- volumes ./ 2
	right_edges = centers .+ volumes ./ 2

	model = let volumes = volumes, left_edges = left_edges, right_edges = right_edges
		@inline function(a,λ)
			return sum(λ[i]*a/(volumes[i])*log((a+right_edges[i])/(a+ left_edges[i])) for i in 1:length(volumes))
		end
	end

	∂model = Function[]
	for i in eachindex(centers)
		∂ = let volumes = volumes, left_edges = left_edges, right_edges = right_edges
			@inline function(a,λ)
				return a/(volumes[i])*log((a+right_edges[i])/(a+ left_edges[i]))
			end
		end
		push!(∂model,∂)
	end

	return generate_model(model,∂model,centers,volumes, parameters,offset)
end







@doc raw"""
	langmuir_inv_const_model(grid::OneDimGrid; offset = nothing)
Create a multi-epitope Langmuir model with `1/K_d = k_a / k_d` as constant domain. Returns `(model,λ,centers,volumes)` where 

* `model` is a [`ModelFunctions`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.ModelFunctions) object.
* `λ` is an initial parameter array (the weights of the grid and the offset if `offset != nothing`). If `offset != nothing`, the last element is the offset parameter `λ[end] = offset`.
* `centers` and `volumes` are the remaining properties of the grid, see [`export_all`](https://antibodypackages.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.export_all)

**Model function**

The following model function and partial derivatives are used:

```math 
\text{model}(a,\lambda) = \lambda_e + \sum_i \lambda_i \left(1+ \frac{1}{a\cdot (u_i-l_i)} \ln\left(\frac{a l_i +1}{a u_i +1}\right)\right)  = \lambda_e + \sum_i \int_{l_i}^{u_i} \frac{\frac{\lambda_i}{u_i - l_i}}{1+\frac{1}{a\cdot k}} \ dk  
```
```math 
\partial_{\lambda_j} \text{model}(a,\lambda) = 1+ \frac{1}{a\cdot (u_i-l_i)} \ln\left(\frac{a l_i +1}{a u_i +1}\right)\quad, \qquad \partial_{\lambda_e} \text{model}(a,\lambda) = 1
```
where ``a`` is the antibody concentration, ``[l_i,u_i]`` are the intervals in the `grid` and ``\lambda_e`` is the offset (if `offset != nothing`).
"""
function langmuir_inv_const_model(grid::AdaptiveDensityApproximation.OneDimGrid ; offset::Union{Nothing,R} = nothing) where R <: Real
	density_grid_tests(grid)

	centers, volumes, parameters = export_all(grid)
	left_edges = centers .- volumes ./ 2
	right_edges = centers .+ volumes ./ 2

	model = let volumes = volumes, left_edges = left_edges, right_edges = right_edges
		@inline function(a,λ)
			return sum(λ[i]*(1+1/(a*volumes[i])*log((a*left_edges[i]+1)/(a*right_edges[i]+1))) for i in 1:length(volumes))
		end
	end

	∂model = Function[]
	for i in eachindex(centers)
		∂ = let volumes = volumes, left_edges = left_edges, right_edges = right_edges
			@inline function(a,λ)
				return 1+1/(a*volumes[i])*log((a*left_edges[i]+1)/(a*right_edges[i]+1))
			end
		end
		push!(∂model,∂)
	end

	return generate_model(model,∂model,centers,volumes, parameters,offset)
end

























# DoseResponseResult for simulations / resulting dose-response curve from a curve fit.
####################################################################################################

export DoseResponseResult

# Extend dose response check for DoseResponseResult constructor.
# Not documented, only for internal use in the constructor.
function dose_response_check(concentrations,responses)
	regular_positive_numbers_check(concentrations)
	regular_positive_numbers_check(responses)

	if length(concentrations) != length(responses)
		throw(DimensionMismatch("Lengths of `concentrations` and `responses` do not match."))
	end
end



"""
	struct DoseResponseResult
Data type to store dose-response result data (e.g. from a dose-response-curve fitting or a simulation).	

**Fields**

* `concentrations`: Antibody concentrations.
* `responses`: Response values.

**Default constructor**

	DoseResponseResult(concentrations,responses)

**Model constructor**

	DoseResponseResult(grid::OneDimGrid,concentrations; 
		offset::Real = 0,
		model::Function = accumulation_model
	)

Calculate a dose-response curve from a K_τ `grid` and a `model` for given `concentrations`. The `offset` value is a global additive shift for all response values.

The available model functions are [`accumulation_model`](@ref), [`accumulation_inv_const_model`](@ref), [`langmuir_model`](@ref) and [`langmuir_inv_const_model`](@ref).
"""
struct DoseResponseResult
	concentrations
	responses

	function DoseResponseResult(concentrations,responses)
		# Only allow proper dose-response data.
		dose_response_check(concentrations,responses)
		return new(concentrations,responses)
	end
end



# Documentation included in struct docstring.
function DoseResponseResult(grid::AdaptiveDensityApproximation.OneDimGrid,concentrations; offset::Real = 0, model::Function = accumulation_model)
	density_grid_tests(grid)
	regular_positive_numbers_check(concentrations)
	model_struct, λ = model(grid,offset = offset)
	model_function = model_struct.model
	responses = [model_function(c,λ) for c in concentrations]
	return DoseResponseResult(concentrations,responses)
end

