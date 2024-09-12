####################################################################################################
# Adaptive fitting methods
####################################################################################################


# Block variations
####################################################################################################

export area_scaled_variation, log_area_scaled_variation

"""
	area_scaled_variation(center, volume, weight, 
		neighbor_centers, neighbor_volumes, neighbor_weights)

Block variation function for [`refine!`](https://translational-pain-research.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.refine!). Variation value based on the difference of the weights, scaled with the area (volume) of the corresponding block.

	mean(@. abs(weight * volume - neighbor_weights * neighbor_volumes))
"""
function area_scaled_variation(center, volume, weight, neighbor_centers, neighbor_volumes, neighbor_weights)
	mean(@. abs(weight * volume - neighbor_weights * neighbor_volumes))
end

"""
	log_area_scaled_variation(center, volume, weight, 
		neighbor_centers, neighbor_volumes, neighbor_weights)

Block variation function for [`refine!`](https://translational-pain-research.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.refine!). Variation value based on the difference of the weight, scaled with the visible area (in a logarithmic plot) of the corresponding block.

	log_volume = (log10(center + volume / 2) - log10(center - volume / 2))
	neighbor_log_volumes = @. (log10(neighbor_centers + neighbor_volumes / 2) - log10(neighbor_centers - neighbor_log_volumes / 2))
	mean(@. abs(weight * log_volume  - neighbor_weights * neighbor_log_volumes)) 
"""
function log_area_scaled_variation(center, volume, weight, neighbor_centers, neighbor_volumes, neighbor_weights)
	log_volume = (log10(center + volume / 2) - log10(center - volume / 2))
	neighbor_log_volumes = @. (log10(neighbor_centers + neighbor_volumes / 2) - log10(neighbor_centers - neighbor_volumes / 2))
	mean(@. abs(weight * log_volume  - neighbor_weights * neighbor_log_volumes)) 
end




















# Adaptive fitting recipe
####################################################################################################


export AdaptiveOptions, AdaptiveResult, adaptive_dose_response_fit




@inline function default_prior_generator(centers, volumes, offset)
	return FittingObjectiveFunctions.log_uniform_prior
end

@inline function default_prior_gradient_generator(centers, volumes, offset)
	return nothing
end

"""
	mutable struct AdaptiveOptions
Data type to define [`adaptive_dose_response_fit`](@ref) options.

**Constructor**

	AdaptiveOptions(keywords...)

The following keywords (with default values) are available:

* `name::AbstractString = "Adaptive optimization"`: The name that is used when `show_progress==true`.
* `show_progress::Bool = true`: Show progress in standard output.
* `iterations::Integer = 1`: Number of refinement iterations.
* `model::Function = accumulation_model`: The model-function that is used for the data-fit. The available model functions are [`accumulation_model`](@ref), [`accumulation_inv_const_model`](@ref), [`langmuir_model`](@ref) and [`langmuir_inv_const_model`](@ref).
* `offset = nothing`: Offset parameter for the model function. If `nothing`, no offset is used.
* `objective::Symbol = :lsq`. The objective function for the data-fit. Available are `:lsq`, `:posterior` and `:log_posterior`.
* `prior_generator::Function = default_prior_generator`: The function that generates the prior. The function must have the signature `(grid_centers,grid_volumes,offset)` and must return a function `λ-> prior(λ)` or `λ-> log_prior(λ)` in case of a `:log_posterior` objective. The `default_prior_generator` generates a uniform prior `λ-> 0` for the log-posterior objective.
* `distribution_derivatives = nothing`: Array of partial derivatives of the logarithmic distributions for the log-posterior objective. See [`log_posterior_gradient`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.log_posterior_gradient).
* `prior_gradient_generator = default_prior_gradient_generator`: The function that generates the log-prior gradient (see  [`log_posterior_gradient`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.log_posterior_gradient)). The function must have the signature `(grid_centers,grid_volumes,offset)` and must return a function `λ-> ∇log_prior(λ)`. The `default_prior_gradient_generator` returns `nothing` which internally corresponds to the uniform prior for the log-posterior objective.
* `block_variation::Function =` [`log_area_scaled_variation`](@ref) and `selection::Function = maximum` are the refinement options of [`refine!`](https://translational-pain-research.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.refine!).
"""
mutable struct AdaptiveOptions
	name::AbstractString
	show_progress::Bool
	iterations::Integer
	model::Function
	offset
	objective::Symbol
	prior_generator::Function
	distribution_derivatives
	prior_gradient_generator
	block_variation::Function
	selection::Function

	function AdaptiveOptions(; 
		name="Adaptive optimization", 
		show_progress=true, 
		iterations = 1,
		model = accumulation_model, 
		offset=nothing, 
		objective = :lsq, 
		prior_generator = default_prior_generator, 
		distribution_derivatives=nothing, 
		prior_gradient_generator = default_prior_gradient_generator, 
		block_variation = log_area_scaled_variation, 
		selection = maximum)

		return new(name, 
			show_progress,
			iterations,
			model,
			offset,
			objective,
			prior_generator, 
			distribution_derivatives,
			prior_gradient_generator, 
			block_variation, 
			selection)
	end
end



"""
	mutable struct AdaptiveResult
Data type used by [`adaptive_dose_response_fit`](@ref) to summarize the results.

The struct has the following fields:

* `result`: The [`DoseResponseResult`](@ref) object corresponding to the fit result.
* `grid`: The grid (with imported weights) corresponding to the fit result.
* `optimizer`: The raw result parameter.
* `objective_value`: The objective-function value of `optimizer`.
* `time`: The elapsed time for the model fit (in seconds).
"""
mutable struct AdaptiveResult
	result
	grid
	optimizer
	objective_value
	time
end







# Create transformed objective function to be minimized (internal).
# offset as extra argument (despite being in options) to pass on changes during iterative calls.
function get_objective(grid::AdaptiveDensityApproximation.OneDimGrid,data::FittingData,offset,options::AdaptiveOptions)

	# Get model and initial parameters.
		model, initial_parameters, centers,volumes = options.model(grid, offset = offset)
		prior = options.prior_generator(centers, volumes, offset)
		prior_gradient = options.prior_gradient_generator(centers, volumes,offset)

	# Get objective function.
	if options.objective == :lsq 
		objective_function = lsq_objective(data,model)
		objective_gradient! = lsq_gradient(data,model)

	elseif options.objective == :posterior
		obj = posterior_objective(data,model,prior)
		# Maximizing posterior <=> minimizing -posterior
		objective_function = @inline function(λ) -obj(λ) end
		objective_gradient! = nothing

	elseif options.objective == :log_posterior
		obj = log_posterior_objective(data,model,prior)
		# Maximizing posterior <=> minimizing -posterior
		objective_function = @inline function(λ) -obj(λ) end

		if isnothing(options.distribution_derivatives)
			objective_gradient! = nothing
		else
			# Sign flip in objective function must also be applied to the gradient.
			grad! = log_posterior_gradient(data,model,options.distribution_derivatives, prior_gradient)
			objective_gradient! = @inline function(gradient_vector,λ) return -grad!(gradient_vector,λ) end
		end
	end

	return objective_function, initial_parameters, objective_gradient!
end









# Import result into grid and create DoseResponseResult object.
function apply_result!(grid::AdaptiveDensityApproximation.OneDimGrid,data::FittingData, parameters,options::AdaptiveOptions)
	# Without offset, all parameters belong to the grid. With offset, parameters[1:end-1] belong to the grid and parameters[end] is the offset.
	if isnothing(options.offset)
		import_weights!(grid,parameters)
		dose_response_result = DoseResponseResult(grid,data.independent, model = options.model)
	else
		import_weights!(grid,parameters[1:end-1])
		dose_response_result = DoseResponseResult(grid, data.independent, offset = parameters[end], model = options.model)
	end
	return grid, dose_response_result, parameters
end








"""
	adaptive_dose_response_fit(initial_grid::OneDimGrid, 
		data::FittingData, 
		minimizer::Function; 
		options::AdaptiveOptions=AdaptiveOptions()
	)

Fit dose-response `data` (with adaptive grid refinements depending on the `options`) and return an [`AdaptiveResult`](@ref) object. The initial gird is not mutated. For the [`AdaptiveResult`](@ref) object a copy of the gird is created.


**Minimizer function**

* The sign of posterior and log-posterior objectives is flipped for consistency reasons. A minimizer needs to be used for all objectives (`:lsq`, `:posterior`, `:log_posterior`).
* `minimizer`: The function that minimizes the objective function. It needs to be specified by the user and must have the signature `(objective_function, objective_gradient!,parameters)`.
* The `objective_function` always has the signature `(parameters)`.
* The `objective_gradient!` can be `nothing`. Otherwise it must have the signature `(gradient_vector, parameter)`. It must mutate the `gradient_vector` and return the mutated `gradient_vector`.
* `parameters` is the initial parameter array to start minimization from.

**Gradients for minimization**

Whether `objective_gradient!` is `nothing` or a proper function depends on the specified `options` (see [`AdaptiveOptions`](@ref)).

* The `:lsq` objective always provides analytical gradients.
* The `:posterior` objective never provides analytical gradients. 
* The`:log_posterior` objective only provides analytical gradients if `distribution_derivatives != nothing`.

"""
function adaptive_dose_response_fit(initial_grid::AdaptiveDensityApproximation.OneDimGrid, data::FittingData, minimizer::Function; options::AdaptiveOptions=AdaptiveOptions())
	
	# Prevent mutation of user-specified grid.
	grid = deepcopy(initial_grid)

	# Declare all variables outside of loop scope to access result after the loop has completed.
	total_time = 0
	offset = options.offset
	objective_function, parameters, objective_gradient! = get_objective(grid, data,offset, options)
	grid, dr_result, parameters = apply_result!(grid,data,parameters,options)


	for i in 1:options.iterations
		objective_function, parameters, objective_gradient! = get_objective(grid, data,offset, options)

		start_time = time()
		parameters = minimizer(objective_function,objective_gradient!,parameters)
		total_time += time()-start_time
		
		grid, dr_result, parameters = apply_result!(grid,data,parameters,options)

		if i < options.iterations
			refine!(grid,block_variation = options.block_variation, selection = options.selection, split_weights = true)
		end
		
		if options.show_progress
			minutes = Int64(round(total_time ÷ 60))
			seconds = Int64(round(total_time % 60))
			print("$(options.name) | completed $(i) of $(options.iterations) | passed time $(minutes):$(seconds) \r")
		end
	end

	return AdaptiveResult(dr_result,grid,parameters,objective_function(parameters),total_time)
end