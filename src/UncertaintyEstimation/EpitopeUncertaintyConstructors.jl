####################################################################################################
# Additional EpitopeUncertainty constructors
####################################################################################################


export EpitopeUncertainty





# EpitopeUncertainty from bin-wise shifting.
# Documentation in struct docstring.
function EpitopeUncertainty(data::FittingData,grid::AdaptiveDensityApproximation.OneDimGrid, bins = collect(1:length(grid)); levels = collect(0.1:0.1:1), options::AdaptiveOptions = AdaptiveOptions(), steps::Integer = 10^4, bisections::Integer = 10^2)

	level_check(levels, non_zero = true)
	sorted_levels = sort_levels(levels)

	# Use the methods of AdaptiveFitting.jl to construct the objective.
	# The methods of AdaptiveFitting.jl return an objective function that is sign-flipped to be a minimization objective.
	minimization_objective, parameters = get_objective(grid,data,options.offset, options)[1:2]
	objective = λ -> -minimization_objective(λ)
	lower, upper = initialize_bounds(sorted_levels, parameters, !isnothing(options.offset))
	
	# Skip empty bins (avoid infinite while loop in find_largest_increase).
	nonempty_bins = [bin for bin in bins if !isempty(bin)]

	for (completed_bins,bin) in enumerate(nonempty_bins)
		samples, objective_values, thresholds = shift_data(objective,parameters, bin,sorted_levels,!(options.objective == :posterior),steps,bisections)
		fill_boundaries!(lower,upper,thresholds,samples,objective_values,bin)
		print("Uncertainty estimation | completed $completed_bins of $(length(bins))\r")
	end

	# Offset parameters are not estimated (copy offset value for all levels if applicable).
	if isnothing(options.offset)
		lower_offset = nothing
		upper_offset = nothing
	else
		lower_offset = [options.offset for i in 1:length(sorted_levels)]
		upper_offset = [options.offset for i in 1:length(sorted_levels)]
	end

	EpitopeUncertainty(sorted_levels,lower,upper,lower_offset,upper_offset)
end






# EpitopeUncertainty from samples.
# Documentation in struct docstring.
#Different order of axes for samples [parameter index, sample index] and bounds [level index, parameter index], because bound indices should match offset indices (level index first).
function EpitopeUncertainty(samples::AbstractMatrix{T}; levels = collect(0.1:0.1:1) , offset::Bool = false) where T <: Real

	level_check(levels)

	if size(samples,1) == 1
		# Otherwise, using size(samples,1)-1 leads to errors if offset != nothing.
		throw(DimensionMismatch("Grid parameter length must be larger than 1!"))
	end

	sorted_levels = sort_levels(levels)

	lower_quantiles = @. 0.5 - sorted_levels/2
	upper_quantiles = @. 0.5 + sorted_levels/2

	lower, upper = initialize_bounds(sorted_levels, samples[:,1], offset)

	if offset
		for level_ind in eachindex(sorted_levels)
			lower[level_ind,:] .= [quantile(samples[parameter_ind,:], lower_quantiles[level_ind]) for parameter_ind in 1:size(samples,1)-1]
			upper[level_ind,:] .= [quantile(samples[parameter_ind,:], upper_quantiles[level_ind]) for parameter_ind in 1:size(samples,1)-1]
		end

		lower_offsets = [quantile(samples[end,:], q) for q in lower_quantiles]
		upper_offsets = [quantile(samples[end,:], q) for q in upper_quantiles]

		return EpitopeUncertainty(sorted_levels,lower,upper,lower_offsets,upper_offsets)
	else
		for level_ind in eachindex(sorted_levels)
			lower[level_ind,:] .= [quantile(samples[parameter_ind,:], lower_quantiles[level_ind]) for parameter_ind in 1:size(samples,1)]
			upper[level_ind,:] .= [quantile(samples[parameter_ind,:], upper_quantiles[level_ind]) for parameter_ind in 1:size(samples,1)]
		end

		return EpitopeUncertainty(sorted_levels,lower,upper)
	end
end


# Dispatch version for samples as array of arrays.
function EpitopeUncertainty(samples; args...)
	return EpitopeUncertainty(hcat(samples...); args...)
end



