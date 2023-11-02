####################################################################################################
# Additional DoseResponseUncertainty constructors
####################################################################################################


export DoseResponseUncertainty


# DoseResponseUncertainty construction from EpitopeUncertainty struct.
# Documentation in struct docstring.
# No offset keyword necessary, offsets are contained in the EpitopeUncertainty object.
function DoseResponseUncertainty(grid::AdaptiveDensityApproximation.OneDimGrid,eu::EpitopeUncertainty,concentrations::AbstractVector{T}; bins = [collect(1:length(grid))], model::Function = accumulation_model) where T <: Real

	# Avoid mutation of original grid during response calculation.
	temp_grid = deepcopy(grid)

	# Use exemplary parameters from EpitopeUncertainty to determine type of responses.
	import_weights!(temp_grid, eu.lower[1,:])
	matrix_type = promote_type(Float64,eltype(DoseResponseResult(temp_grid,concentrations,model = model).responses), eltype(concentrations))

	# If eu.upper_offset exists, the maximal response needs to be calculated using the maximal offset.
	# No need to calculate min_response, as lower bound is 0 anyway.
	if isnothing(eu.upper_offset)
		max_response = maximum(DoseResponseResult(import_weights!(temp_grid, eu.upper[1,:]),concentrations).responses)
	else
		max_response = maximum(DoseResponseResult(import_weights!(temp_grid, eu.upper[1,:]),concentrations, offset = maximum(eu.upper_offset)).responses)
	end

	# Initialize with largest/smallest possible values, as point-wise minima/maxima will be taken later on.
	lower = initialize_bounds(eu.levels,fill(matrix_type(max_response), length(concentrations)),false)[1]
	upper = initialize_bounds(eu.levels,zeros(matrix_type, length(concentrations)),false)[2]
	

	# Best parameters that are varied one bin at a time (without mutation using setindex).
	result_parameters = export_weights(grid)

	for i in eachindex(eu.levels)
		for bin in bins
			# lower.
			import_weights!(temp_grid, setindex(result_parameters,eu.lower[i,bin],bin))
			offset = isnothing(eu.lower_offset) ? 0 : eu.lower_offset[i]
			lower[i,:] .= min.(DoseResponseResult(temp_grid,concentrations,model = model, offset = offset).responses,lower[i,:])

			# upper.
			import_weights!(temp_grid, setindex(result_parameters,eu.upper[i,bin],bin))
			offset = isnothing(eu.upper_offset) ? 0 : eu.upper_offset[i]
			upper[i,:] .= max.(DoseResponseResult(temp_grid,concentrations,model = model, offset = offset).responses,upper[i,:])
		end
	end
	return DoseResponseUncertainty(eu.levels,concentrations,lower,upper)
end