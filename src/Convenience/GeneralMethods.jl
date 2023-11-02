####################################################################################################
# General convenience methods
####################################################################################################


export LogRange, peak_detection




"""
    LogRange(start::Real,stop::Real,n::Integer, base::Real = 10.0)
Return a vector with `n` values logarithmically distributed between `start` and `stop`. The logarithm `base` can be changed with the last, optional argument.
"""
function LogRange(start::Real,stop::Real,n::Integer, base::Real = 10.0)
    if base <= 0 
        throw(DomainError(base, "`base` must be larger than 0."))
    end
    if start <= 0 
        throw(DomainError(start, "`start` must be larger than 0."))
    end
    if stop <= 0 
        throw(DomainError(stop, "`stop` must be larger than 0."))
    end
    if n <= 0 
        throw(DomainError(n, "Number of intermediate points must be larger than 0."))
    end
	
    return base .^ LinRange(log(base,start),log(base,stop),n)
end





# Return groups of indices (a group is defined by adjacent indices).
function index_groups(indices::Vector{T}) where T <: Integer

	groups = Vector{Int64}[]
	if isempty(indices)
		return groups
	end

	current_group = Int64[indices[1]]

	for i in 2:length(indices)
		# Index difference larger than 1 defines a new group.
		if abs(current_group[end]-indices[i]) > 1
			push!(groups,current_group)
			current_group = Int64[indices[i]]
		else
			push!(current_group,indices[i])
		end
	end

	# Add last group (if not empty) to groups.
	if !isempty(current_group)
		push!(groups,current_group)
	end

	return groups
end




"""
	peak_detection(grid::OneDimGrid, relative_threshold::Real = 0.1; volume_normalization::Symbol = :log, fill::Bool = true)
Return `(peak_group_indices,peak_group_domains)` w.r.t. the `relativ_ threshold`.

* `peak_group_indices` the index-vectors of the peaks. The indices correspond to exported grid centers/volumes/weights (see [`export_weights`](@ref) and [`export_all`](@ref)).
* `peak_group_domains` contains the intervals covered by the respective peaks.
* The cutoff threshold is determined by `relative_threshold * largest_weight`.
* If `fill == true`, the gaps between the peaks are added to `peak_group_indices` and `peak_group_domains`.
* `volume_normalization` normalizes the weights of the grid (without mutation). `:none` uses the raw weights. `:linear` divides the weights by the block volume and `:log` divides the weight by the block volume in a logarithmic scale.
"""
function peak_detection(grid::AdaptiveDensityApproximation.OneDimGrid, relative_threshold::Real = 0.1; volume_normalization::Symbol = :log, fill::Bool = true)

	if relative_threshold <= 0 || relative_threshold > 1
		throw(DomainError("Relative threshold must be in (0,1]."))
	end

	centers, volumes, weights = export_all(grid)

	if volume_normalization == :linear
		weights = weights ./ volumes
	elseif volume_normalization == :log
		log_volumes = [log(centers[i]+volumes[i]/2) - log(centers[i] - volumes[i]/2) for i in eachindex(centers)]
		weights = weights ./ log_volumes
	end

	absolute_threshold = maximum(weights) * relative_threshold
	above_threshold = findall(x-> x >= absolute_threshold, weights)
	
	peak_groups = index_groups(above_threshold)


	if fill
		remaining_indices = findall(x-> x < absolute_threshold, weights)
		remaining_groups = index_groups(remaining_indices)
		# Prevent pushing empty groups (see below).
		if !isempty(remaining_groups)
			push!(peak_groups, remaining_groups...)
		end
	end

	if isempty(peak_groups)
		return peak_groups, Vector{Float64}[]
	end
	

	peak_domains = Vector{promote_type(eltype(centers), eltype(volumes), Float64)}[]
	for peak in peak_groups
		# Calculate peak intervals from centers and volumes, using the first and the last index.
		push!(peak_domains,[centers[peak[1]]-0.5*volumes[peak[1]], centers[peak[end]]+0.5*volumes[peak[end]]])
	end
	
	return peak_groups, peak_domains
end



