####################################################################################################
# Internal methods for the uncertainty estimation.
####################################################################################################


# General methods.
####################################################################################################

# Levels must be probabilities.
function level_check(levels; non_zero::Bool = false)

	# Test indexablity before individual element tests are performed (otherwise an unspecific getindex method error is thrown).
	if !(applicable(eachindex,levels))
		throw(ArgumentError("Indexable collection expected, but got $(typeof(levels))."))
	end

	for i in eachindex(levels)
		if !(0<= levels[i] <= 1)
			throw(DomainError("Levels must be within `[0,1]`!"))
		end
		if non_zero && iszero(levels[i])
			throw(DomainError("Levels must not be zero if keyword `non_zero` is set to `true`!"))
		end
	end
end



# Create sorted array of unique levels.
function sort_levels(l)
	return sort(unique(l))
end

# Dispatch version for single level passed as number.
function sort_levels(l::Real)
	return [l]
end


# Relative threshold values (w.r.t. to most probable value p_max) from the levels, depending on whether logarithmic densities are used.
function get_thresholds(levels, p_max, log_density::Bool)
	if log_density
		# p_max is expected to be the density value (i.e. already logarithmic) 
		# Multiplication becomes addition for logarithmic densities.
		return [log(l) + p_max for l in levels]
	else
		return [l*p_max for l in levels]
	end
end


# Construct matrices for the lower and upper bound estimates.
# Levels and parameters determine the size (order: level, parameter index).
# Fill the matrices with the parameters for each level.
function initialize_bounds(levels, parameters, offset::Bool = false)
	if offset
		# With offset, parameters[end] does not belong to the grid -> length(parameters)-1 and parameters[1:end-1].
		lower = Matrix{promote_type(eltype(parameters),Float64)}(undef, length(levels), length(parameters)-1)
		for i in 1:length(levels)
			lower[i,:] = parameters[1:end-1]
		end
	else
		lower = Matrix{promote_type(eltype(parameters),Float64)}(undef, length(levels), length(parameters))
		for i in 1:length(levels)
			lower[i,:] = parameters
		end
	end
	upper = deepcopy(lower)
	return lower,upper
end






















# EpitopeUncertainty from bin-wise shifting.
# Methods are tested in runtests.
####################################################################################################



# Standard interval bisection algorithm. Returns the last interval.
function interval_bisection(f,initial_interval,iterations)
	f_values = f.(initial_interval)

	# Sign change of function is a prerequisite.
	if (f_values[1] > 0 && f_values[2] > 0) || (f_values[1]<0 && f_values[2]<0)
		throw(ArgumentError("Interval must satisfy f(I[1]) > 0 and f(I[2]) < 0 or vice versa. Got f($(initial_interval[1]))=$(f_values[1]) and f($(initial_interval[2])) = $(f_values[2])."))
	end

	# If zero is one of the interval boundaries, return it as trivial interval immediately (no need for additional calculations).
	if iszero(f_values[1])
		return [initial_interval[1], initial_interval[1]]
	elseif iszero(f_values[2])
		return [initial_interval[2],initial_interval[2]]
	end
	
	# Only a 2-element array is needed. Also, sort the interval [a,b] s.t. f(a) > 0 > f(b) to standardize the iteration algorithm.
	reduced_interval = Vector{promote_type(eltype(initial_interval), Float64)}(undef,2)
	if f_values[1] > f_values[2]
		reduced_interval[1] = initial_interval[1]
		reduced_interval[2] = initial_interval[2]
	else
		reduced_interval[1] = initial_interval[2]
		reduced_interval[2] = initial_interval[1]
	end


	for i in 1:iterations
		center = mean(reduced_interval)
		if f(center) > 0
			reduced_interval[1] = center
		elseif iszero(f(center)) # Finding the zero should end the iteration immediately -> return trivial interval.
			return [center,center]
		else
			reduced_interval[2] = center
		end
	end
	return reduced_interval
end



# Approximate the largest/smallest value s.t. f(x) < threshold.
# Stepping out to cross threshold by multiplying x_init with 2 followed by an interval bisection to approximate threshold point.
function find_largest_increase(f, x::Real, threshold::Real, iterations::Integer)
	
	if iszero(x)
		@warn("x=0 was set to x=1 to prevent an infinite loop.")
		x = 1
	end

	# Stepping out. While condition is skipped if already f(previous_x) < threshold
	previous_x = x
	while f(x) > threshold
		previous_x = x
		x *= 2
	end

	# Without stepping out, there is no interval with sign change for an interval bisection.
	if x == previous_x
		return x
	else
		# Approaching the threshold.
		return mean(interval_bisection(x-> f(x)-threshold,[previous_x,x],iterations))
	end
end

# setindex without mutation (convenience function to increase readability)
function setindex(A,x,ind)
	return setindex!(deepcopy(A),x,ind)
end



# index_shift_generator(f,[a_1,...,a_n],1) generates a function of the shift parameter g(s) = f([a_1 + s,a_2,...,a_n]), etc..
function index_shift_generator(f::Function, argument::AbstractArray,index)
	return function(s::Real)
		return f(setindex(argument,s .+ argument[index],index))
	end
end


# shift_samples([a_1,...,a_n],S,1) creates a vector of samples [[max(a_1 + s,0),...,a_n] for s in S], etc..
function shift_samples(parameters, shifts, index)
	samples = Matrix{promote_type(eltype(parameters),Float64, eltype(shifts[1]))}(undef, length(shifts)+1, length(parameters))
	for i in 1:size(samples,1)-1
		samples[i,:] = parameters
		samples[i,index] = max.(samples[i,index] .+ shifts[i],0)
	end
	samples[end,:] = parameters
	return samples
end



# shift_data(f, A, bin, levels, log_density, steps, bisections)
	# Calculate `thresholds` from levels (different for log_density = true and log_density = false)
	# Use smallest threshold to determine the largest shift `max_shift` (n=bisections interval bisections)
	# Determine largest negative shift `min_shift`, s.t. all values in bin are shifted below 0.

	# Create shift `samples` using `n=steps` shifts from 0 to max_shift and from 0 to min_shifts.
	# Calculate `density_values` by applying f.(samples)
	# Recalculate `thresholds` from best density value.

	# Return `samples`, `density_values`, `thresholds`
function shift_data(objective::Function, parameters, bin, sorted_levels, log_density::Bool, steps::Integer, bisections::Integer)
	thresholds = get_thresholds(sorted_levels, objective(parameters), log_density)

	f = index_shift_generator(objective,parameters,bin)
	max_increase = find_largest_increase(f,1,thresholds[1],bisections)
	samples = shift_samples(parameters,vcat(LinRange(0,max_increase,steps)...,LinRange(0,-maximum(parameters[bin]),steps)...),bin)
	objective_values = [objective(samples[i,:]) for i in 1:size(samples,1)]
	thresholds = get_thresholds(sorted_levels,maximum(objective_values),log_density)

	return samples, objective_values, thresholds
end






# For each threshold, fill the bins of the boundaries with the min/max of the samples that are admissible (determined by objective_values and the threshold).
function fill_boundaries!(lower, upper, thresholds, samples, objective_values, bin)
	for i in eachindex(thresholds)
		index_list = findall(x-> x >= thresholds[i], objective_values)
		lower[i,bin] = [minimum(samples[index_list,b]) for b in bin]
		upper[i,bin] = [maximum(samples[index_list,b]) for b in bin]
	end
	return lower, upper
end



# Dispatch version for single index bin.
function fill_boundaries!(lower, upper, thresholds, samples, objective_values, bin::Integer)
	for i in eachindex(thresholds)
		index_list = findall(x-> x >= thresholds[i], objective_values)
		lower[i,bin] = minimum(samples[index_list,bin]) 
		upper[i,bin] = maximum(samples[index_list,bin]) 
	end
	return lower, upper
end