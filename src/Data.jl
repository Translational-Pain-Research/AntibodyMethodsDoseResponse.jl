####################################################################################################
# Methods for measurement data.
####################################################################################################






# Test and normalize data
####################################################################################################
export dose_response_check, normalize_data, normalize_data!

# Internal function to test that a collection consists of positive, finite, real values.
function regular_positive_numbers_check(numbers)
	# Test indexablity before individual element tests are performed (otherwise an unspecific getindex method error is thrown).
	if !(applicable(eachindex,numbers))
		throw(ArgumentError("Indexable collection expected, but got $(typeof(numbers))."))
	end
	for ind in eachindex(numbers)
		if !(typeof(numbers[ind]) <: Real)
			throw(ArgumentError("Real numbers expected, but got $(typeof(numbers[ind])) at index $ind."))
		elseif isnan(numbers[ind])
			throw(ArgumentError("Numbers must not be `NaN`. Got $(numbers[ind]) at index $ind."))
		elseif isinf(numbers[ind])
			throw(ArgumentError("Numbers must not be `±Inf`. Got $(numbers[ind]) at index $ind."))
		elseif numbers[ind] < 0
			throw(DomainError(numbers[ind], "Numbers must not be negative. Got $(numbers[ind]) at index $ind."))
		end
	end
end



"""
	dose_response_check(fitting_data::FittingData)
Test if the `FittingData` object satisfies properties for dose-response data:

* only real numbers
* numbers must be positive
* no `NaN` or `Inf` in the data set
"""
function dose_response_check(fitting_data::FittingData)
	regular_positive_numbers_check(fitting_data.independent)
	regular_positive_numbers_check(fitting_data.dependent)
end





"""
	normalize_data!(fitting_data::FittingData; offset::Real = 0, reference::Union{Nothing,Real} = nothing)
Normalize `fitting_data` object by mutation.

* `reference = nothing`: The reference signal to normalize the responses to. If `reference = nothing`, the maximal response is used a reference point.
* `offset = 0`: Signal offset to be subtracted (applies to both the responses and the reference point).

The normalization of the `responses` and `errors` are:

	new_responses = (responses - offset) / (reference - offset)
	errors =  errors / (reference - offset)
"""
function normalize_data!(fitting_data::FittingData; offset::Real = 0, reference::Union{Nothing,Real} = nothing)
	dose_response_check(fitting_data)

	# Offset tests.
	if offset < 0
		throw(DomainError(offset,"The offset must be a positive number."))
	elseif isnan(offset)
		throw(DomainError(offset,"The offset must not be NaN."))
	end


	# Reference tests.
	if !isnothing(reference)
		if reference <= 0
			throw(DomainError(reference,"The reference must be a positive number."))
		elseif isinf(reference)
			throw(DomainError(reference,"The reference must not be Inf."))
		elseif isnan(reference)
			throw(DomainError(reference,"The reference must not be NaN."))
		end
	end

	responses = fitting_data.dependent

	if isnothing(reference)
		reference = maximum(responses)
	end

	if offset > minimum(responses) 
		@warn("Offset to large, reduced the offset to avoid negative responses!")
		offset = minimum(responses)
	end

	if offset >= reference
		throw(DomainError(offset,"The offset value must be smaller than the reference value."))
	end

	if !iszero(offset)
		responses = responses .- offset
		# Assuming the reference point also contains the offset component.
		reference -= offset
	end
	

	fitting_data.dependent = responses ./ reference
	fitting_data.errors = fitting_data.errors ./ reference

	FittingObjectiveFunctions.constructor_checks_fitting_data(fitting_data.independent,fitting_data.dependent,fitting_data.errors,fitting_data.distributions)
	return fitting_data
end



"""
	normalize_data!(fitting_data::FittingData; offset::Real = 0, reference::Union{Nothing,Real} = nothing)
Return `FittingData` object with normalized response data and errors.

* `reference = nothing`: The reference signal to normalize the responses to. If `reference = nothing`, the maximal response is used a reference point.
* `offset = 0`: Signal offset to be subtracted (applies to both the responses and the reference point).

The normalization of the `responses` and `errors` are:

	new_responses = (responses - offset) / (reference - offset)
	errors =  errors / (reference - offset)
"""
function normalize_data(fitting_data::FittingData; args...)
	new_data = deepcopy(fitting_data)
	return normalize_data!(new_data; args...)
end



















# Depletion correction.
####################################################################################################

export  scale_bound, simple_depletion_correction, simple_depletion_correction!



"""
	scale_bound(fitting_data::FittingData)
Get the upper bound for the scale factor `β` between the responses `r_i` and the (initial) antibody concentrations `a_i`, s.t.  `a_i - β*r_i ≥ 0`.
"""
function scale_bound(fitting_data::FittingData)
	# Improper data, e.g. negative values, could lead to an infinite while loop.
	dose_response_check(fitting_data)
	scale =  minimum(fitting_data.independent ./ fitting_data.dependent)

	# Division may lead to floating point errors. Decrease scale until a_i - β r_i ≥ 0
	while minimum(fitting_data.independent .- scale .* fitting_data.dependent ) <0
		scale -= eps(typeof(scale))
		# If floating point error leads to negative scale values, return 0.
		if scale < 0
			return 0
		end
	end
	
	return scale
end





# Internal function to get the correction values (additive shifts).
function get_depletion_correction(fitting_data::FittingData,scale::Real)
	dose_response_check(fitting_data)
	
	# Proper scale factors are positive finite numbers.
	if isnan(scale)
		throw(ArgumentError("Scale factor must not be `NaN`."))
	elseif  isinf(scale)
		throw(ArgumentError("Scale factor must not be `±Inf`."))
	elseif  scale < 0
		throw(DomainError(scale,"Scale must not be negative."))
	end

	corrections = - scale .* fitting_data.dependent 
	
	# Corrected concentrations must not become negative, otherwise the scale factor was too large.
	for c in fitting_data.independent .+ corrections
		if c < 0
			throw(DomainError(scale,"Scale factor to large, resulting in negative concentrations."))
		end
	end
	
	return corrections
end





# Non-mutating depletion corrections.

"""
	simple_depletion_correction(fitting_data::FittingData,scale::Real)
Return depletion-corrected `FittingData` object

The concentrations are corrected to `concentration - scale * response`.
"""
function simple_depletion_correction(fitting_data::FittingData,scale::Real)
	corrections = get_depletion_correction(fitting_data,scale)
	return FittingData(fitting_data.independent .+ corrections,fitting_data.dependent,fitting_data.errors, distributions = fitting_data.distributions)
end



"""
	simple_depletion_correction(fitting_data::FittingData)
Return depletion-corrected `FittingData` object using the largest possible scale factor.

The concentrations are corrected to `concentration - scale * response`.
"""
function simple_depletion_correction(fitting_data::FittingData)
	return simple_depletion_correction(fitting_data,scale_bound(fitting_data))
end





# Non-mutating depletion corrections.

"""
	simple_depletion_correction!(fitting_data::FittingData,scale::Real)
Depletion-correct `fitting_data` by mutation.

The concentrations are corrected to `concentration - scale * response`.
"""
function simple_depletion_correction!(fitting_data::FittingData,scale::Real)
	corrections = get_depletion_correction(fitting_data,scale)
	fitting_data.independent = fitting_data.independent .+ corrections
	FittingObjectiveFunctions.constructor_checks_fitting_data(fitting_data.independent,fitting_data.dependent,fitting_data.errors,fitting_data.distributions)
	return fitting_data
end



"""
	simple_depletion_correction!(fitting_data::FittingData)
Depletion-correct `fitting_data` by mutation, using the largest possible scale factor.

The concentrations are corrected to `concentration - scale * response`.
"""
function simple_depletion_correction!(fitting_data::FittingData)
	return simple_depletion_correction!(fitting_data,scale_bound(fitting_data))
end

