# Testing the details of some EpitopeUncertainty constructors is difficult. Thus, testing the auxiliary functions is helpful.
@testset "Auxiliary functions" begin

	@testset "interval_bisection" begin
		interval_bisection = AntibodyMethodsDoseResponse.interval_bisection
		
		# Test that opposite signs are required.
		@test_throws ArgumentError interval_bisection(x->x ,[1,2], 10^3)
		@test_throws ArgumentError interval_bisection(x->x ,[-1,-2], 10^3)

		# Hitting the zero (should immediately return the zero as point-interval).
		@test interval_bisection(x-> x, [0,2], 10) == [0,0]
		@test interval_bisection(x-> x, [2,0], 10) == [0,0]
		@test interval_bisection(x-> x, [-1,1], 10) == [0,0]

		# Test approximation of zero by interval-halving property.
		for i in [10,20,30]
			interval = interval_bisection(x-> x,[-1,2],i)
			@test abs(interval[2]-interval[1]) ≈ 3/2^i
		end
	end










	# find_largest_increase(f,x_init,threshold,iterations) approximates the largest/smallest value s.t. f(x) < threshold.
	# Stepping out to cross threshold by multiplying x_init with 2 followed by an interval bisection to approximate threshold point.
	@testset "find_largest_increase" begin
		find_largest_increase = AntibodyMethodsDoseResponse.find_largest_increase

		# Sign of initial point determines direction -> 0 will be automatically replaced and a warning should be thrown.
		@test_warn "x=0" find_largest_increase(x-> -x, 0,1,10^3)

		# If initial point already yields function value below the threshold the initial point should be returned.
		@test find_largest_increase(x-> -x, 2, 0,10^3) == 2

		# Test approximation of threshold by function property (linearity).
		# Hitting the exact value (multiple of 2 and 1/2 by combination of stepping out by doubling and stepping in by interval bisection).
		@test find_largest_increase(x-> -x,1,-64,10^3) == 64
		# Approximate value (no multiple of 2 and 1/2).
		@test find_largest_increase(x-> -x,1,-exp(2),10^3) ≈ exp(2)
	end






	# setindex without mutation.
	@testset "setindex" begin
		setindex = AntibodyMethodsDoseResponse.setindex

		original_array = collect(1:5)
		returned_array = setindex(original_array,[20,30],[2,3])
		@test original_array == collect(1:5)
		@test returned_array == [1,20,30,4,5]
	end




	# index_shift_generator(f,[a_1,...,a_n],1) generates a function of the shift parameter g(s) = f([a_1 + s,a_2,...,a_n]), etc..
	@testset "index_shift_generator" begin
		index_shift_generator = AntibodyMethodsDoseResponse.index_shift_generator

		# Summing the squares from 1:10 allows to determine at which index the shift was added (non-linearity (2+1)^2 + 3^2 != 2^2 + (3+1)^2 etc.).
		F(A) = sum(x->x^2,A)
		A = collect(1:10)

		# 1-dim index.
		f = index_shift_generator(F,A,1)
		@test f(1) == 4+ F(A[2:end]) # (1+1)^2
		@test f(11) == 12^2 + F(A[2:end]) # (1+11)^2

		# 2-dim index.
		f = index_shift_generator(F,A,[1,3])
		@test f(1) == 4+ 2^2 + 4^2 +  F(A[4:end]) # (1+1)^2 + 2^2 + (3+1)^2
		@test f(11) == 12^2 + 2^2 + 14^2 + F(A[4:end]) # (1+11)^2 + 2^2 + (3+11)^2
	end











	# shift_samples([a_1,...,a_n],S,1) creates a vector of samples [[max(a_1 + s,0),...,a_n] for s in S], etc..
	@testset "shift_samples" begin
		shift_samples = AntibodyMethodsDoseResponse.shift_samples

		original_array = ones(10)
		shifts = collect(-5:10)
		
		# 1-index bin.
		samples = shift_samples(original_array,shifts,1)

		# Original array should not be mutated.
		@test original_array == ones(10)
		
		for i in eachindex(shifts)
			# i-th sample should use the i-th shift.
			@test samples[i,:] == vcat(max(1+shifts[i],0), ones(9)...)
		end
		# shift_samples should always include the original array (trivial shift).
		@test samples[end,:] == original_array

		# 2-index bin.
		shifts = [[i,i+2] for i in -5:10]
		samples = shift_samples(original_array,shifts,[1,3])
		for i in eachindex(shifts)
			# i-th sample should use the i-th shift.
			@test samples[i,:] == vcat(max(1+shifts[i][1],0), 1, max(1+shifts[i][2],0) , ones(7)...)
		end
		# shift_samples should always include the original array (trivial shift).
		@test samples[end,:] == original_array
	end












	# shift_data(f, A, bin, levels, log_density, steps, bisections)
		# Calculate `thresholds` from levels (different for log_density = true and log_density = false)
		# Use smallest threshold to determine the largest shift `max_shift` (n=bisections interval bisections)
		# Determine largest negative shift `min_shift`, s.t. all values in bin are shifted below 0.

		# Create shift `samples` using `n=steps` shifts from 0 to max_shift and from 0 to min_shifts.
		# Calculate `density_values` by applying f.(samples)
		# Recalculate `thresholds` from best density value.

		# Return `samples`, `density_values`, `thresholds`
	@testset "shift_data" begin
		shift_data = AntibodyMethodsDoseResponse.shift_data
		
		@testset "log_density = true" begin
			# Dummy function to analyze properties.
			lp = X -> sum(-x for x in X)
			λ = ones(5)
			# exp(-64) since the logarithm is taken if log_density = true.
			levels = [exp(-64),0.5,0.7]

			# 1-index bin (using Sets because concatenation of shift samples no longer ordered).
			samples, objective_values, thresholds = shift_data(lp,λ,1,levels,true,2,10^3)
			# 65 because shift is 64 and initial parameter value is 1.
			@test Set(samples[:,1]) == Set([0,1,65]) # steps = 2 leads to the shifts -1, 0, +64
			# No changes outside the bin.
			@test Set(samples[:,2:end]) == Set([1])
			# Only the first element was changed (cf. dummy function: sum(-x), λ = ones(5)).
			@test Set(objective_values) == Set([0,-1,-65] .- 4)
			# Thresholds should be calculated from the best function value.
			@test Set(thresholds) == Set(log.(levels) .+ lp([0,1,1,1,1]) )

			# 2-index bin.
			samples, objective_values, thresholds = shift_data(lp,λ,[1,5],levels,true,2,10^3)
			# 33 because shift is 32=64/2 and initial parameter value is 1.
			@test Set(samples[:,1]) == Set([0,1,33])
			@test Set(samples[:,5]) == Set([0,1,33])
			# No changes outside the bin.
			@test Set(samples[:,2:end-1]) == Set([1])
			# Only the first and the last element were changed (cf. dummy function: sum(-x), λ = ones(5)).
			@test Set(objective_values) == Set([0,-2,-66] .- 3)
			# Thresholds are calculated from best function value.
			@test Set(thresholds) == Set(log.(levels) .+ lp([0,1,1,1,0]) )
		end
		



		@testset "log_density = false" begin
			# Dummy function to analyze properties.
			# +105 to ensure positivity and allow for easy shift calculations (-5+105 = 100)
			p = X -> sum(-x for x in X) + 105 
			λ = ones(5)
			# level 1/10 -> single bin shift 90: (-4-(1+90)+105 = 10).
			levels = [1/10,0.5,0.7]

			# 1-index bin.
			samples, objective_values, thresholds = shift_data(p,λ,1,levels,false,2,10^3)
			# 91 because shift is 90 and initial parameter value is 1.
			@test Set(samples[:,1]) == Set([0,1,91])
			# No changes outside the bin.
			@test Set(samples[:,2:end]) == Set([1])
			# Only the first element was changed (cf. dummy function: sum(-x)+105, λ = ones(5)).
			@test Set(objective_values) == Set([21,20,-91+21] .+ 20*4)
			# Thresholds are calculated from best function value.
			@test Set(thresholds) == Set(levels .* p([0,1,1,1,1]) )

			# 2-index bin.
			samples, objective_values, thresholds = shift_data(p,λ,[1,5],levels,false,2,10^3)
			# 46 because shift is 45=90/2 and initial parameter value is 1.
			@test Set(samples[:,1]) == Set([0,1,46])
			@test Set(samples[:,5]) == Set([0,1,46])
			# No changes outside the bin.
			@test Set(samples[:,2:end-1]) == Set([1])
			# Only the first and the last element were changed (cf. dummy function: sum(-x)+105, λ = ones(5)).
			@test Set(objective_values) == Set([2*21,2*20,2*(-46+21)] .+ 20*3)
			# Thresholds are calculated from best function value.
			@test Set(thresholds) == Set(levels .* p([0,1,1,1,0]) )
		end
	end



	














	# fill_boundaries!(lower, upper, thresholds, samples, objective_values, bin)
		# For each threshold, fill the bins of the boundaries with the min/max of the samples that are admissible (determined by objective_values and the threshold).
	@testset "fill_boundaries!" begin
		fill_boundaries! = AntibodyMethodsDoseResponse.fill_boundaries!

		@testset "1-dim bin" begin
			lower = zeros(4,5)
			upper = zeros(4,5)

			# Manually construct samples, objective_values and threshold for the testing.
			f = x-> -abs(sum(x))
			samples = zeros(21,5) # zeros, s.t. only values in the bin matter.
			samples[:,2] .= collect(-10:10)
			objective_values = [f(samples[i,:]) for i in 1:size(samples,1)]
			thresholds = [-8,-6,-4,-2]
			
			new_lower, new_upper = fill_boundaries!(lower,upper,thresholds,samples, objective_values,2)

			# Test mutation.
			@test new_lower == lower
			@test new_upper == upper

			# Test values in the bin (x-> -abs(sum(x)), thus shifts should equal the thresholds up to the sign).
			@test lower[:,2] == [-8,-6,-4,-2]
			@test upper[:,2] == [8,6,4,2]

			# Test outside of the bin (should be unchanged, i.e. 0).
			@test Set([lower[i,j] for i in 1:4, j in 1:5 if j != 2]) == Set([0])
			@test Set([upper[i,j] for i in 1:4, j in 1:5 if j != 2]) == Set([0])
		end


		@testset "2-dim bin" begin
			lower = zeros(4,5)
			upper = zeros(4,5)

			# Manually construct samples, objective_values and threshold for the testing.
			f = x-> -abs(sum(x))
			samples = zeros(21,5) # zeros, s.t. only values in the bin matter.
			samples[:,2] .= collect(-10:10)
			samples[:,5] .= collect(-10:10)
			objective_values = [f(samples[i,:]) for i in 1:size(samples,1)]
			thresholds = 2 .* [-8,-6,-4,-2]
			
			new_lower, new_upper = fill_boundaries!(lower,upper,thresholds,samples, objective_values,[2,5])

			# Test mutation.
			@test new_lower == lower
			@test new_upper == upper

			# Test values in the bin (x-> -abs(sum(x)), thus shifts should equal the thresholds up to the sign).
			@test lower[:,2] == [-8,-6,-4,-2]
			@test lower[:,5] == [-8,-6,-4,-2]
			@test upper[:,2] == [8,6,4,2]
			@test upper[:,5] == [8,6,4,2]

			# Test outside of the bin (should be unchanged, i.e. 0).
			@test Set([lower[i,j] for i in 1:4, j in 1:5 if j != 2 && j != 5]) == Set([0])
			@test Set([upper[i,j] for i in 1:4, j in 1:5 if j != 2 && j != 5]) == Set([0])
		end
		

	end
end
