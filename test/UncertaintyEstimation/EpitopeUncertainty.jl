@testset "EpitopeUncertainty" begin
	
	@testset "Constructor" begin
		# Arguments must have the same lengths.
		@test_throws DimensionMismatch EpitopeUncertainty([0.7], [1 2; 2 3],[1 2; 2 2])
		@test_throws DimensionMismatch EpitopeUncertainty([0.7,0.9], [2 3],[1 2 ; 2 2])
		@test_throws DimensionMismatch EpitopeUncertainty([0.7,0.9], [1 2; 2 3],[2 2])

		# Sizes of matrices for lower bounds and upper bounds must match.
		@test_throws DimensionMismatch EpitopeUncertainty([0.7], [1 2 3],[1 2])

		# Argument order should be (levels, lower, upper).
		@test EpitopeUncertainty([0.1],[0.2 0.2],[0.3 0.3]).levels == [0.1]
		@test EpitopeUncertainty([0.1],[0.2 0.2],[0.3 0.3]).lower == [0.2 0.2]
		@test EpitopeUncertainty([0.1],[0.2 0.2],[0.3 0.3]).upper == [0.3 0.3]

		# Offset argument tests

		# Length checks include offset argument.
		@test_throws DimensionMismatch EpitopeUncertainty([0.7,0.9], [1 2 ; 2 3],[1 2 ; 2 2],[1],[1,2])
		@test_throws DimensionMismatch EpitopeUncertainty([0.7,0.9], [1 2 ; 2 3],[1 2 ; 2 2],[1,2],[1])

		# Partial offset information is ignored.
		@test isnothing(EpitopeUncertainty([0.1],[0.2 0.2],[0.3 0.3]).lower_offset)
		@test isnothing(EpitopeUncertainty([0.1],[0.2 0.2],[0.3 0.3]).upper_offset)
		@test isnothing(EpitopeUncertainty([0.1],[0.2 0.2],[0.3 0.3],[0.4]).lower_offset)
		@test isnothing(EpitopeUncertainty([0.1],[0.2 0.2],[0.3 0.3],[0.4]).upper_offset)

		# Argument order should be (levels, lower, upper, lower_offset, upper_offset)
		@test EpitopeUncertainty([0.1],[0.2 0.2],[0.3 0.3],[0.4],[0.5]).lower_offset == [0.4]
		@test EpitopeUncertainty([0.1],[0.2 0.2],[0.3 0.3],[0.4],[0.5]).upper_offset == [0.5]
	end






	
	@testset "Bin-wise shift construction" begin
		# Internal functions are tested in detail -> Only property tests here.


		grid = create_grid(LogRange(1e-10,1e-2,10))
		dummy_data = FittingData([1,2,3,4],[1,2,3,4])

		# Levels must be probabilities (0 excluded to avoid infinite search for largest shift).
		@test_throws DomainError EpitopeUncertainty(dummy_data,grid,levels = 0)
		@test_throws DomainError EpitopeUncertainty(dummy_data,grid,levels = -1)
		@test_throws DomainError EpitopeUncertainty(dummy_data,grid,levels = 2)

		# Test that empty bins are skipped.
		@test EpitopeUncertainty(dummy_data,grid,[[]]).lower[1,:] == export_weights(grid)
		@test EpitopeUncertainty(dummy_data,grid,[[]]).upper[1,:] == export_weights(grid)


		# Convenience for the following tests.
		function setindex(A,x,ind)
			B = deepcopy(A)
			return setindex!(B,x,ind)
		end


		# Test all combinations of bin dimensions, objective functions and offset configurations.
		for bins in [collect(1:length(grid)),[[1,2],[3,5]]], objective_type in [:lsq, :posterior, :log_posterior], offset in [nothing,1]

			# Use tested methods (Data.jl and Models.jl) to create data.
			dr_result = DoseResponseResult(grid,LogRange(1e-10,1e-2,10), offset = isnothing(offset) ? 0 : offset)
			data = FittingData(dr_result.concentrations, dr_result.responses)

			# Set options and obtain objective function manually for the tests.
			options = AdaptiveOptions(objective = objective_type, offset = offset)
			minimization_objective, parameters = AntibodyMethodsDoseResponse.get_objective(grid, data,options.offset,options)[1:2]
			objective = λ -> -minimization_objective(λ)

			# Create EpitopeUncertainty.
			levels = [0.1,0.9,1.0,0.2,0.4,0.3,0.5,0.6,0.8,0.7,0.1,0.1] # Unordered with duplicates to test processing of levels.
			e = EpitopeUncertainty(data,grid, bins, levels = levels, options = options, steps = 10^2, bisections = 10)

			# Test that the levels are processed correctly.
			@test e.levels == sort(unique(levels))
			
			# Test-thresholds need different constructors depending on the offset.
			if isnothing(offset)
				thresholds = AntibodyMethodsDoseResponse.get_thresholds(sort(unique(levels)),objective(e.lower[end,:]),!(objective_type == :posterior))
			else
				thresholds = AntibodyMethodsDoseResponse.get_thresholds(sort(unique(levels)),objective(vcat(e.lower[end,:]...,offset)),!(objective_type == :posterior))
				# Test offset fields.
				@test Set(e.lower_offset) == Set([offset])
				@test Set(e.upper_offset) == Set([offset])
			end

			# Test properties of bounds (objective value must be larger than the corresponding threshold).
			for i in eachindex(thresholds)
				# setindex to avoid mutation of parameters from testing different bins.
				@test prod([objective(setindex(parameters,e.lower[i,bin],bin)) >= thresholds[i] for bin in bins])
				@test prod([objective(setindex(parameters,e.upper[i,bin],bin)) >= thresholds[i] for bin in bins])
			end
		end

	end















	@testset "Construction from samples" begin

		# Identical samples in different shapes.
		matrix_samples = hcat([[i,i,i] for i in 1:5]...)
		vector_samples = [[i,i,i] for i in 1:5]

		# Quantiles for testing level=0.5 (0.25 above and 0.25 below the median)
		l = quantile(collect(1:5),0.25)
		u = quantile(collect(1:5),0.75)

		# 1 level = 0.5 (use quantiles) | 2 levels = [0.5,1] (use quantiles for 0.5 and extrema for 1).
		lower_1 = [l l l]
		upper_1 = [u u u]
		lower_2 = [l l l ; 1 1 1]
		upper_2 = [u u u ; 5 5 5]

		# If offset=true, the last sample value is used for the offset bounds -> 1:end-1 for lower and upper.
		lower_reduced_1 = [l l]
		upper_reduced_1 = [u u]
		lower_reduced_2 = [l l ; 1 1]
		upper_reduced_2 = [u u ; 5 5]

		lower_offset_1 = [l]
		upper_offset_1 = [u]
		lower_offset_2 = [l,1]
		upper_offset_2 = [u,5]


		# Levels must be probabilities.
		@test_throws DomainError EpitopeUncertainty(matrix_samples,levels = 2)
		@test_throws DomainError EpitopeUncertainty(matrix_samples,levels=-0.1)
		@test_throws DomainError EpitopeUncertainty(matrix_samples,levels = [0.3,1.1])

		# Grid parameter length must be larger than 1.
		@test_throws DimensionMismatch EpitopeUncertainty([1,2,3,4,5],levels = 0.5)
		@test_throws DimensionMismatch EpitopeUncertainty([1 2 3 4 5],levels = 0.5)

		# Test both shapes for samples.
		for samples in [matrix_samples, vector_samples]


			# Test single levels.

			@test EpitopeUncertainty(samples,levels = 0.5).levels == [0.5]
			@test EpitopeUncertainty(samples,levels = 0.5).lower == lower_1
			@test EpitopeUncertainty(samples,levels = 0.5).upper == upper_1

			@test EpitopeUncertainty(samples,levels = 0.5,offset = true).levels == [0.5]
			@test EpitopeUncertainty(samples,levels = 0.5,offset=true).lower == lower_reduced_1
			@test EpitopeUncertainty(samples,levels = 0.5,offset =true).upper == upper_reduced_1
			@test EpitopeUncertainty(samples,levels = 0.5,offset=true).lower_offset == lower_offset_1
			@test EpitopeUncertainty(samples,levels = 0.5,offset =true).upper_offset == upper_offset_1



			# Test multiple levels.

			@test EpitopeUncertainty(samples,levels = [0.5,1]).levels == [0.5,1]
			@test EpitopeUncertainty(samples,levels = [0.5,1]).lower == lower_2
			@test EpitopeUncertainty(samples,levels = [0.5,1]).upper == upper_2

			@test EpitopeUncertainty(samples,levels = [0.5,1],offset = true).levels == [0.5,1]
			@test EpitopeUncertainty(samples,levels = [0.5,1],offset=true).lower == lower_reduced_2
			@test EpitopeUncertainty(samples,levels = [0.5,1],offset =true).upper == upper_reduced_2
			@test EpitopeUncertainty(samples,levels = [0.5,1],offset=true).lower_offset == lower_offset_2
			@test EpitopeUncertainty(samples,levels = [0.5,1],offset =true).upper_offset == upper_offset_2



			# Test sorting of multiple levels without offset.
			E_non_unique = EpitopeUncertainty(samples,levels = [0.5,0.5,1])
			E_not_sorted = EpitopeUncertainty(samples, levels = [1,0.5])

			@test E_non_unique.levels == E_not_sorted.levels
			@test E_non_unique.lower == E_not_sorted.lower
			@test E_non_unique.upper == E_not_sorted.upper

			# Test sorting of multiple levels with offset.
			E_non_unique = EpitopeUncertainty(samples,levels = [0.5,0.5,1], offset = true)
			E_not_sorted = EpitopeUncertainty(samples, levels = [1,0.5], offset = true)

			@test E_non_unique.levels == E_not_sorted.levels
			@test E_non_unique.lower == E_not_sorted.lower
			@test E_non_unique.upper == E_not_sorted.upper
			@test E_non_unique.lower_offset == E_not_sorted.lower_offset
			@test E_non_unique.upper_offset == E_not_sorted.upper_offset
		end

	end





end