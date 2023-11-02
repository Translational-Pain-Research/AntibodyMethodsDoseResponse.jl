@testset "DoseResponseUncertainty" begin

	@testset "DoseResponseUncertainty constructor" begin
		# Concentrations and responses (lower, upper) must have the same length.
		@test_throws DimensionMismatch DoseResponseUncertainty([0.5],[1,2],[1 2 3],[1 2 3])
		@test_throws DimensionMismatch DoseResponseUncertainty([0.5],[1,2,3],[1 2],[1 2 3])
		@test_throws DimensionMismatch DoseResponseUncertainty([0.5],[1,2,3],[1 2 3],[1 2])

		# Number of levels and number lower/upper bounds must match.
		@test_throws DimensionMismatch DoseResponseUncertainty([0.5],[1,2],[1 2 ; 1 2],[1 2 ; 1 2])
		@test_throws DimensionMismatch DoseResponseUncertainty([0.5,1],[1,2],[1 2],[1 2 ; 1 2])
		@test_throws DimensionMismatch DoseResponseUncertainty([0.5,1],[1,2],[1 2 ; 1 2],[1 2])

		# Argument order should be (levels, concentrations, lower, upper, lower_offset, upper_offset).
		@test DoseResponseUncertainty([0.5],[1,2],[3 4],[5 6]).levels == [0.5]
		@test DoseResponseUncertainty([0.5],[1,2],[3 4],[5 6]).concentrations == [1,2]
		@test DoseResponseUncertainty([0.5],[1,2],[3 4],[5 6]).lower == [3 4]
		@test DoseResponseUncertainty([0.5],[1,2],[3 4],[5 6]).upper == [5 6]
	end










	@testset "DoseResponseUncertainty from EpitopeUncertainty" begin
		@testset "Global uncertainty (default bin)" begin
			# Convenience function to obtain responses from weights.
			function test_response(grid,weights,concentrations,model = accumulation_model, offset = 0)
				temp_grid = deepcopy(grid)
				import_weights!(temp_grid,weights)
				return DoseResponseResult(temp_grid,concentrations,model = model, offset = offset).responses
			end



			# Setup concentrations and the K_τ grid.
			concentrations = LogRange(1e-10,1e-2,10)
			grid = create_grid(LogRange(1e-10,1e-2,5))

			# Default bins (one bin containing all indices).
			for offset in [0,2]
				# Create EpitopeUncertainty from random samples.
				samples = iszero(offset) ? rand(4,10) : rand(5,10) # Additional parameter if offset is used.
				e_uncertainty = EpitopeUncertainty(samples,levels = [0.5,1], offset = !iszero(offset))


				# Setup test conditions.

				# DoseResponseUncertainty with accumulation_model.
				dr_uncertainty = DoseResponseUncertainty(grid,e_uncertainty,concentrations)
				# DoseResponseUncertainty with langmuir_model.
				dr_uncertainty_lang = DoseResponseUncertainty(grid,e_uncertainty,concentrations, model = langmuir_model)

				# Collect data to iterate tests for both models.
				dr_uncertainties = [dr_uncertainty, dr_uncertainty_lang]
				e_uncertainties = [e_uncertainty,e_uncertainty]
				models = [accumulation_model,langmuir_model]
				
				for i in eachindex(models)
					dr_uncertainty = dr_uncertainties[i]
					e_uncertainty = e_uncertainties[i]
					model = models[i]
					
					# Test that levels and concentrations are passed correctly.
					@test dr_uncertainty.levels == e_uncertainty.levels
					@test dr_uncertainty.concentrations == concentrations

					# Using single bin containing all indices evaluates lower/upper bound by plugging in the lower/upper bounds of the EpitopeUncertainty object.
					if isnothing(e_uncertainty.lower_offset)
						@test dr_uncertainty.lower[1, :] == test_response(grid,e_uncertainty.lower[1,:],concentrations, model)
						@test dr_uncertainty.lower[end, :] == test_response(grid,e_uncertainty.lower[end,:],concentrations, model)
						@test dr_uncertainty.upper[1, :] == test_response(grid,e_uncertainty.upper[1,:],concentrations, model)
						@test dr_uncertainty.upper[end, :] == test_response(grid,e_uncertainty.upper[end,:],concentrations, model)
					else
						@test dr_uncertainty.lower[1, :] == test_response(grid,e_uncertainty.lower[1,:],concentrations, model, e_uncertainty.lower_offset[1])
						@test dr_uncertainty.lower[end, :] == test_response(grid,e_uncertainty.lower[end,:],concentrations, model, e_uncertainty.lower_offset[end])
						@test dr_uncertainty.upper[1, :] == test_response(grid,e_uncertainty.upper[1,:],concentrations, model, e_uncertainty.upper_offset[1])
						@test dr_uncertainty.upper[end, :] == test_response(grid,e_uncertainty.upper[end,:],concentrations, model, e_uncertainty.upper_offset[end])
					end

				end
			end
		end




		@testset "Bin-wise uncertainty" begin 

			# Setup data.
			concentrations = LogRange(1e-10,1e-2,10)
			grid = create_grid(LogRange(1e-10,1e-2,5))
			dr = DoseResponseResult(grid,concentrations)
			data = FittingData(dr.concentrations,dr.responses)
			# Bin is set to [1].
			e_uncertainty = EpitopeUncertainty(data,grid,[1],levels = [0.5,1])

			# DoseResponseUncertainty tests all bins and keeps the smallest/largest values for the respective response bounds (lower and upper).
			# To simplify testing, only a single bin [1] is chosen.
			# Then, expected responses can be obtained by importing weights into the grid, where only a single component is replaced (by the corresponding EpitopeUncertainty bound).


			# Setup the test uncertainty.
			dr_uncertainty = DoseResponseUncertainty(grid,e_uncertainty,concentrations, bins = [1])
			
			# Test that levels and concentrations are passed correctly.
			@test dr_uncertainty.levels == e_uncertainty.levels
			@test dr_uncertainty.concentrations == concentrations
			
			# Temporary gird to import altered weights into to obtain responses from.
			# Allows to avoid accidental mutation of the original grid.
			temp_grid = deepcopy(grid)
			# Test both levels [0.5,1]
			for i in 1:2
				temp_grid = import_weights!(temp_grid,setindex!(export_weights(grid),e_uncertainty.lower[i,1],1))
				# See choice of bin comment above for construction of expected lower bound.
				@test dr_uncertainty.lower[i,:] == DoseResponseResult(temp_grid,concentrations).responses

				temp_grid = import_weights!(temp_grid,setindex!(export_weights(grid),e_uncertainty.upper[i,1],1))
				# See choice of bin comment above for construction of expected upper bound.
				@test dr_uncertainty.upper[i,:] == DoseResponseResult(temp_grid,concentrations).responses
			end
		end
	end


end