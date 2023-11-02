@testset "Adaptive fitting" begin


	@testset "AdaptiveOptions" begin
		# Test default values.
		default_options = AdaptiveOptions()

		@test default_options.name == "Adaptive optimization"
		@test default_options.show_progress
		@test default_options.iterations == 1
		@test default_options.model == accumulation_model
		@test isnothing(default_options.offset)
		@test default_options.objective == :lsq
		@test default_options.prior_generator == AntibodyMethodsDoseResponse.default_prior_generator
		@test isnothing(default_options.distribution_derivatives)
		@test default_options.prior_gradient_generator == AntibodyMethodsDoseResponse.default_prior_gradient_generator
		@test default_options.block_variation == log_area_scaled_variation
		@test default_options.selection == maximum


		# Test keyword matching.
		custom_options = AdaptiveOptions(
			name = "name",
			show_progress = false,
			iterations = 10,
			model = sum,
			offset = 1,
			objective = :log_posterior,
			prior_generator = FittingObjectiveFunctions.uniform_prior, # not sensible (no prior generator), only for testing.
			distribution_derivatives = "a", # not sensible, only for testing.
			prior_gradient_generator = "b", # not sensible, only for testing.
			block_variation = area_scaled_variation,
			selection = minimum
		)

		@test custom_options.name == "name"
		@test custom_options.show_progress == false
		@test custom_options.iterations == 10
		@test custom_options.model == sum
		@test custom_options.offset == 1
		@test custom_options.objective == :log_posterior
		@test custom_options.prior_generator == FittingObjectiveFunctions.uniform_prior
		@test custom_options.distribution_derivatives == "a"
		@test custom_options.prior_gradient_generator == "b"
		@test custom_options.block_variation == area_scaled_variation
		@test custom_options.selection == minimum

	end






	@testset "AdaptiveResult" begin
		# Test field names.
		@test hasfield(AdaptiveResult,:result)
		@test hasfield(AdaptiveResult,:grid)
		@test hasfield(AdaptiveResult,:optimizer)
		@test hasfield(AdaptiveResult,:objective_value)
		@test hasfield(AdaptiveResult,:time)
	end
















	# Create dose-response data set with previously tested methods.
	grid = create_grid(LogRange(1e-10,1e-2,40))
	approximate_density!(grid,x-> FittingObjectiveFunctions.normal_distribution(x,1e-6,1e-5), volume_normalization = true)
	dr_result = DoseResponseResult(grid,LogRange(1e-10,1e-2,10))
	data = FittingData(dr_result.concentrations,dr_result.responses)










	@testset "get_objective" begin

		# Test that the offset argument, not the offset option is used.

		# Both offsets are nothing.
		objective, parameters, gradient = DR.get_objective(grid,data,nothing,AdaptiveOptions(objective=:lsq))
		@test length(parameters) == 39 # parameter length equals grid length, thus no offset is used.

		# offset argument is nothing, offset option is not.
		objective, parameters, gradient = DR.get_objective(grid,data,nothing,AdaptiveOptions(objective=:lsq,offset = 1))
		@test length(parameters) == 39 # parameter length equals grid length, thus no offset is used.

		# offset argument is 2, offset option is nothing.
		objective, parameters, gradient = DR.get_objective(grid,data,2,AdaptiveOptions(objective=:lsq))
		@test parameters[end]==2
		# offset argument is 2, offset option is 1.
		objective, parameters, gradient = DR.get_objective(grid,data,2,AdaptiveOptions(objective=:lsq,offset = 1))
		@test parameters[end]==2



		# Setup objective functions manually for testing.
		model, λ = accumulation_model(grid)
		lsq = lsq_objective(data,model)
		∇lsq = lsq_gradient(data,model)
		p = posterior_objective(data,model)
		lp = log_posterior_objective(data,model)
		∇lp = log_posterior_gradient(data,model,(x,m,σ)-> 2*(m-x)/σ)

		# Test lsq objective.
		objective, parameters, gradient = DR.get_objective(grid,data,nothing,AdaptiveOptions(objective=:lsq))
		@test λ == parameters
		@test prod(objective(a*λ) == lsq(a*λ) for a in 1:10)
		@test prod(gradient(a*λ,a*λ) == ∇lsq(a*λ,a*λ) for a in 1:10)

		# Test posterior objective.
		# prior_generator should take (centers, volumes, offset) and return a generated function (in this case uniform_prior regardless of the arguments).
		unifrom_prior_generator(centers,volumes,offset) = FittingObjectiveFunctions.uniform_prior
		objective, parameters, gradient = DR.get_objective(grid,data,nothing,AdaptiveOptions(objective=:posterior, prior_generator = unifrom_prior_generator))
		@test λ == parameters
		@test prod(objective(a*λ) == -p(a*λ) for a in 1:10)
		@test isnothing(gradient)
		
		# Test log_posterior objective without gradient.
		# AdaptiveOptions by default uses the correct log_uniform_prior. 
		objective, parameters, gradient = DR.get_objective(grid,data,nothing,AdaptiveOptions(objective=:log_posterior))
		@test λ == parameters
		@test prod(objective(a*λ) == -lp(a*λ) for a in 1:10)
		@test isnothing(gradient)

		# Test log_posterior objective with gradient.
		objective, parameters, gradient = DR.get_objective(grid,data,nothing,AdaptiveOptions(objective=:log_posterior,distribution_derivatives = (x,m,σ)-> 2*(m-x)/σ))
		@test λ == parameters
		@test prod(objective(a*λ) == -lp(a*λ) for a in 1:10)
		@test prod(gradient(a*λ,a*λ) == -∇lp(a*λ,a*λ) for a in 1:10)



		# Test if the model-option is used correctly (using langmuir_model here):
		model, λ = langmuir_model(grid)
		lsq = lsq_objective(data,model)
		∇lsq = lsq_gradient(data,model)

		objective, parameters, gradient = DR.get_objective(grid,data,nothing,AdaptiveOptions(objective=:lsq,model = langmuir_model))
		@test λ == parameters
		@test prod(objective(a*λ) == lsq(a*λ) for a in 1:10)
		@test prod(gradient(a*λ,a*λ) == ∇lsq(a*λ,a*λ) for a in 1:10)
	end






















	@testset "apply_result" begin

		# parameters without offset (length equals grid length).
		test_parameters_no = ones(length(grid))
		# parameters with offset (grid length + 1).
		test_parameters_o = vcat(ones(length(grid))...,5)


		@testset "Offset application and errors" begin
			# Avoid accidental mutation of grid.
			test_grid = deepcopy(grid)

			# Parameters with offset should fail (offset = nothing is the default of AdaptiveOptions)
			@test_throws DimensionMismatch DR.apply_result!(test_grid,data,test_parameters_o, AdaptiveOptions())
			# Parameters without offset should work.
			returned_grid, dr_result, parameters = DR.apply_result!(test_grid,data,test_parameters_no, AdaptiveOptions())
			# Without offset, the last parameter element should be the last element of the grid weights (1).
			@test parameters[end] == 1

			# Parameters without offset should fail (offset=10 explicitly called in AdaptiveOptions).
			@test_throws DimensionMismatch DR.apply_result!(test_grid,data,test_parameters_no, AdaptiveOptions(offset = 10))
			# Parameters with offset should work.
			returned_grid, dr_result, parameters = DR.apply_result!(test_grid,data, test_parameters_o, AdaptiveOptions(offset = 10))
			# With offset, the last parameter element should be the offset from the input parameters (5) not from the options (10).
			@test parameters[end] == 5
		end
		


		# Dummy grid to obtain models from.
		temp_grid = deepcopy(grid)
		import_weights!(temp_grid, test_parameters_no)
		default_model = accumulation_model(temp_grid)[1].model
		explicit_model = langmuir_model(temp_grid)[1].model



		@testset "Application of results without offset" begin
			for (i,model) in enumerate([default_model,explicit_model])
				
				# Test mutation of grid without mutating the original gird.
				# Inside loop to reset the grid (otherwise mutation is not tested properly)
				test_grid = deepcopy(grid)

				# Explicit model needs to be defined in the options.
				if i == 1
					returned_grid, dr_result, parameters = DR.apply_result!(test_grid,data, test_parameters_no, AdaptiveOptions())
				else
					returned_grid, dr_result, parameters = DR.apply_result!(test_grid,data, test_parameters_no, AdaptiveOptions(model = langmuir_model))
				end

				# Test mutation of grid / returned grid.
				@test returned_grid == test_grid
				
				# Test if results are applied correctly.
				@test dr_result.concentrations == data.independent
				# Define responses manually for the test.
				@test dr_result.responses == [model(c,test_parameters_no) for c in data.independent]
				# Grid weights are applied correctly 
				@test export_weights(returned_grid) == test_parameters_no
				# Parameters should be returned without mutation.
				@test parameters == test_parameters_no
			end
		end




		@testset "Application of results with offset" begin
			for (i,model) in enumerate([default_model,explicit_model])
				
				# Test mutation of grid without mutating the original gird.
				# Inside loop to reset the grid (otherwise mutation is not tested properly)
				test_grid = deepcopy(grid)

				# Explicit model needs to be defined in the options.
				if i == 1
					returned_grid, dr_result, parameters = DR.apply_result!(test_grid,data, test_parameters_o, AdaptiveOptions(offset = 10))
				else
					returned_grid, dr_result, parameters = DR.apply_result!(test_grid,data, test_parameters_o, AdaptiveOptions(offset = 10, model = langmuir_model))
				end

				# Test mutation of grid / returned grid.
				@test returned_grid == test_grid
				
				# Test if results are applied correctly.
				@test dr_result.concentrations == data.independent
				# Define responses manually for the test. Again offset argument (5) should be used, not the options offset (10).
				@test dr_result.responses == [model(c,test_parameters_no)+5 for c in data.independent]
				# Grid weights are applied correctly (only first 1:end-1 elements because of the offset).
				@test export_weights(returned_grid) == test_parameters_o[1:end-1]
				# Parameters should be returned without mutation.
				@test parameters == test_parameters_o
			end
		end

	end























	@testset "adaptive_dose_response_fit" begin
		# Dummy function for minimization. Just halve all parameters.
		function minimizer(f,∇f,λ)
			return λ ./ 2
		end


		@testset "Grid" begin
			# Test iterations and immutability of input grid (using the effect of the dummy minimization, halving all parameters).
			# Iterations are tested by the power of 2 and immutability is tested by using the sum of the original gird divided by 2^n.
			@test sum(adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 1)).grid) ≈ sum(grid)/2
			@test sum(adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 4)).grid) ≈ sum(grid)/(2^4)


			# Test that there is no redundant refinement after the last optimization.
			# 1 optimization iteration = no refinement (grid length remains 39).
			@test length(adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 1)).grid) == 39
			# 1 optimization iteration = 1 refinement step (grid length increases by 1).
			@test length(adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 2)).grid) == 40
		end

		

		@testset "DoseResponseResult" begin
			# Test default model.
			result = adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 2))
			# Responses created from the result-grid should match the result-responses.
			@test (result.result).responses == DoseResponseResult(result.grid,data.independent).responses

			# Test explicit model.
			result = adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 2,model = langmuir_model))
			# Responses created from the result-grid should match the result-responses.
			@test (result.result).responses == DoseResponseResult(result.grid,data.independent,model = langmuir_model).responses


			# The dummy minimizer function also divides the offset by 2.
			# This also tests the effect of the offset option.
			result = adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 2, offset = 2))
			# Responses created from the result-grid should match the result-responses. offset = 1 since minimizer halved the original offset (2).
			@test (result.result).responses == DoseResponseResult(result.grid,data.independent, offset = 1).responses
		end



		# Optimizer as in the returned parameter array that optimizes the objective function.
		@testset "Optimizer and objective" begin
			# Dummy prior generator to check if the prior option is used.
			test_prior_generator(centers,volumes, offset) = λ-> sum(λ)


			result = adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(objective = :log_posterior, prior_generator = test_prior_generator))

			# Use previously tested get_objective to get the objective function.
			objective = DR.get_objective(result.grid,data,nothing,AdaptiveOptions(objective = :log_posterior,prior_generator = test_prior_generator))[1]

			# Returned objective value should match the result of applying the objective to the returned optimizer.
			# This tests the internally used objective function, the returned optimizer and the returned objective value.
			@test objective(result.optimizer) == result.objective_value



			# Use offset = 1 to test that result.optimizer is the optimizing parameter and not just the grid weights.
			result = adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(objective = :log_posterior, offset = 1,prior_generator = test_prior_generator))

			# Again, use previously tested get_objective to get the objective function.
			objective = DR.get_objective(result.grid,data,1,AdaptiveOptions(objective = :log_posterior,prior_generator = test_prior_generator))[1]

			# As before returned objective value should match the result of applying the objective to the returned optimizer.
			@test objective(result.optimizer) == result.objective_value
		end


		@testset "Time" begin
			# Internal time should be enclosed by 0 and the external time (started before and stopped after internal time measurement).
			start_time = time()
			result = adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 2))
			maximal_time = time() - start_time
			@test 0 < result.time < maximal_time
		end

		
		@testset "Derivative options" begin

			# Test derivatives by checking if errors are thrown for nonsensical derivatives.
			∂Φ(x,m,σ) = log(-1) # Should throw DomainError
			∇p(grad,λ) = Int64(0.3) # Should throw InexactError

			# Define minimizer that actually calls the gradient (otherwise no error will be thrown).
			gradient_minimizer(f,∇f,λ) = ∇f(λ,λ)

			# Not setting the prior_gradient generator leads to a log-likelihood objective. Since distribution_derivatives are passed the gradient is still calculated, only using ∂Φ.
			@test_throws DomainError adaptive_dose_response_fit(grid,data,gradient_minimizer, options = AdaptiveOptions(objective = :log_posterior,distribution_derivatives = ∂Φ))

			# Set up the prior_gradient_generator with ∇p, which should throw a MethodError.
			test_gradient_generator(centers,volumes,offset) = ∇p
			# Use non-failing distribution_derivatives to check if ∇p fails (InexactError thrown).
			@test_throws InexactError adaptive_dose_response_fit(grid,data,gradient_minimizer, options = AdaptiveOptions(objective = :log_posterior,prior_gradient_generator = test_gradient_generator, distribution_derivatives = (x,m,σ)-> 1))
		end



		@testset "Refine options" begin
			# Test refine options (block_variation and selection) by checking if errors are thrown.
			error_variation(args...) = log(-1) # Should throw DomainError
			error_selection(variations) = Int64(0.3) # Should throw InexactError

			# block_variation is used (DomainError).
			@test_throws DomainError adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 2, block_variation = error_variation))
			# selection is used (InexactError).
			@test_throws InexactError adaptive_dose_response_fit(grid,data,minimizer, options = AdaptiveOptions(iterations = 2, selection = error_selection))
		end


	end
end