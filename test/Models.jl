@testset "Models" begin

	@testset "grid-tests for all models" begin
		# Setup to test all models in a loop.
		models = [langmuir_model,accumulation_model,langmuir_inv_const_model,accumulation_inv_const_model]
		names = ["langmuir_model","accumulation_model","langmuir_inv_const_model","accumulation_inv_const_model"]

		for (model,name) in zip(models,names)
			@testset "grid-tests" begin
				# Test all dispatch versions here.

				# K_τ must not be negative.
				@test_throws DomainError langmuir_model(create_grid([-1,2,3]))

				# K_τ must not be NaN.
				@test_throws ArgumentError langmuir_model(create_grid([NaN,2,3]))
				@test_throws ArgumentError langmuir_model(create_grid([1,2,NaN]))

				# K_τ must not be ±Inf.
				@test_throws ArgumentError langmuir_model(create_grid([-Inf,2,3]))
				@test_throws ArgumentError langmuir_model(create_grid([1,2,Inf]))

				# Weights must be real numbers.
				@test_throws ArgumentError langmuir_model(create_grid([1,2,3],initial_weight = [1,2]))
			end
		end
	end





	@testset "function value test for all models" begin

		# Here λ = [1,1] and c_i = [1.5,2.5], i.e. l_i = [1,2] and r_i = [2,3] for all tests.

		@testset "langmuir_model" begin
			model_1, λ_1 = langmuir_model(create_grid([1,2,3]))

			# Function values given by ∑_i λ_i*a/(r_i-l_i) * ln((a+r_i)/(a+l_i)).
			for a in 1:5
				@test model_1.model(a,λ_1) == a*log((a+2)/(a+1))+ a*log((a+3)/(a+2))
			end
		end

		@testset "accumulation_model" begin
			model_1, λ_1 = accumulation_model(create_grid([1,2,3]))

			# Function values given by ∑_i λ_i*(1-exp(-a/c_i)).
			for a in 1:5
				@test model_1.model(a,λ_1) == (1-exp(-a/1.5)) + (1-exp(-a/2.5))
			end
		end

		@testset "langmuir_inv_const_model" begin
			model_1, λ_1 = langmuir_inv_const_model(create_grid([1,2,3]))

			# Function values given by ∑_i λ_i* ( 1+1/(a*(r_i-l_i)) * log((a*l_i + 1) / (a*r_i +1)) ).
			for a in 1:5
				@test model_1.model(a,λ_1) == (1+1/a * log((a*1 + 1)/(a*2 + 1))) + (1+1/a * log((a*2 + 1)/(a*3 + 1)))
			end
		end

		@testset "langmuir_inv_const_model" begin
			model_1, λ_1 = accumulation_inv_const_model(create_grid([1,2,3]))

			# Function values given by ∑_i λ_i* ( 1+1/(a*(r_i-l_i)) * (exp(-a*r-i) - exp(-a*l_i)) ).
			for a in 1:5
				@test model_1.model(a,λ_1) == (1+1/a * (exp(-a*2) - exp(-a*1))) + (1+1/a * (exp(-a*3) - exp(-a*2)))
			end
		end
	end




















	@testset "DoseResponseResult constructor" begin
		# Data must contain only real numbers.
		@test_throws ArgumentError DoseResponseResult([1,2,3+3im],[1,2,3])
		@test_throws ArgumentError DoseResponseResult([1,2,3],[1,2,3+3im])

		# Data must consist of positive numbers.
		@test_throws DomainError DoseResponseResult([1,-2,3],[1,2,3])
		@test_throws DomainError DoseResponseResult([1,2,3],[1,-2,3])

		# Data must not contain NaN or Inf.
		@test_throws ArgumentError DoseResponseResult([1,NaN,3],[1,2,3])
		@test_throws ArgumentError DoseResponseResult([1,Inf,3],[1,2,3])
		@test_throws ArgumentError DoseResponseResult([1,2,3],[1,NaN,3])
		@test_throws ArgumentError DoseResponseResult([1,2,3],[1,Inf,3])

		# Dimensions of data must match.
		@test_throws DimensionMismatch DoseResponseResult([1,2,3],[1,2])

	end



	@testset "DoseResponseResult grid tests" begin
		# Test if K_τ grid check is called (using negative number).
		@test_throws DomainError DoseResponseResult(create_grid([-1,2,3]),[1,2,3])

		# Test if concentrations are checked with regular_positive_numbers_check (using negative number).
		@test_throws DomainError DoseResponseResult(create_grid([1,2,3]),[-1,1])

		# Test return type.
		@test typeof(DoseResponseResult(create_grid([1,2,3]),[1,2,3])) <: DoseResponseResult

		# Test that concentrations are passed correctly:
		@test DoseResponseResult(create_grid([1,2,3]),[1,2,3,4]).concentrations == [1,2,3,4] 
	end



	@testset "DoseResponseResult model selection" begin
		grid = create_grid([1,2,3])
		weights = export_weights(grid)
		concentrations = [1,2,3,4]

		model= accumulation_model(grid)[1].model
		# accumulation_model is the default
		@test DoseResponseResult(grid,concentrations,).responses == [model(c,weights) for c in concentrations]
		@test DoseResponseResult(grid,concentrations,offset = 1).responses == [model(c,weights) + 1 for c in concentrations]

		model= langmuir_model(grid)[1].model
		@test DoseResponseResult(grid,concentrations,model = langmuir_model).responses == [model(c,weights) for c in concentrations]
		@test DoseResponseResult(grid,concentrations,model = langmuir_model,offset = 1).responses == [model(c,weights) + 1 for c in concentrations]

	end
end