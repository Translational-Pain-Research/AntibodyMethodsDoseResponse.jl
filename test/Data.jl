@testset "Measurement data" begin

	@testset "Dose-response data-check function" begin
		# Data must contain only real numbers.
		@test_throws ArgumentError dose_response_check(FittingData([1,2,3+3im],[1,2,3]))
		@test_throws ArgumentError dose_response_check(FittingData([1,2,3],[1,2,3+3im]))

		# Data must consist of positive numbers.
		@test_throws DomainError dose_response_check(FittingData([1,-2,3],[1,2,3]))
		@test_throws DomainError dose_response_check(FittingData([1,2,3],[1,-2,3]))

		# Data must not contain NaN or Inf.
		@test_throws ArgumentError dose_response_check(FittingData([1,NaN,3],[1,2,3]))
		@test_throws ArgumentError dose_response_check(FittingData([1,Inf,3],[1,2,3]))
		@test_throws ArgumentError dose_response_check(FittingData([1,2,3],[1,NaN,3]))
		@test_throws ArgumentError dose_response_check(FittingData([1,2,3],[1,Inf,3]))
	end





	@testset "Data normalization" begin
		# Test if dose_response_check is called (using negative number).
		@test_throws DomainError normalize_data(FittingData([1,-2,3],[1,2,3]))
		@test_throws DomainError normalize_data!(FittingData([1,-2,3],[1,2,3]))

		# Test referenced range (positive, finite number).
		@test_throws DomainError normalize_data(FittingData([1,2,3],[1,2,3]), reference = -1)
		@test_throws DomainError normalize_data(FittingData([1,2,3],[1,2,3]), reference = 0)
		@test_throws DomainError normalize_data!(FittingData([1,2,3],[1,2,3]), reference = -1)
		@test_throws DomainError normalize_data!(FittingData([1,2,3],[1,2,3]), reference = 0)
		@test_throws DomainError normalize_data(FittingData([1,2,3],[1,2,3]), reference = Inf)
		@test_throws DomainError normalize_data!(FittingData([1,2,3],[1,2,3]), reference = Inf)
		@test_throws DomainError normalize_data(FittingData([1,2,3],[1,2,3]), reference = NaN)
		@test_throws DomainError normalize_data!(FittingData([1,2,3],[1,2,3]), reference = NaN)

		# Test offset range (non-negative, finite number).
		@test_throws DomainError normalize_data(FittingData([1,2,3],[1,2,3]), offset = -1)
		@test_throws DomainError normalize_data!(FittingData([1,2,3],[1,2,3]), offset = -1)
		@test_throws DomainError normalize_data(FittingData([1,2,3],[1,2,3]), offset = NaN)
		@test_throws DomainError normalize_data!(FittingData([1,2,3],[1,2,3]), offset = NaN)

		# Test that the offset is limited by the minimal data response value.
		@test_warn "Offset" normalize_data(FittingData([1,2,3],[1,2,3]), offset = Inf)
		@test_warn "Offset" normalize_data!(FittingData([1,2,3],[1,2,3]), offset = Inf)
		@test_throws DomainError normalize_data(FittingData([1,2,3],[1,2,3]), offset = 0.5, reference = 0.5)
		@test_throws DomainError normalize_data!(FittingData([1,2,3],[1,2,3]), offset = 0.5, reference = 0.5)





		# Check that concentrations are not altered.
		@test normalize_data(FittingData([1,2,3],[1,2,3])).independent == [1,2,3]
		@test normalize_data(FittingData([1,2,3],[1,2,3]), reference = 2).independent == [1,2,3]
		@test normalize_data!(FittingData([1,2,3],[1,2,3])).independent == [1,2,3]
		@test normalize_data!(FittingData([1,2,3],[1,2,3]), reference = 2).independent == [1,2,3]
		# Test if scaling is normalized to 1.
		@test maximum(normalize_data(FittingData([1,2,3],[1,2,3])).dependent) ≈ 1
		@test maximum(normalize_data(FittingData([1,2,3],[1,2,3]), reference = 2).dependent) ≈ 1.5
		@test maximum(normalize_data!(FittingData([1,2,3],[1,2,3])).dependent) ≈ 1
		@test maximum(normalize_data!(FittingData([1,2,3],[1,2,3]), reference = 2).dependent) ≈ 1.5
		# Errors should be multiplied by scaling.
		@test normalize_data(FittingData([1,2,3],[1,2,3],[1,1,1])).errors[1] ≈ 1/3
		@test normalize_data(FittingData([1,2,3],[1,2,3],[1,1,1]), reference = 2).errors[1] ≈ 1/2
		@test normalize_data(FittingData([1,2,3],[1,2,3],[1,1,1]), reference = 2, offset = 0.5).errors[1] ≈ 1/(1.5)
		@test normalize_data!(FittingData([1,2,3],[1,2,3],[1,1,1])).errors[1] ≈ 1/3
		@test normalize_data!(FittingData([1,2,3],[1,2,3],[1,1,1]), reference = 2).errors[1] ≈ 1/2
		@test normalize_data!(FittingData([1,2,3],[1,2,3],[1,1,1]), reference = 2, offset = 0.5).errors[1] ≈ 1/(1.5)


		# Automatic offset correction (data - offset, then normalizing to reference - offset).
		@test maximum(normalize_data(FittingData([1,2,3],[1,2,3]),offset = Inf).dependent) ≈ 1
		@test maximum(normalize_data(FittingData([1,2,3],[1,2,3]),offset = Inf, reference = 2).dependent) ≈ 2
		@test maximum(normalize_data!(FittingData([1,2,3],[1,2,3]),offset = Inf).dependent) ≈ 1
		@test maximum(normalize_data!(FittingData([1,2,3],[1,2,3]),offset = Inf, reference = 2).dependent) ≈ 2
		# Minimal response must be 0 in all cases.
		@test minimum(normalize_data(FittingData([1,2,3],[1,2,3]),offset = Inf).dependent) ≈ 0
		@test minimum(normalize_data(FittingData([1,2,3],[1,2,3]),offset = Inf, reference = 2).dependent) ≈ 0
		@test minimum(normalize_data!(FittingData([1,2,3],[1,2,3]),offset = Inf).dependent) ≈ 0
		@test minimum(normalize_data!(FittingData([1,2,3],[1,2,3]),offset = Inf, reference = 2).dependent) ≈ 0


		# Fixed offset correction (data - offset, then normalizing to reference - offset).
		@test maximum(normalize_data(FittingData([1,2,3],[1,2,3]),offset = 0.5).dependent) ≈ 1
		@test maximum(normalize_data(FittingData([1,2,3],[1,2,3]),offset = 0.5, reference = 2).dependent) ≈ 2.5/(1.5)
		@test maximum(normalize_data!(FittingData([1,2,3],[1,2,3]),offset = 0.5).dependent) ≈ 1
		@test maximum(normalize_data!(FittingData([1,2,3],[1,2,3]),offset = 0.5, reference = 2).dependent) ≈ 2.5/(1.5)
		# Minimal response must be 0 in all cases.
		@test minimum(normalize_data(FittingData([1,2,3],[1,2,3]),offset = 0.5).dependent) ≈ 0.5/2.5
		@test minimum(normalize_data(FittingData([1,2,3],[1,2,3]),offset = 0.5, reference = 2).dependent) ≈ 0.5/1.5
		@test minimum(normalize_data!(FittingData([1,2,3],[1,2,3]),offset = 0.5).dependent) ≈ 0.5/2.5
		@test minimum(normalize_data!(FittingData([1,2,3],[1,2,3]),offset = 0.5, reference = 2).dependent) ≈ 0.5/1.5



		# Mutability tests
		fitting_data = FittingData([1,2,3],[1,2,3])
		normalize_data(fitting_data, offset = true)
		
		# fitting_data should not be mutated.
		@test fitting_data.independent == [1,2,3]
		@test fitting_data.dependent == [1,2,3]
		@test fitting_data.errors == ones(3)

		returned_data = normalize_data!(fitting_data,offset = true)
		# fitting_data should be mutated.
		@test fitting_data.independent == [1,2,3]
		@test fitting_data.dependent == [0,0.5,1]
		@test fitting_data.errors == 0.5 .* ones(3)
		# The mutated array should be returned.
		@test fitting_data == returned_data
	end




	@testset "Upper bound for scale factor" begin
		# Test if dose_response_check is called (using negative number).
		@test_throws DomainError scale_bound(FittingData([1,-2,3],[1,2,3]))
		# Bound must always be larger than 0.
		@test scale_bound(FittingData([1,2,3],[2,3,4])) >= 0
		# Constructed such that upper bound must be 1/0.5 = 2.
		@test scale_bound(FittingData([1,2,3],[0.5,0.5,0.5])) == 2
	end

















@testset "Simple depletion correction" begin
		# Test if dose_response_check is called (using negative number).
		@test_throws DomainError simple_depletion_correction(FittingData([1,-2,3],[1,2,3]),1)
		@test_throws DomainError simple_depletion_correction!(FittingData([1,-2,3],[1,2,3]),1)

		# Provided scale factor must not be Nan, ±Inf, or negative.
		@test_throws ArgumentError simple_depletion_correction(FittingData([1,2,3],[1,2,3]),NaN)
		@test_throws ArgumentError simple_depletion_correction!(FittingData([1,2,3],[1,2,3]),NaN)
		@test_throws ArgumentError simple_depletion_correction(FittingData([1,2,3],[1,2,3]),Inf)
		@test_throws ArgumentError simple_depletion_correction!(FittingData([1,2,3],[1,2,3]),Inf)
		@test_throws DomainError simple_depletion_correction(FittingData([1,2,3],[1,2,3]),-1)
		@test_throws DomainError simple_depletion_correction!(FittingData([1,2,3],[1,2,3]),-1)

		# Provided scale factor is too large.
		@test_throws DomainError simple_depletion_correction(FittingData([1,2,3],[1,2,3]),2)
		@test_throws DomainError simple_depletion_correction!(FittingData([1,2,3],[1,2,3]),2)

		
		fitting_data = FittingData([1,2,3],[1,2,3])

		# Correction is calculated by c -α * r, where c is the concentration, α is the scale factor and r is the response.
		@test simple_depletion_correction(FittingData([1,2,3],[1,2,3]),1).independent == zeros(3)
		# Only the concentrations are changed.
		@test simple_depletion_correction(FittingData([1,2,3],[1,2,3]),1).dependent == [1,2,3]
		@test simple_depletion_correction(FittingData([1,2,3],[1,2,3]),1).errors == ones(3)

		# Automatic scale should be 1.
		@test simple_depletion_correction(FittingData([1,2,3],[1,2,3])).independent == zeros(3)
		@test simple_depletion_correction(FittingData([1,2,3],[1,2,3])).dependent == [1,2,3]
		@test simple_depletion_correction(FittingData([1,2,3],[1,2,3])).errors == ones(3)


		# fitting_data should not be mutated by simple_depletion_correction.
		@test fitting_data.independent == [1,2,3]
		@test fitting_data.dependent == [1,2,3]
		@test fitting_data.errors == ones(3)

		# simple_depletion_correction! should mutate fitting_data.
		returned_data = simple_depletion_correction!(fitting_data,1)

		@test fitting_data.independent == zeros(3)
		@test fitting_data.dependent == [1,2,3]
		@test fitting_data.errors == ones(3)
		# The mutated FittingData object should be returned.
		@test returned_data == fitting_data


		# simple_depletion_correction! with automatic scale should also mutate fitting_data.
		fitting_data = FittingData([1,2,3],[1,2,3])
		returned_data = simple_depletion_correction!(fitting_data,1)

		@test fitting_data.independent == zeros(3)
		@test fitting_data.dependent == [1,2,3]
		@test fitting_data.errors == ones(3)
		# The mutated FittingData object should be returned.
		@test returned_data == fitting_data
	end
end