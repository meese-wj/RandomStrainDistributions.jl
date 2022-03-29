using RandomStrainDistributions
using Test

names(RandomStrainDistributions)

@testset "PhysicalVectors.jl" begin

    @testset "In-Place Operations" begin
        A = Vector2D(1, 2)
        @test equal(A, A)
        
        add!(A,A)
        ans = Vector2D(2, 4)
        @test equal(A, ans)
        
        A = Vector2D(1, 2)
        B = Vector2D(3., 4.)
        subtract!(A, B)
        @test equal(A, Vector2D(-2., -2.))
        
        A = Vector2D(1, 2)
        multiply!(A, 2.)
        @test equal(A, Vector2D(2., 4.))
        
        A = Vector2D(2, 4)
        divide!(A, 2.)
        @test equal(A, Vector2D(1., 2.))
        
        A = Vector2D(-100., 0.)
        normalize!(A)
        @test equal(A, Vector2D(-1., 0.))
    end

    @testset "Normal Operations" begin
        A = Vector2D(1., 2.)
        B = Vector2D(3., 4.)
        @test equal(A + B, Vector2D(4., 6.))

        @test A ⋅ B == 11 && dot(A, B) == A ⋅ B
        @test equal(2 * A, Vector2D(2., 4.)) && equal(2 * A, A * 2)
        @test equal(A / 2, Vector2D(0.5, 1.))
        @test magnitude2(B) == 25. && magnitude(B) == 5.
        @test equal(normalize(Vector2D(0., 20.)), Vector2D(0., 1.))
    end
end

@testset "RandomStrainDistributions.jl" begin
    
end
