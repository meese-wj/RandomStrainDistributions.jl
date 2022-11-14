using Test
using BenchmarkTools

@testset "PhysicalVectors.jl" begin

    @testset "Normal Operations" begin
        A = Vector2D(1., 2.)
        B = Vector2D(3., 4.)
        @test isequal(A + B, Vector2D(4., 6.))

        @test A ⋅ B == 11 && dot(A, B) == A ⋅ B
        @test isequal(2 * A, Vector2D(2., 4.)) && isequal(2 * A, A * 2)
        @test isequal(A / 2, Vector2D(0.5, 1.))
        @test magnitude2(B) == 25. && magnitude(B) == 5.
        @test isequal(normalize(Vector2D(0., 20.)), Vector2D(0., 1.))
    end

    @testset "Memory allocations" begin
        
        A = Vector2D(10., 23)
        B = Vector2D(3., 42)

        bm = @benchmark $(A + B)
        @test bm.allocs == zero(Int)
        
        bm = @benchmark $(A - B)
        @test bm.allocs == zero(Int)
        
        bm = @benchmark $(rand() * A)
        @test bm.allocs == zero(Int)
        
        bm = @benchmark $(A / rand())
        @test bm.allocs == zero(Int)
        
        bm = @benchmark $(A ⋅ B)
        @test bm.allocs == zero(Int)
        
        bm = @benchmark unit($A)
        @test bm.allocs == zero(Int)
        
        bm = @benchmark normalize($A)
        @test bm.allocs == zero(Int)
        
        bm = @benchmark magnitude($A)
        @test bm.allocs == zero(Int)
        
        bm = @benchmark magnitude2($A)
        @test bm.allocs == zero(Int)
        
        bm = @benchmark subtract_PBC($A, $B, 12)
        @test bm.allocs == zero(Int)
        
        bm = @benchmark subtract_PBC($A, $B, 12, 16)
        @test bm.allocs == zero(Int)

    end
end