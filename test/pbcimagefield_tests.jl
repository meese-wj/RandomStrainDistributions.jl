using Test
using BenchmarkTools
using RandomStrainDistributions

@testset "PBCImageFields Tests" begin
    
    @testset "Convergence without allocations" begin
        _Lx = 128
        _Ly = _Lx
        _r0 = Vector2D( _Lx ÷ 2 + 0.5, _Ly ÷ 2 + 0.5 )

        for _b0 ∈ tetragonal_burgers_vectors
            _dis = Dislocation2D( _b0, _r0 )
            bm = @benchmark PBCField( $b1g_shear, $(_r0 + Vector2D( 0.5, 0. )), $_dis, $_Lx, $_Ly )
            @test bm.allocs == zero(bm.allocs)
            
            bm = @benchmark PBCField( $b2g_shear, $(_r0 + Vector2D( 0.5, 0.5 )), $_dis, $_Lx, $_Ly )
            @test bm.allocs == zero(bm.allocs)
        end
    end

end