using Test
using BenchmarkTools
using RandomStrainDistributions

function PBC_lattice(ϕfunc::F, _dis::Dislocation2D{T}, _Lx, _Ly, tolerance::T = sqrt(eps())) where {F, T}
    val::T = zero(T)
    for (xdx, ydx) ∈ Iterators.product( UnitRange(1, _Lx), UnitRange(1, _Ly) )
        val += PBCField_value(ϕfunc, Vector2D( T(xdx), T(ydx) ), _dis, _Lx, _Ly, tolerance)
    end
    return val
end

@info "Testing PBCImageFields"
@testset "PBCImageFields Tests" begin
    
    println("  Testing convergence without allocations")
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

    # TODO: Fix this test
    # println("  Testing allocations in a lattice sweep")
    # @elapsed @testset "Lattice sweep without allocations" begin
    #     _Lx = 64
    #     _Ly = _Lx
    #     _r0 = Vector2D( _Lx ÷ 2 + 0.5, _Ly ÷ 2 + 0.5 )

    #     for _b0 ∈ tetragonal_burgers_vectors
    #         _dis = Dislocation2D( _b0, _r0 )

    #         tm = @timed PBC_lattice(b1g_shear, _dis, _Lx, _Ly)
    #         tm = @timed PBC_lattice(b1g_shear, _dis, _Lx, _Ly)
    #         # @test tm.bytes == zero(tm.bytes)
            
    #         tm = @timed PBC_lattice(b2g_shear, _dis, _Lx, _Ly)
    #         tm = @timed PBC_lattice(b2g_shear, _dis, _Lx, _Ly)
    #         # @test tm.bytes == zero(tm.bytes)
    #     end
    # end
    # println()

end