cd(@__DIR__)
using Pkg
Pkg.activate(".")

using RandomStrainDistributions
using Profile
using PProf


const b0 = tetragonal_burgers_vectors[1]
const b1 = -1 * b0

strain_func(r::Vector2D{T}, r0) where T = b2g_shear(r - r0, b0) 
strain_dis(r::Vector2D, dis::Dislocation) = b2g_shear( r, dis )

rposition = Vector2D( 3., 3. )

Lx = 128; Ly = Lx
r0 = Vector2D(Lx ÷ 2 + 0.5, Ly ÷ 3 + 0.5)
r1 = Vector2D(Lx ÷ 8 + 0.5, Ly ÷ 8 + 0.5)

sq_idx = 80

val = PBCField(strain_func, rposition, r0, Lx, Ly)

Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1 PBCField(strain_func, rposition, r0, Lx, Ly)
results = Profile.Allocs.fetch()
PProf.Allocs.pprof(results; from_c = false)

while true end


