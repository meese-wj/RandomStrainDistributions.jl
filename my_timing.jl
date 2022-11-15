using RandomStrainDistributions


const b0 = tetragonal_burgers_vectors[1]
const b1 = -1 * b0

strain_funcT(r::Vector2D{T}, r0s) where T = b2g_shear(r - r0s[1], b0) + b2g_shear(r - r0s[2], b1) 

rposition = Vector2D( 3., 3. )

Lx = 128; Ly = Lx
r0 = Vector2D(Lx ÷ 2 + 0.5, Ly ÷ 3 + 0.5)
r1 = Vector2D(Lx ÷ 8 + 0.5, Ly ÷ 8 + 0.5)

sq_idx = 80

@time _left_edge(strain_funcT, sq_idx, rposition, (r0, r1), Lx, Ly)

