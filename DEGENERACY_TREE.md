# Degenerate corners of the compactified parameter space (Stage 1, Q2)

Parameters: e ∈ [0, 1], l = sin I ∈ [0, 1], ω ∈ [0, π/2], circle radius s ∈ [0, ∞]. The rows below are the places where
a plain Rolle certificate cannot work, because the polynomial vanishes identically or a root has multiplicity > 10.
Each needs its own chart.

| # | Corner | Physical meaning | Slack (numerical) | Status |
|---|---|---|---|---|
| D1 | e = 0, l = 0 (all s) | concentric coplanar circles (excluded) | max 8 nearby | chart C1 (e = ε cos θ, l² = ε sin θ) |
| D1a | D1 with s = 1 | coinciding conics (excluded) | max 8 (30,000 samples) | **resolved**: cone charts B1, B2, B3± in (e, θ = l²/e, τ) (gen_bu.py) |
| D1b | inside D1a: s = 1, l = e, ω = 0 at leading order | special near-coincident configuration | – | **resolved**: level-2 cones L2E, L2P, L2A±, L2T± |
| D2b | e → 0, s ~ e, l → 1, ω → 0 (inside C1p) | circular orbit at 90° with a tiny circle | – | **resolved**: cones D2E, D2P, D2K, D2T± in (e, p, 1 − q, τ − ½) |
| D2 | e = 0, s = 0 | circle shrunk to the focus of a circular orbit | max 8 | chart C1p (t = e²λ, l ≥ 0.4) |
| D2a | D1 ∩ D2 | | | chart C1pp (e = ε'' l², t = ε''² λ'') |
| D3 | e = 1, μ = q/R = 0 | parabola tip vs a huge circle (q ≪ R ≪ a) | max 10 (tight); on the face itself P = −D²B² has 12 real roots (6-fold at aphelion + three double roots of B²) for l ≥ l*(ω) | **resolved**: engine KQ (runs KQ1, KQ2, KQK1, KQK2, KQB) away from the face; on the face neighbourhood ε, μ ≤ 4·10⁻⁴ the weighted blow-up (τ, ε, μ) = (1, 2, 2) with leading part 16c²y²[4m² − c²(E + 2y²)²] (≤ 4 near zeros): charts FE, FM (near), FIP, FIM (gap), run KQF (outer ≤ 6) (C3_CORE.md §6) |
| D4 | D3 with P_z = 1 (aphelion direction on the circle's axis) | C3 core | P̂0 ≤ 8 real roots near aphelion + ≤ 2 away | **resolved**: C3 cones (near) + Φ± (intermediate) + KA (away) + level-2 zone charts (C3_CORE.md) |
| D4' | D3 with l ≈ 1, ω ∈ [3π/20, 9π/20] | the orbit crosses the circle's axis inside the aphelion arc (σ = 0) | – | **resolved**: engine KQ, run KQB (C3_CORE.md §5) |
| S | symmetric faces ω = 0, ω = π/2 (pitchfork surfaces, butterfly points) and the t → 1 six-fold curves | clusters of real roots on the unit circle | N_mult = 10 attained | handled by Rolle runs with J = 6, the winding certificate and the deep fallback, with no special chart; in the large-circle KQ runs the boxes are split in t first |
| D5 | e = 1 (parabola), μ finite | regular: f-chart face (parabola–circle ≤ 8 + 2 aphelion) | gap √2 | f-chart run F0 (t ≤ 0.5) and engine KQ (t ≥ 0.5) |

## Branch normal form near e = 0 (nearest-point branch g₋)

h₋ := g₋'/2 = e sin u (1 − e cos u) − s(ρ²)'/(2ρ)
     ≈ e(1 − s) sin u − (e²/2) sin 2u + (l²/2) sin 2(u + ω)  + higher order.

This has trigonometric degree 2, so at most 4 zeros per branch and N ≤ 8 near e = 0 at leading order. It degenerates
(vanishes identically) only when s = 1, l = e and ω = 0 at that order (D1b). At e = 0, l > 0 the count is exactly 8
(the 8 critical points of two inclined concentric circles); they are non-degenerate, with Hessian determinant
±4s² sin²I ≠ 0 (four minima/maxima and four saddles; corrected by audit R2-4).
