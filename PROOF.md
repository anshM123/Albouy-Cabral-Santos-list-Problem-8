# Problem 8, second question (Albouy–Cabral–Santos 2012): circle and ellipse — at most 10 stationary points

**Theorem** (computer-assisted). Let a circle of radius R and an ellipse of eccentricity e ∈ (0, 1) in ℝ³ have a common focus, located
at the centre of the circle (inclination and orientation arbitrary). Then the squared distance
d²(u, v) = |X(u) − C(v)|² between a point of the ellipse and a point of the circle has at most **10** stationary
points on the torus. The bound is attained.

This is the second question of Problem 8 (G. F. Gronchi) in A. Albouy, H. E. Cabral, A. A. Santos, *Some problems on
the classical n-body problem*, Celest. Mech. Dyn. Astr. 113 (2012) 369–375: *"If only one of the eccentricities is 0,
can this upper bound [12] be reduced to 10?"* Answer: **yes**. The canonical statement and the excluded cases are in
STATEMENT.md. The cases with infinitely many stationary points (coinciding conics, concentric coplanar circles) have
e = 0 and are outside the question.

**Sharpness.** The configuration below has exactly 10 stationary points. Three independent counts agree:
- Pol has 10 real roots, by 60-digit root finding;
- the branch derivatives have 6 + 4 sign changes;
- a 2-D Newton search on the torus finds 10 points (code/verify_mp.py, code/ce_core.py).

The configuration is e = 0.673, sin I = 0.992, ω = 0.338, R/a = 0.7772.

<!-- STATUS -->
**Certification status (2026-10-23 18:30): complete.** Every run family in RESULTS.md is complete (`code/check_all.py`: ALL COMPLETE). Every recorded split tree replays exactly, and every open or failed box is covered by a continuation run or by a documented chart. Sampled leaves of every run file re-certify with the final code (verify/), and the final audit (audit/AUDIT_15.md) re-certifies every tau- and disc-winding leaf. 60-digit spot checks (verify/) cover sampled leaves of every whole-circle family and C3 chart, and the disc-winding leaves of KQF and KQF2. The other outer-count leaves of KQF and KQF2, and the D2 charts, are covered by the audit's 100-digit spot checks.
<!-- /STATUS -->

---

## 1. Reduction to a trigonometric polynomial (STATEMENT.md; audited in AUDIT_1 §2)

**Coordinates.** Semi-major axis a = 1 (scaling), circle C(v) = s(cos v, sin v, 0) with s = R/a, and
ellipse X(u) = (cos u − e)P + β sin u Q, where β = √(1 − e²).

**Symmetry reduction.** ω ↦ −ω, ω ↦ ω + π and I ↦ π − I reduce the parameters to l = sin I ∈ [0, 1] and
ω ∈ [0, π/2]. Then p = l sin ω and q = l cos ω.

**Branch functions.** For fixed u the stationary v are the nearest and farthest points of the circle, so the
stationary points are the critical points of g∓(u) = r² + s² ∓ 2sρ. Squaring gives

  Pol(u) = A²ρ² − s²((ρ²)')² = h₋h₊,  h∓ = ρ g∓'(u),

a trigonometric polynomial of degree ≤ 6, where A = (r²)', r = 1 − e cos u, and ρ² = r² − Z².

**Counting lemma.** Pol ≢ 0 for e > 0. The number N of stationary points is at most the number of real roots of Pol
on the circle counted with multiplicity, including the special points with A = (ρ²)' = 0 and with ρ = 0.

**Positive normalisations.** Every chart below certifies a function of the form c·(Pol ∘ φ). Here c > 0 is smooth
and φ is a diffeomorphism of the circle, or of the circle minus one point together with a separate count there. Such
a function has the same real roots with the same multiplicities as Pol, so the counting lemma applies to it. In detail:
- compactification t = s²/(1+s²): f_t = (1 − t)Pol;
- the f-chart (true anomaly): P_f = F₋F₊ with F∓ = (D³σ/(μ(1+e)))·g∓'(f) in the R-normalisation (C3_CORE.md §1),
  hence P_f = c(f)·Pol(u(f)) with c > 0;
- blow-up charts: division by e³, e²θ, e²η, e³ν, ν¹⁰ or ν¹⁰η⁶ (all > 0 on the physical domain);
- tan-half-angle form: Q(τ) = P(f)(1 + τ²)⁶ with τ = tan((f − π)/2).

## 2. The certificate (Rolle runs)

**Lemma (Rolle).** Let F be C^∞ on a closed interval I. If F^(j) has no zero on I, then F has at most j zeros on I,
counted with multiplicity.

**How a box is certified.** For a parameter box, the circle (or an interval) is covered by overlapping arcs. On every
arc the engine certifies a strict sign of F^(j) for some j ≤ 6, uniformly over the box (AUDIT_1 §4).
- The enclosures use the chart's coefficient polynomials, evaluated at the box centre, plus a mean-value bound for
  the parameter variation and a Taylor bound in the anomaly. All floating-point errors are padded; see
  FLOAT_RIGOUR.md.
- **Second-order variant** (the a-chart bulk run from its 3826th initial box on; engine KQ: away arcs, winding and
  near count; C3 chart Ca). The parameter variation of
  each derivative is also bounded by Σ_v r_v|g_v| + Σ_{v≤w} q_vw|H_vw| + Σ_k|T_k|R_k. Here:
  - g and H are the gradient and Hessian at the box centre, exact by the chain rule;
  - q_vv = r_v²/2 and q_vw = r_v r_w;
  - R_k bounds the third-order Taylor remainder of coefficient k, from absolute-value polynomials of the third
    derivatives over the box, plus the evaluation errors.

  The engine uses the smaller of this bound and the mean-value bound; both are valid. Code: code/nb_cert2.py
  (a- and f-chart; driver option bb_run.py --order2), code/c3/make_c3_2.py (polynomial charts); generators
  code/gen_nb2a.py, code/gen_nb2.py and code/c3/gen_c3_2.py;
  FLOAT_RIGOUR.md §5.
- A run of consecutive arcs on which one F^(j) keeps one strict sign contains at most j zeros.
- `dp_count` / `dp_linear` take the minimum, over partitions into such runs, of the sum of the j's. That minimum
  is a certified upper bound for the number of zeros, for every parameter in the box.
- A box is a **leaf** when this bound is ≤ 11 for a count over the whole circle (then the true count is ≤ 10 by the
  parity lemma below), or ≤ the stated partial target for a count over part of the circle.

**Lemma (parity).** For a real trigonometric polynomial F ≢ 0, the number of real zeros on the circle counted with
multiplicity is even. Zeros of odd multiplicity are exactly the sign changes, and a periodic function changes sign an
even number of times.

**Lemma (winding certificate).** Let F(u) = Σ_{|k|≤6} c_k e^{iku} be real (so c_{−k} = conj(c_k)), F ≢ 0, and let
Q(z) = Σ c_k z^{k+6}.
- Q is self-inversive: z¹²·conj(Q(1/z̄)) = Q(z). Hence its 12 roots on the Riemann sphere (a degree drop puts roots
  at ∞, matched by roots at 0) are symmetric under z ↦ 1/z̄, with multiplicity.
- The real zeros of F correspond, with multiplicity, to the roots of Q on |z| = 1.
- If Q has no root on |z| = ρ < 1 and W ≥ 1 roots in |z| < ρ, then it also has W roots in |z| > 1/ρ. So F has at most
  12 − 2W ≤ 10 real zeros, counted with multiplicity.

**How the winding certificate is applied (code/wind_cert.py).** With L = ln ρ:

  Q(ρe^{iθ}) = ρ⁶e^{6iθ}·G(θ),  G = Σ cosh(kL)(a_k cos kθ + b_k sin kθ) + i·Σ sinh(kL)(a_k sin kθ − b_k cos kθ).

- On every arc of an overlapping cover of the circle, the engine certifies a strict half-plane condition, uniformly
  over the box: Re G > 0, Re G < 0, Im G > 0 or Im G < 0. The enclosures are those of §2, and G is affine in t.
- Then G ≠ 0 on |z| = ρ throughout the box, and the winding number is constant on the box.
- The winding number is evaluated at the box centre from G at the arc boundaries. Consecutive boundaries lie in one
  labelled arc, so each argument step is < π in modulus.
- W = 6 + wind(G) is the number of roots in |z| < ρ (argument principle). The choice of ρ, from the roots at the box
  centre, affects only success, never validity.

**Lemma (tau-winding certificate; C3_CORE.md §8).** Let Q(τ) = P(f)(1 + τ²)⁶ with τ = tan((f − π)/2). Then the number of real zeros of P on the circle, counted with multiplicity, is 12 minus the number of non-real zeros of Q. So if Q has no zero on a circle in the open upper half-plane for any parameter of the box, and has a zero inside it at the box centre, P has at most 10 real zeros on the whole box. The same holds with Q replaced by:
- its reverse G(t_p) = t_p¹²Q(−1/t_p) (chart KBPe), since −1/t_p is non-real when t_p is;
- the chart polynomials Q(νy)/ν¹⁰ (C3 chart Ca) or Q(sy)/s⁶ (face charts FE, FM), since ν, s > 0 at every physical
  point.

Engine KQ (option --tw) and these charts use this certificate with leaf code '7'. The contour is stored, and it is
re-certified by code/bb_replay.py (KQ) and code/recert_tw_c3.py (c3 charts).

**Lemma (disc-winding outer count; C3_CORE.md §9).** With Q as above, suppose that for every parameter of the box Q has
no zero on the circle |τ| = T, and that Q has W zeros inside it at the box centre. Then P has at most 12 − W real zeros,
counted with multiplicity, on the outer region {|τ| ≥ T} ∪ {f = 0}, at every point of the box. The outer runs KQF and
KQF2 (option --disc) use this with T = T_X and W ≥ 6, as leaf code '8' with outer bound 6. The contour is fixed by the
family's T_X, which check_all.py pins.

**Why both certificates are used.** The winding certificate does not resolve real roots at all. So it does not
degrade at root clusters on the unit circle (pitchforks and six-fold points), which is exactly where Rolle runs need
tiny boxes. The bb_run driver tries Rolle runs first, then the winding certificate, then deep Rolle runs on small
boxes; engine KQ tries the winding certificate first, then its counts (§4).

The driver splits boxes deterministically and records the full split tree, one character per box:
- winding leaves are coded '1'..'6' (W);
- their radii L are stored in order, so bb_replay.py can re-certify them.

bb_replay.py replays every tree and checks that the leaves tile the initial box exactly.

**Run families.** A run family is a base run on a grid of initial boxes, plus possibly continuation runs on the
open boxes of the previous level (an initial box whose time cap expired, or a failed leaf). A continuation run may use a
different split scale; for example, it may split the circle-size coordinate t first in the large-circle KQ runs, or
the ratio K first in the C3 chart Ca. The continuation levels of KQK1, KQK2 and KQB near the face e = 1, μ = 0 use
the split rule 'face' of code/splitrule.py (recorded in the run header and replayed exactly). The split rule affects
only efficiency.

`code/check_all.py` replays every file with its own recorded scale. It checks:
- that the base run's initial boxes are exactly its grid;
- that every open box of the base run is the root of a complete record of the continuation run;
- that the continuation run itself has no open box.

## 3. The atlas (coverage lemma)

All parameter points with e ∈ (0, 1), l ∈ [0, 1], ω ∈ [0, π/2], s ∈ (0, ∞) are covered as follows. The runs are in
code/launch_runs.sh, and every boundary below overlaps its neighbour.

| region | chart / engine | run(s) | target |
|---|---|---|---|
| e ∈ [0.199, 0.83617] (η = asin e ∈ [0.2003, 0.99025]) | a-chart (eccentric anomaly), all l, ω, t | bulk (its initial boxes with η ≥ 0.99025 are excluded) | ≤ 11 (parity) |
| e ∈ [0.835, 0.9501], t ≤ 0.5 | f-chart | F1 | ≤ 11 (parity) |
| e ∈ [0.835, 0.9501], t ∈ [0.5, 0.9975] | engine KQ (below) | KQ0 | ≤ 10 (winding), ≤ 11 (parity) |
| e ∈ [0.835, 0.9501], t ≥ 0.9975 | engine KQ | KQ0T | ≤ 10 (winding), ≤ 11 (parity) |
| e ≤ 0.2, l ≥ 0.4 | C1p: t = e²τ/(1−τ), τ ∈ [0,1] (all s > 0) | C1p | ≤ 10 |
| e ≤ 0.2, l < 0.4, l² ≤ 2e | C1: e = ε cos θ, l² = ε sin θ (ε ≤ 0.26), minus the D1a block | C1 | ≤ 10 |
| e ≤ 0.2, l < 0.4, l² > 2e | C1pp: e = ε''l², t = ε''²τ/(1−τ) | C1pp | ≤ 10 |
| D1a block: ε ≤ 0.065, θ ≤ 0.1846, \|t − ½\| ≤ 0.125 | cones B1, B2, B3± in (e, θ = l²/e, τ = t − ½), cube e ≤ 0.07, θ ≤ 0.2, \|τ\| ≤ 0.13 (D1_CORNER below) | B1, B2, B3p, B3m | ≤ 10 |
| D1b block of B1 | level-2 cones L2E, L2P, L2A±, L2T± | L2* | ≤ 10 |
| failed leaves of C1p (D2b corner: e ≤ 6.107e-6, l ≥ 0.997986, ω ≤ 6.711e-4, \|τ − ½\| ≤ 2.14e-3) | cones D2E, D2P, D2K, D2T± on the cube e, p, 1 − q, \|τ − ½\| ≤ 5e-3, plus D2E2 (D2E extended to τ'/e ∈ [1.04, 1.28]); coverage rules below | D2*, D2E2 | ≤ 10 (D2E/P/K), ≤ 11 (parity; D2T±, D2E2) |
| e ≥ 0.9499, t ≤ 0.5 | f-chart | F0 | ≤ 11 (parity) |
| e ≥ 0.9499, t ∈ [0.5, 0.9975] | engine KQ: the winding certificate of the f-chart polynomial, or else the tan-half-angle polynomial in ε = 1 − e on \|f − π\| ≤ 0.9 (monomial basis around aphelion) + f-chart arcs away (code/kq_run.py) | KQ1 | ≤ 10 (winding), ≤ 11 (parity) |
| 0.9499 ≤ e ≤ 0.99, t ≥ 0.9975 | engine KQ | KQ2 | ≤ 11 (parity) |
| e ≥ 0.99, t ≥ 0.997456 (⊇ μ ≤ 0.0505), outside l ≥ 0.95 & ω ≥ 3π/20, minus the face neighbourhood N | engine KQ | KQK1, KQK2 | ≤ 11 (parity) |
| e ≥ 0.99, t ≥ 0.997456, l ≥ 0.95, ω ∈ [3π/20, 9π/20] (the band), minus N | engine KQ | KQB | ≤ 11 (parity) |
| face neighbourhood N: ε = 1 − e ≤ 4·10⁻⁴, μ ≤ 4·10⁻⁴, outside the C3 cell and the corner K_lw = {l ≥ 0.88, ω ≥ 1.29} | weighted blow-up (τ, ε, μ) = (1, 2, 2): near FE, FM; gap FIP, FIM; outer KQF (C3_CORE.md §6) | c3/FE, FM, FIP, FIM, KQF | 5 + 0 + 6 ≤ 11 (parity) |
| N ∩ K_lw with ε, μ ≤ SF2² (SF2 = 0.010005) | the same blow-up with s ≤ SF2; outer count from \|τ\| ≥ 3·SF2 | c3/FE, FM, FIP, FIM, KQF2 | 5 + 0 + 6 ≤ 11 (parity) |
| N ∩ K_lw, the rest (SF2 < max(√ε, √μ) ≤ 0.02) | engine KQ | KQS | ≤ 11 (parity) |
| e ≥ 0.99, μ ≤ 0.0505, l ≥ 0.95, ω ≥ 9π/20 (C3 cell) | near: C3 cones Ca, CK, CE, Cm (+ level 2); intermediate: PhiP, PhiM; away: KA | c3/*, KA | 9 + 0 + 2 ≤ 11 (parity) |

**Why the e ≤ 0.2 split is exhaustive.** If l ≥ 0.4, C1p applies. Otherwise either l² ≤ 2e or l² > 2e. In the
first case ε = √(e² + l⁴) < √(0.04 + 0.0256) < 0.26, so C1 applies. In the second, ε'' = e/l² < 0.5 ≤ 0.51 with
l < 0.4 < 0.52, so C1pp applies.

**Why the D2 charts cover the failed leaves of C1p.** Every D2 chart function is the C1p function divided by a
positive factor, so every chart counts the same zeros. The 659 failed leaves of C1p lie in the cube
e, p, κ = 1 − q, |τ − ½| ≤ 5e-3. Each point of the cube lies in the cone chart of its dominant coordinate, with all
ratios ≤ 1 ≤ 1.05.

D2E, D2P, D2K and D2Tm have no failed or open box. D2Tp has failed and open boxes (bound 12), all around
P = p/τ' ≈ 0, 2K + E ≈ 2. There its leading term s²·[4Ks² + (8K − 8)c² − 4Ec³] has a 4-fold zero at u = π. Each of
these boxes is harmless for one of two reasons (both checked in code/check_all.py):
1. Its smallest physical e exceeds the largest e of every failed leaf of C1p (6.11e-6). Then it contains no point of
   those leaves.
2. Its image in the D2E variables (ν_E = Eν, P/E, K/E, τ'/e = 1/E) lies in the domain of D2E where τ'/e ≤ 1.05, and
   in the domain of D2E2 where τ'/e ≥ 1.04. D2E2 is D2E with its affine variable τ'/e running over [1.04, 1.28], and
   it has no failed box. This case covers the one remaining box,
   (ν, E, P, K) ∈ [0, 0.12625] × [0.7875, 1.05] × [0, 0.2625] × [0.2625, 0.525].

**Why the e ≥ 0.95 split is exhaustive.** In the f-chart, t = κ²/(1+κ²) with κ = R/q and μ = q/R = √((1−t)/t). F0 and
KQ1 cover t ≤ 0.9975 for all e ∈ [0.9499, 1). For t ≥ 0.9975, KQ2 covers e ≤ 0.99. For e ≥ 0.99 and t ≥ 0.997456
(μ ≤ 0.0505024), KQK1 ∪ KQK2 ∪ KQB ∪ C3 cover everything; the C3 cell is μ ≤ 0.0505, and its sliver
0.0505 < μ ≤ 0.0505024 lies in KQ1 (t ≤ 0.9975). The boxes of KQK1, KQK2 and KQB that lie inside the face
neighbourhood N = {ε ≤ 4·10⁻⁴, μ ≤ 4·10⁻⁴} are pruned (leaf code 'p'; check_all.py re-checks the predicate) and are
covered by the face construction (C3_CORE.md §6): the face charts with KQF, and in the corner K_lw = {l ≥ 0.88, ω ≥ 1.29}
the face charts for s ≤ SF2 with KQF2 plus the whole-circle run KQS for SF2 < s ≤ 0.02. The C3-cell parts are pruned
('o') and covered by the C3 core.

**Why the seam at e ≈ 0.836 is exhaustive.** The bulk grid has η-edges at multiples of (1.25356 − 0.20034)/8; its
initial boxes with η ≥ 0.9902515840381316 (e ≥ 0.836164) are excluded (recorded in the run header) and the bulk
covers e ≤ 0.836164 for all l, ω, t. F1 (t ≤ 0.5), KQ0 (t ∈ [0.5, 0.9975]) and KQ0T (t ≥ 0.9975) cover
e ∈ [0.835, 0.9501] for all l, ω; F0, KQ1, KQ2, ... cover e ≥ 0.9499. Engine KQ uses exact polynomial identities
(Q(τ) = P(f)(1 + τ²)⁶ in ε = 1 − e, the f-chart polynomial, the winding lemma) and is valid for every e ∈ (0, 1)
with t ≥ 0.5; the a-chart bulk run had become slow above e ≈ 0.84, where the f-chart variables are better suited. The corner engine K and the band engine KB (CORNER_E1.md,
C3_CORE.md §5) were superseded by KQ and are not part of the final atlas.

**Coverage tests.** The auditor sampled 800,000 points over the actual run grids and exclusions and found no
uncovered point (AUDIT_2 §coverage). The C3 cell was added later and is checked in AUDIT_3. The final e ≥ 0.95 atlas
(F0, KQ1, KQ2, KQK1, KQK2, KQB, C3 + KA) was checked with 500,000 samples stressing every seam, again with no
uncovered point (AUDIT_6 §b). All seams use identical floating-point literals on both sides.

## 4. The degenerate corners

These are the places where the polynomial vanishes identically, or has clusters of roots beyond J = 6.

**D1a: near-coincident circle and ellipse (e, l → 0, s → 1).**
- In CT coordinates, f/e² = 2 sin u [θ sin 2(u+ω) − 4τ sin u − e sin 2u] + O(2), homogeneous of degree 1. The
  weighted blow-up gives three cone charts.
- The leading term vanishes identically only in direction D1b (l = e, ω = 0, s = 1).
- A second homogeneous blow-up there has the leading term 2 sin u·[2a sin 2u + 2ph cos 2u − 4τ sin u + e·r₁(u)],
  where r₁ contains sin 3u. It never vanishes identically.
- Exact identities and divisibility were checked in AUDIT_2.
- Code: code/gen_bu.py, code/make_engines_bu.py.

**C3: e → 1, circle ≫ pericentre distance, aphelion direction → circle axis.**
- The tan-half-angle polynomial has weighted order exactly 10 under the weights (1, 1, 2, 2, 3): its terms have
  weighted degrees 10, 12, …, 28, and the degree-10 part P̂₀ is the leading term. It has a
  nested degenerate point, which gets a second blow-up.
- C3_CORE.md; code/c3/.

**The near-parabolic region (engine KQ).** As e → 1, D = 1 + e cos f is small near aphelion, so the Fourier-basis
f-chart has values of order ε⁴ against coefficients of order 1. KQ counts the zeros on |f − π| ≤ φ₀ with the polynomial
Q(τ) = P(f)(1 + τ²)⁶, τ = tan((f − π)/2), written in ε = 1 − e. Its monomial basis centred at aphelion has no such
cancellation: the aphelion pair at τ = O(ε²) is an ordinary, well-conditioned double cluster. The f-chart counts the
rest of the circle. The split angle is φ₀ = 0.9. For nearly polar orbits (l ≥ 0.9) whose axis crossing,
|f − π| = π/2 − ω, lies near 0.9, it is φ₀ = 1.25 instead. φ₀ is a fixed function of the box, and any φ₀ gives a
valid certificate, because the two parts overlap. The engine also covers the l ≈ 1 band, where the orbit crosses
the circle's axis inside the aphelion arc, because Q is a polynomial and σ = 0 does no harm.

On each box KQ first tries the winding certificate (§2) of the f-chart polynomial on the whole circle, and uses the
count above only when that fails. The winding certificate never resolves the real roots, so it is effective where up
to 10 real roots cluster but the remaining complex pair stays off the unit circle. That is typical for e ≈ 0.95. As
e → 1 the Fourier cancellation makes it fail, and the count above takes over. Winding leaves have bound 12 − 2W ≤ 10
directly.

**ω = 0 pitchfork surfaces and the t → 1 six-fold curves.** These are clusters of up to 6 roots. They are handled
by Rolle runs with J = 6, the winding certificate and the deep fallback, with no special chart. Near t = 1 (R ≫ q)
the counts depend mainly on the width in t, so the large-circle KQ continuation runs split t first.

**The face e = 1, μ = 0 (the parabolic limit with a point-like pericentre).** There P = −D²B² has up to 12 real
roots (a 6-fold root at aphelion and three double roots of B²), so no box of the KQ runs touching it can be
certified. A weighted blow-up with weights (τ, ε, μ) = (1, 2, 2) resolves it: the leading part
16c²y²[4m² − c²(E + 2y²)²] has at most 4 real zeros near aphelion (c = σ at aphelion > 0.156 outside the C3 cell),
the rest of the circle has at most 6, and the near, gap and outer parts are certified by the charts FE, FM, FIP, FIM
and the run KQF (C3_CORE.md §6).

## 5. Results

See RESULTS.md: per-run counts of initial boxes, leaves, bound histogram, failed leaves and open boxes, and the
status of each family. A family is complete when every failed leaf and open box is covered: by a complete
continuation run, or (C1p) by the D2 corner charts, or (D2Tp) by the D2 coverage rules of §3, or (FE, FM, FIP,
FIM) because it lies outside the needed range s ≤ SF of the face charts (C3_CORE.md §6). The
re-certification and spot-check outputs are in verify/.

## 6. Verification and audit

- **audit/AUDIT_1.md – AUDIT_15.md.** Fifteen rounds of independent adversarial audit of every definition, generated
  formula, engine, coverage claim and run, with the auditor's own scripts and outputs in audit/. Every issue raised
  was fixed. Round 8 found the obstruction at the face e = 1, μ = 0 (R8-1), resolved by the face charts (C3_CORE.md
  §6, audited in round 9); round 10 audited the perihelion count and the split rules (C3_CORE.md §7); round 11 the
  corner repair of the face construction; round 12 the snapshot of 2026-10-23 15:00; round 13 the tau-winding
  certificate (C3_CORE.md §8); round 14 the disc-winding outer count (C3_CORE.md §9) and the final contour choice;
  AUDIT_15 is the final pass over the complete run set.
- **code/bb_replay.py.** Exact replay of every split tree; optional re-certification of leaves.
- **code/verify_mp.py.** 60-digit spot checks of sampled leaves against the original polynomial. It has no map for
  the D2 cone charts and reports 0 points there. Their spot checks are the auditor's 100-digit checks
  (audit/round7/g3_spot_*.out: 0 violations, max N = 8). code/disc_mp_check.py does the same for the disc-winding
  leaves, counting the outer region only.
- **FLOAT_RIGOUR.md.** Floating-point error analysis. The only assumption beyond IEEE-754 arithmetic is that libm
  sin and cos are within 2 ulps for |x| ≤ 64 (all arguments satisfy |x| ≤ 62.2; checked empirically in
  code/libm_check.py).
