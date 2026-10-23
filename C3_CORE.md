# The C3 core and the l ≈ 1 band (e → 1 corner)

Code: `code/c3/` (C3 core), `code/c3_away.py` (away part, run KA) and `code/kq_run.py` (the band, run KQB).

## 1. Setting

R-normalised f-chart: the circle has radius 1, the ellipse has pericentre distance μ = q/R and eccentricity e. True
anomaly f; ζ = P_z cos f + Q_z sin f, σ² = 1 − ζ², B = (σ²)' + 2e(sin f − Q_z ζ), D = 1 + e cos f, k = 2e(1+e).

  P(f) = μ²k² sin²f σ² − D²B² = F₋F₊,  F∓ = μk sin f σ ∓ DB.

F∓ = (D³σ/(μ(1+e)))·g∓'(f) in the R-normalisation (direct computation, see PROOF.md §1; corrected per audit R3-4;
only positivity matters). So the zeros of F∓ with σ > 0 are
exactly the stationary points on the nearest/farthest branch, and P is a positive multiple of the a-chart polynomial
Pol composed with the anomaly change u ↦ f. Roots of P therefore bound N with multiplicity.

**Anomaly variable.** τ = tan((f − π)/2) maps f ∈ (0, 2π) bijectively onto ℝ (f = 0 ↔ τ = ∞). Then
Q(τ) := P(f)(1 + τ²)⁶ is a polynomial of degree ≤ 12 in τ, and its real roots are the roots of P in (0, 2π) with the
same multiplicities.

## 2. The C3 region and the weighted blow-up

The KQ runs for e ≥ 0.99 (KQK1, KQK2, KQB; PROOF.md §3) leave out the cell e ∈ [0.99, 1], l ∈ [0.95, 1],
ω ∈ [9π/20, π/2], μ ∈ [0, 0.0505]. Write
Q_z = l cos ω, κ = 1 − P_z, ε = 1 − e. On the cell:
- Q_z ≤ cos(9π/20) = 0.15643;
- κ ≤ 1 − 0.95 sin(9π/20) = 0.061696;
- ε ≤ 0.01, μ ≤ 0.0505.

The weighted radius is δ = max(Q_z, √κ, √ε, μ^{1/3}) ≤ 0.36963 < NU1 = 0.37.

**Quasi-homogeneity (checked exactly, `c3/c3_sym.py`).** Substitute τ = νy, Q_z = νa, P_z = 1 − ν²K, e = 1 − ν²E,
μ = ν³m. Then

  Q = ν¹⁰ [ P̂₀(2y) + ν² P̂₁ + … + ν¹⁸ P̂₉ ]   (only even powers of ν),

with
- P̂₀(x) = 16m²x²(λ + (x − a)²) − (E + x²/2)² B̂(x)²;
- B̂ = 2E(x − a) − 2λx − 2a²x + 3ax² − x³;
- λ = 2K − a².

## 3. Three regions in τ, for each parameter point p in the cell

| region | τ-range | certificate | target |
|---|---|---|---|
| near | \|τ\| ≤ Y0·δ(p) | level-1 cone charts Ca, CK, CE, Cm (dominant variable Q_z, κ, ε or μ); root variable y = τ/ν ∈ [−3, 3]; G = Q(νy)/ν¹⁰ | ≤ 9 (Ca), ≤ 8 (CK, CE, Cm), ≤ 5 + 4 in the level-2 zone |
| intermediate | Y0·δ(p) ≤ \|τ\| ≤ Y0·NU1 = 1.11 | cone charts PhiP / PhiM (τ = ±3ν dominant): G = Q(±3ν)/ν¹⁰ has no zero | 0 |
| away | \|f − π\| ≥ 2 atan(1.11) | f-chart arcs over the physical cell (bb_run engine KA, c3_away.away_count_t) | ≤ 2 |

Total ≤ 9 + 0 + 2 = 11 real roots of P on the whole circle, counted with multiplicity. P is a real trigonometric
polynomial in f, so this number is even (parity lemma, PROOF.md §2), hence ≤ 10.

**Cones.** Each chart's ratios lie in [0, R1] (R1 = 1.001 ≥ 1), with radial range ν ∈ [0, 0.37]. Every p with δ(p) ≤ 0.37
lies in the chart of its dominant weighted component, with all ratios ≤ 1. For the intermediate charts the condition
δ(p) ≤ |τ|/Y0 = ν means exactly ratios ≤ 1.

**Physical pruning (`c3_run.phys_bounds`).**
- A box is skipped ('u') only if every point has l > 1, i.e. Q_z² > 2κ − κ² everywhere on the box.
- A box is skipped ('o') only if every point is outside the cell, using lower bounds of Q_z, κ, ε, μ.

Empirical motivation: in 20,000 samples of the cell, the near group had ≤ 8 roots (8 occurred 16 times), the
intermediate region had 0 roots, and the far group (pericentre, |τ| ≥ 12.7) had exactly 2 (`c3/c3_roots.py`).

## 4. Nested degeneracy (level 2)

At E = m = 0 in Ca or CK, P̂₀ = −(x⁴/4)·x²·(x² − 3ax + 4K)². That is a 6-fold root at y = 0, plus two double roots
which are real only in the wedge 16K/9 ≤ a² ≤ 2K. Here 6 + 4 = 10 > 8, so these points need a second blow-up.

Near (y, E, m) = 0 the function G vanishes to weighted order 6 for weights (1, 2, 2). This is checked exactly: the η⁶
division in gen_c3. The leading part is y²·[4m² − 2K(E + 2y²)²] up to a positive factor, so there are at most 4 roots
near y = 0. This is the D3 corner structure (F₋ ≤ 3, F₊ ≤ 1) at a smaller scale.

**Level-2 zone** Z = {|y| ≤ Y2 = 0.25, E ≤ E2 = 0.01, m ≤ M2 = 0.01}.
- A level-1 box of Ca or CK lying entirely in E ≤ E2, m ≤ M2 is certified by the outer count only, on
  |y| ∈ [Y2, Y0], with target ZONE_OUTER = 5.
- The inner part |y| ≤ Y2 is certified by level-2 cone charts in (y/Y2, √(E/E2), √(m/M2)):
  - L2·E, L2·m (dominant E or m; root variable Y = y/(Y2η) ∈ [−R2, R2]): target 4;
  - L2·Yp, L2·Ym (dominant y): no zero (root variable mp ∈ [0, R2²]).
- In either case the near count is ≤ 5 + 4 = 9.

## 5. The l ≈ 1 band (engine KQ, run KQB)

For l = 1 (an orbit plane containing the circle's axis), σ vanishes where the orbit crosses the axis, at
|f − π| = π/2 − ω. For ω ≥ 0.67 that point lies inside the aphelion arc |f − π| ≤ 0.9. This does no harm to engine KQ
(code/kq_run.py), because KQ counts the real roots of P itself, not of the branch factors F∓ = μk sin f σ ∓ DB:
- on |f − π| ≤ φ₀ (0.9, or 1.25 when the axis crossing lies near 0.9) it counts the real roots of
  Q(τ) = P(f)(1 + τ²)⁶, a polynomial in τ whose coefficients are
  polynomials in (ε, P_z, Q_z, μ);
- on the rest of the circle it counts the roots of the f-chart polynomial.

Neither involves σ = √σ², so σ = 0 is not a singularity of either count. (On each box KQ first tries the winding
certificate of the f-chart polynomial, which does not involve σ either; PROOF.md §4.) The band e ∈ [0.99, 1], μ ∈ [0, 0.0505],
l ∈ [0.95, 1], ω ∈ [3π/20, 9π/20] (as t-range: t ≥ 0.997456) is certified by run KQB, with total ≤ 11 and hence ≤ 10
by parity.

An earlier band engine (KB, which counted F₋ and F₊ separately near aphelion) is superseded and is not part of the
final atlas.

## 6. The face e = 1, μ = 0 (charts FE, FM, FIP, FIM and run KQF)

**The obstruction (AUDIT_8, R8-1).** At e = 1 and μ = 0, P = −D²B². For l ≥ l*(ω) it has 12 real roots counted with
multiplicity:
- a 6-fold root at aphelion τ = 0;
- three double roots of B² elsewhere.

l*(ω) runs from 0.707 at ω = 0 to 0.9985 at ω = 1.41 (`audit/round8/h4_corner_limit.py`). So a KQ box that touches
e = 1 and t = 1 can never have total ≤ 11, and the runs KQK1, KQK2 and KQB cannot close near this face on their own.
Physical points next to the face have N_mult = 10 (`audit/round8/h3_kqk2_corner.py`).

**Weighted blow-up.** Give weights (τ, ε, μ) = (1, 2, 2), with ε = 1 − e, and let s ≥ 0 be the radius. Then:
- G = Q(sy)/s⁶ is a polynomial: every monomial τ^k ε^a μ^b of Q has k + 2a + 2b ≥ 6. This is checked exactly
  (`c3/gen_c3.py`, face_G asserts that the lowest power of s is 6).
- Its leading part is, exactly (sympy):

  G(0) = 16c²y²[4m² − c²(E + 2y²)²],  with c² = 1 − P_z² (c = σ at aphelion), E = ε/s², m = μ/s².

**Real zeros of the leading part.** For c > 0:
- y = 0 (double);
- 4m² = c²(E + 2y²)² ⟺ 2m = c(E + 2y²), i.e. y² = (2m/c − E)/2, which gives at most one pair.

So there are at most 4 zeros near aphelion. The three double roots of B² stay away from aphelion, which leaves 6 for
the rest of the circle.

Outside the C3 cell (l < 0.95 or ω < 9π/20), P_z < max(0.95, sin(9π/20)) = 0.98769, so c > 0.1564. Then the
leading-order pair has |y| ≤ (R1/c)^{1/2} < 2.53 < Y_F = 3.

**Face neighbourhood and charts** (SF = 0.02, Y_F = 3, T_X = Y_F·SF = 0.06, R1 = 1.001). Let
N = {ε ≤ SF², μ ≤ SF²} and s(p) = max(√ε, √μ) ≤ SF.

| region of τ, for p ∈ N | chart | box, root variable | target |
|---|---|---|---|
| near: \|τ\| ≤ Y_F·s(p) | FE: ε = s², μ = s²m (ε dominant); FM: μ = s², ε = s²E (μ dominant) | (s, ratio, P_z, Q_z) ∈ [0, SF·R1] × [0, R1] × [0, 1]²; root y = τ/s ∈ [−3, 3] | ≤ 5 |
| gap: Y_F·s(p) ≤ ±τ ≤ T_X | FIP (τ > 0), FIM (τ < 0): τ = ±Y_F·σ, ε = σ²E, μ = σ²m | (σ, E, P_z, Q_z); root variable m ∈ [0, R1] | 0 (no zero) |
| outer: \|τ\| ≥ T_X, and the away arc | KQF (`kq_run.py --face`): the near polynomial on T_X ≤ \|τ\| ≤ tan(φ₀/2), both signs, plus f-chart arcs on \|f − π\| ≥ φ₀ | (e, l, ω, t) over e ∈ [0.9996, 1], t ∈ [0.99999983, 1] ⊇ N | ≤ 6 |

The total is ≤ 5 + 0 + 6 = 11 roots of P on the circle, counted with multiplicity, so ≤ 10 by parity (PROOF.md §2).

**Coverage.** Let p ∈ N.
- Near: p lies in FE if ε ≥ μ, with ratio m = μ/ε ≤ 1, and in FM otherwise. Its near τ-range is y = τ/s(p) ∈ [−3, 3].
- Gap: a point τ ∈ [Y_F·s(p), T_X] has σ = τ/Y_F ∈ [s(p), SF], E = ε/σ² ≤ 1 and m = μ/σ² ≤ 1. So it lies in FIP,
  and similarly in FIM for τ < 0.
- Only s ≤ SF and σ ≤ SF are needed. The chart domains run to SF·R1, and `check_all.py` (and `next_level.py`) excuse an
  open or failed box of FE, FM, FIP or FIM whose lower s-bound exceeds SF.
  - FIP has failed leaves (bound 1), all with σ ≥ 0.0200070 > SF.
  - They sit at P_z ≈ 0.986, Q_z ≈ 0.165, i.e. l ≈ 1 and ω ≈ 80.4°, just outside the C3 cell. There a root of Q
    approaches τ = Y_F·σ for σ slightly above SF.

**The corner K (repair).** Near l = 1 and ω just below 9π/20 (c close to its minimum 0.156 outside the C3 cell) the
construction above fails as stated: the gap charts have genuine zeros.
- Mechanism (AUDIT_11, R11-4): the outer member of the near-aphelion pair y² = (2m/c − E)/2 sits at |y| up to
  (R1/c)^{1/2} ≈ 2.5 at leading order. The next-order terms push it further out, so for s close to SF it reaches the
  gap boundary y = 3. This happens from about σ ≈ (c − R1/9)/(3 sin ω) ≈ 0.0153 at the corner of the cell.
- B itself has no zero there: at l = 1, ω = 9π/20 its zeros next to aphelion are at τ = 0.0787 and 0.1584.
- Example: at σ = 0.01877, E = 0.25, P_z = 0.9863, Q_z = 0.1646, G changes sign at m ∈ (0.995, 1] (50 digits), so Q
  has a root at y ≈ 3.004, inside the gap.
- FIP's failed leaves lie at σ ≥ 0.0175. AUDIT_11 finds that the smallest σ with a real gap zero in K outside the C3
  cell is 0.0152, i.e. about 1.5·SF2. For s ≤ SF2 the near zeros stay at |y| ≤ 2.79, the gap is empty, the outer count
  is 6, and the nearest zero beyond the near region is at τ ≈ 0.079 > T_X2.

The repair uses K = {P_z ≥ 0.875, 0.125 ≤ Q_z ≤ 0.25}, which contains every face-chart box with a failed or open leaf.
Its image in (l, ω) lies in K_lw = {l ≥ 0.88, ω ≥ 1.29}, because l ≥ 0.884 and atan2(P_z, Q_z) ≥ 1.2925 on K. Let
SF2 = 0.010005 and T_X2 = 3·SF2 = 0.030015. For p ∈ N with (l, ω) ∈ K_lw:
- if s(p) ≤ SF2:
  - near: FE/FM with s ≤ SF2;
  - gap: FIP/FIM with σ ≤ SF2 (their records with root s ≤ 0.01001 are complete everywhere);
  - outer: the run KQF2 (`kq_run.py --face --faces=0.010005`), the outer count from |τ| ≥ T_X2 over
    ε, μ ≤ SF2², l ≥ 0.88, ω ≥ 1.29, target 6;
  - total ≤ 11;
- if SF2 < s(p) ≤ SF: the run KQS (`kq_run.py --faces=0.010005 --c3prune --splitrule=shell/30`), a whole-circle KQ
  count ≤ 11 over the face neighbourhood restricted to l ≥ 0.88, ω ≥ 1.29. Its 'p' leaves lie in {ε, μ ≤ SF2²} and
  its 'o' leaves in the C3 cell.

For (l, ω) ∉ K_lw, (P_z, Q_z) ∉ K, and the construction above applies unchanged. `check_all.py` therefore excuses:
- a face-chart box inside K with s_lo > SF2;
- a KQF box inside K_lw.

It requires FE, FM, FIP, FIM, KQF, KQF2 and KQS to be complete for any family with 'p' leaves.

**Pruning.**
- 'u': unphysical, P_z² + Q_z² > 1 on the whole box.
- 'o': the whole box lies in the C3 cell (l ≥ 0.95 and ω ≥ 9π/20, from lower bounds of l and ω). For e ≥ 0.99 and
  μ ≤ 0.0505 that cell is covered by the C3 core (§2–§4), and N ⊂ {ε ≤ 0.01, μ ≤ 0.0505}.
- KQF prunes 'o' by the same predicate.
- In the KQ runs (engine KQ), a box inside N is pruned 'p'. Here eps_hi ≤ SF²(1 − 1e-9) and
  mu_hi = √((1 − t_lo)/t_lo) ≤ SF²(1 − 1e-9).
- `check_all.py` re-evaluates each of these predicates on every pruned leaf, with the face radius pinned per family
  (SF2 for KQS and KQF2, SF otherwise; also the outer start T_X = 3·radius of KQF and KQF2). It accepts a family that
  uses 'p' only when FE, FM, FIP, FIM, KQF, KQF2 and KQS are all complete.

**Split rule near the face** (`code/splitrule.py`, `kq_run.py --facesplit`, header field split_rule). It is used for the
continuation levels KQK1 r4+, KQK2 r4+ and KQB r5+, and it replays exactly in `check_all.py` and `bb_replay.py`:
- A box that meets {ε < S2, μ < S2}, with S2 = SF²(1 − 5e-9), is cut at e = 1 − S2 (or else at t = 1/(1 + S2²)). The
  child on the face side is then pruned 'p'.
- Otherwise the ε- and μ-widths are measured relative to s² = max(ε_hi, μ_hi, SF²), the scale of the blow-up, and the
  l- and ω-widths are measured as before. The widest coordinate is bisected.

The variant 'face/K' (`kq_run.py --splitrule=face/K`, used with K = 30 for KQB from level r5 on) divides the ε- and
μ-widths by K, so l and ω are split about K times finer; on the band ridge the counts need l, ω about 30 times finer
than ε, μ. The split rule changes only which boxes are certified, never what a leaf certifies.

## 7. The perihelion count (kq_run.py --peri) and the split angle

**The ridge l = l*(ω).** Near the face, on the ridge where two zeros of B merge, the f-chart arcs need very small
boxes. At ω = 0 the merging zeros meet the zero of B at perihelion: at e = 1, ω = 0,
B = 2 sin f·[(1 − l²) − l² cos f], which has a triple zero at f = 0 when l = 1/√2. Near this ridge the circle carries 10
real zeros (four near aphelion, six near perihelion; e.g. e = 0.9925, l = 0.714, ω = 0.001, μ = 0.0103). There the
Fourier-basis enclosure of the sixth derivative needed boxes about 10⁻³ wide in every relative coordinate.

**Perihelion polynomial.** With tp = tan(f/2) = −1/τ:

  G(tp) = tp¹²·Q(−1/tp) = P(f)(1 + tp²)⁶,

i.e. the coefficients of the aphelion polynomial (chart KBQe) reversed, with signs (−1)^k. This is chart KBPe
(`c3/gen_c3.py`; self-test against the direct formula for P, relative error 1.8e-14). Its real zeros in tp are the
zeros of P with the same multiplicities.

With --peri, engine KQ covers the circle by four overlapping parts and adds their bounds:
- the aphelion polynomial on |f − π| ≤ φ₀;
- the perihelion polynomial on |f| ≤ φ_P = 1.2 (|tp| ≤ tan 0.6 + 10⁻⁹);
- the f-chart arcs [π + φ₀, 2π − φ_P] and [φ_P, π − φ₀], each padded at both ends.

When φ₀ + φ_P ≥ π the two polynomials already cover the circle and the arcs are dropped. A box gets the smaller of
this count and the plain one. On the ridge points above, the perihelion count certifies the total 10 on boxes about
100 times wider per coordinate.

**Split angle.** Engine KQ always tries the box's split angle φ₀ = phi_split(bx) and then the default 0.9. Any split
angle gives a valid count. On the band ridge (l ≈ 0.95, ω ≈ 0.76), φ₀ = 1.25 gave aphelion counts of 13 where
φ₀ = 0.9 certifies.

Re-certification (`bb_replay.py`, DEEP_ALWAYS) tries every combination, so records made with either option re-certify.

## 8. The tau-winding certificate (code/tau_wind.py, leaf code '7')

Near the ridge l ≈ l*(ω), N = 10 is attained, so the whole-circle Rolle count has no slack. There the one non-real
pair of zeros of the aphelion polynomial sits near aphelion, well off the real axis: τ ≈ ±0.12i on the KQB and KQK1
ridges, ±0.16i on the KQ2 band. In the face blow-up this is the imaginary branch y² = −(2m/c + E)/2 of the near cluster.
The tau-winding certificate certifies that this pair is non-real.

**Lemma (counting).** Q(τ) = P(f)(1 + τ²)⁶ is a real polynomial of degree ≤ 12. P has, counted with multiplicity:
- the real zeros of Q, at f ≠ 0;
- a zero of order 12 − deg Q at f = 0 (τ = ∞).

The reverse G(t_p) = t_p¹²Q(−1/t_p), with t_p = tan(f/2), has lowest coefficients equal to the top ones of Q. Hence
N_P = 12 − #{non-real zeros of Q}. If Q has one non-real zero, its conjugate is another, so N_P ≤ 10. A non-real zero
t_p of G gives the non-real zero −1/t_p of Q. The same holds for the chart polynomials in their root variable y:
- Q(νy)/ν¹⁰ on the C3 chart Ca;
- Q(sy)/s⁶ on the face charts FE and FM.

Here ν > 0 and s > 0 at every physical point, so a non-real zero y gives the non-real zero νy (or sy) of Q.

**Certificate.** Take F = Q (chart KBQe), F = G (chart KBPe), or one of these chart polynomials, and a circle
|z − c| = r with Im c > r > 0. Suppose
F(·; x) ≠ 0 on the circle for every x in the box. The box is connected, so the number of zeros inside the circle is
the same for all x; a vanishing leading coefficient only moves zeros to infinity, outside the circle. It therefore
equals the winding number of F(·; x_c) along the circle. If that number is ≥ 1, every point of the box has N_P ≤ 10.

**Enclosure on each arc** (contiguous arcs of the angle, centres z_j, |z − z_j| ≤ D_j). For all x in the box,
F(z; x) lies within R_j of the computed centre value b₀ = F(z_j; x_c), where:
- the deviation along the arc is Σ_{m≥1} |b_m| D_j^m, with b_m the exact Taylor coefficients at z_j of the polynomial
  with the centre coefficients;
- the parameter variation is the engine's bound (c3/make_c3_2.py: the smaller of the mean-value and the second-order
  form, at a complex point). Every |p(z)| is replaced by its disc maximum Σ_m |b_m(p)| D_j^m, and every |z|^k by
  (|z_j| + D_j)^k;
- all floating-point errors are covered by 10⁻¹² of the corresponding sums of absolute terms.

If |b₀| > R_j on every arc, F ≠ 0 on the circle. The argument of F then stays within π/2 of arg b₀ on each arc,
consecutive arcs share an end point, and the winding number is the sum of the principal arguments of b₀ between
consecutive arcs.

**Contour choice.** The contour comes from the zeros of F(·; x_c). Its centre is one of:
- a zero with Im > 0;
- the midpoint of two such zeros;
- the mean of all of them.

Its radius is halfway between two consecutive distances from the centre to the zeros, and below Im c. Candidates are
tried in order of margin, i.e. the distance of the circle from the zeros and from the real axis. A cluster of nearby
zeros is enclosed together, since only W ≥ 1 is needed. The contour uses chart KBPe when the zeros lie far from
aphelion. The choice affects only success, never validity.

**Records and checks.**
- Leaf code '7' means bound 10 for the whole circle. Its contour [which, c_re, c_im, r] is stored in the record field
  'tw', in tree order.
- `check_all.py` counts these leaves and checks that each contour lies in the open upper half-plane.
- `bb_replay.py --recert` re-certifies each stored contour (tau_wind.check_tw). `recert_tw_c3.py` re-certifies every
  '7' leaf of the c3 charts Ca, FE and FM (tau_wind.check_tw_chart). There which = 0, which `check_all.py` enforces.
- In the outer runs KQF and KQF2, a '7' leaf certifies N_P ≤ 10 on its box directly, so those points need no face
  decomposition.
- The final audit (AUDIT_15) re-certified every '7' leaf of the complete run set (175,932) and every '8' leaf
  (4,451).

**Validation.** On 107 certified boxes from the KQB, KQK1, KQ2, KQ0T and KQS levels, 642 random points were checked:
- every point has N ≤ 10 at 60 digits;
- at every point, a 100-digit zero of F lies inside the certified circle.

## 9. The disc-winding outer count (tau_wind.certify_disc, leaf code '8', runs KQF and KQF2 with --disc)

At the face, Q has a 6-fold zero at aphelion τ = 0 and three double zeros of B² away from it (§6). Near the face the
six aphelion zeros stay in a small disc about τ = 0, and the Rolle count of the outer region must resolve the three
near-double zeros. Counting the zeros in the disc instead gives the outer bound directly.

**Lemma (outer count).** Let T > 0. Suppose Q(·; x) has no zero on the circle |τ| = T, and W zeros, counted with
multiplicity, in the open disc |τ| < T. Then P has at most 12 − W real zeros, counted with multiplicity, on the outer
region O_T = {f : |τ(f)| ≥ T} ∪ {f = 0}.

*Proof.* By the counting lemma of §8, the real zeros of P on O_T are:
- the real zeros of Q with |τ| ≥ T, which all have |τ| > T since none lies on the circle;
- a zero of order 12 − deg Q at f = 0.

Q has deg Q zeros in ℂ, counted with multiplicity. W of them lie in the disc, so at most deg Q − W lie outside the
closed disc. The total is at most (deg Q − W) + (12 − deg Q) = 12 − W. ∎

**Certificate.** F = Q (chart KBQe), on the circle with centre 0 and radius T_X, the outer start of the run (T_X = 0.06
for KQF, T_X2 = 0.030015 for KQF2). The arc enclosure is the one of §8. If |b₀| > R_j on every arc, then for every x
in the box Q(·; x) ≠ 0 on the circle. The box is connected, so the number W of zeros in the disc is the same for every
x (a vanishing leading coefficient only moves zeros to infinity, outside the disc), and it equals the winding number of
Q(·; x_c). If W ≥ 6, then every point of the box has at most 6 real zeros of P on O_{T_X}.

This is exactly the outer target of §6. count_outer certifies the same region: T_X ≤ |τ| ≤ tan(φ₀/2) on both signs,
plus the away arc |f − π| ≥ φ₀. So does count_peri, the variant with --peri (§7): its parts cover the same closed
region {|τ| ≥ T_X} ∪ {f = 0}. With the near count ≤ 5 (FE, FM) and the empty gap (FIP, FIM), the total stays ≤ 11, so
≤ 10 by parity.

**Records and checks.**
- Leaf code '8' means outer bound 6. It is valid only in KQF-engine files, and stores no data: its contour is fixed by
  the family's T_X, which `check_all.py` pins to 3 × the face radius.
- `check_all.py` counts '8' leaves as bound 6 and rejects them in any other engine.
- `bb_replay.py --recert` re-certifies them (tau_wind.check_disc with the header T_X).
- `disc_mp_check.py` is an independent spot check. It maps sample points of '8' leaves to physical (e, l, ω, s) and
  builds the original a-chart polynomial at 60 digits. It then counts its zeros on the unit circle whose true anomaly
  has |tan(f/2)| ≤ 1/T_X, i.e. |τ| ≥ T_X, and checks that the count is ≤ 6.

**Validation.** On test runs over random KQF r2 and KQF2 r3 boxes, every box was certified W = 6 at full box size: 57
'8' leaves, 228 sample points, maximum 60-digit outer count 6 (attained), no violation. In the runs, KQF2 r4 (543 boxes)
closed in one second, all '8' leaves.
