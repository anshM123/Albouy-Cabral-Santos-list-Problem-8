# The e → 1 corner (near-parabolic orbits): plan and lemmas

> **Status.** This file records the branch-wise analysis of the e → 1 corner (lemmas N1, N2, A), which was used by an
> earlier corner engine. The final atlas does not use it: the whole region e ≥ 0.9499 is certified by the f-chart run
> F0 (t ≤ 0.5) and by engine KQ (t ≥ 0.5), which counts the zeros of P = F₋F₊ directly (PROOF.md §3–4), together with
> the C3 core (C3_CORE.md). The analysis below is kept as background.

## Chart

True anomaly f. The circle has radius R = 1. The pericentre distance is μ = q/R, and ε = 1 − e.
- r = μ(1+e)/D, D = 1 + e cos f.
- ζ = P_z cos f + Q_z sin f, σ = √(1 − ζ²), ρ = rσ.
- B = (σ²)' + 2e(sin f − Q_z ζ).

The stationary points of g∓ = r² + 1 ∓ 2ρ (∓ = nearest/farthest point of the circle) are the zeros of

    F∓(f) = 2μ e(1+e) sin f · σ  ∓  D · B          (positive multiple of g∓', valid where σ > 0, e < 1).

Squaring gives the polynomial P_f = F₋F₊ used in the f-chart certificate.

Facts at f = π (aphelion):
- B(π) = −2P_zQ_z ε;
- B'(π) = 2P_z² − 2e + 2(e − 1)Q_z², which is −2c² at e = 1, with c² = 1 − P_z²;
- D(π) = ε.

## Target in the corner

The corner is K = {ε ≤ ε0, μ ≤ μ1, c ≥ c0}. In K:
- #zeros(F₋) ≤ 6, so g₋ has at most 3 local minima (local MOIDs);
- #zeros(F₊) ≤ 4, so g₊ has at most 2 local maxima.

Hence the total is at most 10.

### Lemma N1 (aphelion arc, nearest branch)
- **Claim:** F₋''' > 0 on |f − π| ≤ φ0 throughout K, so F₋ has at most 3 zeros there.
- **Why:** near π, D·B ≈ (ε + φ²/2)(−2P_zQ_zε − 2c²φ), whose cubic term −c²φ³ gives −(DB)''' ≈ 6c² > 0. The
  μ-term is O(μ).
- **Status:** not degenerate anywhere in K, so it can be certified directly by interval arithmetic.

### Lemma N2 (aphelion arc, farthest branch)
- **Claim:** F₊ is strictly decreasing on |f − π| ≤ φ0 throughout K, so F₊ has at most 1 zero there.
- **Why:** F₊' = 2μe(1+e)(cos f σ + sin f σ') + (DB)'.
  - The first bracket is about −σ(π) < 0. It is non-degenerate and certifiable.
  - (DB)' ≈ −2c²ε − 3c²φ² − 2eP_zQ_z εφ ≤ 0. Equality holds only at φ = ε = 0, and the cross term is absorbed by
    AM–GM when ε ≤ 24c⁴.
- **Status:** needs an explicit Taylor-remainder (or rescaled interval) proof because of the degenerate point
  (φ, ε) = (0, 0).

### Lemma A (away from aphelion)
- **Claim:** on |f − π| ≥ φ0, D ≥ D_min > 0, and F∓ = ∓DB + O(μ). Each branch has at most 3 zeros there.
- **Why:** B has trigonometric degree 2. One of its four zeros lies within φ0 of π (at φ ≈ −P_zQ_zε/c²).
- **Status:** certify by branch and bound over (e, l, ω, μ) with μ ∈ [0, μ1], μ = 0 included, using Rolle arcs.

### Numerical status (test_lemmas_e1.py, 2026-10-22)
- **Lemma A** (at most 3 zeros per branch away from aphelion) and **total ≤ 10**: 0 violations in 2,000 corner samples
  (ε, μ ∈ [1e-8, 1e-2], c² ≥ 0.05).
- **Lemma N2** (F₊ decreasing): 0 violations with φ0 ≤ 0.15. At φ0 = 0.3 one arc end fails when l ≈ 1.
- **Lemma N1** (F₋''' > 0 on the whole arc): fails at the outer parts of the arc. So the aphelion count for F₋ is left
  to the Rolle-run certificate. It is not degenerate at the corner point, because −(DB)''' ≈ 6c².

### Lemma N2: proof that avoids degenerate monotonicity
1. B' < 0 on |φ| ≤ φ0. This is non-degenerate (B'(π) ≈ −2c²) and certifiable by intervals.
2. B has a single zero φ_B there, with |φ_B| ≤ 2|P_zQ_z| ε / min|B'|.
3. Signs:
   - for φ > max(0, φ_B), both terms of F₊ = μk sin f σ + DB are ≤ 0 and not both zero;
   - for φ < min(0, φ_B), both terms are ≥ 0.
   So every zero lies in the interval between 0 and φ_B, of length O(ε).
4. On that interval:
   - μk(cos f σ + sin f σ') < 0;
   - |D'B| = O(ε²);
   - DB' ≤ −c² ε, because D ≥ ε.
   So F₊' < 0 once ε is below an explicit threshold. Hence F₊ has at most one zero on the aphelion arc.

### Mixed certificate for the whole corner {e ≥ 1 − ε0, μ = q/R ≤ μ1, c ≥ c0}

    N  ≤  [real zeros of P = F₋F₊ on the arcs |φ| ≥ φ0]  +  [zeros of F₋ on |φ| < φ0]  +  [zeros of F₊ on |φ| < φ0]
       ≤        6 (branch and bound; P'' < 0 at each zero of B)  +  3 (Rolle j = 3)  +  1 (Lemma N2)  =  10.

The squared polynomial is used away from aphelion, so σ = 0 there (an inclination of 90°, where the orbit crosses the
circle's axis) does no harm. σ > 0 is needed only on the aphelion arc, which is exactly the condition c ≥ c0.

### C3 ∩ corner (aphelion direction close to the circle's axis, c → 0, together with e → 1 and μ → 0)
At P_z = 1, Q_z = 0 we get B = 2 sin f (cos f + e), which equals 2 sin f · D only at e = 1 (corrected, audit m2).
At e = 1 and μ = 0 this gives P = −4 sin²f · D⁴. That is a 10-fold zero at f = π, plus the double
zero at 0: all 12 roots real in the limit. This is the most degenerate point of the whole problem and needs a dedicated
blow-up; this is now resolved by the C3 core charts (C3_CORE.md). Probes near it found N ≤ 10, with one case of
exactly 10 and gap 1.07.

### Remaining pieces of the e → 1 region
- **μ ≥ μ1** (circle comparable to or smaller than the pericentre distance): use the f-chart certificate. The e = 1
  face is the confocal parabola–circle problem: at most 8 numerically, with a complex-pair margin of about √2 (corrected, audit m1).
- **c → 0** (aphelion direction close to the circle's axis, the C3 sub-corner): needs separate treatment, since σ → 0
  makes the branch functions singular. Use the squared polynomial there.
