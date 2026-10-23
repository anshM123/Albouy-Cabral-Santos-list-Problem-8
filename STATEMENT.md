# Stage 1: circle–ellipse case (Problem 8, Q2), statement and reduction

## Canonical statement (literal)

Albouy–Cabral–Santos 2012, Problem 8 (G. F. Gronchi): "Consider the distance d between two points on two distinct
confocal ellipses in R^3. If both eccentricities are > 0, is 12 the maximum number of stationary points of the function
d^2? If only one of the eccentricities is 0, can this upper bound be reduced to 10?" The cases with infinitely many
stationary points (coinciding conics; concentric coplanar circles) are excluded.

**Target Theorem 1 (Q2).** For a circle and an ellipse with a common focus (the centre of the circle), not concentric
coplanar circles, d^2 has at most 10 stationary points on the torus. The bound is attained: Gronchi 2002, and our
numerics.

## Coordinates

Scale the ellipse to a = 1. The circle C(v) = s (cos v, sin v, 0) has radius s = R/a > 0. The ellipse is
X(u) = (cos u − e) P + β sin u Q, with β = √(1 − e²) and
- P = (cos ω, sin ω cos I, sin ω sin I),
- Q = (−sin ω, cos ω cos I, cos ω sin I).

Then
- d²(u, v) = |X|² + s² − 2s (X_x cos v + X_y sin v), with |X| = r = 1 − e cos u;
- ρ² = X_x² + X_y² = r² − σ z², where σ = sin² I and z = (cos u − e) sin ω + β sin u cos ω.

## Reduction lemma (validated numerically against a 2-D Newton search)

- For fixed u, ∂_v d² = 0 exactly at the nearest and farthest points of the circle, provided ρ(u) > 0.
- So the stationary points are the critical points of g∓(u) = r² + s² ∓ 2sρ(u). Equivalently they are the solutions of
  k(u) = ±s, where k = Aρ/(ρ²)' and A = (r²)' = 2e sin u · r.
- Squaring gives Pol(u) = A²ρ² − s²((ρ²)')², a trigonometric polynomial of degree ≤ 6.
- **Counting.** The number of stationary points is at most the number of roots of Pol in [0, 2π), counted with
  multiplicity. This also holds at the special points:
  - where A = (ρ²)' = 0: two stationary points, and Pol has a root of multiplicity ≥ 2;
  - where ρ = 0: at most two stationary points, and Pol again has a root of multiplicity ≥ 2.
- Write Q(z) = z⁶ Pol with z = e^{iu}. Its roots come in pairs z, 1/z̄. Roots lost to a degree drop sit at 0 and ∞,
  so they lie off the circle.
- **Hence Theorem 1 follows from: Pol ≢ 0 and Q has at least one root off the unit circle.** Equivalently, Pol is
  never real-rooted.

## Parameter domain and symmetries

- e ∈ [0, 1), σ ∈ [0, 1], ω ∈ [0, π/2], s ∈ (0, ∞). The reductions use ω → −ω, ω → ω + π and I → π − I.
- Degenerate point: e = 0, σ = 0, which gives Pol ≡ 0 (concentric coplanar circles). It is excluded.
- Pol = Φ − s²Ψ with Φ = A²ρ² ≥ 0 and Ψ = ((ρ²)')² ≥ 0.

## Sanity tests (2026-10-22, code/explore_ce.py and code/margin_ce.py)

- Reduction check: polynomial roots = sign changes of the two branches = 2-D Newton count, on 8 random cases.
- 4,000 random (e, σ, ω), maximised over s: N = 4 (33 cases), 6 (1,900), 8 (1,872), 10 (195). **Never 12.**
- ~~Minimising the off-circle gap δ = |log|z₀|| among configurations with N = 10: δ ≥ 0.836, a robust margin.~~
  **Corrected (audit M4):** this is false. The auditor found δ = 0.5557 at a genuine N = 10 configuration with
  e = 0.94998 (60-digit check), and δ → 0 as e → 1 at small s (δ ≈ 0.28 at e ≈ 0.99). The proof does not use any
  uniform margin.
- Planar check (σ = 0): Pol = A²(r² − s²), which gives N ≤ 6, matching Gronchi–Baù–Grassi.
- Near σ = 0 with e → 1, the off-circle double root from r = 0 (|log z| = arccosh(1/e)) approaches the circle. There N is
  small, so this corner needs its own certificate.

## Proof plan (historical — superseded by PROOF.md: Rolle-run certificates over an atlas of charts)

1. Symbolic coefficients of Q(z), as polynomials in (e, β, σ, cos ω, sin ω, s²), exact (sympy).
2. Certificate per parameter box B, using Rouché: a disk D with closure inside |z| < 1 such that
   |Q_B(z) − Q_c(z)| < |Q_c(z)| on ∂D for every parameter in B, where Q_c is the polynomial at the box centre and has a
   root in D. This gives a root of Q_B off the circle.
3. Branch and bound over the compactified domain: S = s² = t/(1 − t), and e up to 1. Boundary regimes get separate analytic
   or rescaled certificates: s → 0, s → ∞, e → 1, and the corner σ → 0, e → 1.
4. Independent re-check and audit; later a Lean formalisation of the reduction and the certificate checker.

Kill criteria: an interior configuration with δ → 0, which would mean 12 is reachable; then switch to certifying a
counterexample. Or a boundary regime that cannot be certified after two sessions.
