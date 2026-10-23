# Floating-point rigour of the certificates (audit item M3)

All certificates run in IEEE-754 binary64 with round-to-nearest. Write u = 2⁻⁵³ ≈ 1.11·10⁻¹⁶ and
γ_n = nu/(1 − nu). The standard bound for a computed sum of products is
|fl(Σ c_i Π x) − Σ c_i Π x| ≤ γ_N Σ |c_i Π x|, where N bounds the number of roundings along any term plus the length of
the summation chain (Higham, *Accuracy and Stability of Numerical Algorithms*, §3).

**Assumption (libm).** For |x| ≤ 64, the library sin and cos (MSVC CRT through numba) return a value within 2 ulps of
the exact one, and sqrt is correctly rounded; cosh and sinh (winding weights, |kL| ≤ 90) are within 2 ulps; atan,
atan2 and tan are used only where coarse accuracy suffices (the winding evaluation, which needs each argument step
to within 0.5 in total, and the 10⁻⁹ seam margins). These are the only assumptions beyond IEEE arithmetic
(measured: cosh, sinh ≤ 1.06 ulp on [−90, 90]; atan2 ≤ 0.68 ulp; AUDIT_7). (All sin/cos arguments
have |x| ≤ 62.2, §3. An empirical check, `code/libm_check.py`, found a maximum error of 0.755 ulp over 174,467
arguments in [−64, 64], including engine-form arguments and points next to multiples of π/2.) Where engine KQ
converts its box to the variables of its near part (ε = 1 − e, P_z = l sin ω, Q_z = l cos ω, μ = √((1 − t)/t)), every
bound is rounded outward with `dn`/`up` (two nextafter steps), which also widens the libm results by 2 ulps.

## 1. Coefficient polynomials (all `nbcoef_*` modules)

Generated code evaluates each coefficient c(x) = Σ_i a_i x^{m_i} with float constants a_i:
- the powers are built by chained multiplication (at most 17 per power);
- the monomials are formed by at most 7 further multiplications;
- the terms are summed left to right.

Maximum over all modules (`code/` and `code/c3/`): **265 terms per polynomial, total degree ≤ 40**. So
N ≤ 40 + 1 + 264 = 305 and γ₃₀₅ ≤ 3.4·10⁻¹⁴.

Constants that are not exactly representable carry a relative error ≤ 3u (sympy Float products; audit R2-3). That is
included in the same budget. More precisely, the certified object is the exact chart polynomial, with exact
coefficients. Its mean-value expansion at the box centre needs c(x_c) and ∇c(x_c). Their floating-point evaluation,
including the rounded constants, is off by ≤ γ₃₀₅·Pabs, which crad and the 10⁻¹³·a1 term of the derivative pad
absorb.

So |fl(c(x_c)) − c(x_c)| ≤ 3.4·10⁻¹⁴·Pabs(|x_c|) < 10⁻¹³·Pabs(|x_c|) = `crad`. Pabs is the absolute-value polynomial,
and it is evaluated by the same code with |a_i|.

The atoms at the box centre (sin, cos, sqrt of centre values) have relative error ≤ 2 ulps. That changes c(x_c) by at
most deg·2u·Pabs ≤ 40·2.2·10⁻¹⁶·Pabs, which is also inside `crad`.

## 2. Parameter variation (mean-value form)

For x in the box, |c(x) − c(x_c)| ≤ Σ_v r_v sup_box |∂_v c|, and sup_box |∂_v c| ≤ |∂_v c(x_c)| + w_v.
- **The variation bound w_v.** It is Pabs_v(|a_c| + Δa) − Pabs_v(|a_c|), for the absolute polynomial of ∂c/∂(atom),
  pushed through the chain rule. It is computed as fl(a1 − a0)·(1 + 10⁻¹²) + 10⁻¹³·a1. The second term exceeds
  2γ₃₀₅·a1 + u·|a1 − a0|, so the computed w_v is ≥ the exact bound.
- **The atom radii Δa.** They are padded: Δa = max(hi − a_c, a_c − lo)(1 + 10⁻¹²) + 10⁻¹⁵(|a_c| + 1).
- **The radii r_v.** They are padded: r_v = ½(hi − lo)(1 + 10⁻¹²) + 10⁻¹⁵(|x_c| + 1). This covers |fl(midpoint) − midpoint|
  ≤ ulp and the rounding of the half-width.
- **The bulk a/f charts.** The p, q radii carry an explicit +10⁻¹⁵ (nb_cert.box_coeffs).

## 3. Arc enclosures (nb_cert.arc_signs, make_engines_bu MV4 `signs`, c3 `root_signs`)

For each derivative order j and each arc [u_m − h, u_m + h], the code computes the value s_j at u_m, the radius
rad_j = rd + prad, and the tail rad_j + h·|f^{(j+1)}(u_m)|^+ + ½h²·B_{j+2}. Here B_{j+2} = Σ_k k^{j+2}·cmax_k
bounds sup|f^{(j+2)}| over the arc and the box.

**Trigonometric arguments (corrected per audit R3-1).** The argument is fl(fl(k·u_m) + fl(j·fl(π)/2)), with
k ≤ 6, j ≤ 7 and |u_m| ≤ 3π (the away arcs of engines KQ and KA reach 3π − φ₀), so |argument| ≤ 62.2.
Counting each rounding separately:
- fl(k·u_m) contributes ≤ u·56.6;
- fl(j·fl(π)/2) contributes ≤ u·11;
- the final sum contributes ≤ u·62.2;
- fl(π) − π contributes ≤ 1.3·10⁻¹⁶·j/2;
- the libm error is ≤ 2 ulps (≤ 2.3·10⁻¹⁶ for values ≤ 1).

The total is < 1.6·10⁻¹⁴, so |T̃_k − T_k| ≤ 1.6·10⁻¹⁴·k^j < 2·10⁻¹⁴·k^j. The engines add 2·10⁻¹⁴·Σ_k k^j(|c_k| + cr_k + Σ_v r_v(|dc_{v,k}| + dw_{v,k})). That
quantity dominates the trigonometric error of every sum: s, rd, prad, and the centre gradients g_v.

**Summation.** Each sum has at most 26 terms, so its error is ≤ γ₃₀·(sum of absolute terms). The engines add
10⁻¹²·(mag + rad), where mag is the sum of absolute terms of s, and rad already dominates the absolute terms of rd and
prad.

**Final sums.** lo = s − tail and hi = s + tail, with tail ← tail + 10⁻¹⁵(|s| + tail). That covers the at most 4
roundings in forming the tail and the difference.

**Polynomial root variable (C3 charts, near part of KQ).** Powers y^{k−j} have relative error ≤ γ₁₂ and the falling factorials are exact
integers, so the same 10⁻¹² and 10⁻¹⁵ terms cover them. There is no trigonometric evaluation.

**Winding certificate (wind_cert.comp_bounds).** This is the same scheme for the two components of G, applied to the
coefficient data scaled by the weights cosh(kL) and sinh(kL), with phases 0 and −π/2.
- The weights are computed in floating point with relative error ≤ (|kL| + 3)u (≤ 93u for |L| ≤ 15; corrected per
  audit R4-3), and a common positive factor 1 + 10⁻¹² is applied to all of
  them. Neither changes signs. The weight rounding perturbs a component by ≤ 3u·(its absolute sum), which is inside
  the 10⁻¹²·mag term.
- The Taylor term uses B₂ = Σ|w_k|k²(cmax₀ + cmax₁).
- After certification, the winding is evaluated in floating point at the box centre from the arc-boundary values.
  Each value is checked to lie in its arc's half-plane, and the accumulated angle must lie within 0.5 of a multiple of
  2π. Both checks make the integer exact.

**t-affine trick (charts with an affine variable).** At each endpoint t ∈ {t0, t1}, the combinations (1 − t)·X − t·Y
are padded by 10⁻¹⁵(|X| + |Y|) in cr and dw, and cmax is a convex combination with a (1 + 10⁻¹²) factor.

## 4. Arcs and ranges

**Arcs (ARC_PAD = 10⁻¹²).** Arcs [u_m − h, u_m + h] use h = ½·width + ARC_PAD. Consecutive arcs therefore overlap
(float gaps ≤ 10⁻¹⁵), and the first and last arcs reach beyond 0 and 2π. A zero in an overlap may be counted twice;
that is an over-count only.

**Sub-ranges of engines KQ and KA.** The away range [π + φ₀, 3π − φ₀] is widened by PAD = 10⁻¹² at both ends, which
absorbs fl(π) ≠ π and the endpoint rounding. In KQ the near range is |τ| ≤ tan(φ₀/2) + 10⁻⁹ and the away range starts
at φ₀ − 10⁻⁹, so the two ranges overlap. KA's box variable μ is mapped to t = 1/(1 + μ²) in floating point; adjacent
boxes share the same computed endpoint, so the KA boxes tile [fl(1/(1 + 0.0505²)), 1] without gaps, and the
sub-ulp sliver below that endpoint lies inside KQ1 (t ≤ 0.9975).

## 5. Second-order enclosures (code/nb_cert2.py, code/c3/make_c3_2.py)

These engines bound the parameter variation by min(mean-value bound of §3, second-order bound). The mean-value bound
is computed exactly as in §3. The second-order bound is

  Σ_v r_v |g_v| + Σ_{v≤w} q_vw |H_vw| + Σ_k |T_k| R_k + 10⁻¹²·(absolute sums).

The rounding is accounted for as follows.
- **Centre data.** The value, gradient and Hessian of every coefficient are generated straight-line polynomials
  (`nbcoef_f2`: at most 13 terms, total degree ≤ 8; `nbcoef_a2`: at most 9 terms, total degree ≤ 8;
  `nbcoef_c3_*_2`: at most 160 terms (the face charts FIP, FIM; 49 for Ca, 37 for KBQe, FE and FM), total degree
  ≤ 28). So N ≤ 28 + 1 + 159 = 188 in the bound of §1, and each centre quantity has evaluation error
  ≤ γ₁₈₈·Pabs ≤ 2.1·10⁻¹⁴·Pabs < 10⁻¹³·Pabs.
- **Chain rule.** The f-chart variables are (e, l, ω) with p = l sin ω, q = l cos ω; the a-chart variables are
  (η, l, ω) with in addition e = sin η, b = cos η. Their gradients and Hessians take
  float sin/cos at the centre (2 ulps) and at most 6 further products per term. Both errors are covered by the pads
  10⁻¹³·(Σ Pabs₁)·(1 + l_max)·Σr (gradient) and 10⁻¹³·(Σ Pabs₂ + Σ Pabs₁)·(1 + l_max)²·(Σr)² (Hessian). The C3
  charts need no chain rule, and there the pads are 10⁻¹³·(Σ Pabs₁)·Σr and 10⁻¹³·(Σ Pabs₂)·(Σr)².
- **Remainder R_k.** It is a sum of non-negative products of absolute-value polynomials at |x_c| + r and powers of
  the radii, rounded upward by the factor 1 + 10⁻¹².
- **Sums over the harmonics or the powers of the root variable**, in g, H and Σ|T_k|R_k. They are padded by 10⁻¹² of
  the corresponding sums of absolute terms, taken over both t-parts of the data. This covers the γ₃₀ summation error,
  the rounding of the t-combination (1 − t)X − tY and the products with r_v and q_vw.
- **Trigonometric arguments** (f-chart only). The magnitudes Σ q_vw|H|, |H| and R enter `aab`, so the 2·10⁻¹⁴·k^j
  term covers them as in §3.

After these pads both bounds are valid upper bounds of the exact parameter variation, so their minimum is too. The
containment was also tested directly against exact evaluations at random points of random boxes, with each bound
checked separately: no violation in about 10⁶ tests.

## 6. What is not covered by floating point

The decisions are discrete (sign certified or not). Every inequality used to certify a sign has its float error
dominated by the padding terms above. So each certified sign holds for the exact real-number function on the closed
arc and box.

The only external assumption is the libm accuracy stated at the top. An independent re-check of sampled leaves in
50-digit interval arithmetic (mpmath; see audit/ and `code/verify_mp.py`) confirms the recorded bounds.
