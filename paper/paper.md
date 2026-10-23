---
title: "A circle and a confocal ellipse have at most ten distance-stationary points: a computer-assisted answer to the second question of Problem 8 of Albouy, Cabral and Santos"
author: "Ansh Mishra, Aryan Senthilkumar"
---

## Abstract

<!-- STATUS -->
**Certification status (2026-10-23 15:00): not yet complete.** The certificate covers the whole parameter space except
10161 open boxes of the run families below, whose continuation runs are still in progress. The run families
that RESULTS.md marks complete replay exactly, and every one of their open or failed boxes is covered; sampled leaves
re-certify and pass 60-digit spot checks (verify/).

The open boxes lie in the near-parabolic region e ≥ 0.835 with a circle much larger than the pericentre distance
(μ = q/R small). Mostly they sit along the ridge l ≈ l*(ω), where two zeros of B merge and N = 10 is attained, and
at the seams of the C3 cell. There the Rolle certificates need very small boxes.

Until these runs finish, the bound N ≤ 10 is proved outside these open boxes only. Inside them, 60-digit sampling
(code/open_box_sample.py, 200 random points per family, verify/open_box_sample.out) finds at most 10 everywhere, and
often exactly 10. The runs continue with code/chain.sh (commands in code/launch_runs.sh), and `python code/check_all.py`
regenerates RESULTS.md.

| family | region | state |
|---|---|---|
| KQ0T | e ∈ [0.835, 0.9501], t ≥ 0.9975 (circle much larger than the pericentre distance) | open boxes 396 |
| KQ2 | e ∈ [0.9499, 0.99], t ≥ 0.9975 | open boxes 2791 |
| KQB | e ≥ 0.99, μ ≤ 0.0505, l ≥ 0.95, ω ∈ [3π/20, 9π/20] (the band) | open boxes 1000 |
| KQK1 | e ≥ 0.99, μ ≤ 0.0505, l ≤ 0.95 | open boxes 1071 |
| KQK2 | e ≥ 0.99, μ ≤ 0.0505, l ≥ 0.95, ω ≤ 3π/20 | complete itself; waits for: face-neighbourhood boxes pruned (p) but the face families are not all complete |
| KQF | face neighbourhood ε, μ ≤ 4·10⁻⁴: outer count \|τ\| ≥ 0.06 | open boxes 2292 |
| KQF2 | corner K_lw of the face neighbourhood, ε, μ ≤ SF2²: outer count \|τ\| ≥ 0.03 | open boxes 912 |
| KQS | corner K_lw of the face neighbourhood, SF2 < max(√ε, √μ) ≤ 0.02 | open boxes 877 |
| C3 Ca | C3 cell, chart Ca (Q_z dominant) | open boxes 822 |
<!-- /STATUS -->

Let a circle and an ellipse in $\mathbb R^3$ share a focus located at the centre of the circle. We give a
computer-assisted proof, complete except for the open boxes listed in the status note above, that the
squared distance between a point of the circle and a point of the ellipse has at most **ten** stationary points, and
that ten is attained. This answers the second question of Problem 8 (posed by G. F. Gronchi) in the list of
Albouy, Cabral and Santos [ACS12], which asked whether the general bound 12 can be reduced to 10 when one of the
two confocal conics is a circle.

**Reduction.** The count reduces to the real zeros of a trigonometric polynomial of degree six in one anomaly. Its
coefficients depend on four parameters: eccentricity, inclination, argument of pericentre, and the radius ratio.

**Certification.** On a branch-and-bound subdivision of the compactified four-dimensional parameter space, every
box receives one of two certificates:
- Rolle-run certificates, i.e. certified signs of derivatives on overlapping arcs;
- winding-number certificates, which locate at least one zero of the associated self-inversive polynomial strictly
  inside the unit disk.

**Degenerate corners.** These are resolved by exact weighted blow-up charts, including second-level blow-ups at
nested degenerate points:
- nearly coincident orbits;
- a nearly circular orbit perpendicular to a small circle;
- nearly parabolic orbits whose aphelion direction approaches the circle's axis while the circle becomes large;
- the parabolic face itself, where the pericentre distance becomes negligible compared with the circle and the
  polynomial acquires twelve real zeros in the limit.

**Rigour.** All floating-point errors are bounded explicitly. Every split tree is recorded and replayed
independently, and the computation was audited throughout by an independent checker.

## 1. Introduction

The squared distance $d^2(u,v)=|X(u)-C(v)|^2$ between points of two Keplerian orbits is the basic object behind the
**minimum orbit intersection distance (MOID)**, the standard first filter for impact hazard of asteroids and comets.
- Algorithms for the MOID compute *all* stationary points of $d^2$ on the torus [Gro02, BM19].
- Whether all of them have been found depends on knowing how many there can be.
- Kholshevnikov and Vassiliev [KV99] reduced the problem to a trigonometric polynomial.
- Gronchi [Gro02] obtained the bound 16 for two ellipses and gave examples with 12. In the circle–ellipse case he
  obtained 12, with examples reaching 10.
- Gronchi, Baù and Grassi [GBG23] proved sharper bounds in the planar case, and conjectured 10 for two coplanar
  ellipses.

Problem 8 of [ACS12] asks two questions. First, is 12 the maximum for two ellipses with positive eccentricities?
Second, if one of the eccentricities is zero, can the bound be reduced to 10? The configurations with infinitely
many stationary points (coinciding conics, concentric coplanar circles) are excluded. We answer the second question.

**Theorem 1.** Let a circle of radius $R>0$ and an ellipse with eccentricity $e\in(0,1)$ lie in $\mathbb R^3$ with a
common focus at the centre of the circle. Then $d^2$ has at most ten stationary points. The bound is attained, for
instance by $e=0.673$, $\sin I=0.992$, $\omega=0.338$, $R/a=0.7772$.

The proof is computer-assisted. Section 2 gives the reduction, Section 3 the two certificates, Section 4 the atlas
of charts and the coverage argument, Section 5 the degenerate corners, Section 6 the floating-point analysis, and
Section 7 the computation and its verification.

## 2. Reduction to a trigonometric polynomial

**Coordinates.**
- Scale the ellipse to semi-major axis $a=1$.
- Put the circle in the plane $z=0$: $C(v)=s(\cos v,\sin v,0)$ with $s=R/a$.
- Write the ellipse as $X(u)=(\cos u-e)P+\beta\sin u\,Q$, with $\beta=\sqrt{1-e^2}$ and
  $P=(\cos\omega,\sin\omega\cos I,\sin\omega\sin I)$, $Q=(-\sin\omega,\cos\omega\cos I,\cos\omega\sin I)$.
- The symmetries $\omega\mapsto-\omega$, $\omega\mapsto\omega+\pi$ and $I\mapsto\pi-I$ reduce the parameters to
  $e\in(0,1)$, $l=\sin I\in[0,1]$, $\omega\in[0,\pi/2]$, $s\in(0,\infty)$.
- Put $p=l\sin\omega$ and $q=l\cos\omega$.

**Branch functions.** For fixed $u$, $\partial_v d^2=0$ exactly at the nearest and farthest points of the circle,
provided the horizontal distance $\rho(u)$ of $X(u)$ from the axis is positive. So the stationary points are the
critical points of $g_\mp(u)=r^2+s^2\mp2s\rho$, where $r=1-e\cos u$. Squaring the conditions gives

$$\mathrm{Pol}(u)=A^2\rho^2-s^2\big((\rho^2)'\big)^2=h_-h_+,\qquad h_\mp=\rho\,g_\mp'(u),\quad A=(r^2)'.$$

This is a real trigonometric polynomial of degree at most six. Its coefficients are polynomials in $(e,\beta,p,q,s^2)$.

**Lemma 2.1 (counting).** $\mathrm{Pol}\not\equiv0$ for $e>0$. The number of stationary points is at most the number
of real zeros of $\mathrm{Pol}$ on the circle, counted with multiplicity. This holds also at points where
$A=(\rho^2)'=0$ or $\rho=0$; each such point carries at most two stationary points and a zero of multiplicity at
least two.

**Charts.** Every chart below certifies a function $c\cdot(\mathrm{Pol}\circ\varphi)$, where $c>0$ is smooth and
$\varphi$ is a diffeomorphism of the circle or of a sub-arc. The real zeros and their multiplicities are therefore
those of $\mathrm{Pol}$. The charts used are:
- the compactification $t=s^2/(1+s^2)$, with $f_t=(1-t)\Phi-t\Psi$, $\Phi=A^2\rho^2$, $\Psi=((\rho^2)')^2$;
- the true-anomaly chart (pericentre distance 1). Here $P_f=F_-F_+$, where
  $F_\mp=\mu k\sin f\,\sigma\mp DB$ is a positive multiple of $g_\mp'(f)$ and $\mu=q/R$;
- the tangent-half-angle form $Q(\tau)=P_f(1+\tau^2)^6$, with $\tau=\tan((f-\pi)/2)$;
- division by the positive normalisations of the blow-up charts.

## 3. Two certificates

**Lemma 3.1 (Rolle runs).** If $F^{(j)}$ has no zero on a closed arc, then $F$ has at most $j$ zeros there, counted
with multiplicity.

For a parameter box, the circle is covered by overlapping arcs. On each arc the engine certifies a strict sign of
$F^{(j)}$ for some $j\le6$, uniformly over the box. A run of consecutive arcs sharing a certified sign of $F^{(j)}$
contains at most $j$ zeros. A dynamic program then minimises the total over partitions into runs.

**Lemma 3.2 (parity).** The number of real zeros of a real trigonometric polynomial $F\not\equiv0$, counted with
multiplicity, is even.

*Proof.* Zeros of odd multiplicity are exactly the sign changes of $F$, and a continuous periodic function changes
sign an even number of times. $\square$

A box is therefore accepted when the certified total over the whole circle is at most $11$, because the true count is
even and hence at most $10$. A certified Rolle run also shows that $F^{(j)}\not\equiv0$ on its arc. So $F\not\equiv0$
and the lemma applies.

**Lemma 3.3 (winding).** Write $F(u)=\sum_{|k|\le6}c_ke^{iku}$ with $c_{-k}=\bar c_k$, and
$\mathcal Q(z)=\sum c_kz^{k+6}$.
- $\mathcal Q$ is self-inversive: $z^{12}\overline{\mathcal Q(1/\bar z)}=\mathcal Q(z)$.
- Suppose $\mathcal Q$ has no zero on $|z|=\rho<1$ and has $W\ge1$ zeros in $|z|<\rho$. Then it has $W$ zeros in
  $|z|>1/\rho$ (counting zeros at $\infty$ for a degree drop).
- Hence $F$ has at most $12-2W\le10$ real zeros.

On $|z|=e^L$ we have $\mathcal Q=e^{6L}e^{6i\theta}G(\theta)$, where

$$G(\theta)=\sum_k\cosh(kL)\,(a_k\cos k\theta+b_k\sin k\theta)+i\sum_k\sinh(kL)\,(a_k\sin k\theta-b_k\cos k\theta).$$

**How the winding certificate is applied.**
1. On every arc of an overlapping cover, the engine certifies a half-plane label for $G$, uniformly over the box.
   The label is one of $\operatorname{Re}G>0$, $\operatorname{Re}G<0$, $\operatorname{Im}G>0$, $\operatorname{Im}G<0$.
2. It then evaluates the winding number at the box centre, measuring both endpoint arguments of each arc in that
   arc's half-plane branch.
3. Finally $W=6+\operatorname{wind}(G)$.

The two certificates complement each other:
- **Rolle runs** need the box to resolve every real zero. That is expensive where real zeros cluster: at pitchfork
  bifurcations on the symmetric faces $\omega\in\{0,\pi/2\}$, and at six-fold points where the circle becomes large.
- **The winding certificate** never looks at the real zeros, so it is unaffected by those clusters.

## 4. The atlas

The parameter space is covered by the following charts, with overlapping boundaries.
- **Eccentric-anomaly chart** (with $\eta=\arcsin e$) for $e\in[0.199,0.8362]$ and all $l,\omega,t$.
- **Three charts for $e\le0.2$**:
  - $C_1$: polar in $(e,l^2)$;
  - $C_1'$: $t=e^2\tau/(1-\tau)$ for $l\ge0.4$;
  - $C_1''$: $e=\varepsilon''l^2$, $t=\varepsilon''^2\tau/(1-\tau)$;
  - inside them, the blow-up cones of the corners D1 and D2 (Section 5).
- **Charts for $e\ge0.835$.** These use the true-anomaly variables, with $t=\kappa^2/(1+\kappa^2)$, $\kappa=R/q$ and
  $\mu=q/R=\sqrt{(1-t)/t}$:
  - for $t\le0.5$, the true-anomaly chart on the whole circle;
  - for $t\ge0.5$, the engine KQ (Section 5). It counts the zeros on the aphelion arc $|f-\pi|\le0.9$ through
    $Q(\tau)$, and the rest with the true-anomaly chart;
  - for the cell $e\ge0.99$, $\mu\le0.0505$, $l\ge0.95$, $\omega\ge9\pi/20$: the C3 core charts near aphelion, a
    sign certificate in an intermediate range, and the true-anomaly chart away from aphelion;
  - for the face neighbourhood $1-e\le4\cdot10^{-4}$, $\mu\le4\cdot10^{-4}$ outside that cell: the face charts of
    Section 5 near aphelion, a sign certificate in a gap range, and engine KQ restricted to $|\tau|\ge0.06$.

Every whole-circle count is certified to be at most $11$, hence at most $10$ by Lemma 3.2. In the C3 cell the near,
intermediate and away counts are at most $9$, $0$ and $2$, and on the face neighbourhood they are at most $5$, $0$ and
$6$; in both cases the total is at most $11$.

Coverage of each seam is elementary. For example, for $e\le0.2$ and $l<0.4$: either $l^2\le2e$, and then
$\sqrt{e^2+l^4}<0.26$, which is the $C_1$ range; or $l^2>2e$, so $\varepsilon''<0.5$, which is the $C_1''$ range. The complete table is
in `PROOF.md`.

## 5. Degenerate corners

**D1: nearly coincident orbits** ($e,l\to0$, $s\to1$).
- In the coordinates $(e,\theta=l^2/e,\tau=t-\tfrac12)$,
  $$f/e^2=2\sin u\,[\theta\sin2(u+\omega)-4\tau\sin u-e\sin2u]+O(2).$$
  This is homogeneous of degree one. Three cone charts resolve it, with exact divisibility.
- The leading term vanishes identically only in the direction $l=e$, $\omega=0$, $s=1$.
- A second homogeneous blow-up there has leading term
  $2\sin u\,[2a\sin2u+2p\cos2u-4\tau\sin u+e\,r_1(u)]$, where $r_1$ contains $\sin3u$. It never vanishes identically.

**D2: a nearly circular orbit perpendicular to a tiny circle** ($e\to0$, $s\sim e$, $l\to1$, $\omega\to0$).
- In the chart $C_1'$ the function vanishes identically at the limit point.
- Its linear part is $-8\tau's^2c^2+2\kappa(2s^4+4s^2c^2)-4p\,sc^3-4e\,s^2c^3$, with $s=\sin u$, $c=\cos u$,
  $\kappa=1-q$. Its four functions are independent, and one homogeneous blow-up resolves the corner.
- In the cone where $\tau'$ dominates, the leading term has a four-fold zero at $u=\pi$ along $p=0$,
  $2\kappa+e=2\tau'$. That cone chart cannot certify boxes there. The boxes it leaves open are covered in one of two
  ways:
  - their eccentricity is larger than that of every leaf the corner charts must cover;
  - or they lie in the cone where $e$ dominates, extended to $\tau'/e\le1.28$, which certifies them with the winding
    certificate.

**C3: the e → 1 corner with the aphelion direction on the circle's axis.**
- In the tangent-half-angle form, substitute $\tau=\nu y$, $Q_z=\nu a$, $1-P_z=\nu^2K$, $1-e=\nu^2E$, $\mu=\nu^3m$.
  Then
  $$Q=\nu^{10}\big[\hat P_0(2y)+O(\nu^2)\big],\quad \hat P_0(x)=16m^2x^2(\lambda+(x-a)^2)-(E+\tfrac{x^2}2)^2\hat B(x)^2,$$
  with $\hat B(x)=2E(x-a)-4Kx+3ax^2-x^3$ and $\lambda=2K-a^2$.
- The near-aphelion roots are counted in the four parameter cones, a sign certificate covers the intermediate
  range, and an arc count handles the region away from aphelion.
- At $E=m=0$, $\hat P_0$ has a six-fold zero at $y=0$. A second blow-up of $(y,E,m)$ with weights $(1,2,2)$ shows
  that at most four roots lie near $y=0$.

**The parabolic face** ($e=1$, $\mu=0$, general direction). Here $P=-D^2B^2$, and for $l$ above a threshold
$l^*(\omega)$ (from $0.707$ at $\omega=0$ to $0.9985$ at $\omega=1.41$) it has twelve real zeros with multiplicity: a
six-fold zero at aphelion and three double zeros of $B^2$. No box touching the face can therefore be certified with
total $\le11$ in the original variables.
- With weights $(1,2,2)$ for $(\tau,\varepsilon,\mu)$ and radius $s$, $G=Q(sy)/s^6$ is a polynomial, with leading part
  $$G_0=16c^2y^2\big[4m^2-c^2(E+2y^2)^2\big],\qquad c^2=1-P_z^2,\ E=\varepsilon/s^2,\ m=\mu/s^2 .$$
- For $c>0$ its real zeros are $y=0$ (double) and at most one pair $y^2=(2m/c-E)/2$: at most four near aphelion,
  against six from the rest of the circle.
- Two cone charts count the zeros on $|y|\le3$ (at most $5$), two charts certify that $Q$ has no zero for
  $3s\le|\tau|\le0.06$, and engine KQ counts $|\tau|\ge0.06$ (at most $6$). Outside the C3 cell, $c>0.156$.

**The near-parabolic region (engine KQ).** As $e\to1$ the factor $D=1+e\cos f$ becomes small near aphelion. In the
Fourier basis the true-anomaly polynomial then takes values of order $\varepsilon^4$ ($\varepsilon=1-e$) there, while
its coefficients are of order one. No enclosure of moderate width resolves its zeros near aphelion.
- KQ writes $Q(\tau)$, for $|\tau|\le\tan(0.45)$, as a polynomial in $\tau$ in the monomial basis centred at aphelion.
  Its coefficients are polynomials in $(\varepsilon,P_z,Q_z,\mu)$.
- In this form the pair of zeros near aphelion, at $\tau=O(\varepsilon^2)$, is an ordinary, well-conditioned cluster.
- The rest of the circle is counted in the true-anomaly chart, and the two ranges overlap.
- $Q$ is polynomial in the parameters, so the points where the orbit crosses the circle's axis ($\sigma=0$) are not
  singular. The same engine therefore also covers the band where the orbit plane nearly contains the axis.
- On each box KQ first tries the winding certificate of the true-anomaly polynomial, and counts as above only when
  that fails. The winding certificate does not resolve the real zeros. It is therefore effective near $e\approx0.95$,
  where up to ten real zeros cluster while the remaining complex pair stays off the unit circle.

## 6. Floating-point rigour

The certificates run in IEEE binary64.
- **Coefficient polynomials.** At most 265 terms each, of degree at most 40. They are evaluated with error at most
  $\gamma_{305}P_{\mathrm{abs}}<10^{-13}P_{\mathrm{abs}}$.
- **Parameter variation.** It is bounded by the mean-value form, using absolute-value polynomials for the variation
  of the gradients. The radii are padded for the float centre.
- **Second-order bound.** In the near-parabolic engine and in the hardest C3 chart, the engine also bounds the
  variation by the exact gradient and Hessian at the box centre plus a third-order Taylor remainder, and uses the
  smaller of the two bounds. Near $e\approx0.95$ the absolute-value term of the mean-value form exceeds the true
  second-order term by factors of 25 to 500, so this saves several levels of subdivision.
- **Trigonometric arguments.** Arguments are smaller than 63 in modulus. Their rounding, and the error of the
  library sine and cosine (assumed within two ulps), are covered by an explicit padding of $2\cdot10^{-14}k^j$ per
  term.
- **Arcs and final sums.** Arcs overlap by $10^{-12}$, and the final sums carry a relative padding of $10^{-15}$.
- **Engine KQ.** Its box is converted to $(\varepsilon,P_z,Q_z,\mu)$, with every bound rounded outward.

Details are in `FLOAT_RIGOUR.md`.

## 7. The computation and its verification

**Recording.** Every run records, for each initial box, the complete split tree, one character per box: split, leaf
with its bound and method, pruned (unphysical), or failed. Winding leaves also store their test radius. A replay tool
re-derives every leaf exactly. Failed leaves and open boxes are covered by continuation runs, or by the D2 corner
charts in one case, or lie outside the range that a face chart needs. A continuation run may split the coordinate that matters most first: the circle size in the
large-circle runs, one ratio in a C3 chart, and the blow-up scale near the parabolic face. The split rule affects only efficiency.

**Independent checks.**
- 60- and 100-digit spot checks of sampled leaves against the original polynomial.
- Re-certification of leaves, including every winding leaf, with the final code.
- Twelve rounds of independent auditing of every derivation, generated formula, engine, coverage claim and run. Among
  other things, the audits found:
  - a coverage seam;
  - an unrecorded-leaf problem;
  - rounding gaps in the trigonometric arguments;
  - a branch error in the winding evaluation;
  - the twelve-root obstruction at the parabolic face, which led to the face charts.

  All of these were repaired; the final pass (AUDIT_12) re-certified every leaf of the smaller run families and
  large samples of the others.

The per-run statistics are in `RESULTS.md`.

## References

- [ACS12] A. Albouy, H. E. Cabral, A. A. Santos, Some problems on the classical n-body problem, *Celest. Mech. Dyn.
  Astr.* 113 (2012) 369–375.
- [BM19] R. V. Baluev, D. V. Mikryukov, Fast error-controlling MOID computation for confocal elliptic orbits,
  *Astron. Comput.* 27 (2019) 11–22.
- [GBG23] G. F. Gronchi, G. Baù, C. Grassi, Revisiting the computation of the critical points of the Keplerian
  distance, *Celest. Mech. Dyn. Astr.* 135 (2023).
- [Gro02] G. F. Gronchi, On the stationary points of the squared distance between two ellipses with a common focus,
  *SIAM J. Sci. Comput.* 24 (2002) 61–80.
- [KV99] K. V. Kholshevnikov, N. N. Vassiliev, On the distance function between two Keplerian elliptic orbits,
  *Celest. Mech. Dyn. Astr.* 75 (1999) 75–83.
