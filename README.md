# Albouy–Cabral–Santos list, Problem 8 (second question): a circle and a confocal ellipse

**Authors:** Ansh Mishra, Aryan Senthilkumar

## The question

Problem 8 of A. Albouy, H. E. Cabral and A. A. Santos, *Some problems on the classical n-body problem*,
Celest. Mech. Dyn. Astr. 113 (2012) 369–375 (posed by G. F. Gronchi), asks:

> Consider the distance d between two points on two distinct confocal ellipses in ℝ³. If both eccentricities are > 0,
> is 12 the maximum number of stationary points of the function d²? **If only one of the eccentricities is 0, can this
> upper bound be reduced to 10?**

## The result

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

**Theorem (computer-assisted; the certification is not yet complete, see the status above).** For a circle and an
ellipse (eccentricity e ∈ (0, 1)) in ℝ³ with a common focus at the centre of the circle, the squared distance between
their points has at most **10** stationary points. The bound is attained, for
example at e = 0.673, sin I = 0.992, ω = 0.338, R/a = 0.7772.

So the answer to the second question is **yes**, once the open boxes in the status above are certified. The excluded configurations with infinitely many stationary points
(coinciding conics, concentric coplanar circles) have e = 0 and are outside the question.

## Prior work

- Kholshevnikov and Vassiliev (1999) reduced the problem to a trigonometric polynomial of degree 8.
- Gronchi (2002) proved the bound 16 for two ellipses, with examples reaching 12. In the circle–ellipse case he proved
  12, with examples reaching 10.
- Gronchi, Baù and Grassi (2023) proved sharper bounds in the planar case and conjectured 10 for two coplanar
  ellipses.

Before this work the circle–ellipse bound in ℝ³ stood at 12.

## How the proof works

The proof is computer-assisted:
1. **Reduction.** The count reduces to the real roots of a degree-6 trigonometric polynomial in one anomaly, depending
   on 4 parameters (STATEMENT.md, PROOF.md §1).
2. **Certificates.** On a branch-and-bound subdivision of the compactified parameter space, each box gets one of two
   certificates (PROOF.md §2):
   - Rolle runs, i.e. certified signs of derivatives on arcs;
   - a winding-number certificate for the associated self-inversive polynomial.

   A parity lemma (the number of real zeros of a real trigonometric polynomial is even) lets a certified total of 11
   count as 10.
3. **Degenerate corners.** Weighted blow-up charts handle the corners, including second-level blow-ups (PROOF.md §4,
   C3_CORE.md, DEGENERACY_TREE.md). The corners are:
   - near-coincident orbits;
   - a nearly circular orbit perpendicular to a tiny circle;
   - nearly parabolic orbits pointing along the circle's axis;
   - the parabolic face e = 1 with a pericentre much smaller than the circle (q ≪ R), where the polynomial has 12 real
     roots in the limit; a blow-up with weights (1, 2, 2) in (anomaly, 1 − e, q/R) separates them (C3_CORE.md §6).
4. **Near-parabolic engine.** For e ≥ 0.835 and circles not much smaller than the pericentre distance, engine KQ counts
   the zeros near aphelion in a tangent-half-angle polynomial written in ε = 1 − e, which avoids the cancellation of
   the Fourier basis as e → 1.
5. **Rigour and verification.**
   - Explicit floating-point error bounds (FLOAT_RIGOUR.md).
   - Recorded split trees with exact replay.
   - Re-certification and high-precision spot checks.
   - Twelve rounds of independent audit (audit/); AUDIT_12 is the final pass over the finished run families.

## Where to read

| file | contents |
|---|---|
| `PROOF.md` | the proof: reduction, certificates, atlas and coverage lemma, degenerate corners |
| `STATEMENT.md` | the canonical statement, the reduction lemma, sanity tests |
| `C3_CORE.md` | the e → 1 corner with the aphelion direction on the circle's axis (C3 core), the l ≈ 1 band, and the face e = 1, q/R = 0 |
| `DEGENERACY_TREE.md` | the list of degenerate corners and how each is resolved |
| `CORNER_E1.md` | background analysis of the e → 1 corner (not used by the final certificate) |
| `FLOAT_RIGOUR.md` | floating-point error analysis |
| `RESULTS.md` | per-run results: initial boxes, leaves, bound histograms, failures, replay status |
| `audit/` | independent adversarial audits (AUDIT_1–AUDIT_12) with their check scripts and outputs |
| `paper/` | the write-up |
| `code/` | generators, generated engines, drivers, verification tools, and the recorded runs (`runs/`, `c3/runs/`) |

## Reproducing

You need Python 3.12 with numpy, numba, sympy and mpmath. The generated coefficient modules (`nbcoef_*.py`) and
engines (`nbeng_*.py`) are included. To regenerate them, run from `code/`:

```bash
python gen_nb.py a && python gen_nb.py f && python gen_nb.py fe && python gen_nb.py fm5 && python gen_nb.py fm8
python gen_nb2.py && python gen_nb2a.py            # second-order f- and a-chart data (nb_cert2)
python gen_chart.py C1 && python gen_chart.py C1p && python gen_chart.py C1pp && python make_engines.py C1 C1p C1pp
python gen_bu.py B1 B2 B3p B3m L2E L2P L2Ap L2Am L2Tp L2Tm D2E D2E2 D2P D2K D2Tp D2Tm
python make_engines_bu.py B1 B2 B3p B3m L2E L2P L2Ap L2Am L2Tp L2Tm D2E D2E2 D2P D2K D2Tp D2Tm
cd c3 && python gen_c3.py Ca CK CE Cm PhiP PhiM L2aE L2am L2aYp L2aYm L2KE L2Km L2KYp L2KYm KBQe && python make_c3.py Ca CK CE Cm PhiP PhiM L2aE L2am L2aYp L2aYm L2KE L2Km L2KYp L2KYm KBQe && cd ..
cd c3 && python gen_c3.py FE FM FIP FIM KBPe && python make_c3.py FE FM FIP FIM KBPe && cd ..   # face charts (C3_CORE.md §6), perihelion polynomial (§7)
cd c3 && python gen_c3_2.py KBQe KBPe Ca FE FM FIP FIM && python make_c3_2.py KBQe KBPe Ca FE FM FIP FIM && cd ..   # second-order engines
```

To rerun the certificates and check them, run from `code/`:

```bash
PY=python sh launch_runs.sh                       # all certificate runs (many CPU-hours; see the comments in the
                                                  # script for the continuation levels, chain.sh)
python check_all.py                               # replay every recorded run, check coverage, write RESULTS.md
python bb_replay.py runs/bulk.json --recert 500   # re-certify sampled leaves of one run
python verify_mp.py runs/bulk.json --n 500        # 60-digit spot checks of sampled leaves
```

The recorded runs store, for every initial box, the complete split tree with one character per box. So `check_all.py`
and `bb_replay.py` can re-derive every leaf without rerunning the search.

Some runs were produced by earlier versions of an engine, each a valid certificate. The bulk run switched from the
mean-value enclosure to the second-order one after 3,825 initial boxes (and its initial boxes with e ≥ 0.836 were
handed to the f-chart and KQ runs F1, KQ0, KQ0T), and the KQ runs went through three versions; their continuation
levels near the face e = 1, q/R = 0 use the split rule of `splitrule.py`. Also, the C3 chart engines other than Ca and KBQe were generated before `make_c3.py` started returning the
minimum bound over its refinement rounds instead of the last one. The current engines are never weaker, arc by arc,
than the versions that produced the records, so `bb_replay.py --recert` re-certifies the recorded leaves with the
current code.

## Citation

See `CITATION.cff`.

## License

MIT (see `LICENSE`).
