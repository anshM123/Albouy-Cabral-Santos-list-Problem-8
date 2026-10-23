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
**Certification status (2026-10-23 18:30): complete.** Every run family in RESULTS.md is complete (`code/check_all.py`: ALL COMPLETE). Every recorded split tree replays exactly, and every open or failed box is covered by a continuation run or by a documented chart. Sampled leaves of every run file re-certify with the final code (verify/), and the final audit (audit/AUDIT_15.md) re-certifies every tau- and disc-winding leaf. 60-digit spot checks (verify/) cover sampled leaves of every whole-circle family and C3 chart, and the disc-winding leaves of KQF and KQF2. The other outer-count leaves of KQF and KQF2, and the D2 charts, are covered by the audit's 100-digit spot checks.
<!-- /STATUS -->

**Theorem (computer-assisted).** For a circle and an
ellipse (eccentricity e ∈ (0, 1)) in ℝ³ with a common focus at the centre of the circle, the squared distance between
their points has at most **10** stationary points. The bound is attained, for
example at e = 0.673, sin I = 0.992, ω = 0.338, R/a = 0.7772.

So the answer to the second question is **yes**. The excluded configurations with infinitely many stationary points
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
2. **Certificates.** On a branch-and-bound subdivision of the compactified parameter space, each box gets one of these
   certificates (PROOF.md §2):
   - Rolle runs, i.e. certified signs of derivatives on arcs;
   - a winding-number certificate for the associated self-inversive polynomial;
   - a tau-winding certificate: a circle in the upper half-plane that encloses a zero of the aphelion polynomial for
     every point of the box. This forces a non-real pair of zeros, hence at most 10 real ones (C3_CORE.md §8);
   - near the parabolic face, a disc-winding count: six zeros certified inside a small disc about aphelion leave at
     most 6 real zeros outside it (C3_CORE.md §9).

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
   - Fifteen rounds of independent audit (audit/); AUDIT_15 is the final pass over the complete run set.

## Where to read

| file | contents |
|---|---|
| `PROOF.md` | the proof: reduction, certificates, atlas and coverage lemma, degenerate corners |
| `STATEMENT.md` | the canonical statement, the reduction lemma, sanity tests |
| `C3_CORE.md` | the e → 1 corner with the aphelion direction on the circle's axis (C3 core), the l ≈ 1 band, the face e = 1, q/R = 0, and the tau-winding and disc-winding certificates |
| `DEGENERACY_TREE.md` | the list of degenerate corners and how each is resolved |
| `CORNER_E1.md` | background analysis of the e → 1 corner (not used by the final certificate) |
| `FLOAT_RIGOUR.md` | floating-point error analysis |
| `RESULTS.md` | per-run results: initial boxes, leaves, bound histograms, failures, replay status |
| `audit/` | independent adversarial audits (AUDIT_1–AUDIT_15) with their check scripts and outputs |
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
python recert_tw_c3.py Ca FE FM                   # re-certify every tau-winding leaf of the c3 charts
python disc_mp_check.py runs/KQF.r2.json          # 60-digit outer counts at disc-winding leaves
PY=python sh final_verify.sh                      # all of the above for every run file (outputs in verify/)
```

The recorded runs store, for every initial box, the complete split tree with one character per box. So `check_all.py`
and `bb_replay.py` can re-derive every leaf without rerunning the search.

Some runs were produced by earlier versions of an engine, each a valid certificate. The bulk run switched from the
mean-value enclosure to the second-order one after 3,825 initial boxes (and its initial boxes with e ≥ 0.836 were
handed to the f-chart and KQ runs F1, KQ0, KQ0T), and the KQ runs went through three versions; their continuation
levels near the face e = 1, q/R = 0 use the split rule of `splitrule.py`. The last continuation levels of the KQ runs
and of chart Ca also use the tau-winding certificate (leaf code '7', `tau_wind.py`). Those of the outer runs KQF and
KQF2 also use the disc-winding outer count (leaf code '8'); see C3_CORE.md §8–§9. Also, the C3 chart engines other than Ca and KBQe were generated before `make_c3.py` started returning the
minimum bound over its refinement rounds instead of the last one. The current engines are never weaker, arc by arc,
than the versions that produced the records, so `bb_replay.py --recert` re-certifies the recorded leaves with the
current code.

## Citation

See `CITATION.cff`.

## License

MIT (see `LICENSE`).
