# AUDIT 14: the disc-winding outer count (leaf code '8'), the new TW contour choice, and the KQF2 r3 record replacement

Independent audit, 2026-10-23, 17:41–17:50. At most 2 Python processes. Scripts and outputs are in `audit/round14/`.

## Verdict

**All three changes are sound.** The lemma is correct, and its region O_{T_X} is exactly what the §6 outer count
certifies. The code implements it, and every '8' leaf I sampled re-certifies and agrees with my own 60-digit count.
The new TW contour choice cannot affect validity. The KQF2 r3 record replacement leaves a consistent, complete family.

One small hardening is recommended (R14-1); it is not a soundness issue.

| # | item | status |
|---|---|---|
| (a) | lemma, O_T versus count_outer / count_peri, degree drop at f = 0 | OK |
| (b) | code paths (tau_wind, kq_run, check_all, bb_replay) | OK |
| (c) | independent re-certification of '8' leaves (KQF.r2, KQF.r3, KQF2.r3, KQF2.r4) | OK |
| (d) | anything more for check_all | not needed for soundness; R14-1 (minor) |
| 2 | TW contour candidates | OK (validity unaffected) |
| 3 | KQF2 r3 record replacement | OK (family replays complete) |

## (a) The lemma: OK

**Setting.** Q(τ) = P(f)(1 + τ²)⁶, τ = tan((f − π)/2) (chart KBQe), of formal degree 12. T > 0, and
O_T = {f : |τ(f)| ≥ T} ∪ {f = 0}. Suppose Q(·; x) has no zero on |τ| = T and W zeros (with multiplicity) in |τ| < T.

**Proof.** Take a real zero of P in O_T with f ≠ 0. It corresponds to a real zero of Q with |τ| ≥ T, and since
|τ| = T is excluded, with |τ| > T. These are among the deg Q − W zeros of Q outside the closed disc. A zero of P at
f = 0 has order 12 − deg Q. The count of P on O_T is therefore ≤ (deg Q − W) + (12 − deg Q) = 12 − W.

The disc may contain real or non-real zeros; both reduce the outside count. With W ≥ DISC_W = 6, the outer count is
≤ 6.

**Degree drops.**
- At either end the identity is the one verified in AUDIT_13 (a), at 41 points including deg Q = 10.
- f = 0 lies in O_T: it is the point τ = ∞, inside the away arc |f − π| ≥ φ₀.
- f = π (τ = 0) is inside the disc, as it should be: it belongs to the near region.

**The region is the one the §6 outer count certifies.**
- count_outer: KBQe on T_X ≤ |τ| ≤ T_HI = tan(PHI_0/2), both signs, plus the away arc |f − π| ≥ PHI_0 − 1e-9.
- count_peri with t_lo = T_X:
  - KBQe on T_X ≤ |τ| ≤ tan(φ₀/2) + 1e-9;
  - the f-chart arcs [π + φ₀ − 1e-9, 2π − peri] and [peri − 1e-9, π − φ₀];
  - KBPe on |tp| ≤ tan(peri/2) + 1e-9 (around f = 0);
  - coverage audited in AUDIT_10 h16.
- Both cover exactly {|τ| ≥ T_X} ∪ {f = 0} = O_{T_X}, closed at |τ| = T_X.

So '8' certifies the same quantity as the codes 'a'..'g' of a KQF-engine file. The face total near + gap + outer stays
≤ 5 + 0 + 6 = 11, which gives ≤ 10 by parity. The near and gap parts (FE/FM, FIP/FIM, up to |τ| = T_X) are unchanged.

**Homotopy and enclosure.** The same audited contour_ok (AUDIT_13 c, d) is used with centre 0 and radius T_X. The arc
discs cover the circle, and the argument steps are exact. W is constant over the box.

## (b) Code paths: OK

**tau_wind.**
- `certify_disc(bx, tx)` calls `contour_ok(*engine(0).box_coeffs2(*params(bx)), 0j, tx, ...)`. That is chart KBQe only,
  centre 0, radius tx, as the lemma requires. It accepts when ok and W ≥ DISC_W = 6.
- `check_disc` is the same call with rounds 10 / cap 65,536 against 8 / 16,384. The refinement is deterministic, so a
  certified leaf always re-certifies.

**kq_run.**
- `--disc` is rejected without `--face`.
- In face mode, '8' is tried after the C3-cell 'o' prune, with the module T_X. That T_X is set from `--faces` (3·radius)
  and passed to the workers.
- The record's hist key for '8' is 'd'.
- The header records disc = True (KQF.r2, KQF.r3, KQF2.r3 and KQF2.r4 all show it).

**check_all.** `replay` maps '8' to hist[6] only for engine KQF. In any other engine it falls through to "bad code",
so a whole-circle KQ family can never accept an outer-count certificate. The header T_X of every KQF-engine file is
pinned to 3·radius (R11-1).

**bb_replay.** Method 4 calls `check_disc(bx, H['T_X'])`, and check_all ties that header value to the family.

## (c) Independent re-certification of '8' leaves: OK

j4_disc_leaves.py, on the closed files:

| file | T_X | '8' leaves | sampled (× 5 points) | W at every point | max outer real count | min distance of a zero to \|τ\| = T_X (relative) | enclosure checks, max ratio |
|---|---|---|---|---|---|---|---|
| KQF.r2 | 0.06 | 2,375 | 30 | 6 | 6 | 0.364 | 2,400, 0.99974 |
| KQF.r3 | 0.06 | 661 | 20 | 6 | 6 | 0.595 | 1,600, 0.998 |
| KQF2.r3 | 0.030015 | 872 | 30 | 6 | 6 | 0.312 | 2,400, 0.996 |
| KQF2.r4 | 0.030015 | 543 | 30 | 6 | 6 | 0.444 | 2,400, 0.995 |

For every sampled leaf:
- `check_disc` with the header T_X re-certifies.
- With my own 60-digit zeros of Q (my P, FFT, degree-12 roots, τ = −i(z+1)/(z−1)), W is exactly 6 at the centre and
  at 4 random physical points (constant, as the homotopy says).
- The outer real count is at most 6, and never above 12 − W.

The enclosure test checks |F(w; x) − b0| ≤ R on the disc contour (c = 0, r = T_X), at box corners and random points
and at w on and inside the arc discs. There were 0 violations. Ratios close to 1 are expected: on the disc boundary
|b_1|D ≈ dev. The largest, re-run to 12 digits, is 0.999737015, a slack of 2.6e-4 relative, far above float level.

The per-record '8' counts equal the record hists ('d'). The full replay (g2, j5) has 0 problems, including the
leaf-code histogram cross-check.

## (d) Does check_all need more? Not for soundness

An '8' leaf is fixed by (box, T_X). T_X is pinned per file, and bb_replay re-certifies it. Tying '8' to the header
disc flag would only catch accidental code collisions; I don't think it's necessary.

**R14-1 (minor, hardening).** The meaning of '8' (outer ≤ 6) depends on the code constant tau_wind.DISC_W, which must
equal LIMIT_OUTER = 6. If DISC_W were lowered, check_all would still count an '8' leaf as 6 while bb_replay's
check_disc certified only ≤ 12 − DISC_W. Add `assert tau_wind.DISC_W == 6` in check_all and bb_replay (or pass W ≥ 6
explicitly to check_disc).

## 2. TW contour candidates: OK (validity unaffected)

`candidates` now returns centres at each upper zero, each midpoint of two upper zeros, and the mean of all of them,
with radius r = (lo + hi)/2 where hi = min(next zero distance, Im c). So 0 < r < Im c always.

certify_tw still requires c.imag > r > 0, and check_all requires c_im > r > 0 on every stored contour. A contour that
encloses a cluster of zeros only needs W ≥ 1. Any enclosed zero is non-real because the closed disc lies in
Im > c_im − r > 0. The candidate choice affects only success.

## 3. KQF2 r3 record replacement: OK

- All 222 root boxes of the superseded records file (`runs/superseded/KQF2.r3.json.records.with_minw_fails.jsonl`)
  reappear among the 912 roots of the current `KQF2.r3.json.records.jsonl`. None is missing, so the replaced root has
  a fresh record.
- My full replay of KQF2 (base + r1..r4, 2,248 records, 44,488 leaves: 1,415 '8' and 2,069 '7') finds every
  continuation root owed by the previous level, every record consistent, and nothing owed at the end. **KQF2 is
  complete.**
- The same replay of KQF (base + r1..r3, 3,334 records, 346,293 leaves: 3,036 '8' and 8,698 '7') finds only the 28
  owed boxes inside K_lw, which are excused. **KQF is complete.**
- So the replacement is fine for coverage: check_all's owed logic applies to the final file as written. The superseded
  file is not read by any tool.

## My tools (for the final pass)

- **g2:** '8' gives bound 6 in kqf families only; hist 'disc'; the record-hist key is 'd'.
- **g1:** '8' leaves are re-certified with check_disc and the header T_X, plus a 100-digit outer count ≤ 6.
- **g3:** '8' gives an outer count ≤ 6.
- **g9:** a new kind 'disc' uses check_disc with the header T_X.
- **Final-pass scripts:** `--kinds reg,deep,tw,disc`.
- The round-13 versions are kept in `round14/*_r13.py`.
