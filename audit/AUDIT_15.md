# AUDIT 15: final pass over the complete run set (2026-10-23, 18:30 state)

Independent audit, 18:31–18:50, on 8 cores. Scripts are in `audit/round15/`; outputs are in `audit/round15/out/`
(the tools come from `audit/round7`, `round13` and `round14`).

## Verdict

**The certificate is complete, and I found no soundness concern.**
- **Completeness.** check_all.py, run by me, reports **ALL COMPLETE** (51/51 families). My own independent replay of
  every family (g2) agrees: every tree replays exactly, the codes and bounds are within the family limits, the pruning
  predicates are exact, every owed box is covered, and every family's last level ends with nothing owed.
- **Every** tau-winding leaf (175,932, including all 6,218 C3-chart leaves of Ca r8–r10) and **every** disc-winding
  leaf (4,451) of the new files re-certifies.
- **Samples** of every other leaf kind of every new file re-certify, and 100-digit spot checks with my own formulas
  find no violation.

Five documentation issues remain (R15-1 to R15-5). One of them, R15-1, is a false statement in the paper's face
paragraph (the corner K is not described) and should be fixed before the "not yet complete" qualifiers are dropped.

## 1. Completeness

**check_all.py** (unchanged project code, RESULTS.md redirected to `out/RESULTS.md`) reports
**"All families complete: YES"**, 51 families, all `complete`. In particular:
- **KQK2** (base + r1..r6) is now complete. Its only blocker, the face dependency, is gone: FE, FM, FIP, FIM, KQF,
  KQF2 and KQS are all complete.
- **KQF** has base + r1..r3. The empty levels r4–r14 moved to `runs/superseded/KQF_empty_levels` have 0 records each,
  and `KQF_r4_boxes.json` is empty. Nothing was lost. chain.sh now parses only the summary line
  `open boxes N ->` (line 18).
- **The last-level box lists** are all empty ({"boxes": []}): KQ0T_r6, KQ2_r9, KQB_r11, KQF_r4, KQF2_r5, KQK1_r13,
  KQK2_r7, KQS_r6, KQ1_r2, KQ0_r1, Ca_r11.

**My independent replay (g2, with every round-11 to 14 extension)** reports
**"ALL FAMILIES REPLAYED, COMPLETE AND WITHIN THEIR BOUNDS: YES"**, with every family OK, open 0 and problems 0. It
checks:
- per-file pins of FACE_S and T_X, and the base domains of KQF2 and KQS;
- each record's own leaf-code histogram;
- '7' contours in the open upper half-plane, with which = 0 on the C3 charts;
- '8' only in KQF-engine files;
- the p-dependency;
- the K / K_lw excuses;
- the D2 rules.

Family totals:

| family | files | leaves |
|---|---|---|
| bulk | | 9,545,140 |
| KQ0T | 6 | 496,001 |
| KQ2 | 9 | 449,344 |
| KQB | 11 | 230,218 |
| KQK1 | 13 | 408,863 |
| KQK2 | 7 | 127,845 |
| KQF | 4 | 346,293 |
| KQF2 | 5 | 44,488 |
| KQS | 6 | 54,897 |
| C3 Ca | 11 | 2,874,961 |

The KQF2 r3 record replacement (round 14) passes.

## 2. Re-certification of the new files (since 15:00)

The new files are: KQ0T.r4–r5, KQ2.r5–r8, KQB.r7–r10, KQF.r2–r3, KQF2.r3–r4, KQK1.r7–r12, KQS.r2–r5 and Ca.r8–r10.

**Every tau-winding and disc-winding leaf (j6), 0 failures.** Each tree was replayed with its header scale and split
rule, then check_tw, check_disc (with the header T_X) and check_tw_chart were run:

| file | '7' leaves | '8' leaves |
|---|---|---|
| KQ0T.r4 / r5 | 8,050 / 3,233 | – |
| KQ2.r5 / r6 / r7 / r8 | 42,565 / 11,796 / 727 / 163 | – |
| KQB.r7 / r8 / r9 / r10 | 21,590 / 20,309 / 16,454 / 3,681 | – |
| KQS.r2 / r3 / r4 / r5 | 7,980 / 9,168 / 2,916 / 138 | – |
| KQK1.r7 / r8 / r9 / r10 / r11 / r12 | 3,747 / 2,211 / 1,901 / 1,370 / 782 / 166 | – |
| KQF.r2 / r3 | 8,698 / 0 | 2,375 / 661 |
| KQF2.r3 / r4 | 2,069 / 0 | 872 / 543 |
| Ca.r8 / r9 / r10 (check_tw_chart) | 2,650 / 3,427 / 141 | – |

That is 175,932 '7' leaves (33 of them on the G chart: KQ2.r5 3, KQK1.r7–r9 30) and 4,451 '8' leaves.
`verify/recert_tw_c3.out` agrees: 6,218 Ca leaves, 0 problems.

**Sampled regular, deep and winding leaves (g9, parity-aware, final engines), 0 mismatches.** Up to 400 regular, every
deep leaf, and up to 400 winding leaves per file (all where fewer), on all 24 new KQ/KQF files. The winding leaves
need, and get, the nb_cert2 enclosure the runs used.

**g1 on the 24 KQ/KQF files.** It replays each file, and for each file:
- re-certifies 400 winding leaves, 100 count leaves, 100 tau-winding leaves and 100 disc-winding leaves;
- runs 100-digit spot checks (15 leaves × 2 points per kind);
- finds 0 winding, TW or disc problems and 0 spot violations.

The disc-winding outer counts at 100 digits are at most 6 (KQF.r2, KQF.r3, KQF2.r3, KQF2.r4).

g1's strict criterion (today's first-success total ≤ the recorded bound) flags 105 count leaves. All re-certify with a
total ≤ 11, which is the requirement for a whole-circle leaf. The flags come from option order: g1 initialises kq_run
without the run's `--peri`, so it tries the plain count first. **i13** checked 20 of them, including the two largest
gaps (6 → 8 and 7 → 8):
- for every box, the smallest certified total over all of today's options equals the recorded bound;
- the 60-digit real-root counts at 9 points never exceed it.

As in AUDIT_12 i11, this is not a soundness issue.

**g3, 100-digit spot checks with my own chart maps.** 60 leaves × 3 points on each of the 24 files: 0 violations, no
odd root count. The largest count is 10 (whole circle, including TW leaves of KQF files); KQF/KQF2 outer counts reach
at most 6.

**g4 on Ca.r8–r10.** Each file replays exactly:
- 90 points at 100 digits, 0 violations;
- 150 TW leaves re-certified (all 141 in r10), with the zeros inside the contour checked at 100 digits (N ≤ 10);
- 150 regular and pruned leaves re-certified with c3_run, 0 problems.

**j1, my own 60-digit zeros inside the stored TW contours plus the arc enclosure test.**
- 15 leaves of each of KQ2.r5, KQB.r8, KQK1.r7, KQS.r3, KQ0T.r4, KQF.r2 and Ca.r8, plus all 33 G-chart leaves: 0
  problems.
- Each leaf has 1 or 2 zeros inside its contour (clusters are enclosed together under the new contour rule), constant
  over the box.
- The engine's polynomial agrees with mine to float-evaluation level (≤ 8.4e-12).
- 8,280 enclosure checks, largest |F − b₀|/R = 0.996.

The '8' leaves of KQF.r2–r3 and KQF2.r3–r4 (files unchanged since 17:39) were checked independently in AUDIT_14 (j4).

## 3. The authors' verification outputs (verify/): nothing unexpected

- **recert_\*.out:** every run file, 300 leaves (fewer where a file has fewer), problems 0, replay problems 0. This
  includes KQB.r9 and r10 and every KQF/KQF2/KQS level.
- **mp_\*.out** (verify_mp, 60 digits): violations 0 everywhere, including Ca r0–r10 (Ca r5 and r6 finished after my
  first look).
- **disc_mp_\*.out:** KQF.r2 2,375, KQF.r3 661, KQF2.r3 872 and KQF2.r4 543 disc leaves; 400 points each; largest
  outer count 6; 0 violations. Earlier levels have no '8' leaves.
- **recert_tw_c3.out:** 6,218 of 6,218 Ca tau-winding leaves, 0 problems.

## 4. Documentation review

I read:
- C3_CORE.md §8 (with the new contour-choice paragraph) and §9;
- PROOF.md (the TW and disc-winding lemmas);
- paper/paper.md (abstract, Lemmas 3.4 and 3.5, atlas, face paragraph, §6, §7);
- FLOAT_RIGOUR.md §5a;
- README.md;
- code/launch_runs.sh (continuation options, lines 84–95);
- code/final_verify.sh (recert_tw_c3.py, disc_mp_check.py).

The lemmas, the contour enclosure and its float analysis are stated correctly:
- §5a's argument-step margin π − 2·10⁻⁶ follows from the 10⁻¹² slack in R_j.
- C3_CORE §9 matches AUDIT_14.
- The round-13/14 hardenings are in place: check_all requires which = 0 on C3 contours; `assert DISC_W == 6` is in
  check_all and bb_replay; TW_CHARTS = (Ca, FE, FM) with the positive-scaling comment.

**R15-1 (major, documentation: paper.md).**
- *Problem.* The paper never describes the corner-K repair of the face construction. §4 and the face paragraph in §5
  say that "two charts certify that Q has no zero for 3s ≤ |τ| ≤ 0.06" and that "engine KQ counts |τ| ≥ 0.06". That
  is false inside K = {P_z ≥ 0.875, 0.125 ≤ Q_z ≤ 0.25}: FIP has genuine gap zeros there from σ ≈ 0.0153 (AUDIT_11
  h22).
- *What the proof actually does* (PROOF.md table rows 147–149, C3_CORE §6): inside K_lw it uses the face construction
  only for s ≤ SF2 = 0.010005, with the gap up to |τ| = 0.030015 and the outer count from there (KQF2). The shell
  SF2 < s ≤ 0.02 is covered by whole-circle counts (KQS).
- *Fix.* Add this to the face paragraph and to the atlas list of §4. Without it the paper's own description of the
  proof is incorrect.

**R15-2 (minor, documentation).** The TW lemma (PROOF.md, paper Lemma 3.4, C3_CORE §8 "Lemma") is stated for Q only.
The runs also use:
- the reverse G(t_p) = t_p¹² Q(−1/t_p) (33 leaves; C3_CORE §8 mentions KBPe only under "Contour choice");
- the C3 chart polynomial G(y) = Q(νy)/ν¹⁰ (6,218 Ca leaves; TW_CHARTS also allows FE and FM, G(y) = Q(sy)/s⁶).

The chart case needs one sentence: ν, s > 0 at every physical point (ε = ν²E > 0), so a non-real zero y gives the
non-real zero τ = νy (or sy) of Q. Add it, and name recert_tw_c3.py as its re-certification.

**R15-3 (minor, documentation).**
- paper §7 says "Re-certification of leaves, including every winding leaf, with the final code". That is not what was
  done. bb_replay re-certifies 300 sampled leaves per file. Every winding leaf was re-certified only for KQ0, KQ1(+r1),
  KQK2(+r1..r6) and the smaller families (AUDIT_12); bulk had 20,000 of 1.8M.
- What is complete is the re-certification of every tau-winding and disc-winding leaf (AUDIT_15, recert_tw_c3).
- Also update "Twelve rounds" to fifteen, "the final pass (AUDIT_12)" to AUDIT_15, and PROOF.md §6 "AUDIT_1 –
  AUDIT_12".

**R15-4 (minor, documentation).** C3_CORE §9 (and the paper's face paragraph) cite count_outer as the region the '8'
leaves replace. KQF2, and KQF from r1, use `--peri`, i.e. count_peri with t_lo = T_X. Its region is the same
O_{T_X} = {|τ| ≥ T_X} ∪ {f = 0} (AUDIT_10 h16, AUDIT_14 a). Name both.

**R15-5 (minor, documentation).** The README status says that sampled leaves of every run file "pass 60-digit spot
checks (verify/)". final_verify.sh runs verify_mp on every family except KQF and KQF2: their outer-count leaves have
no project 60-digit check, and only their '8' leaves are checked (disc_mp_check). The audit's 100-digit outer counts
(g1, g3: AUDIT_12 and 15) cover them. Either say so, or add an outer-count mode to verify_mp.

## 5. What "verified" means

- **Certificate tree:** every split tree of every family is replayed exactly by two independent replays (check_all
  and my g2), and coverage is complete.
- **Re-certified in full with the final engines:**
  - every leaf of the 22 smaller families (AUDIT_12);
  - every winding leaf of KQ0, KQ1 and KQK2 (AUDIT_12);
  - every tau-winding and disc-winding leaf (AUDIT_15).
- **Sampled:** the other leaves (bulk, the large KQ families, the C3 charts), with 0 mismatches in every sample.
- **Independent of the engines:** the 100-digit spot checks with my own formulas. They are sampled and find no
  violation.
- **Engines:** they are shared between the runs and the re-certification. Their soundness rests on the code review and
  the independent tests of AUDIT_1–14 (including the enclosure tests of the TW and disc contours), not on formal
  verification.

With R15-1 fixed in the paper text, I see no reason to keep the "not yet complete" qualifiers.

## R15 fixes verified (docs, read at about 19:05)

| # | where | status |
|---|---|---|
| R15-1 | paper §4 atlas (K_lω = {l ≥ 0.88, ω ≥ 1.29}: face construction only for s ≤ 0.010005, KQ from \|τ\| ≥ 0.030015, the shell 0.010005 < s ≤ 0.02 by whole-circle KQ counts); paper §5 "The corner K" bullet | **closed** |
| R15-2 | paper Lemma 3.4, PROOF.md tau-winding lemma, C3_CORE §8 (counting lemma and certificate) | **closed** |
| R15-3 | paper §7, README (fifteen rounds, AUDIT_15), PROOF.md §6 (AUDIT_1 – AUDIT_15, rounds 12–15 listed) | **closed** |
| R15-4 | C3_CORE §9 | **closed** |
| R15-5 | status blocks of README, PROOF.md and paper | **closed** |

**R15-1.** The §5 bullet gives:
- the mechanism (the outer near zero pushed beyond |y| = 3; genuine gap zeros from σ ≈ 0.015);
- K = {P_z ≥ 0.875, 0.125 ≤ Q_z ≤ 0.25} ⊂ K_lω;
- the s ≤ SF2 construction (empty gap, near zeros |y| ≤ 2.79, outer 6);
- the KQS shell.

This matches PROOF.md rows 147–149, C3_CORE §6 and AUDIT_11. A wording nit, optional: the preceding bullet's
"3s ≤ |τ| ≤ 0.06" could say "outside the corner K (next bullet)". The paragraph is correct as a whole, because the
next bullet states the exception.

**R15-2.** All three places now cover:
- the reverse G (−1/t_p is non-real when t_p is);
- the chart polynomials Q(νy)/ν¹⁰ (Ca) and Q(sy)/s⁶ (FE, FM), with ν, s > 0 at physical points.

C3_CORE §8 and PROOF.md name recert_tw_c3.py / check_tw_chart, and the which = 0 rule that check_all enforces.

**R15-3.** The re-certification claims now match what was done: samples of the Rolle, deep and winding leaves of every
file; every tau-winding leaf (175,932, of which 6,218 on C3 charts); every disc-winding leaf (4,451).

**R15-4.** C3_CORE §9 now says count_peri (the --peri variant, §7) covers the same closed region
{|τ| ≥ T_X} ∪ {f = 0}.

**R15-5.** The status blocks now say exactly what the 60-digit project checks cover:
- sampled leaves of every whole-circle family and C3 chart;
- the disc-winding leaves of KQF and KQF2.

The other KQF/KQF2 outer-count leaves and the D2 charts rely on the audit's 100-digit checks. PROOF.md §6 notes that
verify_mp has no D2 map and points to the auditor's D2 spot checks. The final ones are in
audit/round7/final12/chainC.out (AUDIT_12): 180 points per D2 chart, 0 violations, max N = 8.

**Qualifiers.**
- PROOF.md: "Theorem (computer-assisted)" … "Answer: **yes**."
- README: "Theorem (computer-assisted)" … "the answer to the second question is **yes**".
- paper abstract: "We give a computer-assisted proof that …".
- All three status blocks read "complete (2026-10-23 18:30)".

These are consistent with this audit.

**No open issue remains from AUDIT_15.**
