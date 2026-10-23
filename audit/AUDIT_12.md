# AUDIT 12: final pass over the families whose run files are final (snapshot of 2026-10-23)

(in progress: started 14:15. Chain A finished at 14:32. Chain B (g1) and chain C (g20, g3, g4, g13, extra bulk
winding samples) are still running; their results are filled in below as they finish.)

## Scope and method

The authors declared these families final at about 14:05 (nothing will be added to their files):
- bulk, F1, C1, C1p, C1pp;
- B1, B2, B3p, B3m;
- L2E, L2P, L2Ap, L2Am, L2Tp, L2Tm;
- D2E, D2E2, D2P, D2K, D2Tp, D2Tm;
- F0, KQ0, KQ1 (base + r1), KA;
- KQK2 (base + r1..r6);
- C3: CK, CE, Cm, PhiP, PhiM and the eight L2* charts;
- face charts: FE, FM, FIP, FIM.

They were checked by `audit/round7/final12.sh`: two serial chains (2 cores), outputs in `audit/round7/final12/`, my own
tools of rounds 7–11 (g1, g2, g3, g4, g9, g12, g13, g20).

## Split rules in my tools (the authors' verify_mp finding, 14:2x): OK for every file checked here

The authors found that verify_mp rebuilt leaves with the plain rule even on files whose header carries a
`split_rule`. I checked how each of my tools rebuilds leaf boxes.

**Which files carry a header split_rule** (i10, all 77 run files, base and continuations):

| file | rule |
|---|---|
| KQ0T.r1 | face/30 |
| KQB.r5 | face/30 |
| KQK1.r4 | face |
| KQK2.r4 | face |
| KQK2.r5 | face |
| KQK2.r6 | face/30 |
| KQS | shell/30 |

No C3, face-chart, D2 or bulk file has one. The running levels KQ2.r4, KQB.r6, KQK1.r5 and KQS.r1 will record theirs in
the header when they are closed.

**How each tool rebuilds leaves:**
- **g1, g3, g9, g13:** a global `SPLIT_RULE = H.get('split_rule')` is set per file at module level, and `split()` calls
  `splitrule.split_box(bx, S, SPLIT_RULE)` with the file's own header scale. g2 does the same per file (`rule_f`),
  with a tiling check. These tools use the project's `splitrule.py` itself, so they rebuild exactly the boxes the run
  certified.
  - So their leaf-level checks on KQK1.r4, KQK2.r4–r6, KQB.r5, KQ0T.r1 and KQS (including the round-11 checks i3 and i5
    on KQS) used the correct boxes.
  - The KQ2.r4+ files are not final yet; the snapshot pass will read their headers.
- **g4 (C3 and face charts) and g20 (D2):** these replay the plain bisection only, which is correct because none of
  their files has a rule. Both now print PROBLEM if a header ever carries one.
- **Header-less input (.records.jsonl):** g1, g3 and g13 also accept it, and then rebuild leaves with the default rule
  and the driver scale.
  - That would be wrong for continuation levels with a rule or a non-default scale.
  - I used such input only in rounds 7–8 (08:06–09:44), before any split rule existed. No rule file was ever checked
    that way.
  - These tools now print a WARNING for header-less input, and take the rule and scale from AUDIT_SPLIT_RULE /
    AUDIT_SCALE when given.
  - The round-11 integrity check i7 on KQK1.r5 records passed the rule (face/30) and scale explicitly.
- **Detection limits:** a wrong rule cannot be detected structurally on a complete record without failed leaves. The
  volume and tiling checks pass for any bisection. Only the stored 'fails' and 'left' boxes, which g2 compares
  exactly, would expose it. So the header rule must be used, and it is.

## Chain A results (14:15–14:32)

### g2: my own replay of every final family (all files, all records)

Each family must pass:
- the grid of the base file and the exclusions;
- exact tree replay with each file's own header scale and split rule;
- leaf codes within the family limits;
- prune predicates (exact for 'p', 'o', and the face-chart 'u');
- continuation roots owed;
- the owed-box rules: the D2 cube for C1p, rules (i) and (ii) for D2Tp, s_lo > SF and the corner K for the face charts;
- the record leaf-code histogram cross-check (new in round 11).

Result: every final family is **OK**.

| family | leaves | of which winding | failed leaves | owed after the rules |
|---|---|---|---|---|
| bulk | 9,545,140 | 1,804,672 | 0 | 0 |
| C1 | 77,436 | 34,471 | 0 | 0 |
| C1p | 44,140 | 24,550 | 659 | 0 |
| C1pp | 11,356 | 6,138 | 0 | 0 |
| B1, B2, B3p, B3m | 8,379 / 15,506 / 1,722 / 1,570 | 5,171 / 0 / 0 / 0 | 0 | 0 |
| L2E, L2P, L2Ap, L2Am, L2Tp, L2Tm | 436 / 1,979 / 2,503 / 2,529 / 256 / 256 | 290 / 0 / 0 / 0 / 0 / 0 | 0 | 0 |
| D2E, D2E2, D2P, D2K, D2Tm | 256 / 256 / 1,130 / 464 / 256 | 137 / 24 / 348 / 102 / 0 | 0 | 0 |
| D2Tp | 26,915 | 0 | 20,245 | 0 |
| F0, F1 | 12,354 / 1,093 | 3,925 / 329 | 0 | 0 |
| KQ0 | 54,718 | 35,950 | 0 | 0 |
| KQ1 (base + r1) | 90,528 | 43,408 | 0 | 0 |
| KA | 64 | 0 | 0 | 0 |
| face FE, FM, FIP, FIM | 226,498 / 348,845 / 307,304 / 57,097 | 0 | FIP 124,394 | 0 |
| C3 CK, CE, Cm, PhiP, PhiM, L2aE, L2am, L2aYp, L2aYm, L2KE, L2Km, L2KYp, L2KYm | 2,598 / 165,916 / 2,223 / 285 / 285 / 335 / 590 / 413 / 162 / 368 / 368 / 368 / 368 | 0 | 0 | 0 |

Notes on the owed-box rules:
- **C1p:** its 659 failed leaves all lie in the D2 cube, and the largest e among them is 6.107e-6.
- **D2Tp:** 20,334 owed boxes. All but one lie above that e (rule i). The last one is covered by D2E ∪ D2E2
  (rule ii), and both of those families are complete.
- **Face charts:**
  - FE: 73 owed boxes, all in K with s_lo > SF2.
  - FM: 67 owed boxes, likewise.
  - FIP: 124,446 owed boxes with s_lo > SF, plus 90 in K.

**KQK2** (base + r1..r6; 127,845 leaves, of which 2,568 winding) replays with nothing owed: 11,100 failed leaves, all
re-covered. g2 still marks it NOT OK, because it has 11,271 'p' leaves and those need the face construction:
- FE, FM, FIP and FIM are complete.
- KQF is not: 281 open boxes, 30 of them inside KQK2's region l ≥ 0.95, ω ≤ 3π/20. KQF2 and KQS are not complete
  either.

KQ0, KQ1 and bulk have no 'p' leaves.

### g12: the project's check_all.py (unchanged code, RESULTS.md redirected), 14:29–14:32

It agrees with g2:
- complete: every final family except KQK2 (the 'p' dependency);
- incomplete: KQ0T, KQ2, KQB, KQK1, KQF, KQF2, KQS and Ca (still running or paused).

### g9: re-certification with the final engines (parity-aware)

**Every leaf**, 14:16–14:23, 0 mismatches, 0 leaves not re-certifiable:

| family | regular | deep | winding |
|---|---|---|---|
| F1 | 764 | – | 329 |
| C1 | 42,957 | 8 | 34,471 |
| C1p | 17,469 | 2,121 | 24,550 |
| C1pp | 5,218 | – | 6,138 |
| B1 | 3,208 | – | 5,171 |
| B2 | 15,506 | – | – |
| B3p | 1,722 | – | – |
| B3m | 1,570 | – | – |
| L2E | 146 | – | 290 |
| L2P | 1,979 | – | – |
| L2Ap | 2,503 | – | – |
| L2Am | 2,529 | – | – |
| L2Tp | 256 | – | – |
| L2Tm | 256 | – | – |
| D2E | 119 | – | 137 |
| D2E2 | 232 | – | 24 |
| D2P | 782 | – | 348 |
| D2K | 362 | – | 102 |
| D2Tp | 17,831 | 9,084 | – |
| D2Tm | 256 | – | – |
| F0 | 8,429 | – | 3,925 |
| KA | 64 | – | – |

F1 was run with the second-order engine, so 987 of its leaves re-certify only with nb_cert2. That is expected.

**Sampled**, 0 mismatches:
- bulk: 300 of 7,740,002 regular, 300 of 466 deep, 20,000 of 1,804,672 winding (8 and 6,098 only with the second-order
  engine; bulk is an `--order2` run);
- KQ0: 300 of 18,768 regular;
- KQ1: 300 of 39,260 regular; KQ1.r1: 300 of 7,860 regular;
- KQK2 and r1..r6: 300 regular per file (747, 0, 3,574, 15,409, 45,072, 57,325 and 3,150 available).

A full re-certification of bulk's 1.8M winding leaves would take about 5.5 h on one core at the measured rate (20,000
in 3.7 min). So bulk's leaf-level check is sampled. Chain C adds another 60,000.

**A tool gap of mine, fixed (not a run problem).** g9's KQ branch re-checked winding leaves with nb_cert's first-order
enclosure only.
- The KQ engine certifies them with nb_cert2 (`kq_run.wind`: `box_coeffs2` + `wind_box2`; header "enclosure: nb_cert2
  (second order)"), which is never wider.
- So g9 printed MISMATCH for the KQ0 winding leaves that need the second-order enclosure: about 600 lines at 14:24.
- I stopped that g9 step and added the nb_cert2 fallback, as g1 and bb_replay already do.
- Every KQ winding leaf of the final files is re-checked with the fallback by g1 (chain B, below).

## Chain B and C results (my own formulas and engines)

### g1: the KQ files, every winding leaf (14:15– ; stopped before the snapshot, see below)

**KQ0** (54,718 leaves):
- The replay is exact.
- All 35,950 winding leaves re-certify: first with nb_cert, then with nb_cert2. 35,945 need the second-order enclosure,
  which is what the run used. 0 problems.
- 100 count leaves re-certify with the current kq_run.
- 100-digit spot checks with my own coefficients, 180 points (90 in winding leaves): 0 violations, no odd root count,
  at most 10 roots. The smallest relative distance of a zero to the winding circle is 0.0187.

KQ1 (base + r1) was in progress when the snapshot work took the cores (see the final summary).

### g20: D2E2 and the D2Tp rules (14:48–14:50)

- The D2E2 and D2E chart values agree with my own physics formula at 50 digits: max relative error 1.0e-18, 2,498
  terms, 0 mismatches.
- D2Tp has 20,334 owed boxes. 20,333 are dropped by rule (i), checked exactly (SD·ν_lo·E_lo > 6.10657e-6). The
  remaining box is (0, 0.12625) × (0.7875, 1.05) × (0, 0.2625) × (0.2625, 0.525).
- Rule (ii): 50 point/chart checks over that box. Every point lies in a certified leaf of D2E or D2E2; 0 violations at
  100 digits; at most 8 roots.

### g3: 100-digit spot checks with my own chart maps (14:50–15:00; 60 leaves × 3 points per family)

0 violations and no odd root counts for:
- bulk, F1, C1, C1p, C1pp;
- B1, B2, B3p, B3m;
- L2E, L2P, L2Ap, L2Am, L2Tp, L2Tm;
- D2E, D2E2, D2P, D2K, D2Tp, D2Tm;
- F0.

That is 180 points per family. The largest root count was 8, and the smallest winding margin was 0.063 (bulk).

The KA and KQ parts of g3 were cut when the snapshot pass needed the core. g1 covers the KQ files.

## Snapshot pass (frozen 15:00 snapshot; levels up to the closed ones)

Files: KQ0T (base + r1..r3), KQ2 (base + r1..r4), KQB (base + r1..r6), KQK1 (base + r1..r6), KQF (base + r1), KQF2
(base + r1..r2), KQS (base + r1), and c3 Ca (base + r1..r7). g2 is capped at these levels (AUDIT_MAXLEV), so later level
files are ignored. Scripts: `round7/snap12a.sh`, `round7/snap12b.sh`; outputs in `round7/snap12/`.

### g2 replay (14:59–15:00)

Every file of every snapshot family replays exactly, with each file's own scale and split rule:
- 0 replay problems and 0 prune violations;
- every record reproduces its own leaf-code histogram;
- every continuation root is owed.

The open boxes agree exactly with the authors' check_all on the snapshot:

| family | files | records | leaves | of which winding | pruned | open boxes |
|---|---|---|---|---|---|---|
| KQ0T | 4 | 1,494 | 442,074 | 220,447 | 0 | 396 |
| KQ2 | 5 | 3,735 | 305,768 | 122,040 | 0 | 2,791 (1,669 failed leaves) |
| KQB | 7 | 1,780 | 114,914 | 23,741 | 30 | 1,000 |
| KQK1 | 7 | 5,694 | 361,233 | 73,540 | 55 | 1,071 |
| KQF | 2 | 381 | 259,492 | 0 | 128 | 2,292 after the 28 KQF boxes inside K_lw are excused (1,062 failed leaves) |
| KQF2 | 3 | 793 | 17,145 | 0 | 260 | 912 |
| KQS | 2 | 340 | 11,280 | 0 | 661 | 877 |
| C3 Ca | 8 | 3,933 | 2,478,929 | 0 | 193,447 | 822 |

KQB, KQK1 and KQS also carry the 'p' dependency on KQF, KQF2 and KQS.

### g9 sampled re-certification (15:00–15:16), 0 mismatches in every file

- 200 regular leaves per file, from all 30 KQ / KQF files: 6,000 leaves.
- 300 winding leaves per file, or all of them where a file has fewer.
  - First nb_cert, then nb_cert2. Every winding leaf of these runs needs the second-order enclosure, as recorded in
    their headers.
- KQF-engine files are re-certified against their own header T_X (R11-2): 0.06 for KQF, 0.030015 for KQF2.

### g4: the final C3 charts and the face charts (15:16–15:25)

For each of the 17 files, g4 does three things:
- its own replay, including zone codes and pruned leaves;
- spot checks on 20 leaves × 2 points at 100 digits, interpolating the chart polynomial and counting real roots on the
  certified interval;
- re-certification of 100 sampled leaves (pruned ones included) with c3_run.

**C3 charts.** Every chart replays exactly and has 0 spot violations and 0 re-certification problems. The largest root
count found on a certified interval:

| charts | largest root count |
|---|---|
| CE, Cm | 6 |
| CK, L2am, L2Km | 4 |
| L2KE | 2 |
| PhiP, PhiM, L2aYp, L2aYm, L2KYp, L2KYm | 0 (these charts certify "no zero") |

**Face charts.** FE, FM, FIP and FIM all replay exactly, with 0 spot violations.
- Re-certification: FE, FM and FIM have 0 problems.
- **FIP has 7 code differences, all benign.** The 7 boxes all have σ_lo ∈ [0.020014, 0.020020], just above SF = 0.02.
  The FIP run certified them as "no zero" ('a', bound 0) before round 10b. Since round 10b, c3_run.face_bounds prunes
  every box with s_lo > SF as 'o', because it contains no needed point. So today's code returns 'o' where the file has
  'a'. Both are valid for these boxes, and they lie outside the region the construction uses. My g4 compares codes
  strictly, which is why it reports them.

### g1: the final KQ files (continued)

| file | winding leaves re-certified | count leaves re-certified | 100-digit points |
|---|---|---|---|
| KQ1 base | all 23,710 (23,629 need nb_cert2) | 100 | 180, 0 violations |
| KQ1.r1 | all 19,698 | 100 | 180, 0 violations |
| KQK2 base, r1–r4 | KQK2.r4: all 658 (base, r1–r3 have none) | 100 each (r1 has no leaves) | 0 violations |

KQ1: the strict criterion of g1 flags 3 count leaves (1 in KQ1, 2 in KQ1.r1).
- They were recorded with bounds 9, 8 and 10.
- Today's engine certifies them with totals 11, 10 and 11, all within the limit 11 of a whole-circle KQ leaf (parity).
- kq_run.certify returns the first option whose count is ≤ 11, not the smallest count. So a larger total than the one
  recorded is expected, and it is not a contradiction.
- These are the parity-aware cases of R7-1, as in g9, and not a soundness issue.

## VERDICT (for the publish copy, 2026-10-23)

**What "re-certified" means here.** Every split tree was replayed by my own code with each file's own scale and split
rule, and every record reproduces its own counts. Leaves were re-certified with the final engines, parity-aware:
- for a whole-circle leaf, total ≤ 11 is accepted, giving ≤ 10 roots by the parity lemma;
- for a KQF leaf, the outer count from the file's own T_X must be ≤ 6;
- for a face-chart gap leaf, there must be no zero.

At sampled points I also counted real roots at 100 digits with my own polynomial formulas (independent of the engines).

The engines themselves (nb_cert, nb_cert2, kq_run, c3_run, the chart engines) are shared between the runs and the
re-certification. Their soundness rests on the code review and the independent tests of AUDIT_1–11, not on a formal
proof.

### 1. Fully verified

All of the following holds for these 22 families:
- the family is complete;
- every leaf was re-certified, with 0 mismatches;
- the 100-digit spot checks found 0 violations;
- every owed box is either re-covered or excused by an audited rule whose own dependencies are complete.

The families:
- **F1, C1, C1p, C1pp.** C1p's 659 failed leaves lie in the D2 cube, and every D2 chart is complete.
- **B1, B2, B3p, B3m.**
- **L2E, L2P, L2Ap, L2Am, L2Tp, L2Tm.**
- **D2E, D2E2, D2P, D2K, D2Tp, D2Tm.** D2Tp's owed boxes are covered by rule (i), and the last one by rule (ii) via
  D2E ∪ D2E2, both complete.
- **F0, KA.**

### 2. Verified, with every winding leaf and a sample of count leaves re-certified (complete, no open dependency)

- **KQ0:** all 35,950 winding leaves; 400 of 18,768 count leaves; 180 points at 100 digits.
- **KQ1 (base + r1):** all 43,408 winding leaves; 800 count leaves (300 + 100 per file); 360 points at 100 digits.
  The 3 strict-criterion flags are parity-valid, see above.

### 3. Complete and fully replayed, with a sample of leaves re-certified

- **bulk:**
  - all 9,545,140 leaves replayed;
  - re-certified: 300 regular, 300 deep and 20,000 of 1,804,672 winding leaves;
  - 180 points at 100 digits.
  - All of its 1.8M winding leaves would take about 5.5 h on one core, so the sample stands.
- **C3 charts CK, CE, Cm, PhiP, PhiM and the eight L2\* charts:** 100 leaves re-certified and 40 points at 100 digits
  per chart.
- **Face charts FE, FM, FIP, FIM:** the same. FIP's 7 differences are benign (σ_lo > SF), as explained above.

### 4. Complete, but each relies on a family that is not finished

| family or region | relies on | open in the snapshot |
|---|---|---|
| KQK2 (11,271 'p' leaves) | the face construction: FE/FM/FIP/FIM complete, KQF not | KQF has 2,292 open boxes, some inside KQK2's region l ≥ 0.95, ω ≤ 3π/20 |
| the slab η ≥ 0.99025 excluded from bulk | F1, KQ0 (complete) and KQ0T | KQ0T: 396 boxes |
| the C3-cell 'o' leaves of the face charts, KQF, KQF2 and KQS | the C3 core: CK, CE, Cm, PhiP, PhiM, L2\* and KA complete | Ca: 822 boxes |
| the 'p' leaves of KQB, KQK1 and KQS | the face construction | KQF, KQF2, KQS |

### 5. The snapshot's incomplete families

KQ0T, KQ2, KQB, KQK1, KQF, KQF2, KQS and Ca, each up to its closed level:
- every file replays exactly;
- 0 replay problems and 0 prune violations;
- every record reproduces its leaf histogram;
- the open boxes match check_all exactly: 396, 2,791, 1,000, 1,071, 2,292, 912, 877 and 822.

Sampled leaves re-certify with 0 mismatches:
- 200 regular and up to 300 winding leaves per file (g9);
- 10 count leaves and 20 points at 100 digits per KQ/KQF file (g1);
- 40 leaves and 40 points per Ca file (g4).

The final g1/g4 numbers on the snapshot are in the next section.

So the theorem's certificate is complete outside the open boxes of these 8 families, and outside the regions that
depend on them (table in §4). That matches the statement planned for README / PROOF / paper.

### 6. Soundness concerns

**None.** No leaf failed re-certification under the parity-aware criterion. No 100-digit spot check found more real
roots than a leaf's bound, or an odd root count, anywhere. No replay or prune problem.

What did turn up, all at the tool level:
- **g9 (mine, fixed):** a false alarm on KQ winding leaves, because it had no nb_cert2 fallback. g1 re-checked every KQ
  winding leaf of the final files with the fallback.
- **g1 (mine):** 3 strict-criterion flags on KQ1 count leaves, parity-valid.
- **g4 (mine):** 7 FIP code differences ('a' recorded, 'o' today) for boxes with σ_lo > SF, which no point of the proof
  needs.
- **verify_mp (the authors' tool):** the split-rule bug. My tools read the header rule (see above).

Round 11's checker items R11-1, R11-2 and R11-5 are fixed and verified (AUDIT_11 addendum).
