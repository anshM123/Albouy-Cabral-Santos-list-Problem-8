# AUDIT 11: the corner-K repair of the face construction (SF2, KQF2, KQS)

Independent audit, 2026-10-23, 13:05–13:40. Scripts and outputs are in `audit/round11/`. At most 2 Python processes
ran at once, all of them light. The audited code is `kq_run.py` (13:01), `check_all.py` (12:52), `next_level.py`,
`splitrule.py`, `bb_replay.py` (11:15), `launch_runs.sh` / `final_verify.sh` (12:57), and C3_CORE.md §6 and PROOF.md §3.

## Verdict

The repair is sound as a construction:
- The case split (s ≤ SF2: near + gap + KQF2; SF2 < s ≤ SF: KQS; outside K_lw: unchanged) covers every point that the
  'p' leaves hand to the face construction.
- The seams at s = SF2, at τ = T_X2 and at the edges of K and K_lw are closed, with exact margins.
- The numbers support it: no zero enters the gap for s ≤ SF2 anywhere in K outside the C3 cell.

Two checking-tool weaknesses (R11-1, R11-2) should be fixed before the final pass. Both could let a wrongly launched
KQS or KQF2 file pass. No current file is wrong.

There is also an urgent operational problem (R11-3): two continuation chains each for KQK1 and KQF2 are running at the
same time on the same output files.

| # | item | status |
|---|---|---|
| R11-1 | check_all takes the meaning of 'p' (FACE_S) and of the outer count (T_X) from each run's header, and does not pin them per family | **ISSUE (major, checker)** |
| R11-2 | bb_replay re-certifies KQF-engine files with the default T_X = 0.06, which is weaker than KQF2's claim | **ISSUE (major, checker)** |
| R11-3 | duplicate chains: KQK1.r5 and KQF2.r1 each written by two concurrent kq_run processes; corrupted records files | **ISSUE (major, operational; not soundness)** |
| R11-4 | C3_CORE §6: the stated mechanism (a B zero at τ₀ ≈ 0.053, split by μ) is wrong; the conclusion and the SF2 choice hold | ISSUE (minor, doc) |
| R11-5 | check_all: winding digits have no family limit (KQF limit 6); the KQ/KQF 'o' test checks l, ω only | ISSUE (minor, robustness) |
| R11-6 | PROOF.md / C3_CORE.md: the KQS/KQF e-range starts at fl(0.9996), so ε ≤ SF² is not literally covered; stale dependency sentences | ISSUE (minor, doc) |
| R11-7 | next_level does not drop face-chart boxes in K with s_lo > SF2 | note (moot now) |
| 1 | coverage logic and seams (s = SF2, τ = T_X2, edges of K and K_lw) | OK |
| 2 | KQF2 and KQS domains contain what is needed, with the 'p' tolerances | OK |
| 3 | near + gap + outer covers the circle | OK |
| 4 | numerics: gap empty for s ≤ SF2 in K; zeros beyond the near region stay far above τ = 0.03 | OK |
| 5 | code changes (kq_run flags, worker init, header, check_all rules, shell/K, next_level for KQF) | OK |
| 6 | run state: face charts complete under the new rules; KQF, KQF2, KQS replay cleanly so far | OK (runs in progress) |

## 1. Coverage logic (item 1): OK

**Pointwise case analysis.** Take p ∈ N with (l, ω) ∈ K_lw, outside the C3 cell.
- If s(p) ≤ SF2:
  - near: FE/FM at s = s(p);
  - gap: FIP/FIM at σ = |τ|/3 ∈ [s(p), SF2], with E = ε/σ² ≤ 1 and m = μ/σ² ≤ 1 ≤ R1;
  - outer: KQF2 on |τ| ≥ T_X2.
- If SF2 < s(p) ≤ SF: KQS. The C3 cell is pruned 'o' and is covered by the C3 core, since every run here has
  e ≥ 0.9996 and μ ≤ 4.2e-4, inside the core region e ≥ 0.99, μ ≤ 0.0505.

Outside K_lw, (P_z, Q_z) ∉ K. There no face-chart excuse applies, and every such point lies in a KQF box that is not
inside K_lw, so that box must be certified.

**K ⊂ K_lw, exact** (i1 §1):
- min l on K = √(0.875² + 0.125²) = 0.88388 ≥ fl(0.88);
- min ω on K = atan2(0.875, 0.25) = 1.29250 ≥ fl(1.29);
- Q_z ≥ 0.125 gives ω < π/2.

**The excuse predicates are "box ⊂ region".**
- face chart: P_z,lo ≥ 0.875, Q_z,lo ≥ 0.125, Q_z,hi ≤ 0.25, s_lo > SF2;
- KQF: l_lo ≥ 0.88, ω_lo ≥ 1.29.

So a straddling box is never excused. The rules are applied to the final owed set of each family, after the
continuation chaining, so the order of the levels does not matter.

**Seam s = SF2.**
- In the face base files every owed box has s_lo ≥ 0.01001 > SF2 (i1 §5). So FE, FM, FIP and FIM are certified on
  s ≤ 0.01001 for all (P_z, Q_z), and the K excuse removes exactly the four incomplete root records of FE, FM and FIP
  (s ∈ [0.01001, 0.02002], inside K).
- KQS prunes 'p' only when ε_hi, μ_hi ≤ fl(SF2²(1 − 1e-9)) = 1.0010002489989998e-4 < SF2² exactly. So a point with
  s(p) ∈ (SF2·√(1 − 1e-9), SF2] is certified twice, and no point is certified by neither.

**Seam τ = T_X2.**
- The KQF2 header has T_X = fl(0.030015), which equals 3·fl(0.010005) exactly (i1 §4).
- KQF2 counts on the closed intervals [T_X, tan(φ₀/2) + 1e-9] and their negatives (count_peri with t_lo = T_X; the
  four parts overlap, AUDIT_10 h16).
- FIP/FIM are certified up to σ = 0.01001, i.e. τ = 0.03003 > T_X2.

**Edges of K_lw.**
- KQF2 and KQS cover l ∈ [fl(0.88), 1] and ω ∈ [fl(1.29), fl(π/2)], the same floats as the KQF excuse.
- fl(π/2) as the upper ω bound is the convention of every run family and is unchanged here.

## 2. Domains (item 2): OK

Exact, from the float values in the headers (i1 §2–3).

| need | provided | margin |
|---|---|---|
| KQF2: ε ≤ SF2² = 1.00100025e-4 | 1 − fl(0.99989989) = 1.0011e-4 | 9.98e-9 |
| KQF2: μ ≤ SF2² | √((1 − t_lo)/t_lo) = 1.004988e-4 | 3.99e-7 |
| KQF2 ⊇ KQS 'p' leaves | 'p' threshold 1.0010002489989998e-4 ≤ SF2² | exact |
| KQS ⊇ 'p' leaves of KQ0, KQ0T, KQ1, KQ2, KQK1, KQK2, KQB (ε_hi, μ_hi ≤ fl(SF²(1 − 1e-9)) = 3.9999999996e-4) | 1 − fl(0.9996) = 3.99999999999956e-4; μ ≤ 4.1231e-4 | 4e-14 (ε) |

1 − fl(0.9996) < SF² = 4e-4, so KQS (like the original KQF) does not literally contain ε ∈ (3.99999999999956e-4,
4e-4]. Those points are not needed, because only the 'p' leaves (with the 1e-9 tolerance) use the face construction
(R11-6, doc only).

## 3. Near + gap + outer (item 3): OK

- near: |τ| ≤ 3s(p) (FE/FM root variable y ∈ [−3, 3]);
- gap: 3s(p) ≤ |τ| ≤ 3·SF2 (FIP/FIM);
- outer: |τ| ≥ T_X2 (KQF2).

The outer region of KQF2 has two forms, both checked in AUDIT_10:
- count_peri: f-chart arcs from π ± (φ₀ − 1e-9), KBPe on |tp| ≤ tan(peri/2) + 1e-9, and KBQe on the two closed
  τ-intervals;
- count_outer: the away arc from φ₀ − 1e-9, and KBQe on [T_X, T_HI = tan(0.45)].

All seams are closed with overlaps, so a zero on a seam is counted at most twice. The sum stays an upper bound, and
5 + 0 + 6 ≤ 11 gives ≤ 10 by the parity lemma.

## 4. Numerics (item 4): OK

All of this uses my own P = μ²k² sin²f (1 − ζ²) − D²B², at 40–100 digits.

**h21 (a).** I reproduce the gap zero at the point σ = 0.01877, E = 0.25, P_z = 0.9863, Q_z = 0.1646:
- G(m = 0.995) = −5.3e-13 and G(m = 1) = +2.1e-13;
- at m = 0.9975 the root is at τ = 0.056267, just above 3s = 0.05624.

**h21 (b)–(d), random points.**
- 50 points in K_lw with s ≤ SF2: near ≤ 4, gap 0, outer ≤ 6, total ≤ 10.
- 25 points with SF2 < s ≤ 0.02: total ≤ 10.
- 50 points just outside K_lw with s ≤ 0.02 and the old T_X = 0.06: near ≤ 4, gap 0, outer ≤ 6.

**h22 (B), the gap threshold.** σ\* is the first σ at which P(τ = ±3σ) is ≥ 0 somewhere on the grid of (E, m) in
[0, 1] × [0, R1], up to σ = 0.04.

| (l, ω) | in K | σ\*_FIP | σ\*_FIM |
|---|---|---|---|
| l = 1, ω → 9π/20⁻ (worst point outside the C3 cell, c = 0.1564) | yes | 0.01525 | > 0.04 |
| l = 1, ω = 1.40 | yes | 0.01975 | > 0.04 |
| l = 1, ω = 1.37 | yes | 0.0298 | > 0.04 |
| l = 0.999, ω → 9π/20⁻ | yes | 0.01675 | > 0.04 |
| l ≤ 0.99 or ω ≤ 1.33 | | > 0.02 | > 0.04 |

So min σ\* over K outside the cell is 0.0152 > SF2 = 0.010005, a factor of 1.5. This is consistent with FIP being
complete for σ ≤ 0.01001.

**h22 (C), zero positions at the worst corner** (l = 1, ω = 9π/20 − 1e-6; also l = 0.99999, ω = 9π/20 − 0.003; and
ω = 1.38), for s ∈ {0.005, 0.008, SF2}, (E, m) ∈ {(0, 1), (1e-6, 1), (1, 1), (1, 0)}:
- The near pair's outer member has |y| ≤ 2.79, below 3.
- The gap is empty, and the outer count is 6.
- The nearest zero beyond the near region is at τ = 0.0787 at l = 1. That is the double zero of cos(f + ω); it splits
  to 0.0799/0.0808 at l = 0.99999. In every case it lies far above T_X2 = 0.030.

So item 4 holds, with a wide margin. There is no "corner pair" near τ = 0.03–0.05 at s ≤ 0.01 (see R11-4).

## 5. Code changes: OK (plus R11-1, R11-2, R11-5)

**kq_run.**
- `--faces=` sets FACE_S and T_X = 3·FACE_S; `--tx=` overrides T_X; `--c3prune` adds 'o' before the 'p' test in
  regular mode.
- `init()` and the Pool initargs pass FACE_S, T_X and C3PRUNE to the workers.
- The header records FACE_S, T_X, c3prune and minw.
- `--minw` only moves the failed-leaf threshold.
- A failed leaf is still '!' and stays owed.

**check_all.**
- FAMILIES include KQF2 and KQS.
- code_limits: KQF2 has engine KQF, limit 6, and rejects upper case. KQS has engine KQ, limit 11 (parity).
- The 'p' test uses the header FACE_S (see R11-1).
- KQ 'o' is accepted only with header c3prune, using the exact cell test on lower bounds.
- Both excuses as in §1.
- uses_p needs FACE_DEPS minus itself. There is no recursion: only KQS among FACE_DEPS has 'p' leaves.

**splitrule `shell/K`.** On 20,000 random KQS-type boxes, shell/30 and face/30 give identical splits wherever face/30
does not cut, and every shell split tiles its parent (i1 §6). 'face' and 'face/K' are unchanged by inspection: the
only change is the `not rule.startswith('shell/')` guard.

**next_level for KQF.** It drops exactly the owed boxes inside K_lw: 309 − 28 = 281 = KQF_r1_boxes.json (i1 §5).

## 6. Run state (13:35)

**Face charts.** My own replay (g2, i6) of FE, FM, FIP and FIM reports all four complete under the round-11 rules:
- FE: 73 owed, all in K with s_lo > SF2;
- FM: 67 owed, likewise;
- FIP: 124,446 owed with s_lo > SF, plus 90 in K;
- FIM: 0 owed;
- 0 prune violations.

**KQF2 base.** 24 records (21 complete), 75 open boxes, leaves {e: 1, f: 5, g: 906, o: 52}, 0 failed. It replays
cleanly (g1, g2). 8 leaves re-certify with T_X = 0.030015 (g1, g9). At 100 digits the outer count is ≤ 6 (8 points).

**KQS base.** 96 records (84 complete), 244 open boxes, leaves {i: 3, j: 8, k: 236, l: 1977, o: 103, p: 59}, 0
failed. It replays cleanly, with exact 'p' (radius SF2) and 'o' checks. 8 leaves re-certify (6 only with the current
kq_run). At 100 digits the total is ≤ 10 and even.

**KQF.** 309 owed, 28 of them inside K_lw. The r1 run (281 boxes) is in progress.

## Issues

**R11-1 (major, checker).**

*Problem.* check_all takes the meaning of the certificates from the header that each run writes:
- 'p' is tested against the header FACE_S;
- nothing checks the T_X of a KQF-engine file;
- nothing checks the domain of KQF2 or KQS against what the construction needs.

Suppose a KQS continuation were launched without `--faces=0.010005` (for example, a missing EXTRA in chain.sh). It
would prune 'p' up to s = 0.02 inside K_lw, exactly where the face construction is excused, and check_all would accept
it. A KQF2 level launched without `--faces` would certify outer counts from |τ| ≥ 0.06 only, which leaves the zeros
in 0.03 < |τ| < 0.06 uncounted.

*Current state.* All present files are right:
- the base headers (i0): KQF2 has FACE_S 0.010005 and T_X 0.030015; KQS has FACE_S 0.010005 and c3prune; KQF has
  FACE_S 0.02 and T_X 0.06;
- the running r1 command lines carry the flags.

*Fix.* Pin the values by family role in check_all, in every file of the family:
- KQS: FACE_S == SF2 (or use SF2 directly in prune_ok);
- KQF2: T_X == 3·SF2 and FACE_S == SF2;
- KQF: T_X ≤ 0.06;
- every other 'p' user: FACE_S == 0.02;
- the base domains of KQF2 and KQS: ε_max ≥ SF2² and ε_max ≥ fl(SF²(1 − 1e-9)) respectively, the μ counterparts, and
  l ≤ 0.88, ω ≤ 1.29 at the lower ends.

My g2 now does all of this.

**R11-2 (major, checker).**

*Problem.* bb_replay re-certifies KQF-engine files with `kq_run.init(H['scale'], True)`, which leaves FACE_S = 0.02
and T_X = 0.06. So every KQF2 leaf is re-checked against "outer count from |τ| ≥ 0.06", a strictly weaker claim.

*Demonstration (i4).* Take a small KQ box at the authors' example point. kq_run.certify_outer gives:
- T_X = 0.06: (True, 6);
- T_X = 0.030015: (False, 7).

final_verify.sh runs exactly this re-certification on KQF2.

*Fix.* `kq_run.init(H['scale'], True, faces=H.get('FACE_S'), tx=H.get('T_X'))`, together with the role assertion of
R11-1. For KQ files, `init(H['scale'], False, H.get('split_rule'), faces=H.get('FACE_S'),
c3prune=bool(H.get('c3prune')))` for consistency (no effect on certify).

My g1, g3 and g9 had the same hardcoded 0.06; they are fixed (below).

**R11-3 (major, operational; not soundness).** Two chain.sh processes run for each of two families. All four were
started by `sh chain.sh` and each has 8–10 worker processes.

| family | chain (sh PID) | started | kq_run PIDs | flags |
|---|---|---|---|---|
| KQK1 | 56308 | 12:57:52 | 24160 / 53460 | no --minw |
| KQK1 | 45716 | 13:02:06 | 3936 / 26348 | --minw=1e-7 |
| KQF2 | 65212 | 12:58:38 | 43716 / 65876 | no --minw |
| KQF2 | 35064 | 13:02:06 | 66948 / 70960 | --minw=1e-7 |

*Effects.* Both of each pair write runs/KQK1.r5.json or runs/KQF2.r1.json, plus the .records.jsonl and .log files:
- KQK1.r5.json.records.jsonl has 3 roots recorded twice and an interleaved fragment at line 1011.
- KQF2.r1.json.records.jsonl has a fragment at line 53 of 73.
- Whichever process finishes last overwrites the .json.
- Each chain then runs next_level and its own next level, so the race continues.

check_all would flag mismatched continuation roots, so this is fail-safe. But a `--resume` from these records files
may fail or lose records, and it doubles about 18 cores.

*Fix.* Stop one chain of each pair, including its kq_run and workers. Then restart that level cleanly, moving the
corrupted .records.jsonl aside first. I did not touch any process.

**R11-4 (minor, doc).** The mechanism paragraph in C3_CORE §6 says there is a B zero at φ₀ ≈ 2c²/(3P_zQ_z), i.e.
τ₀ ≈ 0.053, split by about ±0.047 or ±0.012. This is not right.

*What the zeros of B actually are* (h22 A, e = 1):

| (l, ω) | zeros of B near aphelion (τ) |
|---|---|
| l = 1, ω = 9π/20 | 0.0787, 0.1584 |
| P_z = 0.9863, Q_z = 0.1646 | 0.0836, 0.1662 |
| l = 0.97, ω = 1.40 | none |

The quadratic truncation drops the cubic term: at l = 1, B = 2cos(φ + ω)(sin ω − sin(φ + ω)) ≈ −φ(φ − c)(φ − 2c).

*The actual mechanism.* The gap zero is the outer member of the near pair 2m = c(E + 2y²). At the next order it is
pushed out by the factor (c − τ sin ω): at l = 1, P = cos²(f + ω)·[μ²k² sin²f − 4D²(sin(f + ω) + sin ω)²]. It reaches
y = 3 at E = 0, m = R1 when σ ≈ (c − R1/9)/(3 sin ω). That predicts 0.0151 at the cell corner; I measured 0.0153.

The conclusion is unaffected. The right justification of SF2 is: SF2 < min σ\* ≈ 0.0152 over K outside the C3 cell.
"Genuine zeros for σ ∈ [0.0175, 0.02]" should read "from σ ≈ 0.0153 at the cell corner".

**R11-5 (minor, robustness).**
- check_all counts a digit leaf '1'..'6' as 12 − 2W in any family, with no family limit. A digit in a KQF-engine family
  (limit 6) would count as up to 10. kq_run's face mode cannot emit digits, and no KQF or KQF2 file has any (i0).
  Reject digits for engine KQF and the C3 charts.
- The KQ/KQF 'o' test checks only l_lo and ω_lo. It is valid because the KQF, KQF2 and KQS domains lie inside the
  C3-core region (ε ≤ 0.01, μ ≤ 0.0505). Add those two bounds to the test (g2 does now).

**R11-6 (minor, doc).**
- PROOF.md rows "N" and "N ∩ K_lw, the rest": the e-range starts at fl(0.9996) (§2). State N with the 1e-9 tolerance,
  or say that only 'p' leaves use it.
- C3_CORE.md line 183 ("FE, FM, FIP, FIM and KQF are all complete") is stale; line 173 has the right list.
- PROOF.md line 176 ("covered by the face charts") should mention KQF2/KQS inside K_lw.

**R11-7 (note).** next_level drops face-chart boxes only for s_lo > SF, not those inside K with s_lo > SF2. This is
moot: FE, FM and FIP need no further level (§6). It is worth adding if a face level is ever re-run.

## My final-pass tools (updated; the round-10 versions are kept as round11/*_r10.*)

**g2.**
- Families KQF2 (kqf) and KQS (kq).
- For KQF, KQF2 and KQS, the header fields (engine, FACE_S, T_X, LIMIT_OUTER, c3prune) are checked in every file, and
  the base domains against launch_runs.sh.
- The 'p' radius comes from the family role: SF2 for KQS, 0.02 otherwise.
- KQ 'o' is allowed only in KQS.
- The exact 'o' test now includes ε ≤ 0.01 and μ ≤ 0.0505.
- Excuses: K (face charts) and K_lw (KQF).
- The p-dependency list is FE, FM, FIP, FIM, KQF, KQF2, KQS minus the family itself.
- My own bug, fixed: the C3-cell literal was the decimal mp.mpf('0.95') instead of the float 0.95 used by the code. It
  produced 128 false 'o' alarms on KQF; also fixed in g1.

**g1.**
- T_X and FACE_S from the header for the KQF recertification and the 100-digit outer count.
- The 'p' radius is SF2 for KQS.
- 'o' is accepted in KQS.

**g3.** The outer count uses the header T_X. The engine map now includes KQS and KQF2.

**g9.** KQF-engine files use `init(face=True, faces=H['FACE_S'], tx=H['T_X'])`.

**final_pass.sh.** KQF2 and KQS are in BB, in kqfiles and in the g9 lists (LF line endings checked).

Light test runs: g2 on KQF, KQF2, KQS and the face charts (i2, i6), g1 (i3) and g9 (i5) on KQF2 and KQS. All ran
cleanly.

## Not checked

- Completion of KQF.r1, KQF2.r1 and KQS.r1 (in progress).
- The C3 core's coverage of the cell (earlier rounds).
- The counts of KQF2 and KQS rest on engine KQ with the outer start t_lo = T_X as a parameter, audited in rounds 9–10.
  Here I re-certified only samples (8 + 6 per run).

## Addendum (13:40–13:55): the authors' fixes, verified

| # | fix | status |
|---|---|---|
| R11-1 | check_all.family pins the face radius per family in every file | OK |
| R11-2 | bb_replay passes the header FACE_S and T_X to kq_run.init | OK |
| R11-5 | `in_c3_cell` for 'o'; winding digits rejected in KQF-engine files | OK |
| R11-3 | duplicate runs killed; the levels resumed from the records files | OK; what the resume kept is intact (i7) |
| R11-4, R11-6 | C3_CORE §6 mechanism paragraph; C3_CORE line 188 and PROOF.md line 177 | OK |

**R11-1 (check_all, 13:40).**
- The pin is SF2 for KQS and KQF2 and 0.02 for every other KQ/KQF family. It is applied per file:
  - a header FACE_S different from the pin is a problem;
  - for a KQF-engine file, a header T_X ≠ 3.0·pin is a problem (3.0·0.010005 == 0.030015 and 3.0·0.02 == 0.06 in
    floats);
  - replay and prune_ok get the pin.
- All 30 current KQ/KQF headers satisfy the pins (i9). Sixteen older KQ files (KQ1, KQ2, KQB and KQK1/KQK2 levels up to r2–r3) have no FACE_S and
  are replayed with the pin 0.02.
- Residual (minor): check_all still does not compare the base domains of KQF2 and KQS with the region they must cover.
  My g2 does, and the current domains are verified in §2.

**R11-2 (bb_replay).**
- KQ files use `init(H['scale'], False, None, None, H['FACE_S'], H['T_X'])`; KQF files the same with face = True.
- Re-certification now tests the run's own claim, and check_all ties the header to the pin.

**R11-5.**
- `in_c3_cell` checks e_lo ≥ 0.99, μ_hi ≤ 0.0505, l_lo ≥ 0.95 and ω_lo ≥ 9π/20. It is used for KQ 'o' (c3prune) and
  KQF 'o'.
- Digits in a KQF-engine file are a problem.

**R11-3.** My process list (about 13:20) showed two kq_run processes per output file. They had different command
lines (with and without `--minw`) and separate worker pools. That predates the authors' kill at 13:27–13:28.

Now, at 13:4x, exactly one run process writes each output: KQ0T.r1, KQ2.r4, KQB.r6, KQF.r1, KQF2.r2, KQK1.r5, KQS.r1
and Ca.r5. The extra KQ0T sh process (92272) is the helper of the dead loop 89128. It only waits for its kq_run (82884),
and the live loop 46100 waits for KQ0T.r1.json.

What the resume kept was checked for integrity (i7). Every record replays exactly and reproduces its own stored leaf-code
histogram ('hist'), nleaf, nwind, fails, left and volume. A line spliced from two records would fail these checks.

| file | records | notes |
|---|---|---|
| KQF2.r1.json (final) | 75 | |
| KQF2.r1 .records.jsonl | 75 | 1 skipped fragment |
| KQK1.r5 .records.jsonl | 1317 (1314 roots) | the 3 duplicated roots are identical complete records; 1 skipped fragment |
| KQF2.r2 .records.jsonl | 15 | |

I added the hist cross-check to g2 for the kq and kqf kinds, as a check in the final pass. It passes on KQF, KQF2 (2
files), KQS, KQ2 (4), KQK2 (7) and KQB (6) (i8).

**R11-4 and R11-6.** The new C3_CORE §6 text matches h22 (outer member of the near pair, B zeros at 0.0787 and 0.1584,
onset about 0.0153). C3_CORE line 188 and PROOF.md line 177 now list KQF2 and KQS.

Nothing in round 11 is open except completion of the running levels. I am standing by for the final pass.
