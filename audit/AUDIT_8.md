# AUDIT 8: two-tier counts, resume of capped records, the KQ split angle; status, and a completeness gap at the corner e = 1, μ = 0

Independent audit, 2026-10-23, 10:18–10:55. Scripts and outputs are in `audit/round8/`. The round-7 final-pass
tools were updated for this round (`DEEP_ALWAYS`; see the end).

## Verdict

- **Round-8 changes: no soundness issue.** The two-tier counts, the resume rule for capped records and the
  box-dependent split angle `phi_split` are all sound.
- **R8-1 (major, completeness, blocks the proof):** the continuation chains **KQK1, KQK2 and KQB cannot finish.** At the
  corner e = 1, μ = 0 (t = 1) the polynomial has exactly 12 real zeros for every (l, ω) with l ≥ l*(ω). No sound zero
  count can certify a box that contains such a point. So no amount of subdivision closes these boxes.
  - No recorded certificate is wrong.
  - The atlas claim "D3 (e = 1, μ = 0) resolved by engine KQ" (DEGENERACY_TREE, PROOF §3–4) is false for l ≥ l*(ω).
  - A dedicated blow-up chart is needed (see R8-1 below).

| id | severity | finding |
|---|---|---|
| R8-1 | **major (completeness; blocks the proof)** | Corner e = 1, t = 1 with l ≥ l*(ω): 12 real zeros of P; boxes containing it never certify. Details below. |
| R8-2 | info | Logging only: `bb_run`'s "certified volume fraction" exceeds 1 after a resume (already noted in R7-11). The resume rule for capped records makes `records` include capped records, as intended. |

## (1) Two-tier counts: OK

**kq_run (10:15).**
- `count_total(bx, rounds, cap, phi0)` is the old `certify` body with (rounds, cap) passed through to both the away count and the near count. The away target is 6, then 0; the near count uses the remaining budget.
- `certify` tries (ROUNDS1, CAP1) = (6, 1024) first, then (16, 8192) for boxes whose normalised widths are all below DEEP_W = 2e-3, or always under DEEP_ALWAYS.
- Every `away_count_t2` and `count` call returns a valid upper bound for any rounds/cap, with min over rounds. Every combination is therefore a valid certificate, and so is the returned minimum.

**c3_run.** `cnt2` gives the same two-tier rule for the zone and non-zone counts. The zone outer budget, ZONE_OUTER − c1, is unchanged.

**Headers** record rounds1, cap1 and deep_w. `bb_replay` sets `kq_run.DEEP_ALWAYS = True` for re-certification.

**Stress test of the current `kq_run.certify` at 100 digits** (`h6`, two-tier counts + `phi_split`): 250 random
boxes over the five KQ domains, widths 1e-5–1e-1, biased to the faces; 184 certified. Each was checked at 8
points: **0 violations**, and N_mult ≤ total everywhere (minimum slack 0, the bound attained).

## (2) Resume keeps capped records: OK

- `--resume` builds one record per box: a dict keyed by box, where a complete record replaces anything and a capped record is kept only if none exists yet. It then skips all those boxes.
- The open stacks of capped records stay owed. `next_level.py` passes them to the next level, and check_all/`g2` require every owed box to be a continuation root, and every continuation root to be owed.
- Own replay of the chains (`h0_g2_chains.out`), with each file's own scale; 0 problems everywhere:

  | family | files | records | open boxes |
  |---|---|---|---|
  | KQ2 | 3 | 1,039 | 340 |
  | KQK1 | 3 | 1,623 | 460 |
  | KQK2 | 3 | 20 | 2,228 |
  | KQB | 3 | 347 | 143 |
  | Ca | 3 | 1,426 | 432 |

  These agree with the authors' check_all numbers.

## (3) Split angle `phi_split`: OK

- **Rule.** φ0 = PHI_BAND = 1.25 when l1 ≥ 0.9, π/2 − ω0 > 0.6 and π/2 − ω1 < 1.05; otherwise φ0 = 0.9. It is a deterministic function of the box.
- **Covering.** The near range is |τ| ≤ tan(φ0/2) + 1e-9, i.e. |f − π| < 2·atan(tan(φ0/2) + 1e-9), which is more than φ0. The away arc starts at φ0 − 1e-9 − PAD. So the two parts overlap for every φ0 ∈ (0, π).
- **Engines.** The away arcs keep |u| ≤ 3π, so the 62.2 argument bound is unchanged. The polynomial near count is valid for any root interval; its float pads scale with |y|^(k−j).
- **Re-certification.** Under DEEP_ALWAYS, `certify` also tries PHI_0. A count that fails for every φ0 returns the minimum of valid bounds, which is still a valid bound.

## (4) Status checks

- **KQ1 is complete:** KQ1.json + KQ1.r1.json; 3,009 records, 90,528 leaves (43,408 winding), 0 fails, 0 owed. My `g2` agrees.
- 100-digit spot checks on KQ1 and KQ1.r1 (`h7`): 270 points, 0 violations, all N_mult even, max 10, minimum winding margin 0.028.

## R8-1: the corner e = 1, μ = 0 is not resolved by KQ when l ≥ l*(ω)

**Why the corner has 12 real zeros.**
- At e = 1 and μ = 0, P = −D²B² with D = 1 + cos f. D² has a 4-fold zero at f = π.
- B(π) = −2P_zQ_z(1 − e) = 0, and B'(π) = −2c² ≠ 0, so B² adds a double zero at π.
- B has trigonometric degree 2, so it has 2 or 4 real zeros. With 4, N_mult(P) = 4 + 2·4 = **12**.
- At l = 1, B = −2cos(f + ω)(sin(f + ω) + sin ω) (checked to 1e-15). Its four real zeros are π, −2ω, π/2 − ω and 3π/2 − ω; the last two are the axis crossings.
- **Range of l.** `h4`/`h5` compute l*(ω), the smallest l with 4 real zeros of B: l*(0) = 0.707, l*(0.47) = 0.909, l*(0.74) = 0.953, l*(1.41) = 0.9985, and l*(ω) → 1 as ω → π/2.
- **Where it bites.**
  - KQK1 (l ≤ 0.95), for ω ≲ 0.74;
  - all of KQK2 (l ≥ 0.95, ω ≤ 3π/20);
  - most of KQB (l ≥ 0.95, ω ∈ [3π/20, 9π/20]).
  - The C3 cell (ω ≥ 9π/20) is covered by its own blow-up.

**Consequence.** A sound count bounds the zeros at every point of the closed box, including the corner, so it gives a total of at least 12. The winding certificate is also impossible there, because all 12 zeros lie on the unit circle, so W = 0.

**Demonstration** (`h8`, current `kq_run.certify` in re-certification mode):
- Boxes containing the corner at (l, ω) = (1, 0.2), (0.999, 0.3), (0.97, 0.6) and (0.9, 0.1) return (False, 12) at widths 1e-3, 1e-5 and 1e-7.
- Boxes of the same size just below e = 1, or just below t = 1, certify with total 10 or 11.
- This matches KQK2.r2: all 2,146 failed leaves have total 12 and sit at e ∈ [0.99998, 1], t ∈ [0.99999992, 1], l ≥ 0.9992.

**Open boxes in the degenerate set** (`h5`). Count of open boxes touching e = 1 and t = 1 whose (l, ω) range meets l ≥ l*(ω):

| family | open boxes | in the degenerate set |
|---|---|---|
| KQK1 | 460 | 12 |
| KQK2 | 2,228 | 2,183 |
| KQB | 142 | 21 |

Further levels (r3, r4, …) cannot close them.

**The theorem itself is not in danger here.**
- Physical points next to the corner have N_mult = 10 (`h3`: four points with μ ≤ 2e-4, at 100 digits).
- At leading order in the scaling τ ~ ν, ε ~ ν², μ ~ ν², the near-aphelion part is 4c²y²ν⁶[4m² − c²(E + y²/2)²]. Its real zeros are y = 0 (double) and at most one pair from E + y²/2 = 2m/c, so there are at most 4 near zeros. The away part has at most 6, three double zeros of B². That gives a total of at most 10.

**What is needed:** a C3-like weighted blow-up of the near-aphelion part at (ε, μ) → 0, with weights (τ, ε, μ) = (1, 2, 2).
- Cone charts dominated by E or by m, plus an intermediate τ-dominant chart, with c = σ(π) > 0 as an ordinary parameter (C3_CORE §2–§4 is the template).
- The away count stays as it is.
- The blow-up region must cover a neighbourhood of {e = 1, μ = 0, l ≥ l*(ω) − margin, ω ≤ 9π/20}.
- Until then, the KQK1/KQK2/KQB chains cannot complete. Their open boxes should be split into the corner part, owed to the new charts, and the rest, which can continue as now.

## Tools updated for the final pass

- `round7/g1_kq_wind.py`, `round7/g9_recert_all.py` and `round7/g4_c3_spot.py` set `DEEP_ALWAYS = True` on `kq_run` and `c3_run`. Records certified with (16, 8192), or with PHI_0, then re-certify.
- `round7/final_pass.sh` is unchanged otherwise. It picks up r1..rN and .gz files. Its `g2` will report the KQK1/KQK2/KQB families incomplete as long as R8-1 is open.

## Not checked

- The authors' timings (two-tier speed-up; split-angle measurements).
- The level-3/4 files: they did not exist as JSON at audit time.
