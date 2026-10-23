"""Status table of the unfinished run families (for README.md / PROOF.md): for every family that check_all.py does not
report complete, its region, open boxes and run files.  Run from code/ after check_all.py (or standalone)."""
import check_all

REGION = {
    'KQ0T': 'e ∈ [0.835, 0.9501], t ≥ 0.9975 (circle much larger than the pericentre distance)',
    'KQ2': 'e ∈ [0.9499, 0.99], t ≥ 0.9975',
    'KQK1': 'e ≥ 0.99, μ ≤ 0.0505, l ≤ 0.95',
    'KQK2': 'e ≥ 0.99, μ ≤ 0.0505, l ≥ 0.95, ω ≤ 3π/20',
    'KQB': 'e ≥ 0.99, μ ≤ 0.0505, l ≥ 0.95, ω ∈ [3π/20, 9π/20] (the band)',
    'KQF': 'face neighbourhood ε, μ ≤ 4·10⁻⁴: outer count \\|τ\\| ≥ 0.06',
    'KQF2': 'corner K_lw of the face neighbourhood, ε, μ ≤ SF2²: outer count \\|τ\\| ≥ 0.03',
    'KQS': 'corner K_lw of the face neighbourhood, SF2 < max(√ε, √μ) ≤ 0.02',
    'C3 Ca': 'C3 cell, chart Ca (Q_z dominant)',
}


def rows():
    out = []
    for name, prefix in check_all.FAMILIES:
        r = check_all.family_cached(prefix)
        if r is None or r['complete']:
            continue
        why = 'open boxes %d' % r['owed'] if r['owed'] else 'complete itself; waits for: ' + '; '.join(r['problems'])
        out.append((name, REGION.get(name, ''), r['owed'], r['files'], why))
    return out


if __name__ == '__main__':
    for name, reg, owed, files, why in rows():
        print('| %s | %s | %s |' % (name, reg, why))
