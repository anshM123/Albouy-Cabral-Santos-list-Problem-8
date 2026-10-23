"""Write the certification status (unfinished families and their open boxes) into README.md, PROOF.md and
paper/paper.md, between the markers <!-- STATUS --> ... <!-- /STATUS --> (replacing an earlier status block).
usage: python apply_status.py "YYYY-MM-DD HH:MM" """
import re
import sys
import status_table

when = sys.argv[1]
rows = status_table.rows()
if not rows:
    body = ('**Certification status (%s): complete.** Every run family in RESULTS.md is complete (exact replay, every '
            'open or failed box covered); see RESULTS.md and verify/.' % when)
else:
    tab = '\n'.join('| %s | %s | %s |' % (n, reg, why) for n, reg, owed, files, why in rows)
    tot = sum(owed for n, reg, owed, files, why in rows)
    body = f"""**Certification status ({when}): not yet complete.** The certificate covers the whole parameter space except
{tot} open boxes of the run families below, whose continuation runs are still in progress. The run families
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
{tab}
"""
block = '<!-- STATUS -->\n' + body.strip() + '\n<!-- /STATUS -->'
for p, anchor in (('../README.md', '## The result\n'), ('../PROOF.md', '---\n\n## 1. Reduction'),
                  ('../paper/paper.md', '## Abstract\n')):
    s = open(p, encoding='utf-8').read()
    if '<!-- STATUS -->' in s:
        s = re.sub(r'<!-- STATUS -->.*?<!-- /STATUS -->', lambda m: block, s, flags=re.S)
    else:
        assert s.count(anchor) == 1, (p, anchor)
        if p.endswith('PROOF.md'):
            s = s.replace(anchor, block + '\n\n' + anchor)
        else:
            s = s.replace(anchor, anchor + '\n' + block + '\n')
    open(p, 'w', encoding='utf-8', newline='\n').write(s)
    print('status written to', p)
