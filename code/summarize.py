"""Aggregate numbers over all run families (for README / paper): run after check_all.py.
Prints, per family and in total: run files, records, certified leaves, pruned boxes, failed leaves
(all covered when the family is complete), and the overall completeness."""
import check_all

tot = dict(files=0, records=0, leaves=0, pruned=0, fails=0)
allok = True
rows = []
for name, prefix in check_all.FAMILIES:
    r = check_all.family(prefix)
    if r is None:
        rows.append((name, None)); allok = False
        continue
    rows.append((name, r))
    allok = allok and r['complete']
    for k in ('files', 'records', 'leaves', 'pruned', 'fails'):
        tot[k] += r[k]
for name, r in rows:
    if r is None:
        print('%-10s missing' % name)
    else:
        print('%-10s files %d records %7d leaves %9d pruned %7d failed %6d complete %s' % (
            name, r['files'], r['records'], r['leaves'], r['pruned'], r['fails'], r['complete']))
print('TOTAL', tot, 'families', len(rows), 'all complete', allok)
