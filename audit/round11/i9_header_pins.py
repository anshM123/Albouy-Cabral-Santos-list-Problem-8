"""AUDIT-11 i9 (light): the new check_all pins (R11-1) against every current KQ / KQF run file header."""
import glob, gzip, json, os
R = r'<repo>\code\runs'
SF, SF2 = 0.02, 0.010005


def hdr(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rb') as fh:
        head = fh.read(6000).decode('latin1')
    i = head.index('"header"'); j = head.index('{', i); depth = 0
    for k in range(j, len(head)):
        depth += head[k] == '{'; depth -= head[k] == '}'
        if depth == 0:
            return json.loads(head[j:k + 1])


rows = []
for p in sorted(glob.glob(os.path.join(R, 'KQ*.json')) + glob.glob(os.path.join(R, 'KQ*.json.gz'))):
    b = os.path.basename(p)
    if '_r' in b and 'boxes' in b:
        continue
    fam = b.split('.')[0]
    try:
        H = hdr(p)
    except Exception as ex:
        rows.append((b, 'header unreadable: %s' % ex)); continue
    eng = H.get('engine'); pin = SF2 if fam in ('KQS', 'KQF2') else SF; pr = []
    if eng in ('KQ', 'KQF'):
        if H.get('FACE_S') is not None and H['FACE_S'] != pin:
            pr.append('FACE_S %r != %r' % (H['FACE_S'], pin))
        if eng == 'KQF' and H.get('T_X') != 3.0 * pin:
            pr.append('T_X %r != %r' % (H.get('T_X'), 3.0 * pin))
    rows.append((b, 'engine %s FACE_S %s T_X %s c3prune %s rule %s -> %s' % (eng, H.get('FACE_S'), H.get('T_X'), H.get('c3prune'),
                                                                        H.get('split_rule'), pr or 'OK')))
for r in rows:
    print('%-16s %s' % r)
