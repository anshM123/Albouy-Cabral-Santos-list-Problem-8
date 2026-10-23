"""AUDIT-11/12 i10: which run files carry a split_rule in their header (all families, base + continuations)."""
import glob, gzip, json, os
CODE = r'<repo>\code'


def hdr(path):
    op = gzip.open if path.endswith('.gz') else open
    for size in (8000, 400000, 4000000):
        with op(path, 'rb') as fh:
            head = fh.read(size).decode('latin1')
        i = head.index('"header"'); j = head.index('{', i); depth = 0
        for k in range(j, len(head)):
            depth += head[k] == '{'; depth -= head[k] == '}'
            if depth == 0:
                return json.loads(head[j:k + 1])
    raise ValueError('header longer than 4 MB')


withrule = []; norule = 0; bad = []
for p in sorted(glob.glob(os.path.join(CODE, 'runs', '*.json*')) + glob.glob(os.path.join(CODE, 'c3', 'runs', '*.json*'))):
    b = os.path.basename(p)
    if 'boxes' in b or b.endswith('.jsonl') or '.recert' in b:
        continue
    try:
        H = hdr(p)
    except Exception as ex:
        bad.append(b); continue
    if H.get('split_rule'):
        withrule.append('%s:%s' % (b, H['split_rule']))
    else:
        norule += 1
print('files with a header split_rule:', withrule)
print('files without:', norule, '; unreadable headers:', bad)
