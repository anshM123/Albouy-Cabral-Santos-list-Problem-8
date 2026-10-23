"""Close the in-progress level of every unfinished run family (snapshot): for each family, find the latest level whose
records file exists but whose JSON does not, write that JSON from the records (driver with TMAX = -1 and --resume:
every record returned so far is kept and every box not yet processed is stored as an open record), then list its open
boxes for the next level (next_level.py).  The flags of each level are the ones it was run with, so that the header
(split rule, scale, face radius, outer start) matches the records.  Run from code/ with all drivers stopped.
usage: python close_levels.py"""
import os
import subprocess
import sys

PY = sys.executable
HP = '1.5707963267948966'
# family -> (kind, scale, {level: flags}); a level missing from the dict uses the flags of the highest listed level
# below it
FAM = {
    'KQ2': ('kq', ['1', '1', HP, '0.005'], {4: ['--splitrule=face/30', '--minw=1e-7']}),
    'KQB': ('kq', ['1', '1', HP, '0.0025'], {5: ['--splitrule=face/30'], 6: ['--splitrule=face/30', '--minw=1e-7']}),
    'KQF': ('kq', ['0.024', '1', HP, '1.02e-5'], {1: ['--face'], 2: ['--face', '--minw=1e-7']}),
    'KQF2': ('kq', ['0.003', '1', HP, '3e-7'], {1: ['--face', '--faces=0.010005', '--peri', '--minw=1e-7']}),
    'KQS': ('kq', ['1', '1', HP, '0.0025'], {1: ['--faces=0.010005', '--c3prune', '--splitrule=shell/30', '--peri',
                                                 '--minw=1e-7']}),
    'KQK1': ('kq', ['1', '1', HP, '0.0025'], {5: ['--splitrule=face/30', '--peri', '--minw=1e-7']}),
    'KQ0T': ('kq', ['1', '1', HP, '0.005'], {1: ['--splitrule=face/30'], 2: ['--splitrule=face/30', '--minw=1e-7']}),
    'Ca': ('c3', ['0.37', '0.25', '1.002', '0.5'], {6: []}),
}


def flags_for(spec, level):
    ks = [k for k in spec if k <= level]
    if not ks:
        raise SystemExit('no flags for level %d' % level)
    return spec[max(ks)]


for name, (kind, scale, spec) in FAM.items():
    d = 'runs' if kind == 'kq' else 'c3/runs'
    lev = None
    for L in range(1, 40):
        if os.path.exists('%s/%s.r%d.json.records.jsonl' % (d, name, L)) and not os.path.exists('%s/%s.r%d.json' % (d, name, L)):
            lev = L
    if lev is None:
        print(name, ': no open level'); continue
    fl = flags_for(spec, lev)
    if kind == 'kq':
        cmd = [PY, 'kq_run.py', '--boxes', 'runs/%s_r%d_boxes.json' % (name, lev), '--scale'] + scale + \
              ['-1', '1', 'runs/%s.r%d.json' % (name, lev)] + fl + ['--resume']
        cwd = '.'
    else:
        cmd = [PY, 'c3_run.py', name, '--boxes', 'runs/%s_r%d_boxes.json' % (name, lev), '--scale'] + scale + \
              ['-1', '1', 'runs/%s.r%d.json' % (name, lev), '--resume']
        cwd = 'c3'
    print(name, 'closing level', lev, ' '.join(cmd[1:]), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True, stdout=open('%s/%s.r%d.close.log' % (d, name, lev), 'w'), stderr=subprocess.STDOUT)
    out = subprocess.run([PY, 'next_level.py', kind, name, str(lev)], capture_output=True, text=True, check=True)
    print('  ', out.stdout.strip(), flush=True)
