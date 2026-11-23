#!/usr/bin/env python3
"""Summarize the direct Kasai component-decoding experiment logs."""
import csv
import re
import sys
from pathlib import Path

run_dir = Path(sys.argv[1])
fields = ['f_m', 'p_D', 'trials', 'Gamma_block_errors', 'Delta_block_errors',
          'Gamma_block_error', 'Delta_block_error', 'Gamma_syndrome_fail',
          'Delta_syndrome_fail', 'N_qubits', 'actual_RQ', 'seed', 'max_iters']
rows = []
for path in sorted(run_dir.glob('fm_*.log')):
    content = path.read_text()
    pairs = dict(re.findall(r'([A-Za-z_]+)=([^\s]+)', content))
    if all(name in pairs for name in fields):
        rows.append({name: pairs[name] for name in fields})
with (run_dir / 'summary.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
print(f'Wrote {len(rows)} points to {run_dir / "summary.csv"}')
