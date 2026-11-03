#!/usr/bin/env python3
"""Plot strongest/second-strongest haplotype support from seven-column *.stat tables."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import PercentFormatter, MaxNLocator

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--dataset', action='append', required=True, metavar='LABEL=PATH')
parser.add_argument('--outdir', type=Path, required=True)
parser.add_argument('--min-total', type=float, default=200,
                    help='Plot windows with total support strictly above this value')
args = parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'axes.labelsize': 10,
                     'axes.titlesize': 11, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 300})
colors = ['#0072B2', '#D55E00', '#009E73', '#CC79A7']
datasets, summary, provenance = [], [], []
for specification in args.dataset:
    label, filename = specification.split('=', 1)
    path = Path(filename)
    windows, assembled, predicted, scores = [], [], [], []
    with path.open(encoding='utf-8') as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) != 7:
                raise ValueError(f'{path}:{line_number}: expected seven columns')
            windows.append(fields[0]); assembled.append(int(fields[1])); predicted.append(int(fields[2]))
            scores.append([float(v) for v in fields[3:7]])
    scores = np.asarray(scores)
    assembled, predicted = np.asarray(assembled), np.asarray(predicted)
    if scores.ndim != 2 or scores.shape[1] != 4 or not np.all(np.isfinite(scores)) or np.any(scores < 0):
        raise ValueError('Support must be finite, nonnegative, and have four columns')
    if len(set(windows)) != len(windows) or np.any((assembled < 1) | (assembled > 4)):
        raise ValueError('Duplicate window IDs or invalid haplotype labels')
    if not np.array_equal(predicted, scores.argmax(axis=1) + 1):
        raise ValueError('Stored maximum-support haplotypes disagree with first-maximum rule')
    total = scores.sum(axis=1)
    ranked = np.sort(scores, axis=1)
    valid = total > args.min_total
    first = np.divide(ranked[:, -1], total, out=np.full(len(total), np.nan), where=total > 0)
    second = np.divide(ranked[:, -2], total, out=np.full(len(total), np.nan), where=total > 0)
    delta = first - second
    tie = ranked[:, -1] == ranked[:, -2]
    slug = re.sub(r'[^A-Za-z0-9]+', '_', label).strip('_')
    with gzip.open(args.outdir / f'{slug}.window_metrics.tsv.gz', 'wt', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['window', 'assembled_haplotype', 'maximum_support_haplotype',
                         'total_support', 'strongest_fraction', 'second_fraction',
                         'normalized_margin', 'exact_tie', 'included'])
        for i, window in enumerate(windows):
            writer.writerow([window, assembled[i], predicted[i], total[i], first[i], second[i],
                             delta[i], int(tie[i]), int(valid[i])])
    for hap in [0, 1, 2, 3, 4]:
        base = np.ones(len(total), dtype=bool) if hap == 0 else assembled == hap
        keep = base & valid
        summary.append({'dataset': label, 'haplotype': 'all' if hap == 0 else f'hap{hap}',
                        'recorded_windows': int(base.sum()), 'included_windows': int(keep.sum()),
                        'low_support_windows': int((base & ~valid).sum()),
                        'exact_ties_included': int((keep & tie).sum()),
                        'mismatched_windows_before_interval_filter': int((keep & (assembled != predicted)).sum()),
                        'median_normalized_margin': float(np.median(delta[keep])) if keep.any() else None,
                        'margin_le_0.05_fraction': float(np.mean(delta[keep] <= .05)) if keep.any() else None,
                        'margin_le_0.1_fraction': float(np.mean(delta[keep] <= .1)) if keep.any() else None})
    provenance.append({'dataset': label, 'input': str(path.resolve()),
                       'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    datasets.append(dict(label=label, slug=slug, total=total, valid=valid, first=first,
                         second=second, delta=delta, assembled=assembled, ties=int((valid & tie).sum())))

def panels(fig, axes, item, letters):
    left, right = axes
    valid = item['valid']; n = int(valid.sum())
    if not n:
        raise ValueError(f"No included windows: {item['label']}")
    image = left.hexbin(item['first'][valid], item['second'][valid], gridsize=65,
                        extent=(.24, 1.01, -.01, .51), mincnt=1, cmap='viridis', norm=LogNorm(),
                        linewidths=0, rasterized=True)
    left.plot([.24, .51], [.24, .51], '--', color='#777777', lw=.8)
    left.set(xlim=(.23, 1.02), ylim=(-.015, .52), xlabel='Strongest haplotype support / total support',
             ylabel='Second-strongest support / total support')
    left.set_title(f"{item['label']}  |  n = {n:,}", loc='left', pad=11)
    bar = fig.colorbar(image, ax=left, fraction=.046, pad=.025)
    bar.set_label('Windows per hexagon')
    for h, color in enumerate(colors, 1):
        values = np.sort(item['delta'][valid & (item['assembled'] == h)])
        if len(values):
            x, counts = np.unique(values, return_counts=True)
            y = np.cumsum(counts) / len(values)
            right.step(np.r_[0, x], np.r_[0, y], where='post', color=color, lw=1.5,
                       label=f'Hap{h} (n={len(values):,})')
    right.set(xlim=(0, 1.015), ylim=(0, 1.02), xlabel='Normalized support margin, (strongest - second) / total',
              ylabel='Cumulative fraction of windows')
    right.set_xlabel('Normalized support margin\n(strongest - second) / total support')
    right.set_title(f"Exact ties: {item['ties']:,} ({item['ties']/n:.3%})", loc='left', pad=11)
    right.legend(loc='upper left', frameon=False, fontsize=8, handlelength=1.8)
    right.grid(axis='y', color='#eeeeee', lw=.7)
    inset = right.inset_axes([.14, .32, .43, .27])
    max_tail = 0
    for h, color in enumerate(colors, 1):
        values = np.sort(item['delta'][valid & (item['assembled'] == h)])
        if len(values):
            x, counts = np.unique(values, return_counts=True)
            y = np.cumsum(counts) / len(values)
            inset.step(np.r_[0, x], np.r_[0, y], where='post', color=color, lw=1)
            max_tail = max(max_tail, np.mean(values <= .1))
    inset.set(xlim=(0, .1), ylim=(0, max(max_tail*1.18, .0001)))
    inset.set_title('Low-margin tail', fontsize=8, pad=3)
    inset.yaxis.set_major_formatter(PercentFormatter(1, decimals=2))
    inset.yaxis.set_major_locator(MaxNLocator(3))
    inset.tick_params(labelsize=7, length=2, pad=2)
    inset.set_xlabel('Normalized margin', fontsize=7, labelpad=1)
    for ax, letter in zip(axes, letters):
        ax.text(-.13, 1.10, letter, transform=ax.transAxes, fontweight='bold', fontsize=13)

def save(fig, name):
    fig.savefig(args.outdir / f'{name}.preview.png', dpi=140, bbox_inches='tight', facecolor='white')
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(args.outdir / f'{name}.{extension}', bbox_inches='tight', facecolor='white')
    plt.close(fig)

for item in datasets:
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), layout='constrained', gridspec_kw={'width_ratios': [1.1, 1]})
    panels(fig, axes, item, 'ab')
    save(fig, item['slug'] + '.support_separation')
fig, axes = plt.subplots(len(datasets), 2, figsize=(10, 3.9*len(datasets)),
                         squeeze=False, layout='constrained', gridspec_kw={'width_ratios': [1.1, 1]})
for row, item in enumerate(datasets):
    panels(fig, axes[row], item, [chr(97 + row*2), chr(98 + row*2)])
save(fig, 'combined.support_separation')
with (args.outdir/'support_summary.tsv').open('w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter='\t')
    writer.writeheader(); writer.writerows(summary)
(args.outdir/'provenance.json').write_text(json.dumps({'min_total_exclusive': args.min_total,
    'normalization': 'divide each of four supports by their sum',
    'stage': 'before candidate selection and interval merging', 'inputs': provenance,
    'numpy': np.__version__, 'matplotlib': matplotlib.__version__}, indent=2), encoding='utf-8')
print(json.dumps([row for row in summary if row['haplotype'] == 'all'], indent=2))
