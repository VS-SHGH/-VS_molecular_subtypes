"""Draw a text-and-data graphical abstract from the reported study results.

This schematic depicts associations, not cell-state transitions or drug effects.
The variance bar is read from the same numerical source as Figure 4h.
"""
from pathlib import Path
import argparse
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root
    source = root / 'outputs/mgs/figure4h_variance_partitioning.csv'
    with source.open(newline='') as f:
        fractions = {row['Component']: float(row['Percent']) for row in csv.DictReader(f)}
    assert abs(sum(fractions.values()) - 100) < 1e-6
    out = root / '02_Submission_Materials/05_Graphical_Abstract'
    out.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
    # Inches equivalent to 183 x 119 mm at export.
    fig = plt.figure(figsize=(7.2047244, 4.6850394), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 183), ylim=(0, 119))
    ax.axis('off')
    ink, muted = '#213342', '#4C5D66'
    colors = ['#34886D', '#338EAF', '#C5605C']
    def text(x, y, s, size=8, color=ink, weight='normal', ha='left', va='top'):
        return ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
                       ha=ha, va=va, linespacing=1.35)
    def box(x, y, w, h, fill, edge='none'):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge, linewidth=0.6))

    text(7, 113, 'Transcriptional variation in vestibular schwannoma', size=11, weight='bold')
    text(7, 106, 'Bulk subtypes, a continuous expression gradient and surgical associations', size=8.5, color=muted)
    box(7, 89, 169, 11, '#F0F4F6')
    text(11, 97, 'Discovery bulk: 38 tumours', size=8.2, weight='bold')
    text(70, 97, 'External bulk: 57 + 31', size=8.2)
    text(124, 97, 'Public single-cell: 15', size=8.2)

    for x, color, label, subtitle in zip(
        (7, 65, 123), colors, ('C1  |  n = 13', 'C2  |  n = 12', 'C3  |  n = 13'),
        ('Proliferative expression', 'ECM-associated expression', 'Immune-associated expression')):
        box(x, 69, 53, 16, '#F8FAFB', '#D8E0E4')
        box(x, 69, 1.3, 16, color)
        text(x+4, 81.5, label, size=9, color=color, weight='bold')
        text(x+4, 75, subtitle, size=7.8)

    text(7, 63.5, 'MGS is closely associated with microenvironmental expression scores', size=9, weight='bold')
    # A stacked bar showing sequential explained variance, not independent shares.
    labels = ('Stromal Score', 'Immune Score', 'Residual')
    bar_colors = ('#9CB5A6', '#7DB4C6', '#DEE3E6')
    x, y, width = 7, 48, 169
    for label, color in zip(labels, bar_colors):
        value = fractions[label]
        w = width * value / 100
        box(x, y, w, 10, color)
        text(x+w/2, y+5.2, f'{label}  {value:.1f}%', size=8.3, ha='center', va='center')
        x += w
    joint = fractions['Stromal Score'] + fractions['Immune Score']
    text(7, 45, f'Joint R² = {joint:.1f}%  |  Stromal entered first, Immune second; increments depend on order.', size=7.5)
    text(7, 40, 'An internal expression association, not a causal partition or a temporal trajectory.', size=7.5, color=muted)

    ax.plot([7, 176], [34.5, 34.5], color='#CAD3D9', linewidth=0.65)
    text(7, 31, 'Clinical associations', size=9, weight='bold')
    text(7, 25, 'C2: harder surgical texture\nC3 adhesion and blood NK-cell percentage:\ncovariate-sensitive, exploratory associations', size=7.8)
    text(99, 31, 'Single-cell interpretation', size=9, weight='bold')
    text(99, 25, 'APP, MIF and SPP1: candidate interactions\nNormalisation-sensitive myeloid scoring\nNo validated C3-high-specific network hub', size=7.8)
    text(7, 6, 'Descriptive findings require independent clinical and mechanistic validation.', size=7.7, weight='bold')

    fig.savefig(out / 'Graphical_Abstract.pdf', metadata={'Title': 'Transcriptional variation in vestibular schwannoma',
                                                        'Author': '', 'Subject': 'Study summary'})
    fig.savefig(out / 'Graphical_Abstract.svg', metadata={'Title': 'Transcriptional variation in vestibular schwannoma'})
    fig.savefig(out / 'Graphical_Abstract.png', dpi=600)
    fig.savefig(out / 'Graphical_Abstract.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)
    print(out)


if __name__ == '__main__':
    main()
