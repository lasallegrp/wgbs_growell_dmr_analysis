#!/usr/bin/env python3
"""
Step 6 job builder: PCA plots (06_dmr_pca.R) for every finished outcome/timepoint.

Layout (GROWELL):  <run_root>/<outcome>/<timepoint>/  e.g. results/dmrichr_output/ppwr/36-38wk/
  "run" = outcome folder (apo, egwg, experimental_group, ppwr); "comparison" = timepoint folder.
  A timepoint folder counts only if it holds sample_info.xlsx (DMRichR input), so other
  folders under run_root (DMRichR_overlap, DMRichR_vs_comethyl, jobs/) are ignored.

Writes, in --out_dir (default <run_root>/jobs/pca/<outcomes joined by __, or all>/):
  pca_manifest.tsv   one row per outcome x timepoint x region set
  run_pca.sbatch     SLURM array job, one task per row

Plots: <outcome>/<timepoint>/DMRs/PCA/<outcome>_<timepoint>_<DMR|background>_PCA.pdf/.xlsx
With --collect, the job also copies every PDF into <out_dir>/pdf/ to flip through.

Usage (from anywhere inside the repo, e.g. results/dmrichr_output):
  python3 scripts/dmrichr/06_make_pca_jobs.py                       # all outcomes, all timepoints
  python3 scripts/dmrichr/06_make_pca_jobs.py --runs ppwr,egwg --regions both --pcs 3 --collect
  python3 scripts/dmrichr/06_make_pca_jobs.py --runs ppwr --comparisons Baseline,36-38wk --label
  sbatch results/dmrichr_output/jobs/pca/all/run_pca.sbatch
"""

import argparse
import os
import shlex
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent


def find_repo_root():
    here = [SCRIPT_DIR, *SCRIPT_DIR.parents]
    for q in here:
        if (q / '.git').exists():
            return q
    for q in here:                                   # no git: first parent holding results/
        if (q / 'results').is_dir():
            return q
    return Path.cwd()


REPO_ROOT = find_repo_root()
PIXI_CANDIDATES = ['pixi.toml', 'reproducibility/dmrichr/pixi/pixi.toml', 'pixi/pixi.toml']
METH = {'dmr': ('DMR', 'DMRs/DMR_individual_smoothed_methylation.txt'),
        'background': ('background', 'DMRs/background_region_individual_smoothed_methylation.txt')}


def absol(p):
    """Absolute as given; relative -> cwd if it exists there, else the repo root."""
    p = Path(p)
    if p.is_absolute():
        return p
    return (Path.cwd() / p).resolve() if (Path.cwd() / p).exists() else REPO_ROOT / p


def find_pixi(arg):
    if arg:
        p = absol(arg)
        if not p.exists():
            sys.exit(f"--pixi_manifest not found: {p}")
        return p
    for c in PIXI_CANDIDATES:
        if (REPO_ROOT / c).exists():
            return REPO_ROOT / c
    sys.exit(f"No pixi.toml under {REPO_ROOT} (tried {', '.join(PIXI_CANDIDATES)}); pass --pixi_manifest")


def iter_comparisons(root, runs, keep):
    """Yield (outcome, timepoint, folder) for every timepoint folder with sample_info.xlsx."""
    if runs is None:
        runs = sorted(d.name for d in root.iterdir() if d.is_dir() and
                      any((c / 'sample_info.xlsx').exists() for c in d.iterdir() if c.is_dir()))
        if not runs:
            sys.exit(f"No <outcome>/<timepoint>/sample_info.xlsx under {root}")
    for run in runs:
        rd = root / run
        if not rd.is_dir():
            sys.exit(f"Outcome folder not found: {rd}")
        for comp in sorted(p.name for p in rd.iterdir() if p.is_dir()):
            cd = rd / comp
            if (keep and comp not in keep) or not (cd / 'sample_info.xlsx').exists():
                continue
            yield run, comp, cd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', default=None, help='comma list of outcomes, e.g. ppwr,egwg (default: all found)')
    ap.add_argument('--run_root', default='results/dmrichr_output')
    ap.add_argument('--comparisons', default=None, help='comma list of timepoints, e.g. Baseline,36-38wk (default: all)')
    ap.add_argument('--out_dir', default=None, help='default <run_root>/jobs/pca/<outcomes joined by __ | all>')
    ap.add_argument('--regions', default='dmr', choices=['dmr', 'background', 'both'],
                    help='dmr = DMRs only; background = all tested regions; both [default dmr]')
    ap.add_argument('--shape', default=None, help='sample_info column for point shape, e.g. Sex')
    ap.add_argument('--colour', default=None, help='sample_info column to colour by (default: test covariate)')
    ap.add_argument('--label', action='store_true', help='label points with sample names')
    ap.add_argument('--show_p', action='store_true', help='show PC group-test p-values in the subtitle')
    ap.add_argument('--pcs', type=int, default=2, choices=[2, 3])
    ap.add_argument('--no_scale', action='store_true', help='do not scale regions to unit variance')
    ap.add_argument('--collect', action='store_true', help='also copy all PDFs into <out_dir>/pdf/')
    ap.add_argument('--force', action='store_true', help='redo timepoints that already have a PCA plot')
    ap.add_argument('--scripts_dir', default=str(SCRIPT_DIR), help='folder holding 06_dmr_pca.R [default: this script\'s folder]')
    ap.add_argument('--pixi_manifest', default=None, help=f'pixi.toml [default: first of {", ".join(PIXI_CANDIDATES)} in the repo]')
    ap.add_argument('--mem', default=None, help='default 16G (dmr) or 64G (background/both)')
    ap.add_argument('--time', default='04:00:00')
    args = ap.parse_args()

    runs = [r.strip() for r in args.runs.split(',') if r.strip()] if args.runs else None
    keep = set(c.strip() for c in args.comparisons.split(',')) if args.comparisons else None
    regions = ['dmr', 'background'] if args.regions == 'both' else [args.regions]
    root = absol(args.run_root)
    if not root.is_dir():
        sys.exit(f"--run_root not found: {root}")
    out = absol(args.out_dir) if args.out_dir else root / 'jobs' / 'pca' / ('__'.join(runs) if runs else 'all')
    rs = absol(args.scripts_dir) / '06_dmr_pca.R'
    if not rs.exists():
        sys.exit(f"R script not found: {rs} (use --scripts_dir)")
    pixi = find_pixi(args.pixi_manifest)

    rows, skipped = [], []
    for run, comp, cd in iter_comparisons(root, runs, keep):
        prefix = f"{run}_{comp}"
        for reg in regions:
            tag, mf = METH[reg]
            if not (cd / mf).exists():
                skipped.append(f"{run}/{comp} [{reg}] (no {Path(mf).name})")
            elif (cd / 'DMRs/PCA' / f"{prefix}_{tag}_PCA.pdf").exists() and not args.force:
                skipped.append(f"{run}/{comp} [{reg}] (done; --force to redo)")
            else:
                rows.append((run, comp, reg, str(cd), prefix))
    if not rows:
        sys.exit("Nothing to do.\n  " + '\n  '.join(skipped))

    (out / 'logs').mkdir(parents=True, exist_ok=True)
    man = out / 'pca_manifest.tsv'
    with open(man, 'w') as f:
        f.write('task\trun\tcomparison\tregions\tcomparison_dir\tprefix\n')
        for i, r in enumerate(rows, 1):
            f.write(f"{i}\t" + '\t'.join(r) + '\n')

    q = lambda p: shlex.quote(str(p))
    opts = f"--pcs {args.pcs}"
    opts += f" --shape {q(args.shape)}" if args.shape else ''
    opts += f" --colour {q(args.colour)}" if args.colour else ''
    opts += ' --label' if args.label else ''
    opts += ' --show_p' if args.show_p else ''
    opts += ' --scale FALSE' if args.no_scale else ''
    mem = args.mem or ('16G' if args.regions == 'dmr' else '64G')
    collect = (f'mkdir -p {q(out / "pdf")}\n'
               f'TAG=$([ "$REG" = dmr ] && echo DMR || echo background)\n'
               f'cp "$CDIR/DMRs/PCA/${{PREFIX}}_${{TAG}}_PCA.pdf" {q(out / "pdf")}/\n'
               ) if args.collect else ''
    (out / 'run_pca.sbatch').write_text(f"""#!/bin/bash
#SBATCH --job-name=dmr_pca
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem={mem}
#SBATCH --time={args.time}
#SBATCH --array=1-{len(rows)}
#SBATCH --output={out}/logs/pca_%A_%a.log
# Written by 06_make_pca_jobs.py: one array task per row of pca_manifest.tsv
set -euo pipefail

LINE=$(awk -F'\\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t' {q(man)})
[ -n "$LINE" ] || {{ echo "No task $SLURM_ARRAY_TASK_ID"; exit 1; }}
IFS=$'\\t' read -r TASK RUN COMP REG CDIR PREFIX <<< "$LINE"
echo "Task $TASK: $RUN / $COMP [$REG]  ($(date))"

eval "$(pixi shell-hook --manifest-path {q(pixi)})"
Rscript {q(rs)} --comparison_dir "$CDIR" --prefix "$PREFIX" --regions "$REG" {opts}
{collect}echo "Done ($(date))"
""")

    print(f"{len(rows)} PCA tasks:")
    for i, r in enumerate(rows, 1):
        print(f"  {i:>2}  {r[0]}/{r[1]}  [{r[2]}]")
    if skipped:
        print("Skipped:\n  " + '\n  '.join(skipped))
    rel = os.path.relpath(out, Path.cwd())
    print(f"\nWrote {rel}/pca_manifest.tsv and run_pca.sbatch\nSubmit with:\n  sbatch {rel}/run_pca.sbatch")


if __name__ == '__main__':
    main()
