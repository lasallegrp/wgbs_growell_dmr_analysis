#!/usr/bin/env python3
"""
Step 8: volcano plots (08_volcano_plot.R) for every finished outcome/timepoint.

Layout (GROWELL):  <run_root>/<outcome>/<timepoint>/  e.g. results/dmrichr_output/ppwr/36-38wk/
  "run" = outcome folder; "comparison" = timepoint folder (must hold sample_info.xlsx).

Lists timepoints that have DMRs/background_annotated.xlsx and writes, in --out_dir
(default <run_root>/jobs/volcano/<outcomes joined by __ | all>/):
  volcano_manifest.tsv   one row per outcome x timepoint: task, run, comparison, background, outdir, prefix
  run_volcano.sbatch     SLURM array job, one task per row

Each plot: <outcome>/<timepoint>/<--out_subdir>/<outcome>_<timepoint>_volcano_full.pdf/.png
--out_subdir defaults to DMRs/volcano so the older <timepoint>/volcano/ plots are left alone;
pass --out_subdir volcano to write into that folder instead.
With --collect, the job also copies every PNG into <out_dir>/png/ so you can flip through them.

Usage (from anywhere inside the repo):
  python3 scripts/dmrichr/08_make_volcano_jobs.py
  python3 scripts/dmrichr/08_make_volcano_jobs.py --runs ppwr,egwg --collect
  python3 scripts/dmrichr/08_make_volcano_jobs.py --runs ppwr --sig_by auto --n_labels 30
  sbatch results/dmrichr_output/jobs/volcano/all/run_volcano.sbatch

--sig_by pvalue (default) colours every region with p < 0.05. If DMRichR called this
timepoint's DMRs by q (>= 100 regions with q < 0.05), the coloured count will be larger
than the DMRs in DMRs.bed; --sig_by auto colours exactly DMRichR's DMRs instead.
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
    ap.add_argument('--out_dir', default=None, help='job folder; default <run_root>/jobs/volcano/<outcomes | all>')
    ap.add_argument('--out_subdir', default='DMRs/volcano', help='plot folder inside each timepoint [default %(default)s]')
    ap.add_argument('--scripts_dir', default=str(SCRIPT_DIR), help="folder holding 08_volcano_plot.R [default: this script's folder]")
    ap.add_argument('--pixi_manifest', default=None, help=f'pixi.toml [default: first of {", ".join(PIXI_CANDIDATES)} in the repo]')
    ap.add_argument('--sig_by', default='pvalue', choices=['pvalue', 'qvalue', 'auto'],
                    help='pvalue (default) = colour regions with p < --cutoff; qvalue = q < --cutoff; '
                         'auto = DMRichR rule (q if >= 100 regions pass q, else p)')
    ap.add_argument('--cutoff', type=float, default=0.05)
    ap.add_argument('--n_labels', type=int, default=20)
    ap.add_argument('--fig_style', default='default', choices=['default', 'paper', 'poster'],
                    help='figure sizing passed to 08_volcano_plot.R [default %(default)s]')
    ap.add_argument('--collect', action='store_true', help='also copy all PNGs into <out_dir>/png/')
    ap.add_argument('--force', action='store_true', help='redo timepoints that already have a volcano plot')
    ap.add_argument('--mem', default='16G')
    ap.add_argument('--time', default='02:00:00')
    args = ap.parse_args()

    runs = [r.strip() for r in args.runs.split(',') if r.strip()] if args.runs else None
    keep = set(c.strip() for c in args.comparisons.split(',')) if args.comparisons else None
    root = absol(args.run_root)
    if not root.is_dir():
        sys.exit(f"--run_root not found: {root}")
    out = absol(args.out_dir) if args.out_dir else root / 'jobs' / 'volcano' / ('__'.join(runs) if runs else 'all')
    rs = absol(args.scripts_dir) / '08_volcano_plot.R'
    if not rs.exists():
        sys.exit(f"R script not found: {rs} (use --scripts_dir)")
    pixi = find_pixi(args.pixi_manifest)

    rows, missing, done = [], [], []
    for run, comp, cd in iter_comparisons(root, runs, keep):
        prefix = f"{run}_{comp}"
        bg = cd / 'DMRs' / 'background_annotated.xlsx'
        if not bg.exists():
            missing.append(f"{run}/{comp}")
            continue
        od = cd / args.out_subdir
        if (od / f"{prefix}_volcano_full.pdf").exists() and not args.force:
            done.append(f"{run}/{comp}")
            continue
        rows.append((run, comp, str(bg), str(od), prefix))
    if not rows:
        msg = "Nothing to do."
        if done:
            msg += f" Already plotted (use --force to redo): {', '.join(done)}"
        if missing:
            msg += f" Not finished: {', '.join(missing)}"
        sys.exit(msg)

    (out / 'logs').mkdir(parents=True, exist_ok=True)
    man = out / 'volcano_manifest.tsv'
    with open(man, 'w') as f:
        f.write('task\trun\tcomparison\tbackground\toutdir\tprefix\n')
        for i, r in enumerate(rows, 1):
            f.write(f"{i}\t" + '\t'.join(r) + '\n')

    q = lambda p: shlex.quote(str(p))
    collect = (f'\nmkdir -p {q(out / "png")}\n'
               f'cp "$OUTDIR/${{PREFIX}}_volcano_full.png" {q(out / "png")}/\n') if args.collect else ''
    sb = f"""#!/bin/bash
#SBATCH --job-name=dmr_volcano
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --array=1-{len(rows)}
#SBATCH --output={out}/logs/volcano_%A_%a.log
# Written by 08_make_volcano_jobs.py: one array task per row of volcano_manifest.tsv
set -euo pipefail

MANIFEST={q(man)}
LINE=$(awk -F'\\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t' "$MANIFEST")
[ -n "$LINE" ] || {{ echo "No manifest row for task $SLURM_ARRAY_TASK_ID"; exit 1; }}
IFS=$'\\t' read -r TASK RUN COMP BACKGROUND OUTDIR PREFIX <<< "$LINE"
echo "Task $TASK: $RUN / $COMP  ($(date))"

eval "$(pixi shell-hook --manifest-path {q(pixi)})"
for p in optparse readxl ggplot2 ggrepel dplyr; do
  Rscript -e "q(status = !requireNamespace('$p', quietly = TRUE))" >/dev/null 2>&1 \\
    || {{ echo "R package $p is missing from the pixi env"; exit 1; }}
done

mkdir -p "$OUTDIR"
Rscript {q(rs)} \\
  --background "$BACKGROUND" \\
  --outdir "$OUTDIR" \\
  --prefix "$PREFIX" \\
  --sig_by {args.sig_by} \\
  --cutoff {args.cutoff} \\
  --n_labels {args.n_labels} \\
  --fig_style {args.fig_style}
{collect}
echo "Done: $OUTDIR/${{PREFIX}}_volcano_full.pdf"
"""
    (out / 'run_volcano.sbatch').write_text(sb)

    print(f"{len(rows)} timepoints to plot:")
    for i, r in enumerate(rows, 1):
        print(f"  {i:>2}  {r[0]}/{r[1]}")
    if done:
        print("Already plotted, skipped (--force to redo): " + ', '.join(done))
    if missing:
        print("Not finished (no background_annotated.xlsx), skipped: " + ', '.join(missing))
    rel = os.path.relpath(out, Path.cwd())
    print(f"\nWrote {rel}/volcano_manifest.tsv and run_volcano.sbatch\nSubmit with:\n  sbatch {rel}/run_volcano.sbatch")


if __name__ == '__main__':
    main()
