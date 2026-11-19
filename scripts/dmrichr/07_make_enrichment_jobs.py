#!/usr/bin/env python3
"""
Step 7: enrichment (07_enrichment_plot.R + helper.R) for every finished outcome/timepoint.

Layout (GROWELL):  <run_root>/<outcome>/<timepoint>/  e.g. results/dmrichr_output/ppwr/36-38wk/
  "run" = outcome folder; "comparison" = timepoint folder (must hold sample_info.xlsx).

Writes, in --out_dir (default <run_root>/jobs/enrichment/<outcomes joined by __ | all>/):
  enrichment_manifest.tsv   one row per outcome x timepoint to run
  run_enrichment.sbatch     SLURM array job, one task per row

Each task: --dmrs DMRs/DMRs_annotated.xlsx --background_regions DMRs/backgroundRegions.bed
-> <outcome>/<timepoint>/<--out_subdir>/ with prefix <outcome>_<timepoint> (e.g. ppwr_36-38wk)
   (GREAT + clusterProfiler + Enrichr; log in enrichment.log)
--out_subdir defaults to DMRs/enrichment so the older <timepoint>/enrichment/ results are left alone;
pass --out_subdir enrichment to write into that folder instead.

Usage (from anywhere inside the repo):
  python3 scripts/dmrichr/07_make_enrichment_jobs.py
  python3 scripts/dmrichr/07_make_enrichment_jobs.py --runs ppwr,egwg --min_dmrs 100
  python3 scripts/dmrichr/07_make_enrichment_jobs.py --runs ppwr --comparisons 36-38wk --force \\
      --great_gene_table /quobyte/lasallegrp/projects/LACOFD/reference/GREATv4.genes.hg38.tsv.gz
  sbatch results/dmrichr_output/jobs/enrichment/all/run_enrichment.sbatch
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



GENE_SETS = '/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/reference/gene_sets_hg38.rds'


def count_bed(f):
    with open(f) as fh:
        return sum(1 for l in fh if l.strip() and not l.startswith(('#', 'track', 'browser')))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', default=None, help='comma list of outcomes, e.g. ppwr,egwg (default: all found)')
    ap.add_argument('--run_root', default='results/dmrichr_output')
    ap.add_argument('--comparisons', default=None, help='comma list of timepoints, e.g. Baseline,36-38wk (default: all)')
    ap.add_argument('--out_dir', default=None, help='job folder; default <run_root>/jobs/enrichment/<outcomes | all>')
    ap.add_argument('--out_subdir', default='DMRs/enrichment', help='results folder inside each timepoint [default %(default)s]')
    ap.add_argument('--min_dmrs', type=int, default=50, help='skip timepoints with fewer DMRs [default 50]')
    ap.add_argument('--gene_set_cache', default=GENE_SETS)
    ap.add_argument('--scripts_dir', default=str(SCRIPT_DIR), help="folder holding 07_enrichment_plot.R and helper.R [default: this script's folder]")
    ap.add_argument('--pixi_manifest', default=None, help=f'pixi.toml [default: first of {", ".join(PIXI_CANDIDATES)} in the repo]')
    ap.add_argument('--great_gene_table', default='', help='local GREATv4.genes.hg38.tsv.gz, if rGREAT cannot download it')
    ap.add_argument('--enrich_args', default='', help='extra options for 07_enrichment_plot.R, quoted, e.g. "--mode presentation"')
    ap.add_argument('--force', action='store_true', help='redo timepoints whose --out_subdir is not empty')
    ap.add_argument('--mem', default='48G')
    ap.add_argument('--cpus', type=int, default=4)
    ap.add_argument('--time', default='24:00:00')
    args = ap.parse_args()

    runs = [r.strip() for r in args.runs.split(',') if r.strip()] if args.runs else None
    keep = set(c.strip() for c in args.comparisons.split(',')) if args.comparisons else None
    root = absol(args.run_root)
    if not root.is_dir():
        sys.exit(f"--run_root not found: {root}")
    out = absol(args.out_dir) if args.out_dir else root / 'jobs' / 'enrichment' / ('__'.join(runs) if runs else 'all')
    sd = absol(args.scripts_dir)
    rs, helper = sd / '07_enrichment_plot.R', sd / 'helper.R'
    for f in [rs, helper]:
        if not f.exists():
            sys.exit(f"Not found: {f} (use --scripts_dir)")
    if not Path(args.gene_set_cache).exists():
        print(f"WARNING: gene-set cache not found: {args.gene_set_cache} (GREAT will be skipped)")
    pixi = find_pixi(args.pixi_manifest)

    rows, skipped = [], []
    for run, comp, cd in iter_comparisons(root, runs, keep):
        dmrs, bg, bed = cd / 'DMRs/DMRs_annotated.xlsx', cd / 'DMRs/backgroundRegions.bed', cd / 'DMRs/DMRs.bed'
        if not (dmrs.exists() and bg.exists() and bed.exists()):
            skipped.append(f"{run}/{comp} (not finished: needs DMRs_annotated.xlsx, backgroundRegions.bed, DMRs.bed)")
            continue
        n = count_bed(bed)
        if n < args.min_dmrs:
            skipped.append(f"{run}/{comp} ({n} DMRs < {args.min_dmrs})")
            continue
        od = cd / args.out_subdir
        if od.exists() and any(od.iterdir()) and not args.force:
            skipped.append(f"{run}/{comp} (done: {args.out_subdir}/ not empty; --force to redo)")
            continue
        rows.append((run, comp, n, str(dmrs), str(bg), str(od), f"{run}_{comp}"))
    if not rows:
        sys.exit("Nothing to do.\n  " + '\n  '.join(skipped))

    (out / 'logs').mkdir(parents=True, exist_ok=True)
    man = out / 'enrichment_manifest.tsv'
    with open(man, 'w') as f:
        f.write('task\trun\tcomparison\tn_dmrs\tdmrs\tbackground\toutdir\tprefix\n')
        for i, r in enumerate(rows, 1):
            f.write(f"{i}\t" + '\t'.join(map(str, r)) + '\n')

    q = lambda p: shlex.quote(str(p))
    (out / 'run_enrichment.sbatch').write_text(f"""#!/bin/bash
#SBATCH --job-name=dmr_enrich
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --cpus-per-task={args.cpus}
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --array=1-{len(rows)}
#SBATCH --output={out}/logs/enrich_%A_%a.log
# Written by 07_make_enrichment_jobs.py: one array task per row of enrichment_manifest.tsv
set -euo pipefail

LINE=$(awk -F'\\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t' {q(man)})
[ -n "$LINE" ] || {{ echo "No task $SLURM_ARRAY_TASK_ID"; exit 1; }}
IFS=$'\\t' read -r TASK RUN COMP N DMRS BG OUT PREFIX <<< "$LINE"
echo "Task $TASK: $RUN / $COMP, $N DMRs  ($(date))"

eval "$(pixi shell-hook --manifest-path {q(pixi)})"
export GREAT_GENE_TABLE={q(args.great_gene_table)}
mkdir -p "$OUT"
Rscript {q(rs)} \\
  --dmrs "$DMRS" \\
  --background_regions "$BG" \\
  --gene_set_cache {q(args.gene_set_cache)} \\
  --helper_file {q(helper)} \\
  --outdir "$OUT" \\
  --prefix "$PREFIX" {args.enrich_args} 2>&1 | tee "$OUT/enrichment.log"
echo "Done: $OUT ($(date))"
""")

    print(f"{len(rows)} timepoints:")
    for i, r in enumerate(rows, 1):
        print(f"  {i:>2}  {r[0]}/{r[1]}  ({r[2]:,} DMRs)")
    if skipped:
        print("Skipped:\n  " + '\n  '.join(skipped))
    rel = os.path.relpath(out, Path.cwd())
    print(f"\nWrote {rel}/enrichment_manifest.tsv and run_enrichment.sbatch\nSubmit with:\n  sbatch {rel}/run_enrichment.sbatch")


if __name__ == '__main__':
    main()
