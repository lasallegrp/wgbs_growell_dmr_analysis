#!/usr/bin/env python3
"""
Step 9: enrichment of what is SHARED between DMR sets (after the overlap step, 05_dmrichr_overlap.R).

GROWELL: overlap folders live in results/dmrichr_output/DMRichR_overlap/<overlap_name>/,
e.g. all_outcomes_gwas-obesity, egwg_ppwr_gwas-obesity. Point --overlap_dir at one of them
(or at DMRichR_overlap/ itself to do all). Workbooks are found at any depth:
  .../<group>/gene/<group>_overlap_gene.xlsx     (sheet Genes)
  .../<group>/coord/<group>_overlap_coord.xlsx   (sheet Loci)
Each set (dir_<set> column) is tied to its DMRichR folder by the nearest manifest.tsv (set_id -> path)
in the workbook's folder or a parent (or --manifest). With no manifest, sets are matched by NAME
to <dmrichr_root>/<outcome>/<timepoint>/: <outcome>_<timepoint>, or just the timepoint inside an
outcome group (within_outcome_across_timepoints/egwg), or just the outcome inside a timepoint group
(within_timepoint_across_outcomes/Baseline); "-" and "_" count as the same (36_38wk = 36-38wk).
The match is printed and saved as resolved_sets.tsv; unmatched sets (e.g. GWAS) get no background.

--mode gene  09_shared_genes.R: genes with a DMR in all sets of the group; background = genes
  testable in ALL DMRichR sets (geneSymbol of DMRs/background_annotated.xlsx).
  07_enrichment_plot.R --genes ... --gene_assignment annotation -> clusterProfiler + Enrichr.
  Output: <group>/gene/shared_gene_enrichment/
--mode coord 09_shared_loci.R: loci where DMRs of all sets overlap; background = territory tested
  in ALL DMRichR sets (intersection of backgroundRegions.bed).
  07_enrichment_plot.R --dmrs ... --background_regions ... --gene_assignment great -> GREAT + gene-based.
  Output: <group>/coord/shared_loci_enrichment/
--mode both (default) = one array task per group per mode.

Sets without a DMRichR background (e.g. a GWAS obesity gene list) are left out of the background
intersection: name them with --no_background_sets, or give the manifest a method column (only
method == dmrichr sets carry a background). They still count toward "shared". Note that the
shared list is then pre-selected for GWAS genes, so read terms in that light.

Usage (from anywhere inside the repo):
  python3 scripts/dmrichr/09_make_shared_jobs.py --overlap_dir results/dmrichr_output/DMRichR_overlap/egwg_ppwr_gwas-obesity
  python3 scripts/dmrichr/09_make_shared_jobs.py --overlap_dir results/dmrichr_output/DMRichR_overlap \\
      --mode coord --concordant --split_direction --no_background_sets gwas_obesity
  sbatch results/dmrichr_output/DMRichR_overlap/egwg_ppwr_gwas-obesity/shared_enrichment/run_shared_enrichment.sbatch
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



GENE_SETS = '/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/reference/gene_sets_hg38.rds'
SUB = {'gene': 'shared_gene_enrichment', 'coord': 'shared_loci_enrichment'}
GREAT_TABLE = '/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/reference/GREATv4.genes.hg38.tsv.gz'

def nearest_manifest(x, top):
    for d in [x.parent, *x.parent.parents]:
        if (d / 'manifest.tsv').exists():
            return d / 'manifest.tsv'
        if d == top:
            break
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--overlap_dir', required=True, help='an overlap folder, or DMRichR_overlap/ for all of them')
    ap.add_argument('--manifest', default=None, help='manifest.tsv to use for every group [default: nearest one]')
    ap.add_argument('--dmrichr_root', default='results/dmrichr_output', help='for name matching when there is no manifest [default %(default)s]')
    ap.add_argument('--groups', default=None, help='comma list of group names (workbook name minus _overlap_<mode>.xlsx) [default: all]')
    ap.add_argument('--mode', default='both', choices=['gene', 'coord', 'both'])
    ap.add_argument('--shared', default='all', choices=['all', 'any2'],
                    help='all = in every set of the group; any2 = in at least two [default all]')
    ap.add_argument('--concordant', action='store_true', help='same direction in every set only')
    ap.add_argument('--split_direction', action='store_true', help='also run shared-hyper and shared-hypo lists')
    ap.add_argument('--no_background_sets', default='', help='comma list of set_ids with no DMRichR background (e.g. GWAS lists)')
    ap.add_argument('--min_items', '--min_genes', type=int, default=20,
                    help='lists with fewer genes/loci are not tested [default 20]')
    ap.add_argument('--gene_set_cache', default=GENE_SETS)
    ap.add_argument('--great_gene_table', default=GREAT_TABLE, help='local GREATv4.genes.hg38.tsv(.gz), needed for GREAT on compute nodes')
    ap.add_argument('--scripts_dir', default=str(SCRIPT_DIR), help="folder holding the 07/09 R scripts and helper.R [default: this script's folder]")
    ap.add_argument('--pixi_manifest', default=None, help=f'pixi.toml [default: first of {", ".join(PIXI_CANDIDATES)} in the repo]')
    ap.add_argument('--enrich_args', default='', help='extra options for 07_enrichment_plot.R, quoted')
    ap.add_argument('--out_tag', default='', help='suffix for output folders, e.g. poster -> shared_gene_enrichment_poster')
    ap.add_argument('--force', action='store_true', help='redo groups that already have results')
    ap.add_argument('--mem', default='32G')
    ap.add_argument('--time', default='12:00:00')
    args = ap.parse_args()

    tag = f"_{args.out_tag}" if args.out_tag else ''

    od = absol(args.overlap_dir)
    if not od.is_dir():
        sys.exit(f"--overlap_dir not found: {od}")
    fixed_man = absol(args.manifest) if args.manifest else None
    if fixed_man and not fixed_man.exists():
        sys.exit(f"--manifest not found: {fixed_man}")
    keep = set(g.strip() for g in args.groups.split(',')) if args.groups else None
    sd = absol(args.scripts_dir)
    prep = {'gene': sd / '09_shared_genes.R', 'coord': sd / '09_shared_loci.R'}
    enrich, helper = sd / '07_enrichment_plot.R', sd / 'helper.R'
    for f in [*prep.values(), enrich, helper]:
        if not f.exists():
            sys.exit(f"Not found: {f} (use --scripts_dir)")
    if not Path(args.gene_set_cache).exists():
        print(f"WARNING: gene-set cache not found: {args.gene_set_cache}")
    pixi = find_pixi(args.pixi_manifest)

    droot = absol(args.dmrichr_root)
    rows, done = [], []
    modes = ['gene', 'coord'] if args.mode == 'both' else [args.mode]
    for mode in modes:
        for x in sorted(od.rglob(f'*_overlap_{mode}.xlsx')):
            if x.name.startswith('~$') or any(p.name.startswith(tuple(SUB.values())) for p in x.parents):
                continue
            group = x.name[:-len(f'_overlap_{mode}.xlsx')]
            if keep and group not in keep:
                continue
            man = fixed_man or nearest_manifest(x, od) or 'auto'
            out = x.parent / (SUB[mode] + tag)
            if (out / 'lists.tsv').exists() and not args.force:
                done.append(f"{group} ({mode})")
                continue
            rows.append((group, mode, str(x), str(out), str(man)))
    if any(r[4] == 'auto' for r in rows) and not droot.is_dir():
        sys.exit(f"No manifest.tsv, and --dmrichr_root not found for name matching: {droot}")
    if not rows:
        sys.exit("Nothing to do." + (f" Already done (--force to redo): {', '.join(done)}" if done else
                                     f" No *_overlap_<gene|coord>.xlsx under {od} for --mode {args.mode}"))

    jd = od / ('shared_enrichment' + tag)
    (jd / 'logs').mkdir(parents=True, exist_ok=True)
    tasks = jd / 'tasks.tsv'
    with open(tasks, 'w') as f:
        f.write('task\tgroup\tmode\toverlap_xlsx\tout_dir\tmanifest\n')
        for i, r in enumerate(rows, 1):
            f.write(f"{i}\t" + '\t'.join(r) + '\n')

    q = lambda p: shlex.quote(str(p))
    prep_opts = f"--shared {args.shared}"
    prep_opts += ' --concordant' if args.concordant else ''
    prep_opts += ' --split_direction' if args.split_direction else ''
    prep_opts += f" --no_background_sets {q(args.no_background_sets)}" if args.no_background_sets else ''
    sb = f"""#!/bin/bash
#SBATCH --job-name=shared_enrich
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=2
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --array=1-{len(rows)}
#SBATCH --output={jd}/logs/shared_%A_%a.log
# Written by 09_make_shared_jobs.py: one array task per overlap group in tasks.tsv
set -uo pipefail

LINE=$(awk -F'\\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t' {q(tasks)})
[ -n "$LINE" ] || {{ echo "No task $SLURM_ARRAY_TASK_ID"; exit 1; }}
IFS=$'\\t' read -r TASK GROUP MODE XLSX OUT MAN <<< "$LINE"
echo "Task $TASK: $GROUP, $MODE overlap  ($(date))"
eval "$(pixi shell-hook --manifest-path {q(pixi)})"
export GREAT_GENE_TABLE={q(args.great_gene_table) if args.great_gene_table else '"${GREAT_GENE_TABLE:-}"'}

# 1. shared genes/loci + background testable in ALL DMRichR sets of the group
if [ "$MAN" = auto ]; then SETS=(--dmrichr_root {q(droot)}); else SETS=(--manifest "$MAN"); fi
mkdir -p "$OUT"
if [ "$MODE" = gene ]; then
  Rscript {q(prep['gene'])} --overlap_xlsx "$XLSX" "${{SETS[@]}}" --out_dir "$OUT" {prep_opts} --min_genes {args.min_items} || exit 1
else
  Rscript {q(prep['coord'])} --overlap_xlsx "$XLSX" "${{SETS[@]}}" --out_dir "$OUT" {prep_opts} --min_items {args.min_items} || exit 1
fi

# 2. enrichment of each list that is big enough
status=0
while IFS=$'\\t' read -r LIST N RUN; do
  [ "$LIST" = "list" ] && continue
  if [ "$RUN" != "TRUE" ]; then echo "$LIST: $N (< {args.min_items}) -> skipped"; continue; fi
  echo "$LIST: $N -> $OUT/$LIST"
  mkdir -p "$OUT/$LIST"
  if [ "$MODE" = gene ]; then
    INPUT=(--genes "$OUT/${{LIST}}_genes.txt" --background_genes "$OUT/background_genes.txt" --gene_assignment annotation)
  else
    INPUT=(--dmrs "$OUT/${{LIST}}_loci.bed" --background_regions "$OUT/background_regions.bed" --gene_assignment great)
  fi
  Rscript {q(enrich)} "${{INPUT[@]}}" \\
    --gene_set_cache {q(args.gene_set_cache)} \\
    --helper_file {q(helper)} \\
    --outdir "$OUT/$LIST" \\
    --prefix "${{GROUP}}_${{MODE}}_${{LIST}}" {args.enrich_args} > "$OUT/$LIST/enrichment.log" 2>&1 \\
    || {{ echo "  FAILED, see $OUT/$LIST/enrichment.log"; tail -n 5 "$OUT/$LIST/enrichment.log"; status=1; }}
done < "$OUT/lists.tsv"
echo "Done $GROUP ($(date))"
exit $status
"""
    (jd / 'run_shared_enrichment.sbatch').write_text(sb)

    print(f"{len(rows)} tasks:")
    for i, r in enumerate(rows, 1):
        print(f"  {i:>2}  {r[0]}  [{r[1]}]  sets: {'matched by name' if r[4] == 'auto' else os.path.relpath(r[4], od)}")
    if done:
        print("Already done, skipped (--force to redo): " + ', '.join(done))
    rel = os.path.relpath(jd, Path.cwd())
    print(f"\nSelection: genes/loci in {'every set' if args.shared == 'all' else '>= 2 sets'}"
          f"{', same direction' if args.concordant else ''}{'; + hyper/hypo lists' if args.split_direction else ''}"
          f"; background = genes/territory testable in all DMRichR sets of the group")
    print(f"Wrote {rel}/tasks.tsv and run_shared_enrichment.sbatch\nSubmit with:\n  sbatch {rel}/run_shared_enrichment.sbatch")


if __name__ == '__main__':
    main()