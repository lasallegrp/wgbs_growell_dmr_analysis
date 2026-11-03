#!/usr/bin/env bash
# ============================================================================
# merge_comethyl_into_dmrichr.sh
#
# Move 2025_DBS_comethyl_George into 2025_DBS_dmrich_George, each folder into
# the matching place in the DMRichR layout (namespaced as comethyl/).
#
#   DRY_RUN=1 (default)  print what would happen, change nothing
#   DRY_RUN=0            do it
#   REWRITE=1            also rewrite absolute old paths in scripts/config
#                        (a .bak copy of every changed file is kept)
#
# Usage:
#   bash merge_comethyl_into_dmrichr.sh                  # dry run
#   DRY_RUN=0 bash merge_comethyl_into_dmrichr.sh        # move
#   DRY_RUN=0 REWRITE=1 bash merge_comethyl_into_dmrichr.sh
# ============================================================================
set -euo pipefail

OLD="${OLD:-/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George}"
NEW="${NEW:-/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George}"
DRY_RUN="${DRY_RUN:-1}"
REWRITE="${REWRITE:-0}"

# source (inside OLD)          -> destination (inside NEW)
MAP=(
  "scripts_comethyl            scripts/comethyl"
  "config_comethyl             config/comethyl"
  "analysis_comethyl/notebooks analysis/notebooks/comethyl"
  "analysis_comethyl/logs      analysis/logs/comethyl"
  "analysis_comethyl/README.md docs/comethyl/analysis_README.md"
  "data_comethyl/processed     data/processed/comethyl"
  "data_comethyl/metadata      data/metadata/comethyl"
  "reproducibility_comethyl    reproducibility/comethyl"
  "comethyl_output             results/comethyl_output"
  "00_std_err_logs_comethyl    results/logs/comethyl"
  "tests_comethyl              tests/comethyl"
  "README.md                   docs/comethyl/README.md"
  "CHANGELOG.md                docs/comethyl/CHANGELOG.md"
  "CITATION.cff                docs/comethyl/CITATION.cff"
)
# dropped: results/ (only .gitkeep), Rplots.pdf, .DS_Store, analysis_comethyl/renv.lock (empty),
#          LICENSE (kept if it differs from NEW/LICENSE -> docs/comethyl/LICENSE)

run() { if [[ "${DRY_RUN}" == "1" ]]; then echo "  [dry] $*"; else echo "  $*"; "$@"; fi; }

echo "== Pre-flight"
[[ -d "${OLD}" ]] || { echo "ERROR: ${OLD} not found"; exit 1; }
[[ -d "${NEW}" ]] || { echo "ERROR: ${NEW} not found"; exit 1; }
if [[ "$(stat -c %d "${OLD}")" != "$(stat -c %d "$(dirname "${NEW}")")" ]]; then
  echo "WARNING: different filesystems; mv will copy (slow). Continue only if intended."
fi
if command -v squeue >/dev/null && squeue -u "${USER}" -h -o "%Z" | grep -q -e "${OLD}" -e "${NEW}"; then
  echo "ERROR: you have running/pending jobs in one of these folders; wait for them first."; exit 1
fi
if [[ -d "${OLD}/.git" ]]; then
  echo "NOTE: ${OLD} is a git repo. Its .git is NOT moved (left in place)."
  echo "      Record its last commit for provenance: $(git -C "${OLD}" rev-parse --short HEAD 2>/dev/null || echo '?')"
  git -C "${OLD}" status --porcelain | grep -q . && \
    echo "      WARNING: uncommitted changes in ${OLD}; commit/push them first."
fi
if [[ -d "${NEW}/.git" ]] && git -C "${NEW}" status --porcelain | grep -q .; then
  echo "NOTE: ${NEW} has uncommitted changes; consider committing before the merge."
fi
du -sh "${OLD}" 2>/dev/null | sed 's/^/  size of OLD: /'

echo "== Destination clashes"
clash=0
for m in "${MAP[@]}"; do
  read -r src dst <<<"${m}"
  if [[ -e "${NEW}/${dst}" ]]; then echo "  EXISTS: ${NEW}/${dst}"; clash=1; fi
done
[[ "${clash}" == "1" ]] && { echo "ERROR: resolve clashes above first"; exit 1; } || echo "  none"

echo "== Move"
for m in "${MAP[@]}"; do
  read -r src dst <<<"${m}"
  if [[ ! -e "${OLD}/${src}" ]]; then echo "  (missing, skipped) ${src}"; continue; fi
  run mkdir -p "$(dirname "${NEW}/${dst}")"
  run mv "${OLD}/${src}" "${NEW}/${dst}"
done
if [[ -f "${OLD}/LICENSE" ]]; then
  if cmp -s "${OLD}/LICENSE" "${NEW}/LICENSE" 2>/dev/null; then echo "  LICENSE identical, not copied"
  else run mkdir -p "${NEW}/docs/comethyl"; run mv "${OLD}/LICENSE" "${NEW}/docs/comethyl/LICENSE"; fi
fi

echo "== .gitignore (results, data, logs, pixi envs stay out of git)"
for pat in "results/comethyl/" "results/logs/" "data/processed/comethyl/" "data/metadata/comethyl/" \
           "analysis/logs/comethyl/" "reproducibility/comethyl/pixi/.pixi/" "*.out" "*.err" \
           "Rplots.pdf" ".DS_Store"; do
  if ! grep -qxF "${pat}" "${NEW}/.gitignore" 2>/dev/null; then
    if [[ "${DRY_RUN}" == "1" ]]; then echo "  [dry] add ${pat}"; else echo "${pat}" >> "${NEW}/.gitignore"; echo "  add ${pat}"; fi
  fi
done

echo "== Hard-coded paths that point at the old project"
# old absolute path -> new absolute path (longest first)
REW=(
  "${OLD}/scripts_comethyl|${NEW}/scripts/comethyl"
  "${OLD}/scripts|${NEW}/scripts/comethyl"
  "${OLD}/config_comethyl|${NEW}/analysis/config/comethyl"
  "${OLD}/config|${NEW}/analysis/config/comethyl"
  "${OLD}/data_comethyl/processed|${NEW}/data/processed/comethyl"
  "${OLD}/data/processed|${NEW}/data/processed/comethyl"
  "${OLD}/data_comethyl/metadata|${NEW}/data/metadata/comethyl"
  "${OLD}/data/metadata|${NEW}/data/metadata/comethyl"
  "${OLD}/reproducibility_comethyl|${NEW}/reproducibility/comethyl"
  "${OLD}/reproducibility|${NEW}/reproducibility/comethyl"
  "${OLD}/comethyl_output|${NEW}/results/comethyl_output"
  "${OLD}/00_std_err_logs_comethyl|${NEW}/results/logs/comethyl"
  "${OLD}/00_std_err_logs|${NEW}/results/logs/comethyl"
  "${OLD}/tests_comethyl|${NEW}/tests/comethyl"
)
SEARCH_DIRS=("${NEW}/scripts/comethyl" "${NEW}/analysis/config/comethyl" "${NEW}/reproducibility/comethyl"
             "${NEW}/tests/comethyl" "${NEW}/scripts/DMRichr")
[[ "${DRY_RUN}" == "1" ]] && SEARCH_DIRS=("${OLD}" "${NEW}/scripts/DMRichr")
files=$(grep -rlI --exclude-dir=.pixi --exclude-dir=.git --exclude='*.bak' -e "${OLD}" \
        "${SEARCH_DIRS[@]}" 2>/dev/null || true)
if [[ -z "${files}" ]]; then
  echo "  none found"
else
  echo "${files}" | sed 's/^/  /'
  echo "  -> $(echo "${files}" | wc -l) file(s) reference ${OLD}"
  if [[ "${REWRITE}" == "1" && "${DRY_RUN}" == "0" ]]; then
    for f in ${files}; do
      cp -p "${f}" "${f}.bak"
      for r in "${REW[@]}"; do
        from="${r%%|*}"; to="${r##*|}"
        sed -i "s#${from}\([/\"' ]\|\$\)#${to}\1#g" "${f}"
      done
      sed -i "s#${OLD}#${NEW}#g" "${f}"                     # anything left: project root
      echo "  rewrote ${f} (backup ${f}.bak)"
    done
  else
    echo "  (run with DRY_RUN=0 REWRITE=1 to rewrite them; review the .bak diffs after)"
  fi
fi

echo "== Remaining in ${OLD}"
if [[ "${DRY_RUN}" == "1" ]]; then echo "  [dry] (shown after a real run)"; else ls -A "${OLD}" | sed 's/^/  /'; fi

cat <<EOF

== Next steps (manual)
  1. Rebuild the comethyl pixi env at its new location (conda envs hard-code their path):
       cd ${NEW}/reproducibility/comethyl/pixi && pixi install
  2. Check relative paths inside scripts/comethyl (e.g. "config/", "comethyl_output/"):
       grep -rnE '"(\\.\\./)*(config|comethyl_output|data|scripts)/' ${NEW}/scripts/comethyl
  3. Run the comethyl smoke test:  bash ${NEW}/tests/comethyl/smoke_test.sh
  4. Commit in ${NEW}:
       git add -A && git commit -m "Import comethyl pipeline from 2025_DBS_comethyl_George@<commit>"
  5. When everything runs, archive the old GitHub repo (README pointer to the merged repo)
     and remove what is left of ${OLD}.
EOF