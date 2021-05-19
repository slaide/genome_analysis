#replace first dot in filenames with underscore for metabat compatibility (no dots in filename!)
for f in *.*.fa; do mv -- "$f" "${f/./_}"; done
