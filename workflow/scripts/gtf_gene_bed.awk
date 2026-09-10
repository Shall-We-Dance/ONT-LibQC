# Merge each gene's annotated extent; parse attributes by key, independent of order.
BEGIN {
    FS = OFS = "\t"
    n = split(primary, names, ",")
    for (i = 1; i <= n; i++) keep[names[i]] = 1
}
($3 == "gene" || $3 == "transcript") && ($1 in keep) {
    gene_id = ""
    nattr = split($9, attributes, ";")
    for (a = 1; a <= nattr; a++) {
        value = attributes[a]
        sub(/^[ \t]+/, "", value)
        if (value ~ /^gene_id[ \t]+/) {
            sub(/^gene_id[ \t]+/, "", value)
            gsub(/"/, "", value)
            sub(/[ \t]+$/, "", value)
            gene_id = value
        }
    }
    if (gene_id == "") next
    key = $1 SUBSEP gene_id SUBSEP $7
    if (!(key in left) || $4 < left[key]) left[key] = $4
    if (!(key in right) || $5 > right[key]) right[key] = $5
    chrom[key] = $1
    gene[key] = gene_id
    strand[key] = $7
}
END {
    for (key in left) print chrom[key], left[key] - 1, right[key], gene[key], 0, strand[key]
}
