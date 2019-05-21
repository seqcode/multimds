library("edgeR")

counts = read.delim("nup60_counts.tsv", row.names="Symbol")
bcv = 0.1
y = DGEList(counts=counts, group=1:2)
et = exactTest(y, dispersion=bcv^2)
et$table$q = p.adjust(et$table$PValue, method="BH")
write.table(et$table, "nup60_edgeR_results.tsv")
