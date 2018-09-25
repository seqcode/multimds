library("edgeR")

args = commandArgs(trailingOnly=TRUE)

mat = read.delim(args[1], row.names="Symbol")
group = factor(c(1,2,2))
dge = DGEList(counts=mat, group=group)
nf = calcNormFactors(dge)
design = model.matrix(~group)
disp = estimateDisp(nf, design)
fit = glmQLFit(disp, design)
qlf = glmQLFTest(fit, coef=2)

write.table(qlf$table, args[2])
