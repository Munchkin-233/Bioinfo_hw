library(GEOquery)
gds  <- getGEO("GDS2901", destdir = ".")   # downloads and caches GDS2901.soft.gz
eset <- GDS2eSet(gds, do.log2 = FALSE)     # ExpressionSet: exprs(eset) is probes × samples
pData(eset)                                # genotype / agent / time per sample


