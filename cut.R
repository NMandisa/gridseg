X = as.matrix(read.table("C:/Users/Cássio/workspace/darbseg/resources/Y2.csv"))

sel = X[,1] <= 34 & X[,2] <= 250 & X[,2] >= 208

X = X[sel,]

write.table(X, "C:/Users/Cássio/workspace/gridseg/Y2cut.csv", quote = F, row.names = F, col.names = F)
