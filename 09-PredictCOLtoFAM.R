cat("Modelling data...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoFAMs"), showWarnings = FALSE, recursive = TRUE)

# Predict each family with the collection
## Now we need phenos and genos_ready files
## reinitialise order of IDs
ids<-intersect(rownames(phenos), rownames(genos_ready)) %>% sort
NbID<-length(ids)
cat("Number of common genos/phenos: \n")
print(NbID)
phenos<-phenos[ids,]
genos_pred<-genos_ready[ids,]
nrow(phenos) == nrow(genos_pred)
## simple rrBLUP

## list of IDs in collection and family names
families<-c("FjDe", "FuPi", "FjPL", "GDFj", "GaPi", "GaPL")
NbFAM<-length(families)
WhichCOL<-c(1:NbID)[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
summary(rownames(genos_pred)[WhichCOL] == rownames(phenos)[WhichCOL]) ## 232 with same names

traits=colnames(phenos)
accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA)


for (trait in traits) {
  for (fam in families) {
    WhichVS<-c(1:NbID)[grep(fam, ids)]
    res <- mixed.solve(y=phenos[WhichCOL,trait],Z=genos_pred[WhichCOL,])
    Y_VS_pred<- as.vector(genos_pred[WhichVS,] %*% as.matrix(res$u))
    accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete"))
    rm(res, WhichVS, Y_VS_pred)
  }
}

write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs.txt"), quote=F, sep="\t")

## predict within year 2012 (common to all)

phenos_raw

## compare results with ADI model later  






# Purge obsolete variables
rm()

cat("Data modelled!\n")
