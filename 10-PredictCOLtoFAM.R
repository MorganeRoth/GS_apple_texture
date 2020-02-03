cat("Predicting families with the collection...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoFAMs"), showWarnings = FALSE, recursive = TRUE)

######################################
## LOAD DATA IF YOU START FROM HERE ##
######################################
phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))

## load clusters from DAPC analysis
clusters<-read.table(paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"), h=T)
dim(clusters)
cluster_fams<-read.table(paste0(odir, "/genos_modelled/assignements_families.txt"), h=T)
## useful lists
families<-c("FjDe", "GDFj", "FjPi", "FjPL", "GaPi", "GaPL")
cat("Families:")
print(families)

#############################################
## Predict each family with the collection ##
## Without clustering (model A) #############
#############################################

## reinitialise order of IDs, prepare data and result matrix
ids<-intersect(rownames(phenos), rownames(genos_ready)) %>% sort
NbID<-length(ids)
cat("Number of common genos/phenos: \n")
print(NbID)
phenos<-phenos[ids,]
genos_ready<-genos_ready[ids,]
nrow(phenos) == nrow(genos_ready)
## list of IDs in collection and family names
NbFAM<-length(families)
WhichCOL<-c(1:NbID)[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
summary(rownames(genos_ready)[WhichCOL] == rownames(phenos)[WhichCOL]) ## 259 with same names
traits=colnames(phenos)
accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA, conf_low=NA, conf_high=NA)

for (trait in traits) {
  for (fam in families) {
    WhichVS<-c(1:NbID)[grep(fam, ids)]
    res <- mixed.solve(y=phenos[WhichCOL,trait],Z=genos_ready[WhichCOL,])
    Y_VS_pred<- as.vector(genos_ready[WhichVS,] %*% as.matrix(res$u))
    acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
    cor.test(phenos[WhichVS, trait], Y_VS_pred, conf.level = 0.95)
    conf_int<-cor.test(phenos[WhichVS, trait], Y_VS_pred, conf.level = 0.95)$conf.int[c(1,2)]
    print(c(trait, fam, acc, conf_int))
    accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, acc, conf_int)
    rm(res, WhichVS, Y_VS_pred)
  }
}
write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs_reps.txt"), quote=F, sep="\t")
rm(accuracy)

## plot observed vs predicted for each trait

for (trait in traits){
  pdf(file=paste0(odir, "/predictions/COLLtoFAMs/", trait, ".pdf"), height=6, width=10)
  par(mfrow=c(2,3), mar=rep(4,4,4,2))
  for (fam in families) {
    WhichVS<-c(1:NbID)[grep(fam, ids)]
    res <- mixed.solve(y=phenos[WhichCOL,trait],Z=genos_ready[WhichCOL,])
    Y_VS_pred<- as.vector(genos_ready[WhichVS,] %*% as.matrix(res$u))
    myplot=plot(Y_VS_pred~phenos[WhichVS, trait], xlab="Y-observed", ylab="Y-predicted")
    print(myplot)
    print(title(fam, adj = 0.5, line = 1))
    mylm=lm(Y_VS_pred~phenos[WhichVS, trait])
    abline(col="red", a=coefficients(mylm)[1], b=coefficients(mylm)[2])
  }
  title(trait, outer = TRUE, line=-1, cex = 2)
  dev.off()
}

#############################################
## Predict each family with the collection ##
## With clustering (model B) ################
#############################################

clusters<-matrix(clusters[,"Cluster"], nrow=nrow(clusters), ncol=1, dimnames = list(clusters[,"Name"], "Cluster"))
rownames(clusters)
clusters<-clusters[rownames(genos_ready)[WhichCOL],"Cluster"] %>% as.matrix
rownames(clusters)<-rownames(genos_ready)[WhichCOL]
clusters<-class.ind(clusters)
accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA, conf_low=NA, conf_high=NA)
cluster_fams$cluster<-factor(cluster_fams$cluster, levels=1:6)
rownames(cluster_fams)<-cluster_fams$Name

for (trait in traits) {
  for (fam in families) {
    WhichVS<-c(1:NbID)[grep(fam, ids)]
    cluster_FAM<- cluster_fams[ids[grep(fam,ids)], "cluster" ] %>% class.ind
    print(table(cluster_fams[ids[grep(fam,ids)], "cluster" ] ))
    res <- mixed.solve(y=phenos[WhichCOL,trait],Z=genos_ready[WhichCOL,], X=clusters) ## cluster as fix effect
    Y_VS_pred<- as.vector(genos_ready[WhichVS,] %*% as.matrix(res$u) + cluster_FAM %*% res$beta) ## add effect of cluster
    acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
    conf_int<-cor.test(phenos[WhichVS, trait], Y_VS_pred, conf.level = 0.95)$conf.int[c(1,2)]
    print(c(trait, fam, acc, conf_int))
    accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, acc, conf_int)
    rm(res, WhichVS, Y_VS_pred)
  }
}


write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs_clusters_fix_reps.txt"), quote=F, sep="\t")

rm(accuracy)

cat("Families predicted with scenario 1!\n")
