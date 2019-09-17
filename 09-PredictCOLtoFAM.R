cat("Modelling data...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoFAMs"), showWarnings = FALSE, recursive = TRUE)

## if import, model phenos, mnodel genos scripts are skipped, load data here
id_pheno<-read.table(paste0(odir, "/phenos_modelled/rownames_phenos.txt"))
## problem with duplicated rowname that I do not understand, use id_pheno to replace them (saved with phenos)
phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T, row.names = id_pheno$x %>% as.character())
phenos<-phenos[,-1]
summary(phenos)
head(phenos)
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))
dim(phenos)
dim(genos_ready)
## clusters from DAPC analysis
clusters<-read.table(paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"), h=T)
dim(clusters)
cluster_fams<-read.table(paste0(odir, "/genos_modelled/assignements_families.txt"), h=T)
## useful lists
families<-c("FjDe", "GDFj", "FjPi", "FjPL", "GaPi", "GaPL")
# lapply(families, function(x) grep(x, rownames(genos_ready)[rownames(genos_ready) %in% phenos$Name])) %>% unlist %>% length
cat("Families:")
print(families)
# Predict each family with the collection
## Now we need phenos and genos_ready files
## reinitialise order of IDs
ids<-intersect(rownames(phenos), rownames(genos_ready)) %>% sort
NbID<-length(ids)
cat("Number of common genos/phenos: \n")
print(NbID)
phenos<-phenos[ids,]
genos_ready<-genos_ready[ids,]
nrow(phenos) == nrow(genos_ready)
## simple rrBLUP

## list of IDs in collection and family names
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
NbFAM<-length(families)
WhichCOL<-c(1:NbID)[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
summary(rownames(genos_ready)[WhichCOL] == rownames(phenos)[WhichCOL]) ## 242 with same names

traits=colnames(phenos)
accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA)

for (trait in traits) {
  for (fam in families) {
    WhichVS<-c(1:NbID)[grep(fam, ids)]
    res <- mixed.solve(y=phenos[WhichCOL,trait],Z=genos_ready[WhichCOL,])
    Y_VS_pred<- as.vector(genos_ready[WhichVS,] %*% as.matrix(res$u))
    acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
    print(c(trait, fam, acc))
    accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, acc)
    rm(res, WhichVS, Y_VS_pred)
  }
}

write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs.txt"), quote=F, sep="\t")
## plot
accuracy<-read.table(paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs.txt"), h=T)
png(file=paste0(odir, "/predictions/COLLtoFAMs/COLLtoFAM_no_clusters.png"), height=500, width=1000)
ggplot(accuracy,aes( y=accuracy, x=trait))+
  geom_boxplot()+
  facet_grid(~FAM)+
  labs(x="Trait", title="Predictions rrBLUP COLL to FAMs", y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))+
  scale_y_continuous(limits = c(-0.5, 1))
dev.off()

rm(accuracy)
## try with clusters (as fix effect)
head(clusters)
clusters<-matrix(clusters[,"Cluster"], nrow=nrow(clusters), ncol=1, dimnames = list(clusters[,"Name"], "Cluster"))
rownames(clusters)
clusters<-clusters[rownames(genos_ready)[WhichCOL],"Cluster"] %>% as.matrix
rownames(clusters)<-rownames(genos_ready)[WhichCOL]
clusters<-class.ind(clusters)

accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA)

cluster_fams$cluster<-factor(cluster_fams$cluster, levels=1:5)
rownames(cluster_fams)<-cluster_fams$Name

for (trait in traits) {
  for (fam in families) {
    WhichVS<-c(1:NbID)[grep(fam, ids)]
    cluster_FAM<- cluster_fams[ids[grep(fam,ids)], "cluster" ] %>% class.ind
    print(table(cluster_fams[ids[grep(fam,ids)], "cluster" ] ))
    res <- mixed.solve(y=phenos[WhichCOL,trait],Z=genos_ready[WhichCOL,], X=clusters) ## cluster as fix effect
    Y_VS_pred<- as.vector(genos_ready[WhichVS,] %*% as.matrix(res$u) + cluster_FAM %*% res$beta) ## add effect of cluster
    acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
    print(c(trait, fam, acc))
    accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, acc)
    rm(res, WhichVS, Y_VS_pred)
  }
}

write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs_clusters_fix.txt"), quote=F, sep="\t")

accuracy$accuracy<- accuracy$accuracy %>% as.character %>% as.numeric
tail(accuracy)

png(file=paste0(odir, "/predictions/COLLtoFAMs/COLLtoFAM_with_clusters.png"), height=500, width=1000)
ggplot(accuracy,aes( y=accuracy, x=trait))+
  geom_boxplot()+
  facet_grid(~FAM)+
  labs(x="Trait", title="Predictions rrBLUP COLL to FAMs with 5 clusters", y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))+
  scale_y_continuous(limits = c(-0.5, 1))
dev.off()

# 2012 not up to date, re-run with new data if necessary
# ## predict within year 2012 (common to all)
# 
# phenos_raw %>% head
# pheno2012<-phenos_raw[phenos_raw$Year==2012,] %>% droplevels
# summary(pheno2012)
# traits=colnames(pheno2012)[-c(1:7)]
# accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA)
# 
# for (trait in traits) {
#   pheno<-pheno2012[which(!(is.na(pheno2012[,trait]))),] %>% droplevels
#   ids<-intersect(levels(pheno$Name), rownames(genos_ready))
#   geno<-genos_ready[ids,]
#   pheno<-pheno[which(pheno$Name %in% ids), ] %>% droplevels
#   pheno<- aggregate(x=pheno[,trait],by=list(pheno$Name),FUN= mean)
#   rownames(pheno)<-pheno$Group.1
#   WhichCOL<-c(1:nrow(geno))[-c(lapply(families, function(x) grep(x, rownames(geno)) ) %>% unlist)]
#   print(length(WhichCOL))
#   for (fam in families) {
#     WhichVS<-c(1:NbID)[grep(fam, ids)]
#     res <- mixed.solve(y=pheno[rownames(geno)[WhichCOL],"x"],Z=geno[WhichCOL,])
#     Y_VS_pred<- as.vector(geno[WhichVS,] %*% as.matrix(res$u))
#     accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, cor(pheno[rownames(geno)[WhichVS], "x"], Y_VS_pred, use="na.or.complete"))
#     rm(res, WhichVS, Y_VS_pred)
#   }
# }
# 
# write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs_2012only.txt"), sep="\t", quote=F)
# accuracy<-read.table(paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs_2012only.txt"), h=T)
# 
# png(file=paste0(odir, "/predictions/COLLtoFAMs/COLLtoFAM_2012only.png"), height=500, width=1000)
# ggplot(accuracy,aes( y=as.numeric(accuracy), x=trait))+
#   geom_boxplot()+
#   facet_grid(~FAM) +
#   labs(x="Trait", title="Predictions 2012 rrBLUP", y="Predictive ability")+
#   theme(axis.text.x=element_text(angle = 45,  hjust = 1))+
#   scale_y_continuous(limits = c(-0.5, 1))
# dev.off()




# Purge obsolete variables
rm()

cat("Data modelled!\n")
