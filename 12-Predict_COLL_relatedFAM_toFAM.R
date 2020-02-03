cat("Predicting with COL and one related family...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoFAMs_add_relFAM"), showWarnings = FALSE, recursive = TRUE)

###############
## Load data ##
###############

phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))
## parents
parents<- read.table(paste0(idir, "/families_parents.txt"), sep="\t", h=T)

##########################################################################
# Scenario 3 
# Predict each family with the collection and a related family
# Model A
##########################################################################

ids<-intersect(rownames(phenos), rownames(genos_ready)) %>% sort
NbID<-length(ids)
cat("Number of common genos/phenos: \n")
print(NbID)
phenos<-phenos[ids,]
genos_pred<-genos_ready[ids,]
cat("genos and phenos have the same ids order:\n")
nrow(phenos) == nrow(genos_pred)
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
NbFAM<-length(families)
WhichCOL<-c(1:NbID)[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
summary(rownames(genos_pred)[WhichCOL] == rownames(phenos)[WhichCOL]) ## 259 with same names
traits=colnames(phenos)
accuracy<-data.frame( fam_predicted=NA,fam_training=NA, trait=NA, acc=NA)
count=0

for (FAM in families) {
    Which_FAM<-grep(FAM, ids)  
    par_sel<-parents[-which(parents[,"Family_name2"]==FAM),]
    fam_par<-parents[which(parents[,"Family_name2"]==FAM),c("Parent1", "Parent2")] %>% as.matrix
    fam_related<-par_sel[which(par_sel[,"Parent1"] %in% fam_par  | par_sel[,"Parent2"] %in% fam_par ),"Family_name2"] %>% as.character
    for (FAM_REL in fam_related) {
      Which_Fam_rel<-grep(FAM_REL,ids)  ## select individuals from the related family
      Which_TR<-c(WhichCOL,Which_Fam_rel ) ## add collection individuals
      for (trait in traits) { 
        count=count+1
        res <- mixed.solve(y=phenos[Which_TR,trait],Z=genos_pred[Which_TR,]) ## estimate marker effects with collection
        Y_VS_pred<- as.vector(genos_pred[Which_FAM,] %*% as.matrix(res$u))
        accuracy[count, ] <-c(FAM, FAM_REL,trait, cor(phenos[Which_FAM, trait], Y_VS_pred))
        print( accuracy[count,])
      }
    }
}
write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/COLLtoFAMs_add_relFAM.txt"), quote=F, sep="\t")

## aggregate results and save 
accuracy<-read.table(paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/COLLtoFAMs_add_relFAM.txt"))
accuracy$acc<-as.numeric(accuracy$acc)
res<-aggregate(accuracy$acc, by=list(accuracy$fam_predicted, accuracy$trait), mean)
res_sd<-aggregate(accuracy$acc, by=list(accuracy$fam_predicted, accuracy$trait), sd)
res<-merge(res, res_sd, by=c("Group.1", "Group.2"))
write.table(res, file=paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/rrBLUPs_aggregated.txt"), quote=F, row.names=F, sep="\t")

##########################################################################
# Scenario 3 
# Predict each family with the collection and a related family
# Model B
##########################################################################
## clusters from DAPC analysis
clusters<-read.table(paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"), h=T)
cluster_fams<-read.table(paste0(odir, "/genos_modelled/assignements_families.txt"), h=T)
head(clusters)
clusters<-matrix(clusters[,"Cluster"], nrow=nrow(clusters), ncol=1, dimnames = list(clusters[,"Name"], "Cluster"))
rownames(clusters)
clusters<-clusters[rownames(genos_pred)[WhichCOL],"Cluster"] %>% as.matrix
rownames(clusters)<-rownames(genos_pred)[WhichCOL]
rownames(cluster_fams)<-cluster_fams$Name
cluster_fams$cluster<-factor(cluster_fams$cluster, levels=1:6)
parents<- read.table(paste0(idir, "/families_parents.txt"), sep="\t", h=T)
accuracy<-data.frame( fam_predicted=NA,fam_training=NA, trait=NA, acc=NA)
count=0
for (FAM in families) {
  Which_FAM<-grep(FAM, ids)  
  par_sel<-parents[-which(parents[,"Family_name2"]==FAM),]
  fam_par<-parents[which(parents[,"Family_name2"]==FAM),c("Parent1", "Parent2")] %>% as.matrix
  fam_related<-par_sel[which(par_sel[,"Parent1"] %in% fam_par  | par_sel[,"Parent2"] %in% fam_par ),"Family_name2"] %>% as.character
  cluster_FAM<-cluster_fams[ids[Which_FAM], "cluster" ] %>% class.ind()
  for (FAM_REL in fam_related) {
    Which_Fam_rel<-grep(FAM_REL,ids)  ## select individuals from the related family
    Which_TR<-c(WhichCOL,Which_Fam_rel ) ## add collection individuals
    cluster_FAMREL<- cluster_fams[ids[Which_Fam_rel] , "cluster" ] %>% class.ind
    rownames(cluster_FAMREL)<- ids[Which_Fam_rel]
    cluster_TRS<-rbind(class.ind(clusters),cluster_FAMREL)
    for (trait in traits) { 
      count=count+1
      res <- mixed.solve(y=phenos[Which_TR,trait],Z=genos_pred[Which_TR,],X=cluster_TRS ) ## estimate marker effects with collection
      Y_VS_pred<- as.vector(genos_pred[Which_FAM,] %*% as.matrix(res$u)+  cluster_FAM %*% res$beta)
      accuracy[count, ] <-c(FAM, FAM_REL,trait, cor(phenos[Which_FAM, trait], Y_VS_pred))
      print( accuracy[count,])
    }
  }
}
write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/COLLtoFAMs_add_relFAM_clusters.txt"))

## aggregate results and save 
accuracy=read.table(paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/COLLtoFAMs_add_relFAM_clusters.txt"))
accuracy$acc<-as.numeric(accuracy$acc)
res<-aggregate(accuracy$acc, by=list(accuracy$fam_predicted, accuracy$trait), mean)
res_sd<-aggregate(accuracy$acc, by=list(accuracy$fam_predicted, accuracy$trait), sd)
res<-merge(res, res_sd, by=c("Group.1", "Group.2"))
colnames(res)<-c("Fam", "Trait", "Mean_acc", "SD_acc")
write.table(res, file=paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/rrBLUPs_clusters_aggregated.txt"), quote=F, row.names=F, sep="\t")
