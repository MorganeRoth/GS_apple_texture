cat("Predicting families using collection and a proportion of the family...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoFAMs_prop"), showWarnings = FALSE, recursive = TRUE)

###############
## Load data ##
###############

phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))

##########################################################################
# Scenario 2 
# Predict each family with the collection and proportion of the family
# Model A
##########################################################################

## reinitialise order of IDs
ids<-intersect(rownames(phenos), rownames(genos_ready)) %>% sort
NbID<-length(ids)
cat("Number of common genos/phenos: \n")
print(NbID)
phenos<-phenos[ids,]
genos_pred<-genos_ready[ids,]
cat("genos and phenos have the same ids order:\n")
nrow(phenos) == nrow(genos_pred)
## list of IDs in collection and family names
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
NbFAM<-length(families)
WhichCOL<-c(1:NbID)[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
summary(rownames(genos_pred)[WhichCOL] == rownames(phenos)[WhichCOL]) ## 259 with same names
traits=colnames(phenos)


## predict nreps times with proportion prop

mypred<-function(nreps, prop){
  accuracy<-data.frame(rep=rep(1:nreps,  length(traits)*NbFAM),
                       FAM=rep(lapply(families, function(x) rep(x, nreps)) %>% unlist, length(traits)), 
                       trait=lapply(traits, function(x) rep(x, NbFAM*nreps)) %>% unlist,
                       prop=prop,
                       accuracy=NA)
  count=0
  for (trait in traits) {
    # trait="Acoustic_Linear_Distance_BLUP"
    for(FAM in families) {
      # FAM="FjDe"
      Which_FAM<-grep(FAM, rownames(genos_pred))
      for (i in 1:nreps){
        count=count+1
        Which_FAM_TRS<-sample(Which_FAM, round(length(Which_FAM))*prop)
        Which_FAM_VS<-Which_FAM[-which(Which_FAM %in% Which_FAM_TRS)]
        print(isTRUE(length(Which_FAM_TRS)+length(Which_FAM_VS) == length(Which_FAM)))
        res <- mixed.solve(y=phenos[c(WhichCOL, Which_FAM_TRS),trait],Z=genos_pred[c(WhichCOL, Which_FAM_TRS),])
        Y_VS_pred<- as.vector(as.matrix(genos_pred[Which_FAM_VS,]) %*% as.matrix(res$u))
        accuracy[count,"accuracy"] <-  cor(phenos[Which_FAM_VS, trait], Y_VS_pred)
        print( cor(phenos[Which_FAM_VS, trait], Y_VS_pred))
        # print(head(accuracy))
      }
    }
  }
  saveRDS(accuracy, file = paste0(odir, "/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_",round(prop, digits = 2),"_prop", nreps, "_reps.rds"))
}
mypred(100,0.3)

##########################################################################
# Scenario 2 
# Predict each family with the collection and proportion of the family
# Model B
##########################################################################

## clusters from DAPC analysis
clusters<-read.table(paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"), h=T)
cluster_fams<-read.table(paste0(odir, "/genos_modelled/assignements_families.txt"), h=T)
clusters<-matrix(clusters[,"Cluster"], nrow=nrow(clusters), ncol=1, dimnames = list(clusters[,"Name"], "Cluster"))
clusters<-clusters[rownames(genos_pred)[WhichCOL],"Cluster"] %>% as.matrix
rownames(clusters)<-rownames(genos_pred)[WhichCOL]
rownames(cluster_fams)<-cluster_fams$Nameclusters<-class.ind(clusters)
cluster_fams$cluster<-factor(cluster_fams$cluster, levels=c(1:6))
cluster_fams$cluster

## predict nreps times with proportion prop

mypred_CL<-function(nreps, prop){
  accuracy<-data.frame(rep=rep(1:nreps,  length(traits)*NbFAM),
                       FAM=rep(lapply(families, function(x) rep(x, nreps)) %>% unlist, length(traits)), 
                       trait=lapply(traits, function(x) rep(x, NbFAM*nreps)) %>% unlist,
                       prop=prop,
                       accuracy=NA)
  count=0
  for (trait in traits) {
    # trait="Acoustic_Linear_Distance_BLUP"
    for(FAM in families) {
      # FAM="FjDe"
      Which_FAM<-grep(FAM, rownames(genos_pred))
      for (i in 1:nreps){
        count=count+1
        Which_FAM_TRS<-sample(Which_FAM, round(length(Which_FAM))*prop)
        Which_FAM_VS<-Which_FAM[-which(Which_FAM %in% Which_FAM_TRS)]
        cluster_FAM<- cluster_fams[ids[grep(FAM,ids)], "cluster" ] %>% class.ind
        rownames(cluster_FAM)<-ids[grep(FAM,ids)]
        cluster_TRS<-rbind(clusters,cluster_FAM[ids[Which_FAM_TRS],])
        res <- mixed.solve(y=phenos[c(WhichCOL, Which_FAM_TRS),trait],Z=genos_pred[c(WhichCOL, Which_FAM_TRS),], X=cluster_TRS) ## cluster as fix effect
        Y_VS_pred<- as.vector(genos_pred[Which_FAM_VS,] %*% as.matrix(res$u) + cluster_FAM[ids[Which_FAM_VS],] %*% res$beta)
        summary(rownames(phenos)[c(WhichCOL, Which_FAM_TRS)]==rownames(genos_pred)[c(WhichCOL, Which_FAM_TRS)]) %>% print
        accuracy[count,"accuracy"] <-  cor(phenos[Which_FAM_VS, trait], Y_VS_pred)
        print( cor(phenos[Which_FAM_VS, trait], Y_VS_pred))
        # print(head(accuracy))
      }
    }
  }
  saveRDS(accuracy, file = paste0(odir, "/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_",round(prop, digits = 2),"_prop", nreps, "_reps_5clusters.rds"))
}
mypred_CL(100,0.3)


cat("Families predicted with scenario 2!\n")
