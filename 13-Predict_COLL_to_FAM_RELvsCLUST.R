cat("Optimising training population with relatedness and clusters...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoFAMs_optim_kinship"), showWarnings = FALSE, recursive = TRUE)

###############
## Load data ##
###############
phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))
A<-read.table(paste0(odir, "/genos_modelled/Ka_Amat.txt"), h=T) ### additive relatinship
clusters<-read.table(odir, "/genos_modelled/assignments_COLL_DAPC.txt", h=T)

## prepare some list and matrices
ids<-intersect(rownames(phenos), rownames(genos_ready)) %>% sort
A<-A[ids,ids]
NbID<-length(ids)
cat("Number of common genos/phenos: \n")
print(NbID)
phenos<-phenos[ids,]
genos_pred<-genos_ready[ids,]
nrow(phenos) == nrow(genos_pred)
## list of IDs in collection and family names
traits=colnames(phenos)
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
NbFAM<-length(families)
WhichCOL<-c(1:NbID)[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
summary(rownames(genos_pred)[WhichCOL] == rownames(phenos)[WhichCOL]) ## 259 with same names
traits=colnames(phenos)
## prepare dataframe to store results
accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA)

####################################################################################################################
##################### PREDICTION OPTIMIZATION BASED ON THREE DIFFERENT PARAMETERS ##################################
####################################################################################################################

##############################
### INCREASING RELATEDNESS ###
### MEAN RELATEDNESS #########
##############################

## run predictions with a TRS starting with 10 IDs increasing in size to the whole collection 
## parameter here is mean relatedness to the predicted family

pred_rel<-function(K,Name) {
  results<-list()
  for (trait in traits) {
    tot_length=length(c(10:length(WhichCOL))) * length(families)
    results[[trait]]<-matrix(NA, nrow=tot_length,ncol=4,
                             dimnames=list(c(1:tot_length), c("size_TRS", "family","mean_rel", "accuracy")))
    counter=0
    for (fam in families) {
      WhichVS<-c(1:NbID)[grep(fam, ids)]
      kin_sel<- apply(K[WhichVS,WhichCOL],2,mean) 
      for (nb_col in c(10:length(WhichCOL))) {
        counter=counter+1
        selection<-order(kin_sel, decreasing = T)[1:nb_col]
        mean_rel<-mean(kin_sel[selection])
        print(mean_rel)
        WhichTRS<- names(kin_sel)[selection]
        res <- mixed.solve(y=phenos[WhichTRS,trait],Z=genos_pred[WhichTRS,])
        Y_VS_pred<- as.vector(genos_pred[WhichVS,] %*% as.matrix(res$u))
        acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
        print(c(trait, fam, nb_col, mean_rel,round(acc, digits=2)))
        results[[trait]][counter,]<-c(nb_col,fam, mean_rel, round(acc, digits=3))
        rm(res, Y_VS_pred,mean_rel)   
      }
      rm(WhichVS)
    }
    names(results[[trait]])<-c("nb_col", "fam", "mean_rel", "acc")
  }
  # return(results)
  saveRDS(results, file=paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/",Name,"_TRS_opt_all_traits.rds"))
}
pred_rel(A, "A.mat")

##############################
### INCREASING RELATEDNESS ###
### MAX RELATEDNESS ##########
##############################

## run predictions with a TRS starting with 10 IDs increasing in size to the whole collection 
## parameter here is max relatedness to the predicted family

pred_rel_max<-function(K,Name) {
  results<-list()
  for (trait in traits) {
    tot_length=length(c(10:length(WhichCOL))) * length(families)
    results[[trait]]<-matrix(NA, nrow=tot_length,ncol=5,
                             dimnames=list(c(1:tot_length), c("size_TRS", "family","rel_added_ID_max","mean_rel_TRS" ,"accuracy")))
    counter=0
    for (fam in families) {
      WhichVS<-c(1:NbID)[grep(fam, ids)]
      kin_sel<- apply(K[WhichVS,WhichCOL],2,max) 
      for (nb_col in c(10:length(WhichCOL))) {
        counter=counter+1
        selection<-order(kin_sel, decreasing = T)[1:nb_col]
        min_rel<-min(kin_sel[selection])
        mean_rel<-mean(kin_sel[selection])
        print(min_rel)
        WhichTRS<- names(kin_sel)[selection]
        res <- mixed.solve(y=phenos[WhichTRS,trait],Z=genos_pred[WhichTRS,])
        Y_VS_pred<- as.vector(genos_pred[WhichVS,] %*% as.matrix(res$u))
        acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
        print(c(trait, fam, nb_col, min_rel,mean_rel,round(acc, digits=2)))
        results[[trait]][counter,]<-c(nb_col,fam, min_rel, mean_rel, round(acc, digits=3))
        rm(res, Y_VS_pred,min_rel)   
      }
      rm(WhichVS)
    }
    names(results[[trait]])<-c("nb_col", "fam", "rel_added_ID_max","mean_rel_TRS", "acc")
  }
  # return(results)
  saveRDS(results, file=paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/",Name,"_TRS_opt_all_traits_MAX.rds"))
}
pred_rel_max(A, "A.mat")

#####################################
### INCREASING NUMBER OF CLUSTERS ###
#####################################

## run predictions with a TRS starting with most related cluster and adding clusters one after the other 
## parameter here is mean relatedness of clusters to the predicted family

rownames(clusters)<-clusters$Name
clusters<-clusters[ids[WhichCOL],] ## put the right order

pred_clust<-function(clusters, K) {
  results<-list()
  for (trait in traits) {
    tot_length=length(levels(as.factor(clusters$Cluster))) * length(families)
    results[[trait]]<-matrix(NA, nrow=tot_length,ncol=5,
                             dimnames=list(c(1:tot_length), c("nb_clusters", "size_TRS", "family","mean_rel", "accuracy")))
    counter=0
    for (fam in families) {
      WhichVS<-c(1:NbID)[grep(fam, ids)]
      kin_sel<- data.frame(mean_rel=apply(K[WhichVS,WhichCOL],2,mean))
      kin_sel$cluster<-clusters[ids[WhichCOL],"Cluster"]
      cluster_sel<- aggregate(kin_sel$mean_rel,by=list(kin_sel$cluster),FUN=mean)
      colnames(cluster_sel)<-c("cluster", "mean_rel")
      cluster_sel<-cluster_sel[order(cluster_sel$mean_rel, decreasing=T),]
      for (nb_clust in 1:length(levels(as.factor(clusters$Cluster))) ) {
        counter=counter+1
        selection<-cluster_sel$cluster[1:nb_clust] %>% as.character
        mean_rel<-mean(cluster_sel[cluster_sel$cluster %in% selection,"mean_rel" ])
        print(mean_rel)
        WhichTRS<- clusters[which(clusters$Cluster %in% selection), "Name"] %>% as.character
        res <- mixed.solve(y=phenos[WhichTRS,trait],Z=genos_pred[WhichTRS,])
        Y_VS_pred<- as.vector(genos_pred[WhichVS,] %*% as.matrix(res$u))
        acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
        print(c(trait, fam, nb_clust,selection, mean_rel,round(acc, digits=2)))
        results[[trait]][counter,]<-c(nb_clust,length(WhichTRS), fam, mean_rel, round(acc, digits=3))
        rm(res, Y_VS_pred,mean_rel)   
      }
      rm(WhichVS)
    }
    names(results[[trait]])<-c("nb_col", "fam", "mean_rel", "acc")
  }
  # return(results)
  saveRDS(results, file=paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/",Name,"_TRS_opt_CLUSTERS_all_traits.rds"))
}
pred_clust(clusters, A)


cat("Prediction optimizations done with 3 first methods!\n")
