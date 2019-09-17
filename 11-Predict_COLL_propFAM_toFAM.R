cat("Predicting families using collection and a proportion of the family...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoFAMs_prop"), showWarnings = FALSE, recursive = TRUE)

## if import, model phenos, mnodel genos scripts are skipped, load data here
id_pheno<-read.table(paste0(odir, "/phenos_modelled/rownames_phenos.txt"))
## problem with duplicated rowname that I do not understand, use id_pheno to replace them (saved with phenos)
phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T, row.names = id_pheno$x %>% as.character())
phenos<-phenos[,-1]
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))

# Predict each family with the collection and proportion of the family
## Now we need phenos and genos_ready files
## reinitialise order of IDs
ids<-intersect(rownames(phenos), rownames(genos_ready)) %>% sort
NbID<-length(ids)
cat("Number of common genos/phenos: \n")
print(NbID)
phenos<-phenos[ids,]
genos_pred<-genos_ready[ids,]
cat("genos and phenos have the same ids order:\n")
nrow(phenos) == nrow(genos_pred)

## clusters from DAPC analysis
clusters<-read.table(paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"), h=T)
cluster_fams<-read.table(paste0(odir, "/genos_modelled/assignements_families.txt"), h=T)
head(clusters)
clusters<-matrix(clusters[,"Cluster"], nrow=nrow(clusters), ncol=1, dimnames = list(clusters[,"Name"], "Cluster"))
rownames(clusters)
clusters<-clusters[rownames(genos_pred)[WhichCOL],"Cluster"] %>% as.matrix
rownames(clusters)<-rownames(genos_pred)[WhichCOL]
rownames(cluster_fams)<-cluster_fams$Name
## simple rrBLUP

## list of IDs in collection and family names
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
NbFAM<-length(families)
WhichCOL<-c(1:NbID)[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
summary(rownames(genos_pred)[WhichCOL] == rownames(phenos)[WhichCOL]) ## 242 with same names
traits=colnames(phenos)


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
        Which_FAM_VS<-Which_FAM[-Which_FAM_TRS]
        print(length(Which_FAM_TRS))
        print(length(Which_FAM_VS))
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

nreps=100
accuracy<-readRDS(paste0(odir, "/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_",round(prop, digits = 2),"_prop", nreps, "_reps.rds"))

write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs_prop/Accuracies.txt"), quote=F, sep="\t")
# accuracy<-read.table(paste0(odir, "/predictions/COLLtoFAMs_prop/Accuracies.txt"), h=T)
summary(accuracy)
head(accuracy)
prop=0.3

png(file=paste0(odir, "/predictions/COLLtoFAMs_prop/COLLtoFAM_prop.png"), height=500, width=2000)
ggplot(accuracy,aes( y=round(accuracy, digits=3), x=trait))+
  geom_boxplot()+
  facet_grid(~FAM)+
  labs(x="Trait", title=paste0("Predictions rrBLUP COLL to FAMs + prop FAM (", prop*100, "%)"), y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1)) 
  # scale_y_continuous(limits = c(0.00, 1.00))
dev.off()

## add clusters
clusters<-class.ind(clusters)
cluster_fams$cluster<-factor(cluster_fams$cluster, levels=c(1:5))
cluster_fams$cluster

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
        accuracy[count,"accuracy"] <-  cor(phenos[Which_FAM_VS, trait], Y_VS_pred)
        print( cor(phenos[Which_FAM_VS, trait], Y_VS_pred))
        # print(head(accuracy))
      }
    }
  }
  saveRDS(accuracy, file = paste0(odir, "/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_",round(prop, digits = 2),"_prop", nreps, "_reps_5clusters.rds"))
}
mypred_CL(100,0.3)

acc<-readRDS( paste0(odir, "/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_",round(prop, digits = 2),"_prop", 100, "_reps_5clusters.rds"))


png(file=paste0(odir, "/predictions/COLLtoFAMs_prop/COLLtoFAM_prop_cluster.png"), height=500, width=2000)
ggplot(acc,aes( y=round(accuracy, digits=3), x=trait))+
  geom_boxplot()+
  facet_grid(~FAM)+
  labs(x="Trait", title=paste0("Predictions rrBLUP COLL to FAMs + prop FAM (", prop*100, "%) + 5 clusters"), y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1)) 
# scale_y_continuous(limits = c(0.00, 1.00))
dev.off()


## predict within year 2012 (common to all)

phenos_raw %>% head
pheno2012<-phenos_raw[phenos_raw$Year==2012,] %>% droplevels
summary(pheno2012)
traits=colnames(pheno2012)[-c(1:7)]
accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA)

for (trait in traits) {
  pheno<-pheno2012[which(!(is.na(pheno2012[,trait]))),] %>% droplevels
  ids<-intersect(levels(pheno$Name), rownames(genos_pred))
  geno<-genos_pred[ids,]
  pheno<-pheno[which(pheno$Name %in% ids), ] %>% droplevels
  pheno<- aggregate(x=pheno[,trait],by=list(pheno$Name),FUN= mean)
  rownames(pheno)<-pheno$Group.1
  WhichCOL<-c(1:nrow(geno))[-c(lapply(families, function(x) grep(x, rownames(geno)) ) %>% unlist)]
  print(length(WhichCOL))
  for (fam in families) {
    WhichVS<-c(1:NbID)[grep(fam, ids)]
    res <- mixed.solve(y=pheno[rownames(geno)[WhichCOL],"x"],Z=geno[WhichCOL,])
    Y_VS_pred<- as.vector(geno[WhichVS,] %*% as.matrix(res$u))
    accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, cor(pheno[rownames(geno)[WhichVS], "x"], Y_VS_pred, use="na.or.complete"))
    rm(res, WhichVS, Y_VS_pred)
  }
}

write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs_2012only.txt"), sep="\t", quote=F)
accuracy<-read.table(paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs_2012only.txt"), h=T)

png(file=paste0(odir, "/predictions/COLLtoFAMs/COLLtoFAM_2012only.png"), height=500, width=1000)
ggplot(accuracy,aes( y=as.numeric(accuracy), x=trait))+
  geom_boxplot()+
  facet_grid(~FAM) +
  labs(x="Trait", title="Predictions 2012 rrBLUP", y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))+
  scale_y_continuous(limits = c(-0.5, 1))
dev.off()

## compare results with ADI model later  






# Purge obsolete variables
rm()

cat("Data modelled!\n")
