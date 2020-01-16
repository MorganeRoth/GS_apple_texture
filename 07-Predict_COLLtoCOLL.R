cat("Cross-validations within the collection...\n")
cat("We use model A and B (with and without clustering...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoCOLL"), showWarnings = FALSE, recursive = TRUE)

#######################################
## IF YOU START FROM HERE, LOAD DATA ##
#######################################
## pheno and geno
phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))
## clusters from DAPC analysis
clusters<-read.table(paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"), h=T)
cluster_fams<-read.table(paste0(odir, "/genos_modelled/assignements_families.txt"), h=T)

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

# List all traits
traits=colnames(phenos)

cat("IDs in pheno and geno file:\n")
summary(rownames(genos_ready)[WhichCOL] == rownames(phenos)[WhichCOL]) %>% print ## 259 with same names

## help to understand 5-fold design available here: https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation
nREPs=100
accuracy<-data.frame(Rep=rep(1:nREPs, length(traits)), trait=lapply(traits, function(x) rep(x, nREPs)) %>% unlist,
                     accuracy_Pearson=NA, accuracy_Spearman)
count=0

######## PREDICT WITHOUT CLUSTER EFFECT #########
for (trait in traits) {
# for (trait in c("PC1","PC2")) {
  for (REP in 1:nREPs) {
    folds <- cut(seq(1,length(WhichCOL)),breaks=5,labels=FALSE) ## create folds
    ids_shuffle<-sample(WhichCOL) ## randomize ids
    fold_acc<-c()
    for(i in 1:5){
      #Segement your data by fold using the which() function 
      WhichVS <- ids_shuffle[which(folds==i) ]
      WhichTS <- ids_shuffle[-which(folds==i) ]
      res <- mixed.solve(y=phenos[WhichTS,trait],Z=genos_ready[WhichTS,])
      Y_VS_pred<- as.vector(genos_ready[WhichVS,] %*% as.matrix(res$u))
      fold_acc_P<-append(fold_acc,cor(phenos[WhichVS, trait], Y_VS_pred, method="pearson"))
      fold_acc_S<-append(fold_acc,cor(phenos[WhichVS, trait], Y_VS_pred, method="spearman"))
      rm(WhichVS, WhichTS)
    }
    count=count+1
    accuracy[count,"accuracy_Pearson"]<-mean(fold_acc_P)
    accuracy[count,"accuracy_Spearman"]<-mean(fold_acc_S)
    print(c(mean(fold_acc_P), mean(fold_acc_S)))
  }
}    
write.table(accuracy, file=paste0(odir, "/predictions/COLLtoCOLL/accuracy_5fold_rrBLUP_pearson_spearman.txt"), quote=F, sep="\t")

## quick boxplot
accuracy<-read.table(paste0(odir, "/predictions/COLLtoCOLL/accuracy_5fold_rrBLUP.txt"))
# png(file=paste0(odir, "/predictions/COLLtoCOLL/COLLtoCOLL_5fold_", nreps, "_reps.png"), height=500, width=800)
ggplot(accuracy,aes( y=accuracy, x=trait))+
  geom_boxplot()+
  labs(x="Trait", title="Predictions rrBLUP 100 CVs 5-fold", y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))+ 
  scale_y_continuous(limits = c(0, 1))
# dev.off()

######## PREDICT WITH CLUSTER EFFECT #########

clusters<-as.data.frame(clusters)
clusters$Cluster<-as.factor(clusters$Cluster)
clusters.mat<-class.ind(clusters$Cluster)
rownames(clusters.mat)<-clusters$Name
summary(rownames(clusters.mat)== rownames(genos_ready)[WhichCOL])
nREPs=100
accuracy<-data.frame(Rep=rep(1:nREPs, length(traits)), trait=lapply(traits, function(x) rep(x, nREPs)) %>% unlist,accuracy=NA)
count=0

accuracy<-data.frame(Rep=rep(1:nREPs, length(traits)), trait=lapply(traits, function(x) rep(x, nREPs)) %>% unlist,accuracy=NA)
for (trait in traits) {
  print(trait)
  for (REP in 1:nREPs) {
    folds <- cut(seq(1,length(WhichCOL)),breaks=5,labels=FALSE) ## create folds
    ids_shuffle<-sample(WhichCOL) ## randomize ids
    fold_acc<-c()
    for(i in 1:5){
      #Segement your data by fold using the which() function 
      WhichVS <- ids_shuffle[which(folds==i) ]
      WhichTS <- ids_shuffle[-which(folds==i) ]
      res <- mixed.solve(y=phenos[WhichTS,trait],Z=genos_ready[WhichTS,], X=clusters.mat[ids[WhichTS],])
      Y_VS_pred<- as.vector(genos_ready[WhichVS,] %*% as.matrix(res$u)  + clusters.mat[ids[WhichVS],] %*% res$beta)
      fold_acc<-append(fold_acc,cor(phenos[WhichVS, trait], Y_VS_pred))
      rm(WhichVS, WhichTS)
    }
    count=count+1
    accuracy[count,"accuracy"]<-mean(fold_acc)
    print(mean(fold_acc))
    
  }
  write.table(accuracy, file=paste0(odir, "/predictions/COLLtoCOLL/rrBLUPs_5fold_6clusters_fix.txt"), quote=F, sep="\t")

  
}

cat("Predictions done for cross-validations!\n")
