cat("Predicting families using collection and a proportion of the family...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoFAMs_prop"), showWarnings = FALSE, recursive = TRUE)

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

## simple rrBLUP

## list of IDs in collection and family names
families<-c("FjDe", "FuPi", "FjPL", "GDFj", "GaPi", "GaPL")
NbFAM<-length(families)
WhichCOL<-c(1:NbID)[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
summary(rownames(genos_pred)[WhichCOL] == rownames(phenos)[WhichCOL]) ## 232 with same names
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
mypred(100,1/3)

# for (trait in traits) {
#   for (fam in families) {
#     WhichFAM<- c(1:NbID)[grep(fam, ids)] 
#     for (rep in 1:nreps) {
#       ## sample part of the family to add to TRS
#       WhichVS<- sample(WhichFAM, round ((1-prop)*length(WhichFAM), digits=0)) %>% sort
#       WhichTS<- c(WhichCOL, setdiff(WhichFAM, WhichVS)) %>% sort
#       res <- mixed.solve(y=phenos[WhichTS,trait],Z=genos_pred[WhichTS,])
#       Y_VS_pred<- as.vector(genos_pred[WhichVS,] %*% as.matrix(res$u))
#       count=count+1
#       accuracy[count,"accuracy" ] <- cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
#       print(c(trait, fam, cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")))
#       rm(res, WhichVS, Y_VS_pred)
#     }
#   }
# }

write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs_prop/Accuracies.txt"), quote=F, sep="\t")
# accuracy<-read.table(paste0(odir, "/predictions/COLLtoFAMs/rrBLUPs.txt"), h=T)
summary(accuracy)
head(accuracy)


png(file=paste0(odir, "/predictions/COLLtoFAMs_prop/COLLtoFAM_prop.png"), height=500, width=2000)
ggplot(accuracy,aes( y=accuracy, x=trait))+
  geom_boxplot()+
  facet_grid(~FAM)+
  labs(x="Trait", title=paste0("Predictions rrBLUP COLL to FAMs + prop FAM (", prop*100, "%)"), y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))
  # scale_y_continuous(limits = c(0, 1))
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