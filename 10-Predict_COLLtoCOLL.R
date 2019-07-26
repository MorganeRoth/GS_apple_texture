cat("Modelling data...\n")

# dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoCOLL"), showWarnings = FALSE, recursive = TRUE)

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

## help on 5-fold design: https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation
nreps=100
accuracy<-data.frame(Rep=rep(1:nreps, length(traits)), trait=lapply(traits, function(x) rep(x, nreps)) %>% unlist,accuracy=NA)
count=0
for (trait in traits) {
  for (rep in 1:nreps) {
    folds <- cut(seq(1,length(WhichCOL)),breaks=5,labels=FALSE) ## create folds
    ids_shuffle<-sample(WhichCOL) ## randomize ids
    fold_acc<-c()
    for(i in 1:5){
      #Segement your data by fold using the which() function 
      WhichVS <- ids_shuffle[which(folds==i) ]
      WhichTS <- ids_shuffle[-which(folds==i) ]
      res <- mixed.solve(y=phenos[WhichTS,trait],Z=genos_pred[WhichTS,])
      Y_VS_pred<- as.vector(genos_pred[WhichVS,] %*% as.matrix(res$u))
      fold_acc<-append(fold_acc,cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete"))
      rm(WhichVS, WhichTS)
    }
    count=count+1
    accuracy[count,"accuracy"]<-mean(fold_acc)
    # print(c(trait, rep, i, count))
  }
}    

write.table(accuracy, file=paste0(odir, "/predictions/COLLtoCOLL/accuracy_5fold_rrBLUP.txt"), quote=F, sep="\t")
head(accuracy)
# plot
accuracy<-read.table(paste0(odir, "/predictions/COLLtoCOLL/accuracy_5fold_rrBLUP.txt"))

png(file=paste0(odir, "/predictions/COLLtoCOLL/COLLtoCOLL_5fold_", nreps, "_reps.png"), height=500, width=800)
ggplot(accuracy,aes( y=accuracy, x=trait))+
  geom_boxplot()+
  labs(x="Trait", title="Predictions rrBLUP 100 CVs 5-fold", y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))+ 
  scale_y_continuous(limits = c(0, 1))
dev.off()

## compare with Monte Carlo
nreps=100
accuracy<-data.frame(Rep=rep(1:nreps, length(traits)), trait=lapply(traits, function(x) rep(x, nreps)) %>% unlist,accuracy=NA)
count=0
for (trait in traits) {
  for (rep in 1:nreps) {
    WhichTS<-sample(WhichCOL, round(length(WhichCOL))*4/5)
    res <- mixed.solve(y=phenos[WhichTS,trait],Z=genos_pred[WhichTS,])
    WhichVS<-setdiff(WhichCOL, WhichTS)
    Y_VS_pred<- as.vector(genos_pred[WhichVS,] %*% as.matrix(res$u))
    acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
    count=count+1
    accuracy[count,"accuracy"]<-acc
    # print(c(trait, rep, i, count))
    rm(WhichVS, WhichTS)
  }
}    

write.table(accuracy, file=paste0(odir, "/predictions/COLLtoCOLL/COLLtoCOLL_MC_", nreps, "_reps.txt"), quote=F, sep="\t")

png(file=paste0(odir, "/predictions/COLLtoCOLL/COLLtoCOLL_MC_", nreps, "_reps.png"), height=500, width=800)
ggplot(accuracy,aes( y=accuracy, x=trait))+
  geom_boxplot()+
  labs(x="Trait", title="Predictions rrBLUP 100 Monte Carlo 1/5", y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))+
  scale_y_continuous(limits = c(0, 1))
dev.off()

## validate within same year
head(phenos_raw)
phenos_col<-phenos_raw[which(phenos_raw[,"Family"]=="Collection"),] %>% droplevels
phenos_col<-phenos_col[order(phenos_col$Name),]
phenos_col$Year<-as.factor(phenos_col$Year)
head(phenos_col)
nbY<-length(levels(phenos_col$Year))
traits<-colnames(phenos_col)[-c(1:7)]
accuracy<- data.frame(Rep=rep(1:nreps, length(traits)*nbY), 
                     Year=lapply(levels(phenos_col$Year), function(x) rep(x, nreps*length(traits)))%>% unlist,
                     trait=lapply(traits, function(x) rep(x, nreps)) %>% unlist %>% rep(. ,nbY),
                     accuracy=NA )
nreps=100
count=0
for (Y in levels(phenos_col$Year)) {
  data<-phenos_col[which(phenos_col$Year==Y),] %>% droplevels
  for (trait in traits) {
    NonMissing<-which(!(is.na(data[,trait]))) 
    pheno<-data[NonMissing,] %>% droplevels 
    pheno<- aggregate(x=pheno[,trait],by=list(pheno$Name),FUN= mean)
    rownames(pheno)<-pheno$Group.1
    ids<-intersect(rownames(genos_pred), rownames(pheno)) %>% sort
    pheno<-pheno[ids,]
    WhichCOL<-which(rownames(genos_pred) %in% ids)
    summary(rownames(genos_pred)[WhichCOL] == rownames(pheno))
    for (rep in 1:nreps) {
      folds <- cut(seq(1,length(WhichCOL)),breaks=5,labels=FALSE) ## create folds
      ids_shuffle<-sample(WhichCOL) ## randomize ids
      fold_acc<-c()
      for(i in 1:5){
        #Segment your data by fold using the which() function 
        WhichVS <- ids_shuffle[which(folds==i) ]
        WhichTS <- ids_shuffle[-which(folds==i) ]
        res <- mixed.solve(y=pheno[rownames(genos_pred)[WhichTS],"x"],Z=genos_pred[WhichTS,])
        Y_VS_pred<- as.vector(genos_pred[WhichVS,] %*% as.matrix(res$u))
        fold_acc<-append(fold_acc,cor(pheno[rownames(genos_pred)[WhichVS],"x"], Y_VS_pred, use="na.or.complete"))
        rm(WhichVS, WhichTS)
      }
      count=count+1
      accuracy[count,"accuracy"]<-mean(fold_acc)
      # print(c(trait, rep, i, count))
    }
  }    
}
write.table(accuracy, file=paste0(odir, "/predictions/COLLtoCOLL/COLLtoCOLL_per_year_5Fold_", nreps, "_reps.txt"),sep="\t", quote=F)


accuracy<-read.table(paste0(odir, "/predictions/COLLtoCOLL/COLLtoCOLL_per_year_5Fold_", nreps, "_reps.png"), sep="\t", h=T)


png(file=paste0(odir, "/predictions/COLLtoCOLL/COLLtoCOLL_per_year_5Fold_", nreps, "_reps.png"), height=500, width=1000)

ggplot(accuracy,aes( y=accuracy, x=trait))+
  geom_boxplot()+
  facet_grid(~Year) +
  labs(x="Trait", title="Predictions rrBLUP 100 5-fold", y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))+
  scale_y_continuous(limits = c(-0.5, 1))
dev.off()


## compare results with ADI model later  






# Purge obsolete variables
rm()

cat("Data modelled!\n")
