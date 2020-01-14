cat("Modelling data...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/FULLtoFULL"), showWarnings = FALSE, recursive = TRUE)

#### Here we test out how much predictions are improved when we predict within the whole pop = collection + families

## load pheno geno
phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)
traits=colnames(phenos)
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))

## Now we need phenos and genos_ready files
## reinitialise order of IDs
ids<-intersect(rownames(phenos), rownames(genos_ready)) %>% sort
NbID<-length(ids)
cat("Number of common genos/phenos: \n")
print(NbID)
phenos<-phenos[ids,]
genos_ready<-genos_ready[ids,]
nrow(phenos) == nrow(genos_ready)
## clusters from DAPC analysis - merge from family and from collection
clusters<-read.table(paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"), h=T)
cluster_fams<-read.table(paste0(odir, "/genos_modelled/assignements_families.txt"), h=T)
colnames(cluster_fams)<-colnames(clusters)
clusters<-rbind(clusters, cluster_fams)
rownames(clusters)<-clusters$Name
clusters<-clusters[rownames(genos_ready),]
summary(clusters$Name == rownames(genos_ready))

## simple rrBLUP

################################
#### 5-fold within whole pop ###
################################

nREPs=100
accuracy<-data.frame(Rep=rep(1:nREPs, length(traits)), trait=lapply(traits, function(x) rep(x, nREPs)) %>% unlist,
                     accuracy_Pearson=NA, accuracy_Spearman=NA)
count=0
for (trait in traits) {
  for (REP in 1:nREPs) {
    folds <- cut(seq(1,nrow(genos_ready)),breaks=5,labels=FALSE) ## create folds
    ids_shuffle<-sample(rownames(genos_ready)) ## randomize ids
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
write.table(accuracy, file=paste0(odir, "/predictions/FULLtoFULL/accuracy_5fold_ALLtoALL_pearson_spearman.txt"), quote=F, sep="\t")
summary(accuracy)
##### compare to within COLL cross-validations
colcol<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoCOLL/accuracy_5fold_rrBLUP.txt", h=T)
CCM<-aggregate(colcol$accuracy, by=list(colcol$trait), mean)
CCsd<-aggregate(colcol$accuracy, by=list(colcol$trait), sd)
CC<-merge(CCM, CCsd, by="Group.1")
CC$TRS<-"COLL"
FF<-read.table( file=paste0(odir, "/predictions/FULLtoFULL/accuracy_5fold_ALLtoALL_pearson_spearman.txt"))
FFM<-aggregate(FF$accuracy_Pearson, by=list(FF$trait), mean)
FFsd<-aggregate(FF$accuracy_Pearson, by=list(FF$trait), sd)
FFF<-merge(FFM, FFsd, by="Group.1")
FFF$TRS<-"COLL_6FAMS"
colnames(FFF)<-colnames(CC)<-c("trait", "accuracy", "sd", "TRS")
CC<-rbind(CC,FFF)

CC$trait<-revalue(CC$trait, c("Acoustic_Linear_Distance_BLUP"="ALD",
                              "Acoustic_Max_Pressure_BLUP"="APMax",
                              "Acoustic_Mean_Pressure_BLUP" ="APMean",
                              "Acoustic_Npeak_BLUP"="ANP",
                              "Area_BLUP"="Area",
                              "Final_Force_BLUP" ="FF",
                              "Force_Linear_Distance_BLUP" ="FLD" ,
                              "Initial_Force_BLUP" = "FI",
                              "Max_Force_BLUP" = "FMax",
                              "Mean_Force_BLUP"= "FMean", 
                              "N_Peak_Force_BLUP" ="FNP",
                              "Young_Module_BLUP" = "YM"))


pdf(file=paste0(odir, "/predictions/FULLtoFULL/comp_cross_valid.pdf"), height=5, width=8)
ggplot(CC,aes( y=accuracy, x=trait, fill=TRS))+
  geom_bar(stat="identity",position=position_dodge()) +
  scale_fill_manual(values=c("grey40", "grey70")) +
  labs(x="Trait", title="Predictions rrBLUP, 5-fold CVs 259 vs. 537 IDs", y="Predictive ability")+
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_errorbar(aes(ymin=accuracy-sd, ymax=accuracy+sd), width=.2,
                position=position_dodge(.9)) 
dev.off()

##########################################################
## Predict each family with the rest of the population ###
##########################################################

NbFAM<-length(families)
WhichCOL<-c(1:NbID)[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA, conf_low=NA, conf_high=NA)

for (trait in traits) {
  # trait=traits[1]
    for (fam in families) {
      # fam=families[1]
      WhichVS<-ids[grep(fam, ids)]
      WhichTRS<-setdiff(ids, WhichVS)
      res <- mixed.solve(y=phenos[WhichTRS,trait],Z=genos_ready[WhichTRS,])
      Y_VS_pred<- as.vector(genos_ready[WhichVS,] %*% as.matrix(res$u))
      acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
      conf_int<-cor.test(phenos[WhichVS, trait], Y_VS_pred, conf.level = 0.95)$conf.int[c(1,2)]
      print(c(trait, fam, acc, conf_int))
      accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, acc, conf_int)
      rm(res, WhichVS, Y_VS_pred)
    } 
  }
write.table(accuracy, file=paste0(odir, "/predictions/FULLtoFULL/accuracy_5fold_ALLtoFAM_pearson_spearman.txt"), quote=F, sep="\t")  

accuracy<-read.table(file=paste0(odir, "/predictions/FULLtoFULL/accuracy_5fold_ALLtoFAM_pearson_spearman.txt"))
head(accuracy)
ggplot(accuracy,aes( y=accuracy, x=trait))+
  geom_boxplot()+
  facet_grid(~FAM)+
  labs(x="Trait", y="Accuracy")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))+
  scale_y_continuous(limits = c(-0.5, 1))
### this approach is not very useful


####################################################################################
#### predict each individual of a family with remaining IDs from same family #######
####################################################################################

## list of IDs in collection and family names
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
NbFAM<-length(families)
accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA, conf_low=NA, conf_high=NA)
pred_obs<-data.frame(FAM=NA, id=NA, trait=NA,pred=NA, obs=NA)
count0=0
for (trait in traits) {
  # trait=traits[1]
  for (fam in families) {
    WhichFAM<-ids[grep(fam, ids)]
    count=0
    pred_obs2=data.frame(pred=NA, obs=NA)
    # fam=families[1]
    for (i in WhichFAM) {
      count=count+1
      count0=count0+1
      WhichTRS<-setdiff(WhichFAM, i)
      res <- mixed.solve(y=phenos[WhichTRS,trait],Z=genos_ready[WhichTRS,])
      Y_VS_pred<- as.vector(genos_ready[i,] %*% as.matrix(res$u))
      pred_obs[count0,]<-c(fam, i, trait, Y_VS_pred ,phenos[i, trait])
      pred_obs2[count,]<-c(Y_VS_pred ,phenos[i, trait])
    }
    acc<-cor(pred_obs2[,"pred"], pred_obs2[,"obs"] , use="na.or.complete")
    conf_int<-cor.test(pred_obs2[,"pred"], pred_obs2[,"obs"], conf.level = 0.95)$conf.int[c(1,2)]
    print(c(trait, fam, acc, conf_int))
    accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, acc, conf_int)
    rm(res,  Y_VS_pred)
    }
}

write.table(pred_obs, file=paste0(odir, "/predictions/FULLtoFULL/rrBLUPs_LOO_FAMs_eachID.txt"), quote=F, sep="\t")
write.table(accuracy, file=paste0(odir, "/predictions/FULLtoFULL/rrBLUPs_LOO_FAMs.txt"), quote=F, sep="\t")
accuracy<-read.table(paste0(odir, "/predictions/FULLtoFULL/rrBLUPs_LOO_FAMs.txt"))
head(pred_obs)

## plot
# png(file=paste0(odir, "/predictions/COLLtoFAMs/COLLtoFAM_no_clusters.png"), height=500, width=1000)

pred_obs$trait<-revalue(pred_obs$trait, c("Acoustic_Linear_Distance_BLUP"="ALD",
                                                      "Acoustic_Max_Pressure_BLUP"="APMax",
                                                      "Acoustic_Mean_Pressure_BLUP" ="APMean",
                                                      "Acoustic_Npeak_BLUP"="ANP",
                                                      "Area_BLUP"="Area",
                                                      "Final_Force_BLUP" ="FF",
                                                      "Force_Linear_Distance_BLUP" ="FLD" ,
                                                      "Initial_Force_BLUP" = "FI",
                                                      "Max_Force_BLUP" = "FMax",
                                                      "Mean_Force_BLUP"= "FMean", 
                                                      "N_Peak_Force_BLUP" ="FNP",
                                                      "Young_Module_BLUP" = "YM"))
gg<-ggplot(pred_obs,aes( y=pred %>% as.numeric(), x=obs%>% as.numeric()))+
  geom_point() + 
  facet_wrap( FAM ~ trait , scales="free")

pdf(file=paste0(odir, "/predictions/FULLtoFULL/LOO_per_FAM_trait.pdf"), height=20, width=20)
print(gg)
dev.off()


##################################################################
#### try to incorporate FjPL to improve predictions of PC2 #######
##################################################################

NbFAM<-length(families)
WhichCOL<-ids[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]
accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA, conf_low=NA, conf_high=NA)

for (trait in traits) {
  # trait=traits[1]
  WhichTRS<-c(WhichCOL, ids[grep("FjPL", ids)])
  for (fam in families[-3]) {
    # fam=families[1]
    WhichVS<-ids[grep(fam, ids)]
    res <- mixed.solve(y=phenos[WhichTRS,trait],Z=genos_ready[WhichTRS,])
    Y_VS_pred<- as.vector(genos_ready[WhichVS,] %*% as.matrix(res$u))
    acc<-cor(phenos[WhichVS, trait], Y_VS_pred, use="na.or.complete")
    conf_int<-cor.test(phenos[WhichVS, trait], Y_VS_pred, conf.level = 0.95)$conf.int[c(1,2)]
    print(c(trait, fam, acc, conf_int))
    accuracy[which(accuracy$FAM==fam & accuracy$trait==trait), ] <- c(fam, trait, acc, conf_int)
    rm(res, WhichVS, Y_VS_pred)
  } 
}
write.table(accuracy, file=paste0(odir, "/predictions/FULLtoFULL/accuracy_5fold_COLL_FjPL_toFAM.txt"), quote=F, sep="\t") 
accuracy<-read.table(file=paste0(odir, "/predictions/FULLtoFULL/accuracy_5fold_COLL_FjPL_toFAM.txt"))
write.table(accuracy[which(accuracy$trait =="PC2"),] , file=paste0(odir, "/predictions/FULLtoFULL/PC2_withFjPL_in_TRS.txt"))
summary(accuracy[which(accuracy$trait =="PC2"),])
summary(accuracy[which(accuracy$trait =="PC2" & accuracy$FAM=="FjPL"),]) ## normal to see that

ggplot(accuracy,aes( y=accuracy%>%as.numeric, x=trait))+
  geom_boxplot()+
  facet_grid(~FAM)+
  labs(x="Trait",  y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1))+
  scale_y_continuous(limits = c(-0.5, 1))

## plot observed vs predicted for one trait

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


## try with clusters (as fix effect)
head(clusters)
summary(clusters)
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

accuracy$accuracy<- accuracy$accuracy %>% as.character %>% as.numeric

png(file=paste0(odir, "/predictions/COLLtoFAMs/COLLtoFAM_with_clusters.png"), height=500, width=1000)
ggplot(accuracy,aes( y=accuracy, x=trait))+
  geom_boxplot()+
  facet_grid(~FAM)+
  labs(x="Trait", title="Predictions rrBLUP COLL to FAMs with 6 clusters", y="Predictive ability")+
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
