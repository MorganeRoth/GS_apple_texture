at("Predicting with COL and one related family...\n")

dir.create(paste0(odir, "/predictions"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(odir, "/predictions/COLLtoFAMs_add_relFAM"), showWarnings = FALSE, recursive = TRUE)

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
summary(rownames(genos_pred)[WhichCOL] == rownames(phenos)[WhichCOL]) ## 217 with same names

traits=colnames(phenos)
# accuracy<-data.frame(FAM=rep(families, length(traits)), trait=lapply(traits, function(x) rep(x, NbFAM)) %>% unlist,accuracy=NA)
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

write.table(accuracy, file=paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/rrBLUPs.txt"), quote=F, sep="\t")
accuracy$acc<-as.numeric(accuracy$acc)

png(file=paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/COLLtoFAMs_add_relFAM.png"), height=500, width=2000)
ggplot(accuracy,aes( y=acc, x=trait))+
  geom_bar(stat="identity") +
  facet_grid(~fam_predicted*fam_training)+
  labs(x="Trait", title="Predictions rrBLUP COLL + related FAM", y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1, size=5)) +
  scale_y_continuous(limits = c(-0.5, 1))
dev.off()

res<-aggregate(accuracy$acc, by=list(accuracy$fam_predicted, accuracy$trait), mean)
res_sd<-aggregate(accuracy$acc, by=list(accuracy$fam_predicted, accuracy$trait), sd)
res<-merge(res, res_sd, by=c("Group.1", "Group.2"))
colnames(res)<-c("Fam", "Trait", "Mean_acc", "SD_acc")

png(file=paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/Mean_SD_barplot.png"), height=500, width=1500)
ggplot(res,aes( y=Mean_acc, x=Trait))+
  geom_bar(stat="identity") +
  facet_grid(~Fam)+
  labs(x="Trait", title="Predictions rrBLUP COLL + related FAM", y="Predictive ability")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(limits = c(-0.5, 1)) +
  geom_errorbar(aes(ymin=Mean_acc-SD_acc, ymax=Mean_acc+SD_acc), width=.2,
              position=position_dodge(.9))
dev.off()

write.table(res, file=paste0(odir, "/predictions/COLLtoFAMs_add_relFAM/rrBLUPs_aggregated.txt"), quote=F, row.names=F, sep="\t")

            
