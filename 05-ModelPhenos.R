cat("Modelling pheno data...\n")

dir.create(paste0(odir, "/phenos_modelled"), showWarnings = FALSE, recursive = TRUE)

#############################
### Get BLUPs from phenos ###
#############################

cat("Getting BLUPs")
#### use lme4 for getting BLUPs - model selection for each trait

mylm_function<-function(trait) {
  mylm1<-lmer(data=phenos, get(trait) ~ (1|Name) + (1|Trial))
  mylm2<-lmer(data=phenos, get(trait) ~ (1|Name) + Trial)
  mylm3<-lmer(data=phenos, get(trait) ~ (1|Name) + Year + Location )
  mylm4<-lmer(data=phenos, get(trait) ~ (1|Name) + Year )
  mylm5<-lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year))
  values<-lapply(1:5, function (x) AIC(logLik(get(paste0("mylm", x))))) %>% unlist
  res<-matrix(values, nrow=5, ncol=1, dimnames=list(1:5,trait ))
}

all_mod<-sapply(traits, function(trait) mylm_function(trait)) %>% as.data.frame()
rownames(all_mod)=1:5
cat("AIC values for each model, see script for models")
cat("for each trait, select the best model")
selection_mod<-data.frame(trait=traits, model=lapply(traits, function(i) which(all_mod[[i]] == min(all_mod[[i]]))[1]) %>% unlist)
selection_mod

mylm1<-function(phenos,trait){lmer(data=phenos, get(trait) ~ (1|Name) + (1|Trial))}
mylm2<-function(phenos,trait){lmer(data=phenos, get(trait) ~ (1|Name) + Trial)}
mylm3<-function(phenos,trait){lmer(data=phenos, get(trait) ~ (1|Name) + Year + Location )}
mylm4<-function(phenos,trait){lmer(data=phenos, get(trait) ~ (1|Name) + Year )}
mylm5<-function(phenos,trait){lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year))}
blups_all_traits<-matrix(NA, nrow=length(Names), ncol=length(traits), dimnames=list(Names, traits))

myfun<-function(trait,number){
  mylm_sel<-get(paste0("mylm", number))(phenos, trait)
  blups<-ranef(mylm_sel)$Name
  print(trait)
  print(c(dim(blups), length(which((rownames(blups)%in% Names)))))
  blups_all_traits[rownames(blups), trait]<<-blups$`(Intercept)` ## use recursive operator
}

for (i in 1:dim(selection_mod)[1] ){ myfun(as.character(selection_mod[i,1]), selection_mod[i,2])} ## somehow mapply not working

png(file=paste0(odir, "/phenos_modelled/BLUPs_all_traits_no_nesting.txt.png"))
par(mfrow=c(3,4))
lapply(traits, function(x) print(hist(blups_all_traits[,x], main="", xlab=x))) %>% print
dev.off()

cat("Histogram of BLUPs plotted")

write.table(blups_all_traits, paste0(odir, "/phenos_modelled/BLUPs_all_traits_no_nesting.txt"), sep="\t", quote=F, row.names=T)

### try incorporating families in the model: genotype is nested within family and within year
### try different models in the same way

phenos$Family<-as.character(phenos$Family)
phenos[which(is.na(phenos$Family)), "Family"]<-"Collection"
phenos$Family<-as.factor(phenos$Family)
mylm_function <- function(trait) {
  rm(values, res)
  mylm1<-lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year))
  mylm2<-lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year/Name))
  mylm3<-lmer(data=phenos, get(trait) ~ (1|Name) +(1|Trial/Name)  )
  values<-lapply(1:3, function (x) AIC(logLik(get(paste0("mylm", x))))) %>% unlist
  # print(values)
  res<-matrix(values, nrow=3, ncol=1, dimnames=list(1:3,trait ))
  }
rm(all_mod)
all_mod<-sapply(traits, function(trait) mylm_function(trait)) %>% as.data.frame()
all_mod

## for each trait, selec the best model, we keep 3 for all
selection_mod<-data.frame(trait=traits, model=lapply(traits, function(i) which(all_mod[[i]] == min(all_mod[[i]]))) %>% unlist)
selection_mod
mylm1<-function(phenos, trait){lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year)+ (1|Family))}
mylm2<-function(phenos, trait){lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year/Name))}
mylm3<-function(phenos, trait){lmer(data=phenos, get(trait) ~ (1|Name) +(1|Trial/Name)  )}
blups_all_traits<-matrix(NA, nrow=length(Names), ncol=length(traits), dimnames=list(Names, traits))
head(blups_all_traits)

print("decided to keep only model 3")

myfun<-function(trait){
  mylm_sel<-lmer(data=phenos, get(trait) ~ (1|Name) + (1|Trial/Name) )
  print(dim(ranef(mylm_sel)$`Name:Trial`)) ## we keep only the Name nested in family because the other terms are GxE
  blups<-ranef(mylm_sel)$`Name`
  print(ranef(mylm_sel)$`Trial`)
  ## need to add the familly effects, because they are genetics!
  ids<-rownames(blups)
  blups_all_traits[which(Names %in% ids), trait]<<-blups$`(Intercept)`[which(ids%in% Names)] ## IMPORTANT use recursive operator
  print(trait)
  rm(blups)
}

for (trait in traits){ myfun(trait)} ## somehow mapply not working

png(file=paste0(odir, "/phenos_modelled/BLUPs_all_traits_year_nested_effects.png"))
par(mfrow=c(3,4))
lapply(traits, function(x) print(hist(blups_all_traits[,x], main="", xlab=x))) %>% print
dev.off()

cat("Histogram of BLUPs plotted")

write.table(blups_all_traits,paste0(odir, "/phenos_modelled/BLUPs_all_traits_year_nested_effects.txt"), sep="\t", quote=F, row.names=T)

################################
### Get LS means from phenos ###
################################

cat("calculating LS-Means with year as fix effect")
Means<-list()
for (i in traits) {
  data=phenos[,c("Name", i, "Year")]
  data<-droplevels(data)
  lm_simple<-lm(get(i)~ Name + Year, data=data)
  lsm<-LSmeans(lm_simple, effect="Name")
  # print(dim(lsm$coef))
  # print(length(levels(data$Name)))
  Means[[i]]<-cbind(lsm$coef[, "estimate"], lsm$grid)
  colnames(Means[[i]])[1]<-"value_LSmean"
}
LS_all<-melt(Means)
head(LS_all)
write.table(LS_all, file=paste0(odir, "/phenos_modelled/LS_means_all_traits_year_fix.txt"), quote=F, row.names=F, sep="\t")

cat("compare LS means and BLUPs")
blups_all_traits<-as.data.frame(blups_all_traits)
blups_all_traits$Name<-rownames(blups_all_traits)
blups_all_traits<-blups_all_traits[,c("Name",sort(traits))]
colnames(blups_all_traits)[-1]<-paste0(colnames(blups_all_traits)[-1],"_BLUP")
LS2<-dcast(LS_all,Name~L1, value.var = "value" )

colnames(LS2)[-1]<-paste0(colnames(LS2)[-1],"_LS")
head(blups_all_traits)
blups_LS<-merge(blups_all_traits, LS2, by="Name")
head(blups_LS)
summary(blups_LS)
cor(blups_LS[,c(2:25)])
cor(x=blups_LS[,c(2:13)], y=blups_LS[,c(14:25)])


png(file=paste0(odir, "/phenos_modelled/corrplot_blups_LS.png"))
par(mfrow=c(1,1))
corrplot(cor(x=blups_LS[,c(2:13)], y=blups_LS[,c(14:25)], use="na.or.complete")) %>% print
dev.off()

print("corrplot LS/Blups plotted")

# Get PC1 and PC2 from phenos
## Need to impute missing data
phenos<-blups_all_traits[,-1]
rownames(phenos)<-blups_all_traits[,1]
nb <- estim_ncpPCA(phenos,method.cv = "Kfold", verbose = FALSE)
res.comp <- imputePCA(phenos, ncp = nb$ncp) 
res.pca <- PCA(res.comp, quanti.sup = 1, quali.sup = 12, ncp = nb$ncp, graph=T)
pca_pheno<-data.frame(PC1=res.pca$ind$coord[,1] , row.names=rownames(res.pca$ind$coord))
## plot with groups
WhichCOL<-c(1:nrow(phenos))[-c(lapply(families, function(x) grep(x,rownames(phenos)) ) %>% unlist)]
groups<-data.frame(Names=rownames(phenos), Group=NA)
groups$Group[WhichCOL]<-"Collection"
groups$Group[-WhichCOL]<-lapply(groups$Name[-WhichCOL] ,function(x) substr(x,1,4)) %>% unlist
png(file=paste0(odir, "/phenos_modelled/PCA_all_traits.png"), width=1200, height=800, res=150)
fviz_pca_ind(res.pca,  col.ind = as.factor(groups$Group), palette=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d')) %>% print
dev.off()
## plot just families and parents
whichParents<-which(rownames(phenos) %in% c("FuMoHo", "Pinova", "CriPin", "RoyGal", "GolDel", "Delear"))
whichFAM<-lapply(families, function(x) grep(x,rownames(phenos)) ) %>% unlist
subgroup=groups[c(whichFAM,whichParents),]
nb <- estim_ncpPCA(phenos[c( whichFAM,whichParents),],method.cv = "Kfold", verbose = FALSE)
res.comp <- imputePCA(phenos[c( whichFAM,whichParents),], ncp = nb$ncp) 
res.pca2 <- PCA(res.comp, quanti.sup = 1, quali.sup = 12, ncp = nb$ncp, graph=T)
png(file=paste0(odir, "/phenos_modelled/PCA_all_traits_FAM_PARENTS.png"), width=1200, height=800, res=150)
fviz_pca_ind(res.pca2,  col.ind = as.factor(subgroup$Group), palette=c('black','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d')) %>% print
dev.off()

## add phenos PC1 and PC2
phenos<-as.data.frame(phenos)
phenos$PC1<-phenos$PC2<-NA
phenos[rownames(pca_pheno),]$PC1<-res.pca$ind$coord[,1]
phenos[rownames(pca_pheno),]$PC2<-res.pca$ind$coord[,2]
cat("checking that PC1 and PC2 are added")
summary(phenos)
cat("phenos file consists in BLUPs + PC1 and PC2 - can be modified in script!")

write.table(phenos, file=paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), sep="\t", quote=F, row.names = T)
# Purge obsolete variables
rm()

cat("Phenos analyzed and ready for predictions!\n")
