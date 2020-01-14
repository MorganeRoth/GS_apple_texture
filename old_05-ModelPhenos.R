cat("Modelling pheno data...\n")

dir.create(paste0(odir, "/phenos_modelled"), showWarnings = FALSE, recursive = TRUE)

#############################
### Get BLUPs from phenos ###
#############################

cat("Getting BLUPs")

### year nested in genotypes
### try different models

phenos=readRDS(paste0(idir, "/phenos_ready_for_pred.rds"))
summary(phenos)
phenos$Family<-as.character(phenos$Family)
phenos[which(is.na(phenos$Family)), "Family"]<-"Collection"
phenos$Family<-as.factor(phenos$Family)
mylm_function <- function(trait) {
  mylm1<-lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year))
  mylm2<-lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year/Name))
  mylm3<-lmer(data=phenos, get(trait) ~ (1|Name) +(1|Trial/Name)  )
  values<-lapply(1:3, function (x) AIC(logLik(get(paste0("mylm", x))))) %>% unlist
  # print(values)
  res<-matrix(values, nrow=3, ncol=1, dimnames=list(1:3,trait ))
  return(res)
}
cat("check if necessary which model doesn't converge for which trait\n")

## This model cannot be used for all traits afterwards
all_mod<-sapply(traits, function(trait) mylm_function(trait)) %>% as.data.frame()
# all_mod
# lmer(data=phenos, get(traits[12]) ~ (1|Name) +(1|Trial/Name)  )
# lmer(data=phenos, get(traits[2]) ~ (1|Name) + (1|Year/Name))
## for each trait, selec the best model
selection_mod<-data.frame(trait=traits, model=lapply(traits, function(i) which(all_mod[[i]] == min(all_mod[[i]]))) %>% unlist)
## keep model 3 for all except trait for which model does not converge (look manually)
Names<-levels(phenos$Name)
blups_all_traits<-matrix(NA, nrow=length(Names), ncol=length(traits), dimnames=list(Names, traits))

myfun<-function(trait){
  mylm_sel<-lmer(data=phenos, get(trait) ~ (1|Name) +(1|Trial/Name))
  print(dim(ranef(mylm_sel)$`Name:Trial`)) ## we keep only the Name nested in family because the other terms are GxE
  blups<-ranef(mylm_sel)$`Name`
  print(ranef(mylm_sel)$`Trial`)
  ## need to add the familly effects, because they are genetics!
  ids<-rownames(blups)
  blups_all_traits[which(Names %in% ids), trait]<<-blups$`(Intercept)`[which(ids%in% Names)] ## IMPORTANT use recursive operator
  print(trait)
  rm(blups)
}

# myfun0<-function(trait){
#   mylm_sel<-lmer(data=phenos, get(trait) ~ (1|Name) +(1|Year/Name))
#   print(dim(ranef(mylm_sel)$`Name:Year`)) ## we keep only the Name nested in family because the other terms are GxE
#   blups<-ranef(mylm_sel)$`Name`
#   print(ranef(mylm_sel)$`Year`)
#   ## need to add the familly effects, because they are genetics!
#   ids<-rownames(blups)
#   blups_all_traits[which(Names %in% ids), trait]<<-blups$`(Intercept)`[which(ids%in% Names)] ## IMPORTANT use recursive operator
#   print(trait)
#   rm(blups)
# }
## IMPORTANT
## final version of phenos: all traits work with the selected model
for (trait in traits){ myfun(trait)} 
# for (trait in traits[2]){ myfun0(trait)} 
png(file=paste0(odir, "/phenos_modelled/BLUPs_all_traits_trial_nested_effects.png"), height=600, width=800 )
par(mfrow=c(3,4))
lapply(traits, function(x) print(hist(blups_all_traits[,x], main="", xlab=x))) %>% print
dev.off()

cat("Histogram of BLUPs plotted")

write.table(blups_all_traits,paste0(odir, "/phenos_modelled/BLUPs_all_traits_trial_nested_effects.txt"), sep="\t", quote=F, row.names=T)

## extract heritabilities- need to check code again

myherit<-function(trait){
  mylm_sel<-lmer(data=phenos, get(trait) ~ (1|Name) +(1|Trial/Name))
  myvar<-VarCorr(mylm_sel) %>% as.data.frame()
  genvar<-myvar[which(myvar$grp=="Name"), "vcov"]
  envvar<-myvar[which(myvar$grp=="Name"), "vcov"]+
    # myvar[which(myvar$grp=="Name:Trial"), "vcov"]+
    # myvar[which(myvar$grp=="Trial"), "vcov"]/length(levels(phenos$Trial))+
    myvar[which(myvar$grp=="Residual"), "vcov"]
  htwo<- genvar/envvar 
  return(htwo)
}

herit_all<-lapply(traits, myherit)
names(herit_all)<-traits
herit_all
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
## Need to impute missing data ? there are here but keep script like that in case

## in case we start from here load data
# phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T, row.names = id_pheno$x %>% as.character())[,-c(15,14)]

phenos<-blups_all_traits
rownames(phenos)<-blups_all_traits[,1]
phenos<-phenos[,-1]
nb <- estim_ncpPCA(phenos,method.cv = "Kfold", verbose = FALSE)
res.comp <- imputePCA(phenos, ncp = nb$ncp) 
# whichFAMS<-lapply(families, function(x)grep(x,rownames(phenos))) %>% unlist
par(mfrow=c(1,1))
# res.pca <- PCA(res.comp,  ind.sup =whichFAMS , ncp = nb$ncp, graph=T, axes=c(1,2))
res.pca <- PCA(res.comp,   ncp = nb$ncp, graph=T, axes=c(1,2))
pca_pheno<-data.frame(PC1=res.pca$ind$coord[,1] , row.names=rownames(res.pca$ind$coord))
## plot with groups
WhichCOL<-c(1:nrow(phenos))[-c(lapply(families, function(x) grep(x,rownames(phenos)) ) %>% unlist)]
groups<-data.frame(Names=rownames(phenos), Group=NA)
groups$Group[WhichCOL]<-"Collection"
groups$Group[-WhichCOL]<-lapply(groups$Name[-WhichCOL] ,function(x) substr(x,1,4)) %>% unlist
# png(file=paste0(odir, "/phenos_modelled/PCA_all_traits.png"), width=1200, height=800, res=150)
pdf(file=paste0(odir, "/phenos_modelled/PCA_all_traits.pdf"), width=8, height=6)
## useful infos on plotting here
# http://www.sthda.com/french/wiki/fviz-pca-visualisation-de-l-analyse-en-composante-principale-logiciel-r-et-analyse-de-donn-es
fviz_pca_ind(res.pca,  habillage= as.factor(groups$Group), 
             # palette=c("black","blue", "turquoise",  "purple", "red","green3", "orange"),
             palette=c("black","red", "lightgreen", "skyblue", "orange", "yellow", "purple"),
             label="none", pointshape = 19) %>% print
dev.off()

## plot graph of variables
res.pca$var$coord
# traits_short=c("IF", "MF","FF", "NPF", "Area", "FLD", "YM", "MF", "ANP", "AMP", "ALD", "AMP")
MYcols<-c(rep(2, 4), rep(3, 8))
pdf(file=paste0(odir, "/phenos_modelled/PCA_variables.pdf"), height=6, width=6)
fviz_pca_var(res.pca,  habillage= as.factor(MYcols),palette=c("darkblue", "red"),label="none")+ 
  geom_text(aes(label=c(paste0("     ",c(2,4,3,1)), "      9", "   12", "    7", "     8", "    11", "           10", "     5", "     6")))+
  theme(legend.position = "none")
dev.off()

## plot just families and parents
whichParents<-which(rownames(phenos) %in% c("FuMoHo", "Pinova", "CriPin", "RoyGal", "GolDel", "Delear"))
whichFAM<-lapply(families, function(x) grep(x,rownames(phenos)) ) %>% unlist
subgroup=groups[c(whichFAM,whichParents),]
nb <- estim_ncpPCA(phenos[c( whichFAM,whichParents),],method.cv = "Kfold", verbose = FALSE)
res.comp <- imputePCA(phenos[c( whichFAM,whichParents),], ncp = nb$ncp) 
res.pca2 <- PCA(res.comp,  ncp = nb$ncp, graph=T)
nn<-c(rep("",(length(c( whichFAM,whichParents))-6)), "PinkLady" ,  "Delear" ,  "FuMoHo" ,  "GolDel" ,  "Pinova" ,  "RoyGal"  )

# png(file=paste0(odir, "/phenos_modelled/PCA_all_traits_FAM_PARENTS.png"), width=1200, height=800, res=150)
pdf(file=paste0(odir, "/phenos_modelled/PCA_all_traits_FAM_PARENTS.pdf"), width=8, height=6)
fviz_pca_ind(res.pca2, habillage = as.factor(subgroup$Group), 
             palette=c("black","red", "lightgreen", "skyblue", "orange", "yellow", "purple"),
             pointshape = 19, label="none") +
             geom_text(aes(label=nn, x=res.pca2$ind$coord[,1], y=res.pca2$ind$coord[,2]), size=5)

dev.off()

## add phenos PC1 and PC2 - perfomed with whole population (progenies and collection)
phenos<-as.data.frame(phenos)
phenos$PC1<-phenos$PC2<-NA
phenos[rownames(pca_pheno),]$PC1<-res.pca$ind$coord[,1]
phenos[rownames(pca_pheno),]$PC2<-res.pca$ind$coord[,2]
cat("checking that PC1 and PC2 are added")
summary(phenos)
cat("phenos file consists in BLUPs + PC1 and PC2 - can be modified in script!")

write.table(phenos, file=paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_nosup_for_pred.txt"), sep="\t", quote=F, row.names = T)

## alternative: calculate PC1 and PC2 only with collection
## do the same but put the progenies as supplementary individuals
## look whether plots look differently

# phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_nosup_for_pred.txt"), h=T)[,-c(13,14)]
head(phenos)
families<-c("FjDe", "GDFj", "FjPi", "FjPL", "GaPi", "GaPL")
whichFAMS<-lapply(families, function(x)grep(x,rownames(phenos))) %>% unlist
res.pca <- PCA(phenos,  ind.sup =whichFAMS , ncp = 5, graph=F)
pca_pheno<-data.frame(PC1=res.pca$ind$coord[,1] , row.names=rownames(res.pca$ind$coord))
## plot with groups
WhichCOL<-c(1:nrow(phenos))[-c(lapply(families, function(x) grep(x,rownames(phenos)) ) %>% unlist)]
groups<-data.frame(Names=rownames(phenos)[whichFAMS], Group=NA)
groups$Group<-lapply(groups$Name ,function(x) substr(x,1,4)) %>% unlist
groups$Group<-revalue(groups$Group, c("FjDe"="red", "FjPi"="lightgreen", "FjPL"="skyblue", "GaPi"="orange",
                                      "GaPL" ="yellow", "GDFj"="purple"))

# png(file=paste0(odir, "/phenos_modelled/PCA_all_traits.png"), width=1200, height=800, res=150)

## useful infos on plotting here
# http://www.sthda.com/french/wiki/fviz-pca-visualisation-de-l-analyse-en-composante-principale-logiciel-r-et-analyse-de-donn-es
fviz_pca_ind(res.pca,label='none', col.ind.sup = as.factor(groups$Group),
             pointshape = 19) %>% print
dev.off()

## create new table and plot with ggplot because function makes problem with labels
res<-data.frame(ID= c(rownames(res.pca$ind$coord), rownames(res.pca$ind.sup$coord)),
                PC1=c(res.pca$ind$coord[,1], res.pca$ind.sup$coord[,1]),
                PC2=c(res.pca$ind$coord[,2], res.pca$ind.sup$coord[,2]))
groups<-matrix(NA,ncol=1, nrow=nrow(res.pca$ind$coord)+nrow(res.pca$ind.sup$coord))
groups[c(1:nrow(res.pca$ind$coord))]<- "Collection"
groups[-c(1:nrow(res.pca$ind$coord))]<-lapply(res$ID[-c(1:nrow(res.pca$ind$coord))] %>% as.character ,function(x) substr(x,1,4)) %>% unlist
res$Population<-as.factor(groups)
# res$groups<-revalue(res$groups, c("Collection"="black", "FjDe"="red", "FjPi"="lightgreen", "FjPL"="skyblue", "GaPi"="orange",
#                                   "GaPL" ="yellow", "GDFj"="purple"))
pdf(file=paste0(odir, "/phenos_modelled/PCA_all_traits_fam_sup.pdf"), width=6, height=4)
ggplot(data=res, aes( x=PC1, y=PC2,color=Population))+
        geom_point()+
        scale_color_manual(values=c("black","red", "lightgreen", "skyblue", "orange", "yellow", "purple"))+
  theme_minimal()+
  geom_abline(slope=0, intercept=0, linetype=2) + 
  geom_vline(xintercept=0, linetype=2) +
  labs(x="PC1 (80.9%)", y="PC2 (11.9%)")
dev.off()
## plot selected individuals
whichParents<-which(rownames(res) %in% c("FuMoHo", "Pinova", "CriPin", "RoyGal", "GolDel", "Delear"))
whichFAM<-lapply(families, function(x) grep(x,rownames(res)) ) %>% unlist
subgroup=res[c(whichFAM,whichParents),]
pdf(file=paste0(odir, "/phenos_modelled/PCA_selec_parents_prog_fam_sup.pdf"), width=6, height=4)
ggplot(data=subgroup, aes( x=PC1, y=PC2,color=Population))+
  geom_point()+
  scale_color_manual(values=c("black","red", "lightgreen", "skyblue", "orange", "yellow", "purple"))+
  theme_minimal()+
  geom_abline(slope=0, intercept=0, linetype=2) + 
  geom_vline(xintercept=0, linetype=2) +
  labs(x="PC1 (80.9%)", y="PC2 (11.9%)")+
  theme(legend.position = "none")+
  geom_text(data=res[whichParents,],
            aes( x=PC1, y=PC2,label=ID))
dev.off()

head(phenos)
phenos[rownames(res),"PC1"]<-res$PC1
phenos[rownames(res),"PC2"]<-res$PC2
## write phenos ready for prediction
write.table(phenos, file=paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), sep="\t", quote=F, row.names = T)
write.table(rownames(phenos), file=paste0(odir, "/phenos_modelled/rownames_phenos.txt"))
## plot graph of variables
res.pca$var$coord
# traits_short=c("IF", "MF","FF", "NPF", "Area", "FLD", "YM", "MF", "ANP", "AMP", "ALD", "AMP")
MYcols<-c(rep(2, 4), rep(3, 8))
pdf(file=paste0(odir, "/phenos_modelled/PCA_variables_with_sup_fams.pdf"), height=6, width=6)
fviz_pca_var(res.pca,  habillage= as.factor(MYcols),palette=c("darkblue", "red"),label="none")+ 
  geom_text(aes(label=c(paste0("     ",c(2,4,3,1)), "      9", "   12", "    7", "     8", "    11", "           10", "     5", "     6")))+
  theme(legend.position = "none")
dev.off()

### histograms for each trait for families and collection for sup data

# phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)

WhichCOL<-c(1:nrow(phenos))[-c(lapply(families, function(x) grep(x,rownames(phenos)) ) %>% unlist)]
rownames(phenos)[WhichCOL]
rownames(phenos)[whichFAM]
traits=colnames(phenos)
pops<- c(list(WhichCOL), lapply(families, function(x) grep(x,rownames(phenos)) ) )
names(pops)=c("COL", families)
pdf(file=paste0(odir, "/phenos_modelled/hist_all_traits_supdata.pdf"), height=20, width=30)
par(mfrow=c(7,14), mar=rep(1.5,4))
for (pop in names(pops)) {
  print(pop)
  for (trait in traits) {
    hist(phenos[pops[[pop]],trait], main=pop, xlab=trait)
  }
}
dev.off()

myvec=rep(NA, nrow(phenos))
for (pop in names(pops)) {
  myvec[pops[[pop]]]=pop
}
phenos$pop=myvec
phenos$id<-rownames(phenos)
phenos_melt<- melt(phenos, id=c("id", "pop"))
phenos_melt$variable<-revalue(phenos_melt$variable, c("Acoustic_Linear_Distance_BLUP"="ALD",
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
# mu <- ddply(phenos_melt, c("variable", "pop"), summarise, grp.mean=mean(value))

pdf(file=paste0(odir, "/phenos_modelled/all_histograms.pdf"), height=8, width=15)
ggplot(data=phenos_melt, aes(x=value)) + 
  geom_histogram(fill="black")+
  # geom_density()+
  facet_grid( phenos_melt$pop ~ phenos_melt$variable , scales="free")+
dev.off()

head(phenos_melt)
mm<-merge(phenos_melt, mu, by=c("pop", "variable"),all.x=T)
head(mm)

pdf(file=paste0(odir, "/phenos_modelled/all_histograms.pdf"), height=8, width=15)
ggplot(data=mm, aes(x=value, color=pop)) + 
  geom_histogram()+
  # geom_density()+
  geom_vline( aes(xintercept=grp.mean, color=pop),linetype="dashed") +
  facet_grid( mm$pop ~ mm$variable , scales="free")+
  labs(color = "Population") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

dev.off()
# Purge obsolete variables
rm()

cat("Phenos analyzed and ready for predictions!\n")
