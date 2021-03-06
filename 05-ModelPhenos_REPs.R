cat("Modelling pheno data...\n")

dir.create(paste0(odir, "/phenos_modelled"), showWarnings = FALSE, recursive = TRUE)

#############################
### Get BLUPs from phenos ###
#############################

cat("Getting BLUPs")


### try different models

phenos=readRDS(paste0(idir, "/phenos_ready_for_pred.rds"))
summary(phenos)
traits=colnames(phenos)[-c(1:4,17)]
mylm_function <- function(trait) {
  # mylm1<-lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year))
  # mylm2<-lmer(data=phenos, get(trait) ~ (1|Name) + (1|Year/Name))
  # mylm3<-lmer(data=phenos, get(trait) ~ (1|Name) +(1|Trial/Name)  )
  mylm1<-lmer(data=phenos, get(trait) ~ (1|Name) +Trial  )
  mylm2<-lmer(data=phenos, get(trait) ~ (1|Name) + Year  )
  values<-lapply(1:2, function (x) AIC(logLik(get(paste0("mylm", x))))) %>% unlist
  print(values)
  res<-matrix(values, nrow=2, ncol=1, dimnames=list(1:2,trait ))
  return(res)
}
cat("check if necessary which model doesn't converge for which trait\n")

# ## This model cannot be used for all traits afterwards
# all_mod<-sapply(traits, function(trait) mylm_function(trait)) %>% as.data.frame()
# # all_mod
# lmer(data=phenos, get(traits[10]) ~ (1|Name) +(1|Trial/Name)  )
# lmer(data=phenos, get(traits[2]) ~ (1|Name) + (1|Trial))
# ## for each trait, selec the best model
# selection_mod<-data.frame(trait=traits, model=lapply(traits, function(i) which(all_mod[[i]] == min(all_mod[[i]]))) %>% unlist)
# ## keep model 3 for all except trait for which model does not converge (look manually)
# Names<-levels(phenos$Name)

##  
Names<-levels(phenos$Name)
blups_all_traits<-matrix(NA, nrow=length(Names), ncol=length(traits), dimnames=list(Names, traits))
par(mfrow=c(3,4))
myfun<-function(trait){
  mylm_sel<-lmer(data=phenos, get(trait) ~ (1|Name) + Trial)
  blups<-ranef(mylm_sel)$`Name`
  # print(fixef(mylm_sel))
  ids<-rownames(blups)
  blups_all_traits[which(Names %in% ids), trait]<<-blups$`(Intercept)`[which(ids%in% Names)] ## IMPORTANT use recursive operator
  print(trait)
  print(hist(residuals(mylm_sel), xlab="", main=trait))
  print(shapiro.test(residuals(mylm_sel)))
  rm(blups)
}

lapply(traits, myfun)
summary(blups_all_traits)

## try with technical replicates to look at the effect of the different apples
## this model is anyway singular for all traits!
phenorep<-readRDS(paste0(idir, "/phenos_with_reps_raw_filtered.rds"))
head(phenorep)
phenorep$Trial<-paste0(phenorep$Location, phenorep$Year)
blups_tech<-matrix(NA, nrow=length(Names), ncol=length(traits), dimnames=list(Names, traits))

myfun<-function(trait){
  print(trait)
  mylm_sel<-lmer(data=phenos, get(trait) ~ (1|Name/Apple) + Trial) ## apple nested within names
  print(VarCorr(mylm_sel,comp="Variance")) ## look at random and fix effects: apple has negligeable effect - ignore for parcimony
  blups<-ranef(mylm_sel)$`Name`
  # print(ranef(mylm_sel)$`Trial`)
  ## need to add the familly effects, because they are genetics!
  ids<-rownames(blups)
  blups_all_traits[which(Names %in% ids), trait]<<-blups$`(Intercept)`[which(ids%in% Names)] ## IMPORTANT use recursive operator
  rm(blups)
}
lapply(traits, myfun) 
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
png(file=paste0(odir, "/phenos_modelled/BLUPs_all_traits_5apples_trial_fix.png"), height=600, width=800 )
par(mfrow=c(3,4))
lapply(traits, function(x) print(hist(blups_all_traits[,x], main="", xlab=x))) %>% print
dev.off()

cat("Histogram of BLUPs plotted")

png(file=paste0(odir, "/phenos_modelled/BLUPs_all_traits_with_5apples.png"), height=600, width=800 )
par(mfrow=c(3,4))
lapply(traits, function(x) print(hist(blups_all_traits[,x], main="", xlab=x))) %>% print
dev.off()


png(file=paste0(odir, "/phenos_modelled/residuals_all_traits_with_5apples.png"), height=600, width=800 )
par(mfrow=c(3,4))
lapply(traits, myfun)
dev.off()

write.table(blups_all_traits,paste0(odir, "/phenos_modelled/BLUPs_all_traits_5apples_trial_fix.txt"), sep="\t", quote=F, row.names=T)

## extract heritabilities- need to check code again
## mean number of reps 

mean_reps<-phenos %>% group_by(Name) %>% tally() %>% summarize(mean = mean(n, na.rm=TRUE))

myherit<-function(trait){
  data=phenos[,c("Name", trait, "Trial")] %>% droplevels
  mean_reps<-data %>% group_by(Name, Trial) %>% tally() %>% summarize(mean = mean(n, na.rm=TRUE))
  mylm_sel<-lmer(data=data, get(trait) ~ (1|Name) +Trial)
  myvar<-VarCorr(mylm_sel,comp=c("Variance")) %>% as.data.frame()
  genvar<-myvar[which(myvar$grp=="Name"), "vcov"]
  envvar<-myvar[which(myvar$grp=="Name"), "vcov"]+
    # myvar[which(myvar$grp=="Name:Trial"), "vcov"]+
    # myvar[which(myvar$grp=="Trial"), "vcov"]/length(levels(phenos$Trial))+
    myvar[which(myvar$grp=="Residual"), "vcov"]/mean_reps
  htwo<- genvar/envvar
  return(htwo)

}

herit_all<-lapply(traits, myherit) %>% as.data.frame() %>% t
rownames(herit_all)<-traits
herit_all
write.table(herit_all,file=paste0(odir, "/phenos_modelled/heritabilities.txt"), quote=F)

myherit_BS<-function(trait){
  data=phenos[,c("Name", trait, "Trial")] %>% droplevels
  mylm_sel3<-lmer(data=data, get(trait) ~ (1|Name) )
  myvar3<-VarCorr(mylm_sel3,comp=c("Variance")) %>% as.data.frame()
  genvar3<-myvar3[which(myvar3$grp=="Name"), "vcov"]
  envvar3<-myvar3[which(myvar3$grp=="Name"), "vcov"]+ myvar3[which(myvar3$grp=="Residual"), "vcov"]
  htwo3<- genvar3/envvar3
  print(htwo3)
  return(htwo3)

}
herit_all2<-lapply(traits, myherit_BS) %>% as.data.frame() %>% t
rownames(herit_all2)<-traits
herit_all2
write.table(herit_all2,file=paste0(odir, "/phenos_modelled/heritabilities_BS_revised.txt"), quote=F)

WhichCOL<-c(1:nrow(phenos))[-c(lapply(families, function(x) grep(x,phenos$Name))  %>% unlist)]
length(WhichCOL)
dim(phenos)
myheritCOLL<-function(trait){
  data=phenos[WhichCOL,c("Name", trait, "Trial")] %>% droplevels
  mean_reps<-data %>% group_by(Name, Trial) %>% tally() %>% summarize(mean = mean(n, na.rm=TRUE))
  mylm_sel<-lmer(data=data, get(trait) ~ (1|Name) +Trial)
  myvar<-VarCorr(mylm_sel,comp=c("Variance")) %>% as.data.frame()
  genvar<-myvar[which(myvar$grp=="Name"), "vcov"]
  envvar<-myvar[which(myvar$grp=="Name"), "vcov"]+
    # myvar[which(myvar$grp=="Name:Trial"), "vcov"]+
    # myvar[which(myvar$grp=="Trial"), "vcov"]/length(levels(phenos$Trial))+
    myvar[which(myvar$grp=="Residual"), "vcov"]/mean_reps
  htwo<- genvar/envvar 
  return(htwo)
}

herit_allCOLL<-lapply(traits, myheritCOLL) %>% as.data.frame() %>% t
rownames(herit_allCOLL)<-traits
herit_allCOLL 
write.table(herit_allCOLL,file=paste0(odir, "/phenos_modelled/heritabilitiesCOLLonly.txt"), quote=F)
rm(WhichCOL)
WhichCOL<-c(1:nrow(phenos))[-c(lapply(families, function(x) grep(x,phenos$Name))  %>% unlist)]
length(WhichCOL)
dim(phenos)
myheritCOLL_BS<-function(trait){
  data=phenos[WhichCOL,c("Name", trait, "Trial")] %>% droplevels
  mylm_sel3<-lmer(data=data, get(trait) ~ (1|Name) )
  myvar3<-VarCorr(mylm_sel3,comp=c("Variance")) %>% as.data.frame()
  genvar3<-myvar3[which(myvar3$grp=="Name"), "vcov"]
  envvar3<-myvar3[which(myvar3$grp=="Name"), "vcov"]+ myvar3[which(myvar3$grp=="Residual"), "vcov"]
  htwo3<- genvar3/envvar3
  print(htwo3)
  return(htwo3)
}


herit_allCOLL_BS<-lapply(traits, myheritCOLL_BS) %>% as.data.frame() %>% t
rownames(herit_allCOLL_BS)<-traits
herit_allCOLL_BS
write.table(herit_allCOLL_BS,file=paste0(odir, "/phenos_modelled/heritabilitiesCOLLonly_BS_revised.txt"), quote=F)
rm(WhichCOL)

myheritFAM<-function(trait, fam){
  whichfam<-c(1:nrow(phenos))[ grep(fam,phenos$Name) ]
  data=phenos[whichfam,c("Name", trait, "Trial")] %>% droplevels
  print(levels(data$Name))
  mean_reps<-data %>% group_by(Name, Trial) %>% tally() %>% summarize(mean = mean(n, na.rm=TRUE))
  mylm_sel<-lmer(data=data, get(trait) ~ (1|Name) +Trial)
  myvar<-VarCorr(mylm_sel,comp=c("Variance")) %>% as.data.frame()
  genvar<-myvar[which(myvar$grp=="Name"), "vcov"]
  envvar<-myvar[which(myvar$grp=="Name"), "vcov"]+
    # myvar[which(myvar$grp=="Name:Trial"), "vcov"]+
    # myvar[which(myvar$grp=="Trial"), "vcov"]/length(levels(phenos$Trial))+
    myvar[which(myvar$grp=="Residual"), "vcov"]/mean_reps
  htwo<- genvar/envvar 
  return(unlist(c(trait, fam,htwo)))
}
myheritFAM(traits[1], "GaPL")
ll<-mapply(function(x,y) myheritFAM(x,y), rep(traits, length(families)), lapply(families, function(x) rep(x,length(traits)))%>% unlist)
ll<-as.data.frame(ll) %>% t 
colnames(ll)<-c("trait", "family", "heritability")
write.table(ll, file=paste0(odir, "/phenos_modelled/heritabilities_families.txt"), quote=F, sep="\t")
## means
mymeans<-data.frame(Mean=lapply( traits, function(x) phenos[,c("Name", x, "Trial")] %>% group_by(Name, Trial) %>% summarize(mean = mean(get(x), na.rm=TRUE))) %>% unlist ,
                   SD=lapply( traits, function(x) phenos[,c("Name", x, "Trial")] %>% group_by(Name, Trial) %>% summarize(sd = sd(get(x), na.rm=TRUE))) %>% unlist)
rownames(mymeans)<-traits
mymeans
allsummaries<-cbind(herit_all, herit_allCOLL, mymeans)
head(allsummaries)
write.table(allsummaries,file=paste0(odir, "/phenos_modelled/all_summaries.txt"), quote=F)
################################
### Get LS means from phenos ###
################################

cat("calculating LS-Means with year as fix effect")
Means<-list()
for (i in traits) {
  data=phenos[,c("Name", i, "Trial")]
  data<-droplevels(data)
  lm_simple<-lm(get(i)~ Name + Trial, data=data)
  lsm<-LSmeans(lm_simple, effect="Name")
  # print(dim(lsm$coef))
  # print(length(levels(data$Name)))
  Means[[i]]<-cbind(lsm$coef[, "estimate"], lsm$grid)
  colnames(Means[[i]])[1]<-"value_LSmean"
}
LS_all<-melt(Means)
head(LS_all)
write.table(LS_all, file=paste0(odir, "/phenos_modelled/LS_means_all_traits_trial_fix.txt"), quote=F, row.names=F, sep="\t")

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

phenos<-blups_all_traits
## in case we start from here load data
# phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)[,-c(13,14)]

# colnames(phenos)=lapply(colnames(phenos), function(x) substr(x, 0,(nchar(x)-5))) %>% unlist
traits=colnames(phenos)
traits
## calculate PC1 and PC2 only with collection
## do the same but put the progenies as supplementary individuals
## look whether plots look differently

# phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_nosup_for_pred.txt"), h=T)[,-c(13,14)]
head(phenos)
families<-c("FjDe", "GDFj", "FjPi", "FjPL", "GaPi", "GaPL")
whichFAMS<-lapply(families, function(x)grep(x,rownames(phenos))) %>% unlist
res.pca <- PCA(phenos,  ind.sup =whichFAMS , ncp = 5, graph=F)
pca_pheno<-data.frame(PC1=res.pca$ind$coord[,1] , row.names=rownames(res.pca$ind$coord))

## loadins to calculate PC1 and PC2 from raw data
loadings<-sweep(res.pca$var$coord,2,sqrt(res.pca$eig[1:ncol(res.pca$var$coord),1]),FUN="/")
## plot graph of variables
res.pca$var$coord
res.pca$var$
MYcols<-c(rep(2, 4), rep(3, 8))
fviz_pca_var(res.pca,  habillage= as.factor(MYcols),palette=c("darkblue", "red"))

vargraph<-fviz_pca_var(res.pca,  habillage= as.factor(MYcols),palette=c("darkblue", "red"),label="none")+ 
  geom_text(size=2.5,aes(label=c(paste0("     ",c(2,3,4,1)), "    9 ", "       7", "    6", "   12", "      11", "           -10", "     5", "    8")))+
  theme(legend.position = "none")+
  labs(title="", x="PC1 (80.5%)", y="PC2 (12.7%)")+
  xlim(c(-1.1, 1.1))

pdf(file=paste0(odir, "/phenos_modelled/PCA_variables_with_sup_fams_paper.pdf"), 
    width= 3.54, height=3.54,fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
par(mar=rep(4,4))
print(vargraph)
dev.off()
tiff(file=paste0(odir, "/phenos_modelled/PCA_variables_with_sup_fams_paper.tiff"), 
     width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, units="in", res=300 )
print(vargraph)
dev.off()
## plot with groups
WhichCOL<-c(1:nrow(phenos))[-c(lapply(families, function(x) grep(x,rownames(phenos)) ) %>% unlist)]
groups<-data.frame(Names=rownames(phenos)[whichFAMS], Group=NA)
groups$Group<-lapply(groups$Name ,function(x) substr(x,1,4)) %>% unlist
groups$Group<-revalue(groups$Group, c("FjDe"="red", "FjPi"="lightgreen", "FjPL"="skyblue", "GaPi"="orange",
                                      "GaPL" ="yellow", "GDFj"="purple"))

# png(file=paste0(odir, "/phenos_modelled/PCA_all_traits.png"), width=1200, height=800, res=150)

## useful infos on plotting here
# http://www.sthda.com/french/wiki/fviz-pca-visualisation-de-l-analyse-en-composante-principale-logiciel-r-et-analyse-de-donn-es
# fviz_pca_ind(res.pca,label='none', col.ind.sup = as.factor(groups$Group),
#              pointshape = 19) %>% print
# dev.off()

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
pdf(file=paste0(odir, "/phenos_modelled/PCA_all_traits_fam_sup_reps.pdf"), width=6, height=4)
ggplot(data=res, aes( x=PC1, y=PC2,color=Population))+
  geom_point()+
  scale_color_manual(values=c("black","red", "lightgreen", "skyblue", "orange", "yellow", "purple"))+
  theme_minimal()+
  geom_abline(slope=0, intercept=0, linetype=2) + 
  geom_vline(xintercept=0, linetype=2) +
  labs(x="PC1 (80.5%)", y="PC2 (12.7%)")
dev.off()


idgraph<-ggplot(data=res[which(res$Population=="Collection"),], aes( x=PC1, y=PC2,color=Population))+
  geom_point()+
  scale_color_manual(values=c("black","red", "lightgreen", "skyblue", "orange", "yellow", "purple"))+
  theme_minimal()+
  geom_abline(slope=0, intercept=0, linetype=2) + 
  geom_vline(xintercept=0, linetype=2) +
  labs(x="PC1 (80.5%)", y="PC2 (12.7%)")+
  theme(legend.position = "none")+
  stat_ellipse(data=res[which(!(res$Population=="Collection")),] , aes( x=PC1, y=PC2,color=Population),
               type="norm", size=1)

pdf(file=paste0(odir, "/phenos_modelled/PCA_COLL_fams_ellipses_paper.pdf"), 
    width= 3.54, height=3.54,fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
print(idgraph)
dev.off()
tiff(file=paste0(odir, "/phenos_modelled/PCA_COLL_fams_ellipses_paper.tiff"), 
     width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, units="in", res=300 )
print(idgraph)
dev.off()
## plot selected individuals
whichParents<-which(rownames(res) %in% c("FuMoHo", "Pinova", "CriPin", "RoyGal", "GolDel", "Delear"))
whichFAM<-lapply(families, function(x) grep(x,rownames(res)) ) %>% unlist
subgroup=res[c(whichFAM,whichParents),]
# pdf(file=paste0(odir, "/phenos_modelled/PCA_selec_parents_prog_fam_sup.pdf"), width=6, height=4)
mydata<-res[whichParents,]
mydata$ID<-as.character(mydata$ID)
mydata$ID[c(1,2,3,4,6)]<-c("Pink Lady","Delearly", "Fuji", "Golden Delicious", "Royal Gala")
# mydata[which(mydata$ID=="FuMoHo"),c("PC1", "PC2")]=c(3.7, 1.55)  
famgraph<-ggplot(data=subgroup, aes( x=PC1, y=PC2,color=Population))+
  geom_point(size=1)+
  scale_color_manual(values=c("black","red", "lightgreen", "skyblue", "orange", "yellow", "purple"))+
  # scale_color_startrek() +
  theme_minimal()+
  geom_abline(slope=0, intercept=0, linetype=2) + 
  geom_vline(xintercept=0, linetype=2) +
  labs(x="PC1 (80.5%)", y="PC2 (12.7%)")+
  theme(legend.position = "none")+
  # theme(legend.position = "bottom",
  #       legend.title = element_blank(),
  #       legend.key.size = unit(0.2,"in"),
  #       legend.box.spacing = unit(0,"in"),
  #       legend.box.margin = margin(0,1.2,0,0,unit="cm"))+
  guides(colour = guide_legend(override.aes = list(size=2), ncol = 7))+
  geom_text(data=mydata,
            aes( x=PC1, y=PC2+0.2,label=ID), fontface = "bold",size=3,show.legend = FALSE)
pdf(file=paste0(odir, "/phenos_modelled/PCA_COLL_fams_dots_paper.pdf"), 
    width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
print(famgraph)
dev.off()
tiff(file=paste0(odir, "/phenos_modelled/PCA_COLL_fams_dots_paper.tiff"), 
     width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, units="in", res=300 )
print(famgraph)
dev.off()
phenos[rownames(res),"PC1"]<-res$PC1
phenos[rownames(res),"PC2"]<-res$PC2
## write phenos ready for prediction
write.table(phenos, file=paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), sep="\t", quote=F, row.names = T)
write.table(rownames(phenos), file=paste0(odir, "/phenos_modelled/rownames_phenos.txt"))


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
mm<-merge(phenos_melt, mu, by=c("pop", "variable"),all.x=T)

pdf(file=paste0(odir, "/phenos_modelled/all_histograms.pdf"), height=8, width=15)
ggplot(data=mm, aes(x=value, color=pop)) + 
  geom_histogram()+
  # geom_density()+
  geom_vline( aes(xintercept=grp.mean, color=pop),linetype="dashed") +
  facet_grid( mm$pop ~ mm$variable , scales="free")+
  labs(color = "Population") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

dev.off()

###############################
##### Save correlations #######
###############################

phenos<-read.table(file=paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)
write.table(cor(phenos), file=paste0(odir, "/phenos_modelled/correlations_all_traits_PC12.txt"), quote=F, sep="\t")

#######################################
##### Outliers for fruit texture ######
#######################################
WhichCOL<-c(1:nrow(phenos))[-c(lapply(families, function(x) grep(x,rownames(phenos)))  %>% unlist)]
coll<-phenos[which(rownames(phenos) %in% rownames(phenos)[WhichCOL]),]
rownames(coll[which(coll$PC1 == max(coll$PC1)),])
coll[which(coll$PC1 == max(coll$PC1)),"PC1"]
coll[which(coll$PC1 == max(coll$PC1)),"PC1"]
rownames(coll[which(coll$PC2 == max(coll$PC2)),])
coll[which(coll$PC2 == max(coll$PC2)),"PC2"]
coll[which(coll$PC2 > 1.5 & coll$PC1 >4),c("PC1","PC2")]

coll[which(coll$PC2 + coll$PC1 == max(coll$PC2 + coll$PC1)),c("PC1","PC2")]

coll[which(coll$PC2 == max(coll$PC2)),"PC1"]
rownames(coll[which(coll$PC2 == min(coll$PC2)),])
coll[which(coll$PC2 == min(coll$PC2)),"PC2"]
coll[which(coll$PC2 == min(coll$PC2)),"PC1"]
rownames(phenos[which(phenos$PC1 == max(phenos$PC1) & rownames(phenos) %in% rownames(phenos)[WhichCOL]),])

# calculate heritability of PCA axes
phenos=readRDS(paste0(idir, "/phenos_ready_for_pred.rds"))
IDs<-phenos$Name
summary(phenos)
rownames(loadings)<-lapply(rownames(loadings), function(x)substr(x, 1, nchar(x)-5))
traits<-colnames(phenos )[5:16]
phenos<-phenos[,rownames(loadings)]
## use the really nice function sweep
phenos_load<-lapply(colnames(loadings)[1:2], function(x) sweep(phenos,2, loadings[,x], "*" ) %>% rowSums ) %>% do.call(cbind,. ) %>% as.data.frame()
colnames(phenos_load)<-c("Dim.1", "Dim.2")
phenos_load$Name<-IDs
myherit_BS<-function(trait){
  data=phenos_load[,c("Name", trait)] %>% droplevels
  mylm_sel3<-lmer(data=data, get(trait) ~ (1|Name) )
  myvar3<-VarCorr(mylm_sel3,comp=c("Variance")) %>% as.data.frame()
  genvar3<-myvar3[which(myvar3$grp=="Name"), "vcov"]
  envvar3<-myvar3[which(myvar3$grp=="Name"), "vcov"]+ myvar3[which(myvar3$grp=="Residual"), "vcov"]
  htwo3<- genvar3/envvar3
  print(htwo3)
  return(htwo3)
  
}
herit_PC<-lapply(c("Dim.1", "Dim.2"), myherit_BS) %>% as.data.frame() %>% t
rownames(herit_PC)<-c("Dim.1", "Dim.2")
herit_PC
write.table(herit_PC,file=paste0(odir, "/phenos_modelled/heritabilities_BS_PC1_PC2.txt"), quote=F)

WhichCOL<-c(1:nrow(phenos_load))[-c(lapply(families, function(x) grep(x,phenos_load$Name))  %>% unlist)]
length(WhichCOL)

myheritCOLL_BS_PC<-function(trait){
  data=phenos_load[WhichCOL,c("Name", trait)] %>% droplevels
  mylm_sel3<-lmer(data=data, get(trait) ~ (1|Name) )
  myvar3<-VarCorr(mylm_sel3,comp=c("Variance")) %>% as.data.frame()
  genvar3<-myvar3[which(myvar3$grp=="Name"), "vcov"]
  envvar3<-myvar3[which(myvar3$grp=="Name"), "vcov"]+ myvar3[which(myvar3$grp=="Residual"), "vcov"]
  htwo3<- genvar3/envvar3
  print(htwo3)
  return(htwo3)
}
herit_COLL_PC<-lapply(c("Dim.1", "Dim.2"), myheritCOLL_BS_PC) %>% as.data.frame() %>% t
rownames(herit_COLL_PC)<-c("Dim.1", "Dim.2")
herit_COLL_PC
write.table(herit_PC,file=paste0(odir, "/phenos_modelled/heritabilities_BS_COLL_PC1_PC2.txt"), quote=F)



# Purge obsolete variables
rm(phenos_melt, phenos_sub,mm, WhichCOL,   whichParents,   whichFAM)

cat("Phenos analyzed and ready for predictions!\n")
