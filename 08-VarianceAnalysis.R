cat("Modelling data...\n")

# Variance Analysis with MM4LMM

dir.create(paste0(odir, "/variance_analysis"), showWarnings = FALSE, recursive = TRUE)

## In a first step, only with additive kinship on raw data
phenos_raw<-readRDS(paste0(idir, "/phenos_with_reps_raw_filtered.rds"))
families<-c("FjDe", "GDFj", "FjPi", "FjPL", "GaPi", "GaPL")
## load kinships
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))
geno_add<-round(genos_ready+1, digits=0)
Kvr<-estimGenRel(geno_add,method="vanraden1", relationships = "additive" )
K.Add<-read.table(paste0(odir, "/genos_modelled/Ka_Vitezica17.txt"), h=T)
K.Add[1:5,1:5]
K.Dom<-read.table(paste0(odir, "/genos_modelled/Kd_Vitezica17.txt"), h=T)
K.AA<-read.table(paste0(odir, "/genos_modelled/K.AA_Vitezica17.txt"), h=T)
K.AD<-read.table(paste0(odir, "/genos_modelled/K.AD_Vitezica17.txt"), h=T)
K.DD<-read.table(paste0(odir, "/genos_modelled/K.DD_Vitezica17.txt"), h=T)
hetero<-read.table(paste0(odir, "/genos_modelled/heterozygosity_mean.txt"), h=T)

# Loop over traits and keep all trials
traits<-colnames(phenos_raw)[-c(1:7)]
all_trials<-levels(phenos_raw$Trial)
NbEnv<-length(all_trials)
# define models to be tested
models<-c("A")
Effects<-c("Add", "Dom", "AA", "AD", "DD", "Error")

## filter out individuals without genotype
cat("intersect between phenos and genos data")
intersect(levels(phenos_raw$Name) , colnames(K.Add)) %>% length
levels(phenos_raw$Name) [which(levels(phenos_raw$Name) %in% setdiff(levels(phenos_raw$Name) , colnames(K.Add)))]
write.table( setdiff(levels(phenos_raw$Name) , colnames(K.Add)), file=paste0(odir, "/phenos_modelled/missing_genos_to_check.txt"), sep="\t", quote=F)
ids<-intersect(phenos_raw$Name , colnames(K.Add)) ## update ids
phenos_raw<-phenos_raw[which(phenos_raw$Name %in% ids),]
phenos_raw<-droplevels(phenos_raw)
phenos_raw<-phenos_raw[order(phenos_raw$Name),]
dim(phenos_raw)
## filter out individuals never phenotyped
filter_order_K<-function(K) {
  NAME<-substitute(K) ## very important step to not evaluate K and use it as a name after
  K<-K[which(rownames(K) %in% ids),which(colnames(K) %in% ids)]
  K<-K[order(rownames(K)), order(colnames(K))]
    assign(deparse(NAME), K, .GlobalEnv)
}

# for ( K in list(K.Add, K.Dom, K.AA, K.AD, K.DD)) {filter_order_K(K)} ## loop not working somehow
filter_order_K(K.Add)
filter_order_K(K.Dom)
filter_order_K(K.AA)
filter_order_K(K.AD)
filter_order_K(K.DD)
for ( K in list(K.Add, K.Dom, K.AA, K.AD, K.DD)) { print(dim(K))}

## reorder ids
ids=sort(ids)
## Distance between individuals
Dist0<-data.frame(names=ids, Het=hetero[ids,])
Dist<-merge(phenos_raw, Dist0, by.x="Name", by.y="names")[,"Het"]
# For using MM4LMM we need to prepare variance and incidence matrices that will be filtered for missing values afterwards

# define useful parameters
NbInd<-nrow(K.Add)
NbObs<-nrow(phenos_raw)
## Variance matrices

## Incidence matrices

Zadd<-Zdom<-ZAA<-ZAD<-ZDD<-class.ind(phenos_raw$Name)
Zres<-diag(1, NbObs)
ZList<-list(Zadd, Zdom,ZAA, ZAD, ZDD, Zres)

## variance matrices
VList=list(K.Add, K.Dom, K.AA, K.AD, K.DD, diag(1, NbObs))
# VList=list(Kvr, K.Dom, K.AA, K.AD, K.DD, diag(1, NbObs))
VList<-lapply(VList, function(x) as.matrix(x))
names(ZList)<-names(VList)<-c("Add","Dom", "AA", "AD", "DD", "Error")

## trial effect (year*location)

trial<-phenos_raw$Trial

# models<-c("A","A_Dist", "AD", "AD_Dist")

## create matrix for results
models="A"
res_mat_years<-matrix(NA, nrow=length(traits)*length(models), ncol=7)
colnames(res_mat_years)<-c( "trait", "model", "mean", "sd", "Add", "Error",  "LogLik")
res_mat_years[,"trait"]<-rep(traits,length(models))
res_mat_years[,"model"]<-unlist(lapply(models, function(x) rep(x, length(traits))))

res<-list()
# Dist<-merge(Het,phenos, by.x="names", by.y="Name", all.y=T)["Het"]

phenos_raw$Year<-as.character(phenos_raw$Year)

for (trait in traits) {
  NonMissing<-which(!(is.na(phenos_raw[,trait])))
  Zloc<-lapply(ZList, function(x) x[NonMissing,])
  Trial<-class.ind(phenos_raw$Trial)[NonMissing,]
  res[[trait]][["A"]]<-MMEst(Y=phenos_raw[NonMissing,trait],ZList= Zloc[c("Add", "Error")], VarList =VList[c("Add", "Error")] , Cofactor=Trial)
  coefs<-res[[trait]][["A"]]$NullModel$Sigma2
  LogLik<-res[[trait]][["A"]]$NullModel$`LogLik (Reml)`
  res_mat_years[which(res_mat_years[,"trait"]==trait & res_mat_years[,"model"]=="A"), c("Add",  "Error")]<-as.numeric(coefs)
  res_mat_years[which(res_mat_years[,"trait"]==trait & res_mat_years[,"model"]=="A"), "LogLik"]<-as.numeric(LogLik)
  res_mat_years[which(res_mat_years[,"trait"]==trait & res_mat_years[,"model"]=="A"), c("mean", "sd")] <- c(mean(phenos_raw[NonMissing,trait]), sd(phenos_raw[NonMissing,trait]))
  print(as.numeric(coefs))
}

write.table(res_mat_years, file=paste0(odir, "/variance_analysis/variance_decomp_with_trial_effect_fix.txt"), quote=F, sep="\t")
# write.table(res_mat_years, file=paste0(odir, "/variance_analysis/variance_decomp_with_trial_effect_fix_Kvanraden.txt"), quote=F, sep="\t")

## plot# Plot histograms with all variance components

res<-as.data.frame(res_mat_years)
melted<-melt(res, measure.vars=c("Error", "Add") ) 
melted$trait<-as.factor(melted$trait)
melted$model<-as.factor(melted$model)
melted$value<-as.numeric(melted$value)
melted$trait<-revalue(melted$trait, c("Acoustic_Linear_Distance"="ALD",
                                      "Acoustic_Max_Pressure"="APMax",
                                      "Acoustic_Mean_Pressure" ="APMean",
                                      "Acoustic_Npeak"="ANP",
                                      "Area"="Area",
                                      "Final_Force" ="FF",
                                      "Force_Linear_Distance" ="FLD" ,
                                      "Initial_Force" = "FI",
                                      "Max_Force" = "FMax",
                                      "Mean_Force"= "FMean", 
                                      "N_Peak_Force" ="FNP",
                                      "Young_Module" = "YM"))

pdf(file=paste0(odir, "/variance_analysis/variance_decomp_additive_only_trial_fix_all_pop_vanraden.pdf"), height=4, width=6)
gg<-ggplot(data=melted, aes(fill=variable, y=value, x=trait)) + 
  geom_bar( stat="identity", position="fill")+
  scale_fill_manual(values=c("grey70", "grey40"),name = "Variance", labels = c("Residual", "Additive"))+ 
  labs(x="Trait", title="", y="Proportion")+
  theme_minimal()+ 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1))
print(gg)
dev.off()

## do the same excluding the families
## now only the collection - with year effect FIX
## create matrix for results
models="A"
res_mat_coll<-matrix(NA, nrow=length(traits)*length(models), ncol=7)
colnames(res_mat_coll)<-c( "trait", "model", "mean", "sd", "Add", "Error",  "LogLik")
res_mat_coll[,"trait"]<-rep(traits,length(models))
res_mat_coll[,"model"]<-unlist(lapply(models, function(x) rep(x, length(traits))))

res<-list()
# Dist<-merge(Het,phenos, by.x="names", by.y="Name", all.y=T)["Het"]

phenos_raw$Year<-as.character(phenos_raw$Year)
whichCOL<-c(1:nrow(phenos_raw))[-c(lapply(families, function(x) grep(x, phenos_raw$Name)) %>% unlist)]

for (trait in traits) {
  NonMissing<-c(which(!(is.na(phenos_raw[,trait]))), whichCOL) %>% sort
  Zloc<-lapply(ZList, function(x) x[NonMissing,])
  Year<-class.ind(phenos_raw$Year)[NonMissing,]
  res[[trait]][["A"]]<-MMEst(Y=phenos_raw[NonMissing,trait],ZList= Zloc[c("Add", "Error")], VarList =VList[c("Add", "Error")] , Cofactor=NULL)
  coefs<-res[[trait]][["A"]]$NullModel$Sigma2
  LogLik<-res[[trait]][["A"]]$NullModel$`LogLik (Reml)`
  res_mat_coll[which(res_mat_coll[,"trait"]==trait & res_mat_coll[,"model"]=="A"), c("Add",  "Error")]<-as.numeric(coefs)
  res_mat_coll[which(res_mat_coll[,"trait"]==trait & res_mat_coll[,"model"]=="A"), "LogLik"]<-as.numeric(LogLik)
  res_mat_coll[which(res_mat_coll[,"trait"]==trait & res_mat_coll[,"model"]=="A"), c("mean", "sd")] <- c(mean(phenos_raw[NonMissing,trait]), sd(phenos_raw[NonMissing,trait]))
  print(as.numeric(coefs))
}

# write.table(res_mat_coll, file=paste0(odir, "/variance_analysis/variance_decomp_coll_only_trial_effect.txt"), quote=F, sep="\t")
write.table(res_mat_coll, file=paste0(odir, "/variance_analysis/variance_decomp_coll_only_trial_effect_Kvanraden.txt"), quote=F, sep="\t")

# Plot histograms with all variance components

res<-as.data.frame(res_mat_coll)
melted<-melt(res, measure.vars=c("Error", "Add") ) 
melted$trait<-as.factor(melted$trait)
melted$model<-as.factor(melted$model)
melted$value<-as.numeric(melted$value)
melted$trait<-revalue(melted$trait, c("Acoustic_Linear_Distance"="ALD",
                                "Acoustic_Max_Pressure"="APMax",
                                "Acoustic_Mean_Pressure" ="APMean",
                                "Acoustic_Npeak"="ANP",
                                "Area"="Area",
                                "Final_Force" ="FF",
                                "Force_Linear_Distance" ="FLD" ,
                                "Initial_Force" = "FI",
                                "Max_Force" = "FMax",
                                "Mean_Force"= "FMean", 
                                "N_Peak_Force" ="FNP",
                                "Young_Module" = "YM"))

pdf(file=paste0(odir, "/variance_analysis/variance_decomp_additive_only_trial_fix.pdf"), height=4, width=6)
gg<-ggplot(data=melted, aes(fill=variable, y=value, x=trait)) + 
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual(values=c("grey70", "grey40"),name = "Variance", labels = c("Residual", "Additive"))+ 
    labs(x="Trait", title="", y="Proportion")+
    theme_minimal()+ 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(angle=45, vjust=1, hjust=1))
print(gg)
dev.off()




## Here for more sophisticated variance decomp allowing to test all models

MM4_model<-function(phenos_raw, trait, model, effects_sub, DistYES_NO) {
  if (DistYES_NO == "NO") {Dist_sel=NULL} else { Dist_sel<-  Dist[NonMissing]}
  res[[trait]][[model]]<-MMEst(Y=phenos_raw[NonMissing,trait],ZList= Zloc[effects_sub], VarList =VList[effects_sub] , Cofactor=Dist_sel)
  coefs<-res[[trait]][[model]]$NullModel$Sigma2
  LogLik<-res[[trait]][[model]]$NullModel$`LogLik (Reml)`
  res_mat[which(res_mat[,"trait"]==trait & res_mat[,"model"]==model), effects_sub]<<-as.numeric(coefs)
  res_mat[which(res_mat[,"trait"]==trait & res_mat[,"model"]==model), "LogLik"]<<-as.numeric(LogLik)
  return(res_mat)
}

## reinitialize models
models=c("A", "A_Dist", "AD", "AD_Dist", "ADI_Dist")
Effects<- list(c("Add","Error"),c("Add","Error"), c("Add","Dom", "Error"),  c("Add","Dom", "Error"),  c("Add","Dom", "AA", "AD", "DD", "Error") )
dist_all<-list("NO", "YES", "NO", "YES", "YES")
names(dist_all)<-names(Effects)<-models

## do in a loop
res_mat<-matrix(NA, nrow=length(traits)*length(models), ncol=9)
colnames(res_mat)<-c( "trait", "model",  "Add","Dom", "AA", "AD", "DD" ,"Error", "LogLik")
res_mat[,"trait"]<-rep(traits,length(models))
res_mat[,"model"]<-unlist(lapply(models, function(x) rep(x, length(traits))))

for (model in models) {
  for (trait in traits ) {
    print(model)
    res<-list()
    print( Effects[[model]])
    print( dist_all[[model]])
    MM4_model(phenos_raw, trait, model, Effects[[model]], dist_all[[model]])
  }
}

## plot results

res<-as.data.frame(res_mat)
melted<-melt(res, measure.vars=c("Error", "Add","Dom", "AA", "AD", "DD") ) 
head(melted)
melted$trait<-as.factor(melted$trait)
melted$model<-as.factor(melted$model)
melted$value<-as.numeric(melted$value)
# melted$variable<-factor(melted$variable, levels=c("Error", "Add","Dom", "AA", "AD", "DD"))
# melted2<-melted %>% arrange (variable)
# head(melted2)
png(file=paste0(odir, "/variance_analysis/variance_decomp_ADI_models.png"), height=500, width=1000)
gg<-ggplot(data=melted, aes(fill=variable, y=value, x=model)) + 
  geom_bar( stat="identity", position="fill")+
  facet_grid(~trait) +
  scale_fill_manual(values=c('#e6ab02','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e'))+
  labs(x="Trait and model", title="", y="Proportion")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1))
print(gg)
dev.off()

write.table(melted, file=paste0(odir, "/variance_analysis/results_var_decomp_all_models_melted.txt"), sep="\t", quote=F )
write.table(res, file=paste0(odir, "/variance_analysis/results_var_decomp_all_models.txt"), sep="\t", quote=F )

## Find best model for each trait
# res<-read.table(paste0(odir, "/variance_analysis/results_var_decomp_all_models.txt"), h=T)

summary(res)
res$LogLik<-as.numeric(res$LogLik)
max_lk<-lapply(traits, function(x) {
  data<-res[res$trait==x, ]
  mymin<-data[which(data$LogLik == max(data$LogLik)),"model"] %>% as.character
  return(mymin)
})
max_lk<-cbind(traits, unlist(max_lk)) ## ADI_Dist is always the best model! Use later in predictions
write.table(max_lk, file=paste0(odir, "/variance_analysis/max_likelihood_each_trait_txt"), quote=F, sep="\t")

## genotypic variance per family

## create matrix for results
models="A"
res_mat_FAM<-matrix(NA, nrow=length(traits)*length(models)*length(families), ncol=6)
colnames(res_mat_FAM)<-c("family", "trait", "model",  "Add", "Error",  "LogLik")
res_mat_FAM[,"trait"]<-rep(traits,length(models)*length(families))
res_mat_FAM[,"model"]<-unlist(lapply(models, function(x) rep(x, length(traits)*length(families))))
res_mat_FAM[,"family"]<-lapply(families, function(x) rep(x, length(traits)*length(models))) %>% unlist

res<-list()
for (fam in families) {
  for (trait in traits) {
    NonMissing<-which(!(is.na(phenos_raw[,trait])) & phenos_raw$Family == fam)
    Zloc<-lapply(ZList, function(x) x[NonMissing,])
    Year<-phenos_raw[NonMissing,"Year"] %>% as.factor() %>% droplevels
    Year<-class.ind(Year)
    res[[trait]][["A"]]<-MMEst(Y=phenos_raw[NonMissing,trait],ZList= Zloc[c("Add", "Error")], VarList =VList[c("Add", "Error")] , Cofactor=Year)
    coefs<-res[[trait]][["A"]]$NullModel$Sigma2
    LogLik<-res[[trait]][["A"]]$NullModel$`LogLik (Reml)`
    res_mat_FAM[which(res_mat_FAM[,"trait"]==trait & res_mat_FAM[,"model"]=="A" & res_mat_FAM[,"family"]==fam) , c("Add",  "Error")]<-as.numeric(coefs)
    res_mat_FAM[which(res_mat_FAM[,"trait"]==trait & res_mat_FAM[,"model"]=="A"& res_mat_FAM[,"family"]==fam), "LogLik"]<-as.numeric(LogLik)
    print(as.numeric(coefs))
  }
    
}

write.table(res_mat_FAM, file=paste0(odir, "/phenos_modelled/variance_per_family_year.txt"), sep="\t", quote=F, row.names=F)


res<-as.data.frame(res_mat_FAM)
melted<-melt(res, measure.vars=c("Error", "Add") ) 
melted$trait<-as.factor(melted$trait)
melted$model<-as.factor(melted$model)
melted$value<-as.numeric(melted$value)
head(melted)

melted$trait<-revalue(melted$trait, c("Acoustic_Linear_Distance"="ALD",
                              "Acoustic_Max_Pressure"="APMax",
                              "Acoustic_Mean_Pressure" ="APMean",
                              "Acoustic_Npeak"="ANP",
                              "Area"="Area",
                              "Final_Force" ="FF",
                              "Force_Linear_Distance" ="FLD" ,
                              "Initial_Force" = "FI",
                              "Max_Force" = "FMax",
                              "Mean_Force"= "FMean", 
                              "N_Peak_Force" ="FNP",
                              "Young_Module" = "YM"))
melted$family<-revalue(melted$family, c("FuPi"="FjPi"))

pdf(file=paste0(odir, "/phenos_modelled/variance_per_family.pdf"), height=6, width=12)
gg<-ggplot(data=melted, aes(fill=variable, y=value, x=family)) + 
  geom_bar( stat="identity", position="fill")+
  facet_grid(~trait)+
  labs(x="Trait and model", title="", y="Proportion")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=8))

print(gg) 

dev.off()

xx<-dcast(phenos_raw, value.var="Acoustic_Max_Pressure", Name~Year, mean)
plot(data=xx[grep("FjPL",xx$Name),], `2012`~ `2013`)
par(mfrow=c(3,4))
traits=colnames(phenos_raw)[-c(1:7)]
lapply(traits, function(x) {
  data<-dcast(phenos_raw, value.var=x, Name~Year, mean)
  plot(data=xx[grep("FjPL",xx$Name),], `2012`~ `2013`, main=x) %>% print
})


# Purge obsolete variables
rm(melted,res_mat, res_mat_years, res )

cat("Variance decomposition done for models A, A_Dist, AD, AD_Dist, ADI_Dist !\n")
