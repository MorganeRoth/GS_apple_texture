cat("Modelling data...\n")

# Variance Analysis with MM4LMM

dir.create(paste0(odir, "/variance_analysis"), showWarnings = FALSE, recursive = TRUE)

## In a first step, only with additive kinship on raw data
phenos_raw %>% summary
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
Dist0<-hetero[ids,]
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
names(ZList)<-names(VList)<-c("Add","Dom", "AA", "AD", "DD", "Error")

## trial effect (year*location)

trial<-phenos_raw$Trial

# models<-c("A","A_Dist", "AD", "AD_Dist")

## create matrix for results
models="A"
res_mat_years<-matrix(NA, nrow=length(traits)*length(models), ncol=5)
colnames(res_mat_years)<-c( "trait", "model",  "Add", "Error",  "LogLik")
res_mat_years[,"trait"]<-rep(traits,length(models))
res_mat_years[,"model"]<-unlist(lapply(models, function(x) rep(x, length(traits))))

res<-list()
# Dist<-merge(Het,phenos, by.x="names", by.y="Name", all.y=T)["Het"]

phenos_raw$Year<-as.character(phenos_raw$Year)

for (trait in traits) {
  NonMissing<-which(!(is.na(phenos_raw[,trait])))
  Zloc<-lapply(ZList, function(x) x[NonMissing,])
  res[[trait]][["A"]]<-MMEst(Y=phenos_raw[NonMissing,trait],ZList= Zloc[c("Add", "Error")], VarList =VList[c("Add", "Error")] , Cofactor=NULL)
  coefs<-res[[trait]][["A"]]$NullModel$Sigma2
  LogLik<-res[[trait]][["A"]]$NullModel$`LogLik (Reml)`
  res_mat_years[which(res_mat_years[,"trait"]==trait & res_mat_years[,"model"]=="A"), c("Add",  "Error")]<-as.numeric(coefs)
  res_mat_years[which(res_mat_years[,"trait"]==trait & res_mat_years[,"model"]=="A"), "LogLik"]<-as.numeric(LogLik)
  print(as.numeric(coefs))
}

## Plot results

# Plot histograms with all variance components

res<-as.data.frame(res_mat_years)
melted<-melt(res, measure.vars=c("Error", "Add") ) 
melted$trait<-as.factor(melted$trait)
melted$model<-as.factor(melted$model)
melted$value<-as.numeric(melted$value)

png(file=paste0(odir, "/variance_analysis/variance_decomp_additive_only.png"), height=400, width=600)
gg<-ggplot(data=melted, aes(fill=variable, y=value, x=trait)) + 
    geom_bar( stat="identity", position="fill")+
    labs(x="Trait and model", title="", y="Proportion")+
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

head(res)
max_lk<-lapply(traits, function(x) {
  data<-res[res$trait==x, ]
  mymin<-data[which(data$LogLik == max(data$LogLik)),"model"] %>% as.character
  return(mymin)
})
max_lk<-cbind(traits, unlist(max_lk)) ## ADI_Dist is always the best model! Use later in predictions
write.table(max_lk, file=paste0(odir, "/variance_analysis/max_likelihood_each_trait_txt"), quote=F, sep="\t")

# Purge obsolete variables
rm(melted,res_mat, res_mat_years, res )

cat("Variance decomposition done for models A, A_Dist, AD, AD_Dist, ADI_Dist !\n")
