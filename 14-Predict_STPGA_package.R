cat("Predicting with STPGA package...\n")

dir.create(paste0(odir, "/predictions/STPGA"), showWarnings = FALSE, recursive = TRUE)

## if import, model phenos, mnodel genos scripts are skipped, load data here
id_pheno<-read.table(paste0(odir, "/phenos_modelled/rownames_phenos.txt"))
## problem with duplicated rowname that I do not understand, use id_pheno to replace them (saved with phenos)
phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T, row.names = id_pheno$x %>% as.character())
phenos<-phenos[,-1]
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))
traits=colnames(phenos)
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
NbFAM<-length(families)


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

## kinship
A<-read.table(paste0(odir, "/genos_modelled/Ka_Amat.txt"),h=T )
A[1:5,1:5]
A<-A[ids,ids]

# clusters<-read.table(paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"), h=T)
# cluster_fams<-read.table(paste0(odir, "/genos_modelled/assignements_families.txt"), h=T)

## assignements
## for running the optimization, either PCA-like coordinates a supplied, or an inverted kinship matrix
## I choose to import the DAPC coordinates for collection and families
col_assi<-read.table(paste0(odir, "/genos_modelled/coordinates_DAPC_COL.txt"),h=T) 
fam_assi<-read.table(paste0(odir, "/genos_modelled/coord_DAPC_families.txt"),h=T) 
assi<-rbind(col_assi,fam_assi)
summary(assi)

WhichCOL<-ids[-c(lapply(families, function(x) grep(x, ids) ) %>% unlist)]


################################################################################
### STEP 1 simplified to avoid crashes: no initialization ######################
################################################################################

repeatgenalg0<-function(ntrain, criteria, G,Test, Candidates){
  ListTrain<-GenAlgForSubsetSelection(P=G,Candidates=Candidates, 
                                      Test=Test,ntoselect=ntrain, 
                                      errorstat=criteria, mc.cores=1)
  return(ListTrain)
}

FAM=families[1]
trait=traits[1]
criteria="CDMEAN"

## very important: P parameter must be a matrix !!!!!!!!!!!!!!

pred_fam_trait0 <-function(FAM, trait, criteria) {
  ThisFAM<-ids[ grep(FAM, ids)]
  WhichCOLFAM<-c(WhichCOL, ThisFAM)%>% sort
  # A_COLFAM<-assi[WhichCOLFAM, ] %>% as.matrix
  G<-assi[WhichCOLFAM,]
  Test<-rownames(G)[grep(FAM, rownames(G))  %>% unlist]
  Candidates<- setdiff(rownames(G), Test)
  phenos_sub<-phenos[WhichCOLFAM,]
  phenos_sub$id<-rownames(phenos_sub)
  to_test<- seq(10,(length(WhichCOL)),10) ## TRS sizes to test increment of 10
  ListTrain<-lapply(to_test, function(ntrain){ 
    print(ntrain)
    repeatgenalg0(ntrain, criteria,as.matrix(G), Test, Candidates)
  } )
  ##predictions by optimized sample
  
  deptestopt<-phenos_sub[phenos_sub$id%in%Test,]
  res<-matrix(NA, ncol=2, nrow=length(to_test), dimnames=list(to_test, c("size", "corr")))
  count=0
  for (x in to_test){
    count=count+1
    deptrainopt<-phenos_sub[phenos_sub$id%in%ListTrain[[count]][[1]],]
    deptrainopt$id<-factor(deptrainopt$id, levels=c(rownames(deptrainopt), rownames(deptestopt)))
    deptestopt$id<-factor(deptestopt$id, levels=c(rownames(deptrainopt), rownames(deptestopt)))
    Ztrain<-model.matrix(~-1+deptrainopt$id)
    Ztest<-model.matrix(~-1+deptestopt$id)
    K=A[c(rownames(deptrainopt), rownames(deptestopt)),c(rownames(deptrainopt), rownames(deptestopt))]
    modelopt<-emmreml(y=deptrainopt[,trait],X=matrix(1, nrow=nrow(deptrainopt), ncol=1), 
                      Z=Ztrain, K=as.matrix(K))
    predictopt<-Ztest%*%modelopt$uhat
    CORR<-cor(predictopt, deptestopt[,trait])
    res[count,]<-c(x, as.numeric(CORR) )
  }
  return(res)
}

res_cor=list()
for (i in traits) {
  for (j in families) {
    # xx<-pred_fam_trait0(i, j)
    res_cor[[i]][[j]]<-pred_fam_trait0(j, i, "PEVMEAN")
   
  }
}
saveRDS(res_cor, file=paste0(odir,  "/predictions/STPGA/opt_PEVMEAN.rds"))
par(mfrow=c(2,3))

par(mfrow=c(2:3), mar=c(rep(5,4)))

myplot <- function(trait){
  lapply(families, function(x) {
    data=res_cor[[trait]][[x]]
    print(summary(data))
    plot(data=data, corr ~ size,
         xlab="size_TRS", ylab="", main=paste(trait, x ), ylim=c(-0.5,1), cex=0.2) %>% print
    axis(side = 2) %>% print
    mtext("correlation", side = 2, line = 3, cex=0.75) %>% print
    # par(new = TRUE)
    # plot(data=data, mean_rel %>% as.character %>% as.numeric ~ size_TRS  %>% as.character %>% as.numeric,
    #      type = "l", xaxt = "n", yaxt = "n",
    #      xlab="", ylab="", main=paste(trait, x ), col="red", ylim=c(-0.09,0.25))  %>% print
    # axis(side = 4) %>% print
    # mtext("mean_relatedness", side = 4, line = 3, cex=0.75, col="red") %>% print
    # abline(a=0, b=0)
  })
}
pdf("trial.pdf")
myplot("PC1") %>% print
dev.off()

res_format<-reshape2::melt(res_cor)
head(res_format)
tail(res_format)
summary(res_format)
res_format<-dcast(res_format,L1+L2+Var1~Var2)
res_format<-res_format[,-3]
colnames(res_format)<-c("trait", "fam", "size", "corr")
head(res_format)
res_CDmean<-res_format
write.table(res_CDmean, file=paste0(odir, "/predictions/STPGA/all_res_optimisation_CDmean.txt"), sep="\t", quote=F)

## compare to opt with relatedness
sel_traits<- traits[c(3,11,13,14)]
par(mfrow=c(2,3))
res_rel<-readRDS("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/20190816/predictions/COLLtoFAMs_optim_kinship/A.mat_TRS_opt_all_traits.rds")
lapply(sel_traits,function(trait) {
  # pdf(file=paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/", trait, "_TRS_opt_size.pdf"), height=5, width=8)
  data=res_rel[[trait]] %>% as.data.frame()
  print(summary(data))
  data$mean_rel<-data$mean_rel %>% as.character %>% as.numeric
  data$accuracy<-data$accuracy %>% as.character %>% as.numeric
  data$size_TRS<-data$size_TRS %>% as.character %>% as.numeric
  res_opt<-res_format[which(res_format$trait==trait),] %>% droplevels
  res_opt$model<-"CDmean"
  data$model<-"Rel"
  data<-data[,-3]
  res_opt<-res_opt[,c(3,2,4,5)]
  colnames(res_opt)<-colnames(data)
  data2<-rbind(res_opt, data)
  gg<-ggplot(data=data2, aes(y= accuracy, x= size_TRS, col=model)) + 
    geom_line(aes(linetype=model))+
    facet_wrap(~family, ncol=3)+
    geom_point()+
    scale_color_manual(values=c("#999999", "#E69F00"))+
    theme_bw()
  print(gg)                  
})
## End(Not run)



################################################################################
### STEP 2 find the optimal training pop with a given number of individuals ####
## CRASHES !!! #################################################################
################################################################################


###instead of using the algorithm directly using a wrapper to 
###implement an for multiple starting points for genetic algorithm.


repeatgenalg<-function(numrepsouter,numrepsinner, ntrain){
  StartingPopulation2=NULL 
  for (i in 1:numrepsouter){
    print("Rep:")
    print(i)
    StartingPopulation<-lapply(1:numrepsinner, function(x){
      GenAlgForSubsetSelection(P=G,Candidates=Candidates, 
                               Test=Test, ntoselect=ntrain, InitPop=StartingPopulation2,
                               mutprob=.5, mutintensity = rpois(1,4),
                               niterations=10,minitbefstop=1, tabumemsize = 2,plotiters=F, 
                               lambda=1e-9,errorstat=model, mc.cores=5)})
    StartingPopulation2<-vector(mode="list", length = numrepsouter*1)
    ij=1
    for (i in 1:numrepsinner){
      for (j in 1:1){
        StartingPopulation2[[ij]]<-StartingPopulation[[i]][[j]]
        ij=ij+1
      }
    }
  }
  ListTrain<-GenAlgForSubsetSelection(P=G,Candidates=Candidates, 
                                      Test=Test,ntoselect=ntrain, InitPop=StartingPopulation2,npop=100, 
                                      mutprob=.5, mutintensity = 1,niterations=100,
                                      minitbefstop=50, tabumemsize = 1,plotiters=T,
                                      lambda=1e-9,errorstat="CDMEAN", mc.cores=5)
  return(ListTrain)
}
### Now repeat this for each value of ntrain from 10, increment of 10 until 210
## do this for each family
traits=colnames(phenos)

pred_fam_trait <-function(FAM, trait) {
  res_cor<-list()
  ThisFAM<-ids[ grep(FAM, ids)]
  WhichCOLFAM<-c(WhichCOL, ThisFAM)%>% sort
  assiCOLFAM<-assi[WhichCOLFAM, ]
  Test<-rownames(assiCOLFAM)[grep(FAM, rownames(assiCOLFAM))  %>% unlist]
  Candidates<- setdiff(rownames(assiFAM), Test)
  phenos_sub<-phenos[WhichCOLFAM,]
  phenos_sub$id<-rownames(phenos_sub)
  to_test<- seq(10,(length(WhichCOL)),10) ## TRS sizes to test increment of 10
  ListTrain<-lapply(to_test, function(ntrain){ 
    print(ntrain)
    repeatgenalg(5,3,ntrain, )
  } )
  ##predictions by optimized sample
  
  deptestopt<-phenos_sub[phenos_sub$id%in%Test,]
  
  res_cor[[trait]][[FAM]]<-matrix(NA, ncol=2, nrow=length(to_test), dimnames=list(to_test, c("size", "corr")))
  lapply(1:length(to_test), function(x){
    deptrainopt<-phenos_sub[phenos_sub$id%in%ListTrain[[x]][[1]],]
    deptrainopt$id<-factor(deptrainopt$id, levels=c(rownames(deptrainopt), rownames(deptestopt)))
    deptestopt$id<-factor(deptestopt$id, levels=c(rownames(deptrainopt), rownames(deptestopt)))
    Ztrain<-model.matrix(~-1+deptrainopt$id)
    Ztest<-model.matrix(~-1+deptestopt$id)
    K=A[c(rownames(deptrainopt), rownames(deptestopt)),c(rownames(deptrainopt), rownames(deptestopt))]
    modelopt<-emmreml(y=deptrainopt[,trait],X=matrix(1, nrow=nrow(deptrainopt), ncol=1), 
                      Z=Ztrain, K=K)
    predictopt<-Ztest%*%modelopt$uhat
    CORR<-cor(predictopt, deptestopt[,trait])
    res_cor[[trait]][[FAM]][x,]<<-c(to_test[x], CORR )
    return(res_cor)
  })
}

res<-lapply(families,function(x) pred_fam_trait(x, "PC1"))

## sanity check here:
summary(rownames(phenos) == ids)

deptestopt<-phenos[rownames(phenos)%in%Test,]

##predictions by optimized sample
deptrainopt<-phenos[(rownames(phenos)%in%ListTrain[[1]]),]

## need a factor for the model
deptrainopt$id<-factor(rownames(deptrainopt), levels=c(rownames(deptrainopt), rownames(deptestopt)))
deptestopt$id<-factor(rownames(deptestopt), levels=c(rownames(deptrainopt), rownames(deptestopt)))

Ztrain<-model.matrix(~-1+deptrainopt$id)
Ztest<-model.matrix(~-1+deptestopt$id)

K=A[c(rownames(deptrainopt), rownames(deptestopt)),c(rownames(deptrainopt), rownames(deptestopt))]

colnames(Ztrain)
modelopt<-emmreml(y=deptrainopt$PC1,X=matrix(1, nrow=nrow(deptrainopt), ncol=1), 
                  Z=Ztrain, K=K)

modelopt %>% attributes

predictopt<-Ztest%*%modelopt$uhat

corvecrs<-c()
for (rep in 1:50){
  ###predictions by a random sample of the same size
  rs<-sample(Candidates, 50)
  
  deptestrs<-phenos[ids%in%Test,]
  
  deptrainrs<-phenos[(ids%in%rs),]
  
  deptrainrs$id<-factor(rownames(deptrainrs), levels=c(rownames(deptrainrs), rownames(deptestrs)))
  deptestrs$id<-factor(rownames(deptestrs), levels=c(rownames(deptrainrs), rownames(deptestrs)))
  
  Ztrain<-model.matrix(~-1+deptrainrs$id)
  Ztest<-model.matrix(~-1+deptestrs$id)
  
  modelrs<-emmreml(y=deptrainrs$PC1,X=matrix(1, nrow=nrow(deptrainrs), ncol=1), 
                   Z=Ztrain, K=K)
  predictrs<-Ztest%*%modelrs$uhat
  corvecrs<-c(corvecrs,cor(predictrs, deptestrs$PC1))
  
}
mean(corvecrs)
cor(predictopt, deptestopt$PC1)

##1 is black, 2 is red: red is training set optimized: all belonging to Cluster 5
plot(assi[,4],assi[,5], col=rownames(assi)%in%ListTrain[[1]]+1,
     pch=2*rownames(assi)%in%Test+1, xlab="Cluster4", ylab="Cluster5")

