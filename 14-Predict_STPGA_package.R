cat("Predicting with STPGA package...\n")

dir.create(paste0(odir, "/predictions/STPGA"), showWarnings = FALSE, recursive = TRUE)

## if import, model phenos, mnodel genos scripts are skipped, load data here
## problem with duplicated rowname that I do not understand, use id_pheno to replace them (saved with phenos)
phenos<-read.table(paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)
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
### STEP 2 find the optimal training pop with a given number of individuals ####
################################################################################


###instead of using the algorithm directly using a wrapper to 
###implement an for multiple starting points for genetic algorithm.


repeatgenalg<-function(numrepsouter,numrepsinner,ntrain, criteria, G,Test, Candidates){
  StartingPopulation2=NULL 
  for (i in 1:numrepsouter){
    print("Rep:")
    print(i)
    StartingPopulation<-lapply(1:numrepsinner, function(x){
      GenAlgForSubsetSelection(P=G,Candidates=Candidates, 
                               Test=Test, ntoselect=ntrain, InitPop=StartingPopulation2,
                               # mutprob=.5, mutintensity = rpois(1,4),
                               niterations=20,minitbefstop=5, tabumemsize = 2,plotiters=F, 
                               lambda=1e-9,errorstat=criteria, mc.cores=8)})
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
                                      lambda=1e-9,errorstat=criteria, mc.cores=5)
  return(ListTrain)
}
### Now repeat this for each value of ntrain from 10, increment of 20 until 259
## do this for each family
## the same TRS is used for all traits

pred_fam_trait2 <-function(FAM, criteria, increment) {
  res_cor<-list()
  ThisFAM<-ids[ grep(FAM, ids)]
  WhichCOLFAM<-c(WhichCOL, ThisFAM)%>% sort
  G<-assi[WhichCOLFAM,]
  assiCOLFAM<-assi[WhichCOLFAM, ]
  Test<-rownames(assiCOLFAM)[grep(FAM, rownames(assiCOLFAM))  %>% unlist]
  Candidates<- setdiff(rownames(G), Test)
  phenos_sub<-phenos[WhichCOLFAM,]
  phenos_sub$id<-rownames(phenos_sub)
  to_test<- seq(10,(length(WhichCOL)),increment) ## TRS sizes to test increment of different sizes
  ListTrain<-lapply(to_test, function(ntrain){ 
    print(ntrain)
    repeatgenalg(5,3, ntrain, criteria, as.matrix(G),Test, Candidates)
  } )
  ##predictions by optimized sample
  
  deptestopt<-phenos_sub[phenos_sub$id%in%Test,]
  for (trait in traits) {
    res_cor[[trait]][[FAM]]<-matrix(NA, ncol=5, nrow=length(to_test), dimnames=list(to_test, c("size", "corr", "mean_rel", "max_rel", "min_rel")))
    lapply(1:length(to_test), function(x){
      deptrainopt<-phenos_sub[phenos_sub$id%in%ListTrain[[x]][[1]],] ## choose solutions with rank 1
      deptrainopt$id<-factor(deptrainopt$id, levels=c(rownames(deptrainopt), rownames(deptestopt)))
      deptestopt$id<-factor(deptestopt$id, levels=c(rownames(deptrainopt), rownames(deptestopt)))
      Ztrain<-model.matrix(~-1+deptrainopt$id)
      Ztest<-model.matrix(~-1+deptestopt$id)
      K=A[c(rownames(deptrainopt), rownames(deptestopt)),c(rownames(deptrainopt), rownames(deptestopt))]
      K=as.matrix(K)
      res <- mixed.solve(y=deptrainopt[,trait],Z=Ztrain,X=matrix(1, nrow=nrow(deptrainopt), ncol=1),K=K  )
      predictopt<- Ztest %*% as.matrix(res$u)
      CORR<-cor(predictopt, deptestopt[,trait])
      mean_rel<- apply(K[Test,ListTrain[[x]][[1]]],2,mean) %>% mean
      max_rel<-apply(K[Test,ListTrain[[x]][[1]]],2,mean) %>% max
      min_rel<-apply(K[Test,ListTrain[[x]][[1]]],2,mean) %>% min
      res_cor[[trait]][[FAM]][x,]<<-c(to_test[x], CORR, mean_rel, max_rel, min_rel )
    })

  }
  saveRDS(res_cor, file=paste0(odir, "/predictions/STPGA/init_STPGA_oneTrs_",paste(FAM, criteria,increment, sep="_" ), ".rds"))
}
res<-lapply(families, function(x) pred_fam_trait2(x,"CDMEAN", 20) )


#################################
## reformat data for 4 traits ###
#################################

criteria="CDMEAN"
increment=20

res<-data.frame(size=NA, corr=NA,mean_rel=NA, max_rel=NA, min_rel=NA, FAM=NA, trait=NA)
for (trait in c("N_Peak_Force_BLUP", "Acoustic_Linear_Distance_BLUP", "PC1", "PC2")){
  assign(trait, res)
  for (FAM in families) {
    data= readRDS(paste0(odir, "/predictions/STPGA/init_STPGA_oneTrs_",paste(FAM, criteria,increment, sep="_" ), ".rds")) %>% as.data.frame()
    selec<-grep(trait, colnames(data))
    data<-data[,selec]
    data$FAM<-rep(FAM, nrow(data))
    data$trait<- trait
    colnames(data)<-c("size", "corr","mean_rel", "max_rel", "min_rel", "FAM", "trait")
    assign(trait, rbind(get(trait), data))
  }
  res<-rbind(get(trait), res)
}
res<-res[-which(is.na(res$size)),]
res<-unique(res)
res$FAM<-as.factor(res$FAM)
res$trait<-as.factor(res$trait)
res$trait<-revalue(res$trait, c("Acoustic_Linear_Distance_BLUP"="ALD","N_Peak_Force_BLUP" ="FNP"))
write.table(res,paste0(odir, "/predictions/STPGA/res_4traits_CDMEAN_20.txt"))

cat("Prediction optimizations done with STPGA method!\n")
