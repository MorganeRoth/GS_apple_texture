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
traits=colnames(phenos)

pred_fam_trait <-function(FAM, trait, criteria, increment) {
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
    # modelopt<-emmreml(y=deptrainopt[,trait],X=matrix(1, nrow=nrow(deptrainopt), ncol=1), 
    # Z=Ztrain, K=K)
    # predictopt<-Ztest%*%modelopt$uhat
    CORR<-cor(predictopt, deptestopt[,trait])
    mean_rel<- apply(K[Test,ListTrain[[x]][[1]]],2,mean) %>% mean
    max_rel<-apply(K[Test,ListTrain[[x]][[1]]],2,mean) %>% max
    min_rel<-apply(K[Test,ListTrain[[x]][[1]]],2,mean) %>% min
    res_cor[[trait]][[FAM]][x,]<<-c(to_test[x], CORR, mean_rel, max_rel, min_rel )
    
    saveRDS(res_cor, file=paste0(odir, "/predictions/STPGA/init_STPGA_",paste(FAM, trait, criteria,increment, sep="_" ), ".rds"))
  })
}


# pred_fam_trait("GaPL", "PC1", "PEVMEAN")
myparam<-matrix(c(lapply(families, function(x) rep(x, length(traits)))%>% unlist,
                  rep(traits, length(families))), nrow=length(traits)*length(families), ncol=2,
                dimnames=list(c(1:(length(traits)*length(families))), c("fam", "trait")))
# myparam<-matrix(c(lapply(families[-c(1:3)], function(x) rep(x, length(traits)))%>% unlist,
#                   rep(traits, length(families[-c(1:3)]))), nrow=length(traits)*length(families[-c(1:3)]), ncol=2,
#                 dimnames=list(c(1:(length(traits)*length(families[-c(1:3)]))), c("fam", "trait")))

res<-mapply(function(x,y) pred_fam_trait(x, y, "CDMEAN", 20), myparam[,"fam"],myparam[,"trait"] )

FAM=families[1]
trait=traits[1]
criteria="CDMEAN"
increment=20
mytrait=read.table(paste0(odir, "/predictions/STPGA/init_STPGA_",paste(FAM, trait, criteria,increment, sep="_" ), ".txt"), h=T)

##############################
##### try another method #####
##### one TRS for all traits #
##############################

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
      # modelopt<-emmreml(y=deptrainopt[,trait],X=matrix(1, nrow=nrow(deptrainopt), ncol=1), 
      # Z=Ztrain, K=K)
      # predictopt<-Ztest%*%modelopt$uhat
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


#############################
## save data for 4 traits ###
## create plots too #########
#############################

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
dim(res)
res<-unique(res)
res$FAM<-as.factor(res$FAM)
res$trait<-as.factor(res$trait)
summary(res)
res$trait<-revalue(res$trait, c("Acoustic_Linear_Distance_BLUP"="ALD","N_Peak_Force_BLUP" ="FNP"))
write.table(res,paste0(odir, "/predictions/STPGA/res_4traits_CDMEAN_20.txt"))


# gg<-ggplot(data=res, aes( x = size, y=corr)) +
#   geom_line(aes(y = corr, color=trait), size=0.7, linetype="dashed") +
#   geom_point(aes(y = corr, color=trait), size=2) +
#   facet_wrap(~FAM) +
#   scale_y_continuous( limits=c(-0.4,1),breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1))+
#   labs(color = "Trait accuracy", y="Accuracy", x="TRS size") +
#   theme( axis.title.y.right = element_text( angle = 90),
#          legend.position = "left")
pdf(file=paste0(odir, "/predictions/STPGA/CDMEAN_opt_4traits.pdf"), height=6, width=10)     
print(gg)
dev.off()


#############################
## save data for all traits #
#############################
traits=colnames(phenos)
criteria="CDMEAN"
increment=20
res<-data.frame(size=NA, corr=NA, FAM=NA, trait=NA, criteria=NA)
for (trait in traits){
  assign(trait, res)
  for (FAM in families) {
    data= readRDS(paste0(odir, "/predictions/STPGA/init_STPGA_",paste(FAM, trait, criteria,increment, sep="_" ), ".rds")) %>% as.data.frame()
    data$FAM<-rep(FAM, nrow(data))
    data$trait<- trait
    data$criteria<-criteria
    colnames(data)<-c("size", "corr", "FAM", "trait", "criteria")
    data=unique(data)
    assign(trait, rbind(get(trait), data))
  }
  res<-rbind(get(trait), res)
}
res<-res[-1,]
write.table(res, file="all_res_CDMean.txt", sep="\t", quote=F)
# ################################################################################
# ### STEP 1 simplified to avoid crashes: no initialization ######################
# ################################################################################
# 
# repeatgenalg0<-function(ntrain, criteria, G,Test, Candidates){
#   ListTrain<-GenAlgForSubsetSelection(P=G,Candidates=Candidates, 
#                                       Test=Test,ntoselect=ntrain, 
#                                       errorstat=criteria, mc.cores=1)
#   return(ListTrain)
# }
# 
# # for trials
# # FAM=families[1]
# # trait=traits[1]
# # criteria="CDMEAN"
# 
# ## very important: P parameter must be a matrix !!!!!!!!!!!!!!
# 
# pred_fam_trait0 <-function(FAM, trait, criteria) {
#   ThisFAM<-ids[ grep(FAM, ids)]
#   WhichCOLFAM<-c(WhichCOL, ThisFAM)%>% sort
#   # A_COLFAM<-assi[WhichCOLFAM, ] %>% as.matrix
#   # G<-A[WhichCOLFAM,] %>% as.matrix
#   G<-assi[WhichCOLFAM,]
#   Test<-rownames(G)[grep(FAM, rownames(G))  %>% unlist]
#   Candidates<- setdiff(rownames(G), Test)
#   phenos_sub<-phenos[WhichCOLFAM,]
#   phenos_sub$id<-rownames(phenos_sub)
#   to_test<- seq(10,(length(WhichCOL)),10) ## TRS sizes to test increment of 10
#   ListTrain<-lapply(to_test, function(ntrain){ 
#     print(ntrain)
#     repeatgenalg0(ntrain, criteria,as.matrix(G), Test, Candidates)
#   } )
#   ##predictions by optimized sample
#   
#   deptestopt<-phenos_sub[phenos_sub$id%in%Test,]
#   res<-matrix(NA, ncol=2, nrow=length(to_test), dimnames=list(to_test, c("size", "corr")))
#   count=0
#   for (x in to_test){
#     count=count+1
#     deptrainopt<-phenos_sub[phenos_sub$id%in%ListTrain[[count]][[1]],]
#     deptrainopt$id<-factor(deptrainopt$id, levels=c(rownames(deptrainopt), rownames(deptestopt)))
#     deptestopt$id<-factor(deptestopt$id, levels=c(rownames(deptrainopt), rownames(deptestopt)))
#     Ztrain<-model.matrix(~-1+deptrainopt$id)
#     Ztest<-model.matrix(~-1+deptestopt$id)
#     K=A[c(rownames(deptrainopt), rownames(deptestopt)),c(rownames(deptrainopt), rownames(deptestopt))] %>% as.matrix()
#     modelopt<-emmreml(y=deptrainopt[,trait],X=matrix(1, nrow=nrow(deptrainopt), ncol=1), 
#                       Z=Ztrain, K=as.matrix(K))
#     predictopt<-Ztest%*%modelopt$uhat
#     CORR<-cor(predictopt, deptestopt[,trait])
#     res[count,]<-c(x, as.numeric(CORR) )
#   }
#   return(res)
# }
# 
# res_cor=list()
# for (i in traits) {
#   for (j in families) {
#     # xx<-pred_fam_trait0(i, j)
#     res_cor[[i]][[j]]<-pred_fam_trait0(j, i, "PEVMEAN")
#     
#   }
# }
# saveRDS(res_cor, file=paste0(odir,  "/predictions/STPGA/opt_PEVMEAN.rds"))
# 
# res_cor=list()
# for (i in traits) {
#   for (j in families) {
#     # xx<-pred_fam_trait0(i, j)
#     res_cor[[i]][[j]]<-pred_fam_trait0(j, i, "CDMEAN")
#     
#   }
# }
# saveRDS(res_cor, file=paste0(odir,  "/predictions/STPGA/opt_CDMEAN.rds"))
# 
# 
# par(mfrow=c(2:3), mar=c(rep(5,4)))
# 
# myplot <- function(trait){
#   lapply(families, function(x) {
#     data=res_cor[[trait]][[x]]
#     print(summary(data))
#     plot(data=data, corr ~ size,
#          xlab="size_TRS", ylab="", main=paste(trait, x ), ylim=c(-0.5,1), cex=0.2) %>% print
#     axis(side = 2) %>% print
#     mtext("correlation", side = 2, line = 3, cex=0.75) %>% print
#     # par(new = TRUE)
#     # plot(data=data, mean_rel %>% as.character %>% as.numeric ~ size_TRS  %>% as.character %>% as.numeric,
#     #      type = "l", xaxt = "n", yaxt = "n",
#     #      xlab="", ylab="", main=paste(trait, x ), col="red", ylim=c(-0.09,0.25))  %>% print
#     # axis(side = 4) %>% print
#     # mtext("mean_relatedness", side = 4, line = 3, cex=0.75, col="red") %>% print
#     # abline(a=0, b=0)
#   })
# }
# 
# res_cor<-readRDS(paste0(odir,  "/predictions/STPGA/opt_PEVMEAN.rds"))
# pdf("trial.pdf", height=6, width = 8)
# par(mfrow=c(2,3))
# myplot("PC1") %>% print
# dev.off()
# 
# res_format<-reshape2::melt(res_cor)
# head(res_format)
# tail(res_format)
# summary(res_format)
# res_format<-dcast(res_format,L1+L2+Var1~Var2)
# res_format<-res_format[,-3]
# colnames(res_format)<-c("trait", "fam", "size", "corr")
# head(res_format)
# res_CDmean<-res_format
# write.table(res_CDmean, file=paste0(odir, "/predictions/STPGA/all_res_optimisation_CDmean.txt"), sep="\t", quote=F)
# 
# ## compare to opt with relatedness
# sel_traits<- traits[c(3,11,13,14)]
# par(mfrow=c(2,3))
# res_rel<-readRDS("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_optim_kinship/A.mat_TRS_opt_all_traits.rds")
# lapply(sel_traits,function(trait) {
#   # pdf(file=paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/", trait, "_TRS_opt_size.pdf"), height=5, width=8)
#   data=res_rel[[trait]] %>% as.data.frame()
#   print(summary(data))
#   data$mean_rel<-data$mean_rel %>% as.character %>% as.numeric
#   data$accuracy<-data$accuracy %>% as.character %>% as.numeric
#   data$size_TRS<-data$size_TRS %>% as.character %>% as.numeric
#   res_opt<-res_format[which(res_format$trait==trait),] %>% droplevels
#   res_opt$model<-"CDmean"
#   data$model<-"Rel"
#   data<-data[,-3]
#   res_opt<-res_opt[,c(3,2,4,5)]
#   colnames(res_opt)<-colnames(data)
#   data2<-rbind(res_opt, data)
#   gg<-ggplot(data=data2, aes(y= accuracy, x= size_TRS, col=model)) + 
#     geom_line(aes(linetype=model))+
#     facet_wrap(~family, ncol=3)+
#     geom_point()+
#     scale_color_manual(values=c("#999999", "#E69F00"))+
#     theme_bw()
#   print(gg)                  
# })
# ## End(Not run)
# 
# ## plot 4 traits at the same time
# res_cor<-readRDS(paste0(odir,  "/predictions/STPGA/opt_PEVMEAN.rds"))
# res2<-res_cor[[sel_traits[1]]] %>% as.data.frame()
# res2$trait=sel_traits[1]
# 
# for (trait in sel_traits[-1] ) {
#   data<-res_cor[[trait]] %>% as.data.frame()
#   print(head(data))
#   data$trait<-trait
#   res2<-rbind(res2,data)
# }
# 
# summary(res2)
# head(res2)
# res2<-res2[,-c(3,5,7,9,11)]
# colnames(res2)<-c("size_TRS", families, "trait")
# res2$trait<-revalue(res2$trait, c("Acoustic_Mean_Pressure_BLUP"="APMax","N_Peak_Force_BLUP" ="FNP"))
# res2$size_TRS<-res2$size_TRS %>% as.character %>% as.numeric
# # res2$mean_rel<-res2$mean_rel%>% as.character %>% as.numeric
# res2<-melt(res2, measure.vars=families)
# summary(res2)
# res2$trait<-as.factor(res2$trait)
# colnames(res2)[c(3,4)]<-c("family", "accuracy")
# gg<-ggplot(data=res2, aes( x = size_TRS)) +
#   geom_line(aes(y = accuracy, color=trait), size=0.7, linetype="dashed") +
#   geom_point(aes(y = accuracy, color=trait), size=2) +
#   facet_wrap(~family) +
#   # geom_line(aes(y = mean_rel*8+0.3)) +
#   # scale_y_continuous(sec.axis = sec_axis(~(.-0.3)/8, name = "Mean relatedness"), limits=c(-0.3,1),breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1,1.2))+  
#   scale_y_continuous(limits=c(-0.45,1),breaks=c(-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2))+  
#   labs(color = "Trait accuracy", y="Accuracy", x="TRS size") +
#   theme( axis.title.y.right = element_text( angle = 90),
#          legend.position = "left")
# 
# pdf(file=paste0(odir, "/predictions/STPGA/PEVmean.pdf"), height=6, width=11)        
# print(gg)
# dev.off()


