
list.dirs(paste0(odir, "/predictions"))

list.files(paste0(odir, "/predictions/COLLtoFAMS"))
list.files(paste0(odir, "/predictions/COLLtoCOLL"))
###############################
### COLL TO COLL ##############
### with and without clusters #
###############################

colcol<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoCOLL/accuracy_5fold_rrBLUP.txt", h=T)
CCM<-aggregate(colcol$accuracy, by=list(colcol$trait), mean)
CCsd<-aggregate(colcol$accuracy, by=list(colcol$trait), sd)
CC<-merge(CCM, CCsd, by="Group.1")
CC$Model<-"No clusters"
colcolcluster<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoCOLL/rrBLUPs_5fold_6clusters_fix.txt", h=T)
CCCM<-aggregate(colcolcluster$accuracy, by=list(colcolcluster$trait), mean)
CCCsd<-aggregate(colcolcluster$accuracy, by=list(colcolcluster$trait), sd)
CCC<-merge(CCCM, CCCsd, by="Group.1")
CCC$Model<-"6 clusters"
colnames(CCC)<-colnames(CC)<-c("trait", "accuracy", "sd", "Model")
CC<-rbind(CC,CCC)
CC[CC$trait=="Acoustic_Lienar_Distance_BLUP","trait"]<-"Acoustic_Linear_Distance_BLUP"
CC$trait<-droplevels(CC$trait)
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
levels(CC$trait)
CC$trait<-factor(CC$trait, levels=c("ALD", "ANP", "APMax", "APMean", "Area",
                                      "FF", "FI", "FLD", "FMax", "FMean", "FNP", 
                                      "YM", "PC1", "PC2"))

pdf(file=paste0(odir, "/predictions/COLtoCOL_modelsAB_newPC12.pdf"), height=5, width=8)
ggplot(CC,aes( y=accuracy, x=trait, fill=Model))+
  geom_bar(stat="identity",position=position_dodge()) +
  scale_fill_manual(values=c("grey40", "grey70")) +
  labs(x="Trait", title="Predictions rrBLUP, 5-fold", y="Predictive ability")+
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_errorbar(aes(ymin=accuracy-sd, ymax=accuracy+sd), width=.2,
                position=position_dodge(.9)) 
dev.off()

####################################
## some summaries on COLL to COLL ##
####################################
write.table(CC, file=paste0(odir, "/predictions/COLLtoCOLL/summaries_all.txt"), quote=F, sep="\t")
CC[which(CC$accuracy==max(CC$accuracy)),]
CC[order(CC$accuracy),]

## test on ranking

test<-read.table(paste0(odir, "/predictions/COLLtoCOLL/for_rank_test_h2_cor.txt"), h=T)
test
wilcox.test(x=test$h2,y=test$cor, paired=T, conf.int=T)

#############################
### COLL TO FAM scenarios ###
### NO CLUSTERING ###########
#############################

scenario1<-read.table(paste0(odir, "/predictions/COLLtoFAMS/rrBLUPS_reps.txt"), h=T)
scenario1
scenario1$scenario<-"Scenario 1"
MYCOL<-colnames(scenario1) ## we will keep this format

scenario2<-readRDS("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_0.3_prop100_reps.rds")
scenario2_mean<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), mean )
scenario2_sd<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), sd )
scenario2_conf_inf=scenario2_mean[,"x"]- scenario2_sd[,"x"]
scenario2_conf_sup=scenario2_mean[,"x"]+ scenario2_sd[,"x"]
scenario2<-cbind(scenario2_mean, scenario2_conf_inf, scenario2_conf_sup)
scenario2$scenario<-"Scenario 2"
colnames(scenario2)<-MYCOL

scenario3<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_add_relFAM/rrBLUPs_aggregated.txt", h=T)
scenario3$conf_inf<-scenario3$Mean_acc - scenario3$SD_acc
scenario3$conf_sup<-scenario3$Mean_acc + scenario3$SD_acc
scenario3$scenario<-"Scenario 3"
scenario3<-scenario3[,-c(4)]
colnames(scenario3)<-MYCOL

all<- ldply(lapply(1:3, function(x) get(paste0("scenario",x))), data.frame)
summary(all)
all$accuracy<-round(all$accuracy, digits=4)
## error in some trait names

all[all$trait=="Acoustic_Lienar_Distance_BLUP","trait"]<-"Acoustic_Linear_Distance_BLUP"
all$trait<-droplevels(all$trait)

families<-sort(levels(all$FAM))
all$FAM<-factor(all$FAM, levels=families[c(1:3,6,4,5)])

all$trait<-revalue(all$trait, c("Acoustic_Linear_Distance_BLUP"="ALD",
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
levels(all$trait)
all$trait<-factor(all$trait, levels=c("ALD", "ANP", "APMax", "APMean", "Area",
                                      "FF", "FI", "FLD", "FMax", "FMean", "FNP", 
                                      "YM", "PC1", "PC2"))
summary(all)
write.table(all, paste0(odir,"/predictions/all_scenarios_accuracies_no_clust.txt"),sep="\t", quote=F )
## plot barplot
# source: http://www.sthda.com/french/wiki/ggplot2-barres-d-erreur-guide-de-demarrage-rapide-logiciel-r-et-visualisation-de-donnees

# bmp(file=paste0(odir, "/predictions/all_scenarios2.bmp"), height = 600, width = 1200, res=100)
# pdf(file=paste0(odir, "/predictions/all_scenarios2.pdf"), height = 6, width = 12,)

gg<-ggplot(all,aes( y=accuracy, x=trait, fill=scenario))+
  geom_bar(stat="identity",position=position_dodge()) +
  facet_wrap(~FAM)+
  labs(x="Trait", y="Accuracy")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1, size=5),
        strip.text.x = element_text(size = 10, face = "bold"),
        legend.position = "none",
        # legend.title = element_blank(),
        # legend.key.size = unit(0.5, "cm"),
        # legend.text=element_text(size = 8),
        axis.title = element_text(size=10)) + 
  # scale_fill_grey(start = 0.8, end = 0.2) +
  # scale_fill_jco()+
  scale_fill_manual(values = viridis(3)) +
  # scale_fill_brewer(palette = "PRGn") +
  # scale_fill_manual(values=wes_palette("FantasticFox1", 3, type = "discrete"))+
  # scale_fill_manual(values=wes_palette("Moonrise1", 3, type = "discrete"))+
  scale_y_continuous(limits = c(-0.55, 1)) +
  geom_errorbar(aes(ymin=conf_low, ymax=conf_high), size=0.1, width=0.4,
                position=position_dodge(.9)) 

pdf(file=paste0(odir, "/predictions/pred_fam_figure4A_paper.pdf"), 
    width= 3.54*2, height=4, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
print(gg)
dev.off()
tiff(file=paste0(odir, "/predictions/pred_fam_figure4A_paper.tiff"), 
     width= 3.54*2, height=4, fonts="Helvetica", pointsize=10, units="in", res=300 )
print(gg)
dev.off()


#############################
### COLL TO FAM scenarios ###
### WITH CLUSTERING #########
#############################

scenario1<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMS/rrBLUPs_clusters_fix_reps.txt", h=T)
scenario1
scenario1$scenario<-"Scenario 1"
MYCOL<-colnames(scenario1) ## we will keep this format

scenario2<-readRDS("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_0.3_prop100_reps_6clusters.rds")
scenario2_mean<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), mean )
scenario2_sd<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), sd )
scenario2_conf_inf=scenario2_mean[,"x"]- scenario2_sd[,"x"]
scenario2_conf_sup=scenario2_mean[,"x"]+ scenario2_sd[,"x"]
scenario2<-cbind(scenario2_mean, scenario2_conf_inf, scenario2_conf_sup)
scenario2$scenario<-"Scenario 2"
colnames(scenario2)<-MYCOL

scenario3<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_add_relFAM/rrBLUPs_clusters_aggregated.txt", h=T)
scenario3$conf_inf<-scenario3$Mean_acc - scenario3$SD_acc
scenario3$conf_sup<-scenario3$Mean_acc + scenario3$SD_acc
scenario3$scenario<-"Scenario 3"
scenario3<-scenario3[,-c(4)]
colnames(scenario3)<-MYCOL

all<- ldply(lapply(1:3, function(x) get(paste0("scenario",x))), data.frame)
summary(all)
all$accuracy<-round(all$accuracy, digits=4)

## error in some trait names

# all[all$trait=="Acoustic_Lienar_Distance_BLUP","trait"]<-"Acoustic_Linear_Distance_BLUP"
# all$trait<-droplevels(all$trait)

families<-sort(levels(all$FAM))
all$FAM<-factor(all$FAM, levels=families[c(1:3,6,4,5)])

all$trait<-revalue(all$trait, c("Acoustic_Linear_Distance_BLUP"="ALD",
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
levels(all$trait)
all$trait<-factor(all$trait, levels=c("ALD", "ANP", "APMax", "APMean", "Area",
                                      "FF", "FI", "FLD", "FMax", "FMean", "FNP", 
                                      "YM", "PC1", "PC2"))
summary(all)
write.table(all, paste0(odir,"/predictions/all_scenarios_accuracies_with_clust.txt"),sep="\t", quote=F )

## plot barplot

# bmp(file=paste0(odir, "/predictions/all_scenarios_clusters2.bmp"), height = 600, width = 1200, res=100)
# pdf(file=paste0(odir, "/predictions/all_scenarios_clusters2.pdf"), height = 6, width = 12)

gg2<-ggplot(all,aes( y=accuracy, x=trait, fill=scenario))+
  geom_bar(stat="identity",position=position_dodge()) +
  facet_wrap(~FAM)+
  labs(x="Trait",y="Accuracy")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1, size=5),
        strip.text.x = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),
        legend.text=element_text(size = 8),
        axis.title = element_text(size=10))   + 
  # scale_fill_grey(start = 0.8, end = 0.2) +
  # scale_fill_jco()+
  scale_fill_manual(values = viridis(3)) +
  # scale_fill_brewer(palette = "PRGn") +
  # scale_fill_manual(values=wes_palette("FantasticFox1", 3, type = "discrete"))+
  # scale_fill_manual(values=wes_palette("Moonrise1", 3, type = "discrete"))+
  scale_y_continuous(limits = c(-0.6, 1)) +
  geom_errorbar(aes(ymin=conf_low, ymax=conf_high), size=0.1, width=0.4,
                position=position_dodge(.9)) 

pdf(file=paste0(odir, "/predictions/pred_fam_figure4B_paper_with_legend.pdf"), 
    width= 3.54*2, height=4, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
print(gg2)
dev.off()
tiff(file=paste0(odir, "/predictions/pred_fam_figure4B_paper_no_legend.tiff"), 
     width= 3.54*2, height=4, fonts="Helvetica", pointsize=10, units="in", res=300 )
gg3<-gg2 + theme(legend.position = "none")
print(gg3)
dev.off()


#############################
### Increasing TRS size #####
### 4 methods - PC1 #########
#############################

## load all results: accuracy and size

clust_res<-read.table( paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/res_4traits_Clusters.txt"), h=T)
clust_res$method<-rep("Clusters", nrow(clust_res))
rel_res<-read.table( paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/res_4traits_Mean relatedness.txt"), h=T)
rel_res$method<-rep("Mean relationship", nrow(rel_res))
rel_max_res<-read.table( paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/res_4traits_Max relatedness.txt"), h=T)
rel_max_res$method<-rep("Max relationship", nrow(rel_max_res))
# pevmean_res<-read.table( paste0(odir, "/predictions/STPGA/res_4traits_PEVMEAN_20.txt"), h=T)
# pevmean_res$method<-rep("PEVmean-opt", nrow(pevmean_res))
# cdmean_res<-read.table( paste0(odir, "/predictions/STPGA/8234snps/res_4traits_CDMEAN_20.txt"), h=T)
cdmean_res<-read.table( paste0(odir, "/predictions/STPGA/res_4traits_CDMEAN_20.txt"), h=T)
cdmean_res<-unique(cdmean_res)
cdmean_res$method<-rep("CDmean-opt", nrow(cdmean_res))
rel_res<-rel_res[,c("size_TRS", "accuracy","family", "trait", "method")]
rel_max_res<-rel_max_res[,c("size_TRS", "accuracy","family", "trait", "method")]
clust_res<-clust_res[,c("size_TRS", "accuracy","family", "trait", "method")]
cdmean_res<-cdmean_res[,c("size", "corr","FAM", "trait", "method")]
colnames(rel_max_res)<-colnames(rel_res)<-colnames(clust_res)<-colnames(cdmean_res)
res<-do.call(rbind, list(rel_res,rel_max_res, clust_res,  cdmean_res))
write.table(res, file= paste0(odir, "/predictions/al_res_opt_4methods.txt"), sep="\t", quote=F, row.names=F )
res$method<-as.factor(res$method)
## keep only one trait

summary(res)
res$method<-factor(res$method,levels=c("Mean relationship","Max relationship","Clusters", "CDmean-opt"), order=T)
sub1<-res[which(res$method %in% c("Mean relationship","Max relationship")),] %>% droplevels
sub2<-res[-which(res$method%in% c("Mean relationship","Max relationship")),] %>% droplevels
# trait="PC2"
plottrait<-function(trait){
  sub<-res[which(res$trait==trait),]
 gg<-ggplot(data=sub, aes( x = size, y=corr, color=method, order=method)) +
    geom_point(data=sub2[which(sub2$trait==trait),], aes(y = corr), size=1) +
    geom_line(data=sub2[which(sub2$trait==trait),],aes(y = corr), size=0.7 ,linetype="dashed") +
    geom_line(data=sub1[which(sub1$trait==trait),],aes(y = corr), size=0.5) +
    facet_wrap(~FAM) +
    scale_y_continuous( limits=c(-0.5,1),breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1))+
    # scale_color_manual(values=wes_palette("Zissou1", n=5))+
    # scale_color_manual(values=c("seagreen4", "turquoise3", "gold2", "firebrick2", "darkorchid3")) +
    scale_color_jco(breaks=c("Mean relationship","Max relationship","Clusters", "CDmean-opt", "PEVmean-opt")) +
    labs(color = "Method", y="Accuracy", x="TRS size") +
   theme( axis.title.y.right = element_text( angle = 90),
          axis.title = element_text(size=10),
          strip.text.x = element_text(size = 8, face = "bold"),
           legend.position ="none")
 pdf(file=paste0(odir, "/predictions/pred_fam_figureS5_paper_", trait, ".pdf"), 
     width= 3.54, height=3, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
 print(gg)
 dev.off()
 tiff(file=paste0(odir, "/predictions/pred_fam_figureS5_paper_", trait, ".tiff"), 
      width= 3.54, height=3, fonts="Helvetica", pointsize=10, units="in", res=300 )
 print(gg)
 dev.off()
}

lapply(levels(res$trait), plottrait)

##################################
#### plots by trait for paper ####
##################################

### first the 3 methods on relatedness

traits=levels(res$trait)
model="Mean relatedness"
name_res="A.mat_TRS_opt_all_traits.rds"
name_sum= "mean_rel_summaries.txt"
name_graph="mean_rel4traits.pdf"
plot_sum_res<-function(model, name_res,name_graph){
  results<-readRDS (paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/", name_res))
  traits<-names(results)
  ## keep only 4 traits and plot results
  res2<-results[[1]] %>% as.data.frame()
  res2$trait<-traits[1]
  for (i in c(11,13,14)) {
    data<-results[[i]] %>% as.data.frame()
    data$trait<-traits[i]
    res2<-rbind(res2,data)
  }
  
  res2$trait<-revalue(res2$trait, c("Acoustic_Linear_Distance_BLUP"="ALD","N_Peak_Force_BLUP" ="FNP"))
  res2$size_TRS<-res2$size_TRS %>% as.character %>% as.numeric
  res2$mean_rel<-res2$mean_rel%>% as.character %>% as.numeric
  res2$accuracy<-res2$accuracy %>% as.character %>% as.numeric
  res2$trait<-as.factor(res2$trait)
  ## ggplot help for 2 axes with different scales here: https://rpubs.com/MarkusLoew/226759
  gg<-ggplot(data=res2, aes( x = size_TRS)) +
    geom_line(aes(y = accuracy, color=trait), size=0.5) +
    facet_wrap(~family) +
    geom_line(aes(y = mean_rel*3), size=0.3) +
    scale_y_continuous(sec.axis = sec_axis(~./3, name = model), limits=c(-0.4,1), breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1))+  
    labs(color = "Trait accuracy", y="Accuracy", x="TRS size") +
    # scale_color_manual(values=viridis(4))+
    scale_color_jco()+
    # scale_color_jama()+
    # scale_color_nejm()+
    theme( axis.title.y.right = element_text( angle = 90),
           axis.text.x=element_text(angle = 45,  hjust = 1, size=5),
           axis.title = element_text(size=8),
           axis.text = element_text(size=7),
           strip.text.x = element_text(size = 8, face = "bold"),
           legend.position = "none")
  pdf(file=paste0(odir, "/predictions/pred_fam_figure5_paper_", model, ".pdf"), 
      width= 3.54, height=2.8, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
  print(gg)
  dev.off()
  tiff(file=paste0(odir, "/predictions/pred_fam_figure5_paper_", model, ".tiff"), 
       width= 3.54, height=2.8, fonts="Helvetica", pointsize=10, units="in", res=300 )
  print(gg)
  dev.off()

}


plot_sum_res("Mean relatedness", "A.mat_TRS_opt_all_traits.rds", "mean_rel4traits.pdf")
plot_sum_res("Clusters relatedness", "A.mat_TRS_opt_CLUSTERS_all_traits.rds", "clusters_rel4traits.pdf")

plot_sum_res_max<-function(model, name_res, name_graph){
  results<-readRDS (paste0(odir, "/predictions/COLLtoFAMs_optim_kinship/", name_res))
  traits<-names(results)
  res2<-results[[1]] %>% as.data.frame()
  res2$trait<-traits[1]
  for (i in c(11,13,14)) {
    data<-results[[i]] %>% as.data.frame()
    data$trait<-traits[i]
    res2<-rbind(res2,data)
  }
  res2$trait<-revalue(res2$trait, c("Acoustic_Linear_Distance_BLUP"="ALD","N_Peak_Force_BLUP" ="FNP"))
  res2$size_TRS<-res2$size_TRS %>% as.character %>% as.numeric
  res2$rel_added_ID_max<-res2$rel_added_ID_max%>% as.character %>% as.numeric
  res2$accuracy<-res2$accuracy %>% as.character %>% as.numeric
  res2$trait<-as.factor(res2$trait)
  ## ggplot help for 2 axes with different scales here: https://rpubs.com/MarkusLoew/226759
  gg<-ggplot(data=res2, aes( x = size_TRS)) +
    geom_line(aes(y = accuracy, color=trait), size=0.5) +
    facet_wrap(~family) +
    geom_line(aes(y =rel_added_ID_max*3), size=0.3) +
    scale_y_continuous(sec.axis = sec_axis(~./3, name = model), limits=c(-0.4,1), breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1))+  
    labs(color = "Trait accuracy", y="Accuracy", x="TRS size") +
    # scale_color_manual(values=viridis(4))+
    scale_color_jco()+
    # scale_color_jama()+
    # scale_color_nejm()+
    theme( axis.title.y.right = element_text( angle = 90),
           axis.text.x=element_text(angle = 45,  hjust = 1, size=5),
         axis.title = element_text(size=8),
         axis.text = element_text(size=7),
         strip.text.x = element_text(size = 8, face = "bold"),
         legend.position = "none")
  pdf(file=paste0(odir, "/predictions/pred_fam_figure5_paper_", model, ".pdf"), 
      width= 3.54, height=2.8, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
  print(gg)
  dev.off()
  tiff(file=paste0(odir, "/predictions/pred_fam_figure5_paper_", model, ".tiff"), 
       width= 3.54, height=2.8, fonts="Helvetica", pointsize=10, units="in", res=300 )
  print(gg)
  dev.off()
}  

plot_sum_res_max("Max relatedness", "A.mat_TRS_opt_all_traits_MAX.rds", "max_rel4traits.pdf")

###### now for the STPGA results

res<-read.table(paste0(odir, "/predictions/STPGA/res_4traits_CDMEAN_20.txt"),h=T)


gg<-ggplot(data=res, aes( x = size)) +
  geom_line(aes(y = corr, color=trait), size=0.5) +
  facet_wrap(~FAM) +
  geom_line(aes(y = mean_rel*3), size=0.3) +
  scale_y_continuous(sec.axis = sec_axis(~./3, name = "Mean relatedness"), limits=c(-0.4,1), breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1))+  
  labs(color = "Trait accuracy", y="Accuracy", x="TRS size") +
  # scale_color_manual(values=viridis(4))+
  scale_color_jco()+
  # scale_color_jama()+
  # scale_color_nejm()+
  theme( axis.title.y.right = element_text( angle = 90),
         axis.text.x=element_text(angle = 45,  hjust = 1, size=5),
         axis.title = element_text(size=8),
         axis.text = element_text(size=7),
         strip.text.x = element_text(size = 8, face = "bold"),
         legend.position = "none")
pdf(file=paste0(odir, "/predictions/pred_fam_figure5D_paper.pdf"), 
    width= 3.54, height=2.8, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
print(gg)
dev.off()
tiff(file=paste0(odir, "/predictions/pred_fam_figure5D_paper.tiff"), 
     width= 3.54, height=2.8, fonts="Helvetica", pointsize=10, units="in", res=300 )
print(gg)
dev.off()
# gg<-ggplot(data=res[which(res$trait==trait),], aes( x = size, y=corr)) +
#   geom_line(data=res[which(res$trait==trait),],aes(y = corr, color=method), size=0.7) +
#   facet_wrap(~FAM) +
#   scale_y_continuous( limits=c(-0.5,1),breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1))+
#   labs(color = "Method", y="Accuracy", x="TRS size") +
#   scale_color_startrek()+
#   # scale_color_manual(values=viridis(5))+
#   theme( axis.title.y.right = element_text( angle = 90),
#          legend.position = "right")



### some stats on these results
methods=c("Mean relationship","Max relationship","Clusters", "CDmean-opt", "PEVmean-opt")
traits=res$trait %>% levels
sum50<-data.frame(trait=NA, fam=NA, method=NA, ratio=NA,acc50=NA, max=NA)

for (trait in traits) {
  for (fam in families) {
    for (method in methods[-2]) {
      data<-res[which(res$trait==trait & res$method==method & res$FAM==fam),] %>% droplevels
      maxacc<-data[,"corr"] %>% as.numeric %>% max 
      # print(maxacc)
      acc50<-data[which(data$size==50), "corr"]%>% unique
      myratio<-acc50/maxacc %>% as.numeric
      print(c(trait, fam,method,  myratio,acc50, maxacc))
      sum50<-rbind(sum50, c(trait, fam,method,  myratio,acc50, maxacc))
      rm(maxacc, acc50, myratio)
    }
  }
}
head(sum50)
summary(sum50)
sum50<-sum50[-1,]
sum50[,c(4:6)] <- apply(sum50[,c(4:6)], 2, function(x)as.numeric(as.character(x)))

write.table(sum50, paste0(odir,"/predictions/comparisons_50IDs_max_4traits.txt") , quote=F, row.names=F, sep="\t" )

methods=c("Individuals","Clusters", "CDmean-opt", "PEVmean-opt")
traits=res$trait %>% levels
maxacc_minsize<-data.frame(trait=NA, fam=NA, method=NA, maxacc_minsize=NA, max=NA)
for (trait in traits) {
  for (fam in families) {
    for (method in methods) {
      data<-res[which(res$trait==trait & res$method==method & res$FAM==fam),] %>% droplevels
      print(head(data))
      maxacc<-data[,"corr"] %>% unique %>% max
      sizemax<-data[which(data$corr==maxacc), "size"] %>% unique %>% min
      maxacc_minsize<-rbind(maxacc_minsize, c(trait, fam,method,  sizemax, maxacc))
    }
  }
}
maxacc_minsize<-maxacc_minsize[-1,]
head(maxacc_minsize)
maxacc_minsize[,c(4:5)] <- apply(maxacc_minsize[,c(4:5)], 2, function(x)as.numeric(as.character(x)))
summary(maxacc_minsize)
write.table(maxacc_minsize, paste0(odir,"/predictions/minsize_maxacc_4traits.txt") , quote=F, row.names=F, sep="\t" )

ind<-maxacc_minsize[which(maxacc_minsize$method=="Individuals" & maxacc_minsize$max>=0.2),] %>% droplevels
hist(ind$maxacc_minsize, nclass=5)
ggplot(data=ind, aes(x=maxacc_minsize))+ 
  geom_histogram(binwidth = 50, color="black", fill="white", boundary=0) +
  scale_x_discrete(name="Minimal size for maximum accuracy",limits=c(10,50,100,150,200,242))

################################################################
####### Find the best TRS in each family for each trait ########
################################################################

res<-read.table(file= paste0(odir, "/predictions/al_res_opt_4methods.txt"), h=T, sep="\t")
res$combi<-paste(res$FAM, res$trait)
selec<-c()
for (i in unique(res$combi)) {
  data<-res[which(res$combi==i),]
  sel<-data[which(data$corr==max(data$corr)), c("size", "combi", "corr", "method")]
  selec<-rbind(selec, sel)
}
selec %>% dim
selec %>% unique 

write.table(selec %>% unique, file=paste0(odir, "/predictions/best_TRS_model_optimization.txt"), sep="\t", quote=F)

###############################
### FULL TO FULL ##############
### cross-valid in whole pop ##
###############################

full<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/FULLtoFULL/accuracy_5fold_ALLtoALL_pearson_spearman.txt", h=T)
fullM<-aggregate(full$accuracy_Pearson, by=list(full$trait), mean)
fullsd<-aggregate(full$accuracy_Pearson, by=list(full$trait), sd)
FF<-merge(fullM, fullsd, by="Group.1")
colnames(FF)<-c("trait", "accuracy", "sd")

FF$trait<-revalue(FF$trait, c("Acoustic_Linear_Distance_BLUP"="ALD",
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
FF$trait<-factor(FF$trait, levels=c("ALD", "ANP", "APMax", "APMean", "Area",
                                    "FF", "FI", "FLD", "FMax", "FMean", "FNP", 
                                    "YM", "PC1", "PC2"))

pdf(file=paste0(odir, "/predictions/COLtoCOL_modelsAB_newPC12.pdf"), height=5, width=8)
ggplot(FF,aes( y=accuracy, x=trait))+
  geom_bar(stat="identity",position=position_dodge()) +
  # scale_fill_manual(values=c("grey40", "grey70")) +
  labs(x="Trait", title="Predictions rrBLUP, 5-fold", y="Predictive ability")+
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_errorbar(aes(ymin=accuracy-sd, ymax=accuracy+sd), width=.2,
                position=position_dodge(.9)) 
dev.off()


################################
### FULL TO FAM ################
### each fam with rest of pop ##
################################

ffam<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/FULLtoFULL/accuracy_5fold_ALLtoFAM_pearson_spearman.txt", h=T)
head(ffam)
ggplot(ffam,aes( y=accuracy, x=trait))+
  facet_wrap(~FAM)+
  geom_bar(stat="identity",position=position_dodge()) +
  # scale_fill_manual(values=c("grey40", "grey70")) +
  labs(x="Trait", title="Predictions rrBLUP, 5-fold", y="Predictive ability")+
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_errorbar(aes(ymin=conf_low, ymax=conf_high), width=.2,
                position=position_dodge(.9)) 
