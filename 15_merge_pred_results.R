
list.dirs(paste0(odir, "/predictions"))

list.files(paste0(odir, "/predictions/COLLtoFAMS"))

###############################
### COLL TO COLL ##############
### with and without clusters #
###############################

colcol<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoCOLL/accuracy_5fold_rrBLUP.txt", h=T)
CCM<-aggregate(colcol$accuracy, by=list(colcol$trait), mean)
CCsd<-aggregate(colcol$accuracy, by=list(colcol$trait), sd)
CC<-merge(CCM, CCsd, by="Group.1")
CC$Model<-"NO_CLUST"
colcolcluster<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoCOLL/rrBLUPs_5fold_5clusters_fix.txt", h=T)
CCCM<-aggregate(colcolcluster$accuracy, by=list(colcolcluster$trait), mean)
CCCsd<-aggregate(colcolcluster$accuracy, by=list(colcolcluster$trait), sd)
CCC<-merge(CCCM, CCCsd, by="Group.1")
CCC$Model<-"5_CLUST"
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

pdf(file=paste0(odir, "/predictions/COLtoCOL_modelsAB.pdf"), height=5, width=8)
ggplot(CC,aes( y=accuracy, x=trait, fill=Model))+
  geom_bar(stat="identity",position=position_dodge()) +
  labs(x="Trait", title="Predictions rrBLUP, 5-fold", y="Correlation")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_errorbar(aes(ymin=accuracy-sd, ymax=accuracy+sd), width=.2,
                position=position_dodge(.9)) 
dev.off()

#############################
### COLL TO FAM scenarios ###
### NO CLUSTERING ###########
#############################

scenario1<-read.table(paste0(odir, "/predictions/COLLtoFAMS/rrBLUPS.txt"), h=T)
scenario1
scenario1$sd<-0
MYCOL<-colnames(scenario1) ## we will keep this format
scenario1$scenario<-"scenario1"

scenario2<-readRDS("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_0.3_prop100_reps.rds")
scenario2_mean<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), mean )
scenario2_sd<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), sd )
scenario2<-merge(scenario2_mean, scenario2_sd, by=c("Group.1", "Group.2"))
scenario2$scenario<-"scenario2"
colnames(scenario2)<-MYCOL

scenario3<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_add_relFAM/rrBLUPs_aggregated.txt", h=T)
colnames(scenario3)<-MYCOL
scenario3$scenario<-"scenario3"

all<- ldply(lapply(1:3, function(x) get(paste0("scenario",x))), data.frame)
summary(all)
all$accuracy<-round(all$accuracy, digits=4)
## error in some trait names

all[all$trait=="Acoustic_Lienar_Distance_BLUP","trait"]<-"Acoustic_Linear_Distance_BLUP"
all$trait<-droplevels(all$trait)
all$FAM<-revalue(all$FAM, c("FuPi"="FjPi", "GDFj"="FjGD"))
families<-sort(levels(all$FAM))
all$FAM<-factor(all$FAM, levels=families)

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
## plot barplot
# source: http://www.sthda.com/french/wiki/ggplot2-barres-d-erreur-guide-de-demarrage-rapide-logiciel-r-et-visualisation-de-donnees


png(file=paste0(odir, "/predictions/all_scenarios.png"), height=500, width=800)
ggplot(all,aes( y=accuracy, x=trait, fill=scenario))+
  geom_bar(stat="identity",position=position_dodge()) +
  facet_wrap(~FAM)+
  labs(x="Trait", title="Predictions rrBLUP", y="Correlation")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(limits = c(-0.5, 1)) +
  geom_errorbar(aes(ymin=accuracy-sd, ymax=accuracy+sd), width=.2,
                position=position_dodge(.9)) 
dev.off()


#############################
### COLL TO FAM scenarios ###
### NO CLUSTERING ###########
#############################

scenario1<-read.table(paste0(odir, "/predictions/COLLtoFAMS/rrBLUPS.txt"), h=T)
scenario1
scenario1$sd<-0
scenario1$scenario<-"Scenario 1"
MYCOL<-colnames(scenario1) ## we will keep this format

scenario2<-readRDS("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_0.3_prop100_reps.rds")
scenario2_mean<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), mean )
scenario2_sd<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), sd )
scenario2<-merge(scenario2_mean, scenario2_sd, by=c("Group.1", "Group.2"))
scenario2$scenario<-"Scenario 2"
colnames(scenario2)<-MYCOL

scenario3<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_add_relFAM/rrBLUPs_aggregated.txt", h=T)
scenario3$scenario<-"Scenario 3"
colnames(scenario3)<-MYCOL

all<- ldply(lapply(1:3, function(x) get(paste0("scenario",x))), data.frame)
summary(all)
all$accuracy<-round(all$accuracy, digits=4)
## error in some trait names

all[all$trait=="Acoustic_Lienar_Distance_BLUP","trait"]<-"Acoustic_Linear_Distance_BLUP"
all$trait<-droplevels(all$trait)
all$FAM<-revalue(all$FAM, c("FuPi"="FjPi"))
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
## plot barplot

bmp(file=paste0(odir, "/predictions/all_scenarios.bmp"), height = 600, width = 1000, res=100)
ggplot(all,aes( y=accuracy, x=trait, fill=scenario))+
  geom_bar(stat="identity",position=position_dodge()) +
  facet_wrap(~FAM)+
  labs(x="Trait", title="Predictions rrBLUP", y="Correlation")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold"),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(-0.5, 1)) +
  geom_errorbar(aes(ymin=accuracy-sd, ymax=accuracy+sd), width=.1,
                position=position_dodge(.9)) 
dev.off()

#############################
### COLL TO FAM scenarios ###
### NO CLUSTERING ###########
#############################

scenario1<-read.table(paste0(odir, "/predictions/COLLtoFAMS/rrBLUPS.txt"), h=T)
scenario1
scenario1$sd<-0
MYCOL<-colnames(scenario1) ## we will keep this format
scenario1$scenario<-"scenario1"

scenario2<-readRDS("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_0.3_prop100_reps.rds")
scenario2_mean<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), mean )
scenario2_sd<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), sd )
scenario2<-merge(scenario2_mean, scenario2_sd, by=c("Group.1", "Group.2"))
scenario2$scenario<-"scenario2"
colnames(scenario2)<-MYCOL

scenario3<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_add_relFAM/COLLtoFAMs_add_relFAM_aggregated.txt", h=T)
colnames(scenario3)<-MYCOL
scenario3$scenario<-"scenario3"

all<- ldply(lapply(1:3, function(x) get(paste0("scenario",x))), data.frame)
summary(all)
all$accuracy<-round(all$accuracy, digits=4)
## error in some trait names

all[all$trait=="Acoustic_Lienar_Distance_BLUP","trait"]<-"Acoustic_Linear_Distance_BLUP"
all$trait<-droplevels(all$trait)
all$FAM<-revalue(all$FAM, c("FuPi"="FjPi", "GDFj"="FjGD"))
families<-sort(levels(all$FAM))
all$FAM<-factor(all$FAM, levels=families)

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
## plot barplot
# source: http://www.sthda.com/french/wiki/ggplot2-barres-d-erreur-guide-de-demarrage-rapide-logiciel-r-et-visualisation-de-donnees


png(file=paste0(odir, "/predictions/all_scenarios.png"), height=500, width=800)
ggplot(all,aes( y=accuracy, x=trait, fill=scenario))+
  geom_bar(stat="identity",position=position_dodge()) +
  facet_wrap(~FAM)+
  labs(x="Trait", title="Predictions rrBLUP", y="Correlation")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(limits = c(-0.5, 1)) +
  geom_errorbar(aes(ymin=accuracy-sd, ymax=accuracy+sd), width=.2,
                position=position_dodge(.9)) 
dev.off()


#############################
### COLL TO FAM scenarios ###
### WITH CLUSTERING #########
#############################

scenario1<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMS/rrBLUPs_clusters_fix.txt", h=T)
scenario1
scenario1$sd<-0
scenario1$scenario<-"Scenario 1"
MYCOL<-colnames(scenario1) ## we will keep this format

scenario2<-readRDS("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_prop/all_traits_pred_progenies_0.3_prop100_reps_5clusters.rds")
scenario2_mean<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), mean )
scenario2_sd<-aggregate(scenario2$accuracy,by=list(scenario2$FAM, scenario2$trait), sd )
scenario2<-merge(scenario2_mean, scenario2_sd, by=c("Group.1", "Group.2"))
scenario2$scenario<-"Scenario 2"
colnames(scenario2)<-MYCOL

scenario3<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/predictions/COLLtoFAMs_add_relFAM/rrBLUPs_clusters_aggregated.txt", h=T)
scenario3$scenario<-"Scenario 3"
colnames(scenario3)<-MYCOL

all<- ldply(lapply(1:3, function(x) get(paste0("scenario",x))), data.frame)
summary(all)
all$accuracy<-round(all$accuracy, digits=4)
## error in some trait names

all[all$trait=="Acoustic_Lienar_Distance_BLUP","trait"]<-"Acoustic_Linear_Distance_BLUP"
all$trait<-droplevels(all$trait)
all$FAM<-revalue(all$FAM, c("FuPi"="FjPi"))
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
## plot barplot

bmp(file=paste0(odir, "/predictions/all_scenarios_clusters.bmp"), height = 600, width = 1000, res=100)
ggplot(all,aes( y=accuracy, x=trait, fill=scenario))+
  geom_bar(stat="identity",position=position_dodge()) +
  facet_wrap(~FAM)+
  labs(x="Trait", title="Predictions rrBLUP with 5 clusters", y="Correlation")+
  theme(axis.text.x=element_text(angle = 45,  hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold"),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(-0.5, 1)) +
  geom_errorbar(aes(ymin=accuracy-sd, ymax=accuracy+sd), width=.1,
                position=position_dodge(.9)) 
dev.off()
