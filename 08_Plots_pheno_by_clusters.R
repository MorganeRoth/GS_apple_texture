
#################################################
### histograms for each trait by cluster ########
#################################################

phenos<-read.table(file=paste0(odir, "/phenos_modelled/BLUPs_PC1_PC2_for_pred.txt"), h=T)
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
rownames(phenos)
clusters<-read.table(paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"), h=T)
cluster_fams<-read.table(paste0(odir, "/genos_modelled/assignements_families.txt"), h=T)
colnames(clusters)<-colnames(cluster_fams)
clusters<-rbind(clusters, cluster_fams)
rownames(clusters)<-clusters$Name
clusters<-clusters[rownames(phenos),]
WhichCOL<-c(1:nrow(phenos))[-c(lapply(families, function(x) grep(x,rownames(phenos)) ) %>% unlist)]
traits=colnames(phenos)
colnames(phenos)<-c("ALD", "APMax", "APMean", "ANP", "Area", "FF","FLD" ,"FI","FMax", "FMean","FNP", "YM", "PC1", "PC2")
phenos<-cbind(phenos, clusters)
head(phenos)

phenos_m<-melt(phenos, id=c("Name","cluster"))
head(phenos_m)
mu <- ddply(phenos_m, c("variable", "cluster"), summarise, grp.mean=mean(value))
phenos_m<-merge(phenos_m, mu, by=c("cluster", "variable"),all.x=T)
pdf(file=paste0(odir, "/phenos_modelled/phenos_by_cluster_histograms.pdf"), width=10, height=6)
ggplot(data=phenos_m, aes(x=value, color=cluster%>% as.factor)) + 
  geom_histogram()+
  # geom_density()+
  geom_vline( aes(xintercept=grp.mean, color=cluster %>% as.factor),linetype="dashed") +
  facet_grid( phenos_m$cluster ~ phenos_m$variable , scales="free")+
  labs(color = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

pdf(file=paste0(odir, "/phenos_modelled/phenos_by_cluster_boxplots.pdf"), width=8, height=6)
ggplot(data=phenos_m, aes(x=cluster %>% as.factor, y=value)) + 
  geom_boxplot()+
  facet_wrap( phenos_m$variable , scales="free")+
  labs(x="Cluster", y="BLUP value")
dev.off()


# # pops<- c(list(WhichCOL), lapply(families, function(x) grep(x,rownames(phenos)) ) )
# # names(pops)=c("COL", families)
# pdf(file=paste0(odir, "/phenos_modelled/hist_all_traits_supdata.pdf"), height=20, width=30)
# par(mfrow=c(7,14), mar=rep(1.5,4))
# for (pop in names(pops)) {
#   print(pop)
#   for (trait in traits) {
#     hist(phenos[pops[[pop]],trait], main=pop, xlab=trait)
#   }
# }
# dev.off()
# 
# myvec=rep(NA, nrow(phenos))
# for (pop in names(pops)) {
#   myvec[pops[[pop]]]=pop
# }
# phenos$pop=myvec
# phenos$id<-rownames(phenos)
# phenos_melt<- melt(phenos, id=c("id", "pop"))
# 
# phenos_melt$variable<-revalue(phenos_melt$variable, c("Acoustic_Linear_Distance_BLUP"="ALD",
#                                                       "Acoustic_Max_Pressure_BLUP"="APMax",
#                                                       "Acoustic_Mean_Pressure_BLUP" ="APMean",
#                                                       "Acoustic_Npeak_BLUP"="ANP",
#                                                       "Area_BLUP"="Area",
#                                                       "Final_Force_BLUP" ="FF",
#                                                       "Force_Linear_Distance_BLUP" ="FLD" ,
#                                                       "Initial_Force_BLUP" = "FI",
#                                                       "Max_Force_BLUP" = "FMax",
#                                                       "Mean_Force_BLUP"= "FMean", 
#                                                       "N_Peak_Force_BLUP" ="FNP",
#                                                       "Young_Module_BLUP" = "YM"))
# # mu <- ddply(phenos_melt, c("variable", "pop"), summarise, grp.mean=mean(value))
# mm<-merge(phenos_melt, mu, by=c("pop", "variable"),all.x=T)
# 
# pdf(file=paste0(odir, "/phenos_modelled/all_histograms.pdf"), height=8, width=15)
# ggplot(data=mm, aes(x=value, color=pop)) + 
#   geom_histogram()+
#   # geom_density()+
#   geom_vline( aes(xintercept=grp.mean, color=pop),linetype="dashed") +
#   facet_grid( mm$pop ~ mm$variable , scales="free")+
#   labs(color = "Population") +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
# 
# dev.off()
