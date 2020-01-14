cat("calculating Fst values between clusters...\n")

## load data
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
WhichCOL=c(1:nrow(genos_round))[-c(lapply(families, function(x) grep(x, rownames(genos_round)))%>%unlist)]
clusters<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/genos_modelled/assignments_COLL_DAPC.txt", h=T)

geno<-genos_ready[clusters$Name,]
dim(clusters)
dim(geno)
pop<-cbind(clusters$Cluster, geno) %>% as.data.frame()
dim(pop)
pop[1:5,1:5]
colnames(pop)[1]="POP"

res<-pairwise.WCfst(pop[,1:100], diploid=T)

data(gtrunchier)
pairwise.WCfst(gtrunchier[,-2],diploid=TRUE)
write.table(res,paste0(odir, "/genos_modelled/Fst_between_clusters.txt"), sep="\t", quote=F)
