cat("Analyzing genos: formatting, metrics and DAPC analysis...\n")

### PLOTS TURNED OFF HERE - YOU CAN REACTIVATE

dir.create(paste0(odir, "/genos_modelled"), showWarnings = FALSE, recursive = TRUE)

#################################################################################
######################## Formatting and producing metrics #######################
#################################################################################
parents<- read.table(paste0(idir, "/families_parents.txt"), sep="\t", h=T)
families<-c("FjDe", "GDFj", "FjPi", "FjPL", "GaPi", "GaPL")

# First we need to round imputed values
genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))
genos_round<-apply(genos_ready,2, function(x) as.numeric(x) %>% round(., digits=0))

rownames(genos_round)<-rownames(genos_ready)
genos_round[genos_round==2]<-1
apply(genos_round[,1:20], 2, function(x) as.factor(x) %>%summary) ## now entire values only
cat("sanity check\n")
genos_round[genos_round==2]<-1
genos_round[genos_round==-2]<- -1
genos_ready[1:6,1:6] ## sanity check

## collection identifiers
WhichCOL=c(1:nrow(genos_round))[-c(lapply(families, function(x) grep(x, rownames(genos_round)))%>%unlist)]

## realized additive relationship matrix (rrBLUP) from imputed data

A <- A.mat(genos_round)
colnames(A)<-rownames(A)<-rownames(genos_ready)
cat("plotting heatmap of additive relationship...\n")
heatmap.2(A, trace="none",key.xlab="Relationship", key.ylab="", key.title = "", 
               keysize = 2, col=inferno(100) )


# pdf(file=paste0(odir, "/genos_modelled/Ka_Amat.pdf"), height=6, width=6)
# 
# pdf(file=paste0(odir, "/genos_modelled/Heatmap_Amat2.pdf"), 
#     width= 3.54*1.5, height=3.54*1.5, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
# tiff(file=paste0(odir, "/genos_modelled/Heatmap_Amat.tiff"), 
#      width= 3.54*1.5, height=3.54*1.5, fonts="Helvetica", pointsize=10, units="in", res=300 )
# heatmap.2(A, trace="none",key.xlab="Relationship", key.ylab="", key.title = "",
#           keysize = 2, col=inferno(100), labRow = FALSE, labCol = FALSE ,margins = c(2,2) ) 
# dev.off()

###### SAVE DATA- DEACTIVATED ######
# write.table(A, file=paste0(odir, "/genos_modelled/Ka_Amat.txt"), row.names=T, quote=F, sep="\t")
# write.table(colnames(A)[myh$colInd], file=paste0(odir, "/genos_modelled/Clustering_Amat.txt"), row.names=T, quote=F, sep="\t")


## Mean relatedness (additive) between families and collection
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
cat("Number of individuals in the collection:\n")
cat(paste(length(WhichCOL))) ## 259 ids
rel<-matrix(NA, nrow=length(families), ncol=1, dimnames=list(families, "mean_rel_collection"))
for (i in families) {
  whichFAM<-c(1:nrow(A))[grep(i, rownames(A)) %>% unlist]
  rel[i,1]<-mean(A[whichFAM,WhichCOL])
}
# save
# write.table(rel,paste0(odir, "/genos_modelled/mean_relatedness_families_coll_Amat.txt"), quote=F, row.names=T, sep="\t")

### relatedness between parents

# A<-read.table(paste0(odir, "/genos_modelled/Ka_Amat.txt") ,sep="\t")
parts<-unique(c(parents$Parent1 %>% as.character(), parents$Parent2 %>% as.character()))
Apar<-A[parts,parts]
# write.table(Apar, file=paste0(odir, "/genos_modelled/parents_Add_rel.txt") ,sep="\t")

## Mean relatedness between individuals in each family and the collection

mean_rel<-function(K, families, Name){
  whichCOL<-c(1:nrow(K))[-c(lapply(families, function(x) grep(x,rownames(K)) ) %>% unlist)]
  which_cat<-lapply(families, function(x) grep(x, rownames(K)))
  names(which_cat)<-families
  rel<-matrix(1, nrow=length(families)+1, ncol=length(families)+1, dimnames=list(c(families,"COLL"),c(families,"COLL")))
  for (i in 1:length(families)) {
    whichI<-grep(families[i], rownames(K))
    rel[i,"COLL"]<- rel["COLL",i] <-mean(K[whichI,whichCOL])
    for (j in c(1:length(families))[-i]) {
      whichJ<-grep(families[j], rownames(K))
      rel[i,j]<-rel[j,i]<-mean(K[whichI,whichJ])
    }
  }
  write.table(rel, file=paste0(odir, "/genos_modelled/mean_rel_", Name, ".txt"), sep="\t", row.names=T)
}

rel_Kadd_rblup<-mean_rel(A,families, "Add_mat_rblup")


###### mean additive relatedness of each parents to the collection
Acol<-A[WhichCOL, WhichCOL]
meanrel<-function(parent) {
  rel<-Acol[-which(colnames(Acol)==parent),parent]
  print(c(parent,mean(rel)))
}
lapply(unique(c(as.character(parents$Parent1), as.character(parents$Parent2))), meanrel)


###############################################################
######################## DAPC analysis ########################
###############################################################

## This analysis is performed on the collection and family offsprings are added as "supplementary" individuals

######### here silenced code because cluster names/numbers will change at every run of DAPC analysis ######

## reinitialized collection IDs just in case
# WhichCOL<-c(1:nrow(genos_round))[-c(lapply(families, function(x) grep(x,rownames(genos_round)) ) %>% unlist)]
# ## create new object with different additive coding (just shifted)
# geno_dapc<- new("genlight", (genos_round+1)[WhichCOL,])
# 
# cat("DAPC 1st step: Choose number of PCA axes to retain (all) et number of clusters\n")
# grp <- find.clusters(geno_dapc, max.n.clust=15)## keep 300 PCA axes and 6 clusters
# # saveRDS(grp, file=paste0(odir, "/genos_modelled/dapc_grp_300pca_6clusters.rds"))
# cat("DAPC 2d step: Choose number of PCA axes to retain (enough) and number of discriminant functions\n")
# dapc1 <- dapc(geno_dapc, grp$grp) ## keep 150 PCA axes and 5 discriminants functions
# # saveRDS(dapc1,file=paste0(odir, "/genos_modelled/DAPC_COLL_150PC_5functions.rds"))
# # write.table(dapc1$ind.coord, file=paste0(odir, "/genos_modelled/coordinates_DAPC_COL.txt"), row.names=T, quote=F, sep="\t")

### load already calculated clusters or run again analysis with above code ##

cat("taking previously calculated clusters")
dapc1<-readRDS(file=paste0(odir, "/genos_modelled/DAPC_COLL_150PC_5functions.rds"))
grp=readRDS(file=paste0(odir, "/genos_modelled/dapc_grp_300pca_6clusters.rds"))
# write.table(dapc1$ind.coord, file=paste0(odir, "/genos_modelled/to_keep/coordinates_DAPC_COL.txt"), row.names=T, quote=F, sep="\t")

## plot eigenvalues
# pdf(file=paste0(odir, "/genos_modelled/eingenvalues_myplot.pdf"), height=4, width=4)
barplot(names.arg = c(1:5), height=dapc1$eig, 
        ylab="F-statistic", xlab="Linear discriminants", main="Eigenvalues")
# dev.off()
### plot coordinates in collection
par(mfrow=c(1,1))
col=c("cyan", "blue", "purple", "green","red", "gold")
PARENTS<-unique(c(parents$Parent1 %>% as.character(), parents$Parent2 %>%as.character()))
df <- data.frame(x = dapc1$ind.coord[,1], y = dapc1$ind.coord[,3])
# identify/ create a vector of names for the individuals in your plot
whichPAR<-which(rownames(dapc1$ind.coord) %in% PARENTS)
noms <- rep(" ", nrow(dapc1$ind.coord))
# noms[whichPAR]<-rownames(dapc1$ind.coord)[whichPAR]
noms[whichPAR]<-c("Pink Lady", "Delearly", "Fuji", "Golden Delicious", "Pinova", "Royal Gala")
df["Pinova",2]<-df["Pinova",2]-0.17
df["RoyGal",2]<-df["RoyGal",2]+0.2
df["Delear",2]<-df["Delear",2]+0.2
df["FuMoHo",2]<-df["FuMoHo",2]-0.7
# pdf(file=paste0(odir, "/genos_modelled/DAPC_150PC_6clusters_PC2_PC4.pdf"), height=6, width=6)
# pdf(file=paste0(odir, "/genos_modelled/DAPC_paperA.pdf"), 
    # width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
# tiff(file=paste0(odir, "/genos_modelled/DAPC_paperA.tiff"), 
     # width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, units="in", res=300 )
scatter(dapc1,  xax=1, yax=3, grp=dapc1$grp,ratio.pca=0.5, bg="white", 
        pch=20, cell=0, cstar=0, col=col, solid=.4, cex=3, 
        clab=0, mstree=T, scree.da=FALSE, posi.leg="bottomright",
        leg=TRUE, txt.leg=paste("Cluster",1:6), xlim=c(-8,15), ylim=c(-1,1))
### online help for adding labels http://lists.r-forge.r-project.org/pipermail/adegenet-forum/2014-September/000959.html
# change graphical parameter to subsequently overlay the labels without
# drawing a new plot
par(new=TRUE)
# make a data frame of the dapc coordinates used in scatter
# use the text function to add labels to the positions given by the coordinates you used in plot
s.label(dfxy = df, xax=1, yax=2, label=c(noms),
        clabel=0.8, # change the size of the labels
        boxes=F, # if points are spaced wide enough, can use TRUE to add boxes around the labels
        grid=FALSE, addaxes=F,  xlim=c(-8,15), ylim=c(-1,1))# do not draw lines or axes in addition to the labels
# dev.off()

### format cluster assignements
clusters<- cbind(rownames(genos_round)[WhichCOL],dapc1$assign) ## keep, will be used in prediction
colnames(clusters)<-c("Name", "Cluster")

parents$Parent1<-as.character(parents$Parent1)
parents$Parent2<-as.character(parents$Parent2)
clusters[which(clusters[,"Name"] %in% c(parents$Parent1, parents$Parent2)),]

## all parents except Fuji in cluster 1

## assign offsprings to clusters and add families to the plot

sup_id<-c(1:nrow(genos_round))[-WhichCOL]
sup_x<-new("genlight", (genos_round+1)[sup_id,]) 
dapc_pred <- dapc(geno_dapc,dapc1$grp,n.pca=150,n.da=5)
pred.sup <- predict.dapc(dapc_pred, newdata=sup_x)
# saveRDS(dapc_pred, file=paste0(odir, "/genos_modelled/dapc_pred.rds"))
# saveRDS(pred.sup, file=paste0(odir, "/genos_modelled/pred.sup.rds"))
# write.table(pred.sup$ind.scores, file=paste0(odir, "/genos_modelled/coord_DAPC_families.txt"),sep="\t", quote=F, row.names=T)

cluster_fams<-data.frame(Name=rownames(genos_round)[sup_id], cluster= as.factor(as.character(pred.sup$assign)))
summary(cluster_fams)
# col <- rainbow(length(levels(dapc1$grp)))
dapc_pred<-ueadRDS(paste0(odir, "/genos_modelled/dapc_pred.rds"))
pred.sup<-readRDS(paste0(odir, "/genos_modelled/pred.sup.rds"))
col.points <- transp(col[as.integer(dapc1$grp)],.7)
# pdf(file=paste0(odir, "/genos_modelled/mapping_fam_DAPC.pdf"), height=6, width=6)
# scatter(dapc_pred, col=col, bg="white", scree.da=0, pch="",
#         cstar=0, clab=0, legend=TRUE, posi.leg="bottomright", xlim=c(-7,17), ylim=c(-5,4)) %>% print
# pdf(file=paste0(odir, "/genos_modelled/DAPC_paperB.pdf"), 
    # width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
# tiff(file=paste0(odir, "/genos_modelled/DAPC_paperB.tiff"), 
     # width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, units="in", res=300 )

scatter(dapc_pred, col=col, bg="white", scree.da=0, pch="",xax=1, yax=3,
        cstar=0, clab=0,  legend=F, xlim=c(-9,15), ylim=c(-1,1)) %>% print
par(xpd=TRUE)
col.sup <- col[as.integer(pred.sup$assign)]
points(pred.sup$ind.scores[,1], pred.sup$ind.scores[,3], pch=20,cex=1.5,
       col=transp(col.sup,.3), xlim=c(-8,15), ylim=c(-1,1))

# dev.off()

# export assignements in COL
# write.table(clusters, file=paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"),sep="\t", quote=F, row.names=F)
# write.table(cluster_fams, file=paste0(odir, "/genos_modelled/assignements_families.txt"), sep="\t", quote=F, row.names=F)

## look at distribution of clusters

clusters<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/genos_modelled/assignments_COLL_DAPC.txt", h=T)
clusters_fam<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/genos_modelled/assignements_families.txt", h=T)
clusters_fam$group=substr(clusters_fam$Name, 1,4)
clusters$group="Collection"
colnames(clusters_fam)<-colnames(clusters)
all_clusters=rbind(clusters, clusters_fam)
all_clusters$Cluster<-as.factor(all_clusters$Cluster)
all_clusters$group<-factor(all_clusters$group, levels=levels(all_clusters$group %>% as.factor) %>% rev)

p=ggplot(data=all_clusters, aes(x=Cluster, fill=group) )
# pdf(file=paste0(odir, "/genos_modelled/clusters_assignements_histo1.pdf"), height=6, width=6)
pp=p+geom_bar(stat="count",position = "stack",colour="black", size=0.2) +
  # scale_fill_manual(values=c("grey60",'#993404','#d95f0e', '#fe9929','#fec44f','#fee391','#ffffd4'))+
  # scale_fill_manual(values=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494'))+
  scale_fill_manual(values=c('#f7f7f7','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525'))+
  # scale_fill_manual(values=c('#034e7b','#045a8d','#2b8cbe','#74a9cf','#a6bddb','#d0d1e6',"grey60"))+
  labs(x="Cluster", y="Count", fill="Population")+
  theme_minimal()

# pdf(file=paste0(odir, "/genos_modelled/DAPC_paperC.pdf"),
    # width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, colormodel = "cmyk", )
# tiff(file=paste0(odir, "/genos_modelled/DAPC_paperC.tiff"), 
     # width= 3.54, height=3.54, fonts="Helvetica", pointsize=10, units="in", res=300 )
print(pp)
# dev.off()


# Purge obsolete variables
rm(genos_add, Coding.Add, Coding.Dom, Z, genos_round, genos_dom, genos_df)

# Garbage collection
gc()

cat("Genotyped analyzed!\n")
