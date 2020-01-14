cat("Analyzing genos...\n")

dir.create(paste0(odir, "/genos_modelled"), showWarnings = FALSE, recursive = TRUE)
# Here we compute different kinship and make some plots
parents<- read.table(paste0(idir, "/families_parents.txt"), sep="\t", h=T)
# First we need to round imputed values
# genos_ready=readRDS(paste0(idir, "/genos_imputed_for_pred.rds"))
genos_round<-apply(genos_ready,2, function(x) as.numeric(x) %>% round(., digits=0))
summary(genos_round["FjPL_001",] == genos_round["FjPi_001",])
rownames(genos_round)<-rownames(genos_ready)
genos_round[genos_round==2]<-1
apply(genos_round[,1:20], 2, function(x) as.factor(x) %>%summary) ## now entire values only
cat("sanity check\n")
genos_round[genos_round==2]<-1
genos_round[genos_round==-2]<- -1
genos_ready[1:6,1:6] ## sanity check

## collection identifiers
WhichCOL=c(1:nrow(genos_round))[-c(lapply(families, function(x) grep(x, rownames(genos_round)))%>%unlist)]
rownames(genos_round)[WhichCOL]
## Code from Tristan and I made in June 2019
## Additive Vitezica 2017
## Check a bit the genotypic data
ID <- rownames(genos_round)
NbMarker <- ncol(genos_round)
NbID <- nrow(genos_round) 
dim(genos_round)
genos_round[1:10,1:10]

#### Get summary per marker

## Get contingency table per snp
SnpSummary <- lapply(1:NbMarker, function(snp){
  table(genos_round[,snp])
})

## Any non informative snps ?
InfSnp <- map_lgl(SnpSummary,~ ((.x[1]==NbID) | (.x[3]==NbID))  ) %>% 
  reduce(c)
cat("How many uninformative SNPs (TRUE)?\n")

print(table(InfSnp) )

## Get probs
Probs <- SnpSummary %>% map(., ~ .x/NbID) %>% reduce(rbind)
cat("Head of Probs")
head(Probs)
colnames(Probs)<-c("BB", "AB", "AA") ## rename for downstream coding
## Additive coding

Coding.Add <- Probs %>% 
  as.data.frame %>% 
  mutate(hAA = AB+2*BB, hAB = -(1-AB-2*BB), hBB = -(2-AB-2*BB)) %>% 
  select(hAA,hAB,hBB)

# head(Coding.Add)
## replace genotypes by their probabilities in genos_round
Geno.Add <- sapply(1:NbMarker, function(snp){
  New <- rep(Coding.Add[['hAA']][snp],NbID)
  New[as.character(genos_round[,snp])=='0'] <- Coding.Add[['hAB']][snp]
  New[as.character(genos_round[,snp])=='-1'] <- Coding.Add[['hBB']][snp]
  return(New)
})

## Get the dominance coding
Coding.Dom <- Probs %>% 
  as.data.frame %>% 
  mutate(Den = AA + BB -(AA-BB)**2,
         hAA = 2*AB*BB/Den, 
         hAB = 4*AA*BB/Den, 
         hBB = 2*AA*AB/Den) %>% 
  select(hAA,hAB,hBB)
# head(Coding.Dom)
Geno.Dom <- sapply(1:NbMarker, function(snp){
  New <- rep(Coding.Dom[['hAA']][snp],NbID)
  New[as.character(genos_round[,snp])=='0'] <- Coding.Dom[['hAB']][snp]
  New[as.character(genos_round[,snp])=='-1'] <- Coding.Dom[['hBB']][snp]
  return(New)
})

#### Kinship matrices

MakeK <- function(MatGeno){
  Mat <- tcrossprod(MatGeno) 
  Mat <- NbID*Mat/sum(diag(Mat))
  colnames(Mat) <- rownames(Mat) <- row.names(genos_round)
  return(Mat)
}

K.Add <- MakeK(Geno.Add)
dim(K.Add)
K.Add[300:310,300:310]
K.Dom <- MakeK(Geno.Dom)
summary(colnames(K.Dom)==rownames(K.Dom))
png(file=paste0(odir, "/genos_modelled/Ka_Vitezica17.png"), height=1500, width=1500)
heatmap(K.Add, symm=T) %>% print
dev.off()

png(file=paste0(odir, "/genos_modelled/Kd_Vitezica17.png"), height=1500, width=1500)
heatmap(K.Dom, symm=T) %>% print
dev.off()

write.table(K.Add, file=paste0(odir, "/genos_modelled/Ka_Vitezica17.txt"), row.names=T, quote=F, sep="\t")
write.table(K.Dom, file=paste0(odir, "/genos_modelled/Kd_Vitezica17.txt"), row.names=T, quote=F, sep="\t")

summary(diag(K.Dom))

## Epistasis: Hadamard product between all kinship matrices

K.AA<-hadamard.prod(K.Add,K.Add)
K.AA[1:5,1:5]

K.AD<-hadamard.prod(K.Add,K.Dom)
K.AD[1:5,1:5]

K.DD<-hadamard.prod(K.Dom,K.Dom)
K.DD[1:5,1:5]
for (epi in c("K.AA", "K.AD", "K.DD")) {
  png(file=paste0(odir, "/genos_modelled/", epi, "_Vitezica17.png"), height=1500, width=1500)
  heatmap(get(epi), symm=T) %>% print
  dev.off()
  write.table(get(epi), file=paste0(odir, "/genos_modelled/", epi, "_Vitezica17.txt"), row.names=T, quote=F, sep="\t")
}

# GetClust <- heatmap(K.Dom,symm = TRUE,keep.dendro = TRUE)$colInd

## realized additive relationship matrix (rrBLUP) from imputed data

A <- A.mat(genos_round)
colnames(A)<-rownames(A)<-rownames(genos_ready)
heatmap.2(A, trace="none",key.xlab="Relationship", key.ylab="", key.title = "", 
               keysize = 2, col=inferno(100) ) 
pdf(file=paste0(odir, "/genos_modelled/Ka_Amat.pdf"), height=6, width=6)
## viridis: inferno, viridis magma, plasma, cividis
myh<-heatmap.2(A, trace="none",key.xlab="Relationship", key.ylab="", key.title = "", 
               keysize = 2, col=inferno(100), labRow = FALSE, labCol = FALSE ) 

print(myh)
dev.off()
write.table(A, file=paste0(odir, "/genos_modelled/Ka_Amat.txt"), row.names=T, quote=F, sep="\t")
write.table(colnames(A)[myh$colInd], file=paste0(odir, "/genos_modelled/Clustering_Amat.txt"), row.names=T, quote=F, sep="\t")
## Means relatedness (additive) between families and collection
families<-c("FjDe", "FjPi", "FjPL", "GDFj", "GaPi", "GaPL")
WhichCOL<-c(1:nrow(K.Add))[-c(lapply(families, function(x) grep(x, rownames(K.Add)) ) %>% unlist)]
length(WhichCOL) ## 233 ids
rel<-matrix(NA, nrow=length(families), ncol=1, dimnames=list(families, "mean_rel_collection"))
for (i in families) {
  whichFAM<-c(1:nrow(K.Add))[grep(i, rownames(K.Add)) %>% unlist]
  rel[i,1]<-mean(K.Add[whichFAM,WhichCOL])
}
write.table(rel,paste0(odir, "/genos_modelled/mean_relatedness_families_coll.txt"), quote=F, row.names=T, sep="\t")

## Calculate inbreeding in each genotype: level of homozygosity (or count heterozygosity, easier)
## For this just make a mean of the dominance coding for each individual
## recode for dominance (0 become 1 and -1 and 1 become 0)

genos_df<-as.data.frame(genos_round) ## important: need a data frame for this funciton to work

myfun2<- function (data) {
  data<-recode(data, "1"="AA", "-1" ="aa", "0" ="1")
  data<-recode(data, "AA" ="0", "aa"="0")
}
## go parallel
no_cores<-detectCores()-10
no_cores

cl<- makeCluster(no_cores)
clusterExport(cl, "myfun2" )
clusterExport(cl, "recode" )

genos_dom<-parLapply(cl,genos_df , myfun2)
stopCluster(cl)

genos_dom<-as.data.frame(genos_dom)
genos_dom<-t(genos_dom)
genos_dom[1:5,1:5]
genos_dom<-apply(genos_dom, 2, as.numeric)
# sapply(1:10, function(x) mean(genos_dom[,x], na.rm=T))
hetero<-sapply(1:ncol(genos_dom), function(x) mean(genos_dom[,x], na.rm=T))
hetero<-data.frame(Het=hetero)
rownames(hetero)<-rownames(genos_round)
hetero
### replace bad name of  "Magr\xe8" to "Magr"
## Find its neighbour
# grep("MagGol",rownames(hetero)) ## row number 132
# rownames(hetero)[133]="Magr" 

## which individuals have low heterozygosity/ high inbreeding?

rownames(hetero)[which(hetero$Het == min(hetero$Het))] 
rownames(hetero)[which(hetero$Het == max(hetero$Het))] 

## create one dimension heatmap
# hetero$names<-rownames(hetero)
# hetero<-hetero[order(hetero$Het),]
# ggplot(hetero[1:nrow(hetero)/2,], aes(y = reorder(names, Het), x = 1, fill = Het)) + geom_tile() + theme(axis.text.x = element_text(angle=90))
# ggplot(hetero[nrow(hetero)/2+1:nrow(hetero),], aes(y = reorder(names, Het), x = 1, fill = Het)) + geom_tile() + theme(axis.text.x = element_text(angle=90))
head(hetero, n=20)
tail(hetero, n=20)
## save hetero

write.table(hetero, file=paste0(odir, "/genos_modelled/heterozygosity_mean.txt"), sep="\t", quote=F, row.names=T)

## create file for heterozygosity and relatedness between parents

parents<- read.table(paste0(idir, "/families_parents.txt"), sep="\t", h=T)
parents$relatednessKadd<- mapply(function(x,y) K.Add[x,y],as.character(parents$Parent1), as.character(parents$Parent2 ))
parents$relatednessAmat<- mapply(function(x,y) A[x,y],as.character(parents$Parent1), as.character(parents$Parent2 ))
parents$het1<- sapply(as.character(parents$Parent1),function(x) hetero[x,])
parents$het2<- sapply(as.character(parents$Parent2),function(x) hetero[x,])
parents$mean_Het<-sapply(parents$Family_name2, function(x) hetero[grep(x,rownames(hetero)), "Het"] %>% mean)

K.Add["FuMoHo", "Delear"]
write.table(parents, file=paste0(odir, "/genos_modelled/parents_het_relat.txt"), sep="\t", quote=F, row.names=F)

## create IBS for optimisation of training population -problem with this function
# Calculate the number of cores
# no_cores <- detectCores() - 5 
# Create structural matrix calc0 to call x and y 
# Z<-genos_round
# myIBS<- function(x) {
#   mysum=summary(Z[calc0[x,1],]==Z[calc0[x,2],]) 
#   if (!("FALSE" %in% names(mysum))) {
#     ibs<-as.numeric(mysum[["TRUE"]])/as.numeric(mysum[["TRUE"]])
#   } else   {
#     ibs<-as.numeric(mysum[["TRUE"]])/(as.numeric(mysum[["TRUE"]]) + as.numeric(mysum[["FALSE"]]))
#   }
# }
# calc0<-matrix(c(as.character(sort(rep(rownames(Z), length(rownames(Z))))),   
#                 as.character(rep(rownames(Z), length(rownames(Z))))), 
#               ncol=2, nrow= length(rownames(Z))^2, byrow=F)
# 
# start_time <- Sys.time()
# cl <- makeCluster(no_cores)
# clusterExport(cl, "calc0")
# clusterExport(cl, "Z")
# IBS_values<-parLapply(cl, 1:nrow(calc0), myIBS)
 # function(x) summary(Z[calc0[x,1],]==Z[calc0[x,2],])[3])

# stopCluster(cl)
# end_time <- Sys.time()
# end_time - start_time
# IBS<-matrix(c(unlist(IBS_values)), nrow=length(rownames(Z)), ncol=length(rownames(Z)), byrow=T,
#              dimnames=list(rownames(Z),rownames(Z))) 

IBSnew<-matrix(NA, nrow=nrow(A), ncol=nrow(A), dimnames=list(rownames(A), colnames(A)))

for (i in 1:nrow(IBSnew)){
  for (j in 1:ncol(IBSnew)) {
    mysum=summary(Z[i,]==Z[j,]) 
    if (!("FALSE" %in% names(mysum))) {
      ibs<-as.numeric(mysum[["TRUE"]])/as.numeric(mysum[["TRUE"]])
    } else   {
      ibs<-as.numeric(mysum[["TRUE"]])/(as.numeric(mysum[["TRUE"]]) + as.numeric(mysum[["FALSE"]]))
    }
    IBSnew[i,j]<-ibs
    print(c(i,j,ibs))
  }
}
summary(colnames(Z)==rownames(Z))

write.table(IBSnew, file=paste0(odir, "/genos_modelled/IBS.txt"), sep="\t", quote=F, row.names=T)
png(file=paste0(odir,"/genos_modelled/IBS_heatmap.png"), height=2000, width=2000)
print(heatmap(IBSnew))
dev.off()

pdf(file=paste0(odir, "/genos_modelled/IBS_heatmap2.pdf"), height=6, width=6)
## viridis: inferno, viridis magma, plasma, cividis
myh<-heatmap.2(IBSnew, trace="none",key.xlab="IBS", key.ylab="", key.title = "", 
               keysize = 2, col=viridis(100), labRow = FALSE, labCol = FALSE ) 
print(myh)
dev.off()
write.table(colnames(IBSnew)[myh$colInd], file=paste0(odir, "/genos_modelled/Clustering_IBS.txt"), row.names=T, quote=F, sep="\t")

#write.table(IBS2, file=paste0(odir, "/genos_modelled/IBS2.txt"), sep="\t", quote=F, row.names=T)

### check that there are no clones
sell=list()
for (x in  1:nrow(IBSnew)) {
  sel<-names(which(IBSnew[x, ]==1))
  if (length(sel) >1) {
    sell=append(sell, list(sel))
    print(sel)
  } else {
    next
  }
} 
sell

# output mean genetic relationship of each family to collection and other families
# list positions of identifiers for each category

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

rel_IBSnew<-mean_rel(IBSnew,families,"IBSnew")
rel_Kadd<-mean_rel(K.Add,families, "Kadd_Vitezica17")
rel_Kadd_rblup<-mean_rel(A,families, "Add_mat_rblup")

## Look at structure with MDS, PCA and DAPC
## use only structure
rownames(genos_round)
WhichCOL<-c(1:nrow(genos_round))[-c(lapply(families, function(x) grep(x,rownames(genos_round)) ) %>% unlist)]
groups<-data.frame(Names=rownames(genos_round), Group=NA)
groups$Group[WhichCOL]<-"Collection"
groups$Group[-WhichCOL]<-lapply(groups$Name[-WhichCOL] ,function(x) substr(x,1,4)) %>% unlist
res.pca <- PCA(genos_round, graph = F, ind.sup=c(1:nrow(genos_round))[-WhichCOL])
res.pca2 <- PCA(genos_round, graph = F, ind.sup=c(1:nrow(genos_round))[-WhichCOL])
png(file=paste0(odir, "/genos_modelled/PCA_collection.png"), width=1200, height=800, res=150)
fviz_pca_ind(res.pca)
fviz_pca_ind(res.pca,  col.ind.sup="blue", label="none" )
dev.off()

###############################################################
######################## DAPC analysis ########################
###############################################################

######### here silenced code for running analysis again (cluster names will change!) ######
# WhichCOL<-c(1:nrow(genos_round))[-c(lapply(families, function(x) grep(x,rownames(genos_round)) ) %>% unlist)]
# geno_dapc<- new("genlight", (genos_round+1)[WhichCOL,])
# # table(genos_round+1)
# cat("DAPC 1st step: Choose number of PCA axes to retain (all) et number of clusters\n")
# grp <- find.clusters(geno_dapc, max.n.clust=15)## keep 300 PCA axes and 6 clusters
# saveRDS(grp, file=paste0(odir, "/genos_modelled/dapc_grp_300pca_6clusters.rds"))
# cat("DAPC 2d step: Choose number of PCA axes to retain (enough) and number of discriminant functions\n")
# dapc1 <- dapc(geno_dapc, grp$grp) ## keep 150 PCA axes and 5 discriminants functions
# saveRDS(dapc1,file=paste0(odir, "/genos_modelled/DAPC_COLL_150PC_5functions.rds"))
# write.table(dapc1$ind.coord, file=paste0(odir, "/genos_modelled/coordinates_DAPC_COL.txt"), row.names=T, quote=F, sep="\t")

### load already calculated clusters or run again analysis ##

cat("taking previously calculated clusters")

# dapc1<-readRDS(file=paste0(odir, "/genos_modelled/DAPC_COLL_150PC_5functions.rds"))
# grp=readRDS(file=paste0(odir, "/genos_modelled/dapc_grp_300pca_6clusters.rds"))
# write.table(dapc1$ind.coord, file=paste0(odir, "/genos_modelled/to_keep/coordinates_DAPC_COL.txt"), row.names=T, quote=F, sep="\t")

## plot eigenvalues
pdf(file=paste0(odir, "/genos_modelled/eingenvalues_myplot.pdf"), height=4, width=4)
barplot(names.arg = c(1:5), height=dapc1$eig, 
        ylab="F-statistic", xlab="Linear discriminants", main="Eigenvalues")
dev.off()
### plot coordinates in collection
par(mfrow=c(1,1))
pdf(file=paste0(odir, "/genos_modelled/DAPC_150PC_6clusters_PC2_PC4.pdf"), height=6, width=6)
scatter(dapc1,  xax=2, yax=3, grp=dapc1$grp,ratio.pca=0.5, bg="white", 
        pch=20, cell=0, cstar=0, col=rainbow(6), solid=.4, cex=3, 
        clab=0, mstree=TRUE, scree.da=FALSE, posi.leg="bottomright",
        leg=TRUE, txt.leg=paste("Cluster",1:6), xlim=c(-10,10), ylim=c(-1,1))


### online help for adding labels http://lists.r-forge.r-project.org/pipermail/adegenet-forum/2014-September/000959.html
# change graphical parameter to subsequently overlay the labels without
# drawing a new plot
par(new=TRUE)
# make a data frame of the dapc coordinates used in scatter
PARENTS<-unique(c(parents$Parent1 %>% as.character(), parents$Parent2 %>%as.character()))
df <- data.frame(x = dapc1$ind.coord[,2], y = dapc1$ind.coord[,3])
# identify/ create a vector of names for the individuals in your plot
whichPAR<-which(rownames(dapc1$ind.coord) %in% PARENTS)
noms <- rep(" ", nrow(dapc1$ind.coord))
noms[whichPAR]<-rownames(dapc1$ind.coord)[whichPAR]
# df["CriPin",2]<-df["CriPin",2]+0.1
# df["RoyGal",2]<-df["RoyGal",2]-0.1
# use the text function to add labels to the positions given by the coordinates you used in plot
s.label(dfxy = df, xax=1, yax=2, label=c(noms),
        clabel=1, # change the size of the labels
        boxes=F, # if points are spaced wide enough, can use TRUE to add boxes around the labels
        grid=FALSE, addaxes=F, xlim=c(-10,10), ylim=c(-1,1))# do not draw lines or axes in addition to the labels
dev.off()

### format cluster assignements
clusters<- cbind(rownames(genos_round)[WhichCOL],dapc1$assign) ## keep, will be used in prediction
colnames(clusters)<-c("Name", "Cluster")

parents$Parent1<-as.character(parents$Parent1)
parents$Parent2<-as.character(parents$Parent2)
clusters[which(clusters[,"Name"] %in% c(parents$Parent1, parents$Parent2)),]

## all parents except Fuji in cluster 1

## add families

sup_id<-c(1:nrow(genos_round))[-WhichCOL]
sup_x<-new("genlight", (genos_round+1)[sup_id,]) 
dapc_pred <- dapc(geno_dapc,dapc1$grp,n.pca=150,n.da=5)


pred.sup <- predict.dapc(dapc_pred, newdata=sup_x)

# write.table(pred.sup$ind.scores, file=paste0(odir, "/genos_modelled/to_keep/coord_DAPC_families.txt"),sep="\t", quote=F, row.names=T)

cluster_fams<-data.frame(Name=rownames(genos_round)[sup_id], cluster= as.factor(as.character(pred.sup$assign)))
summary(cluster_fams)
col <- rainbow(length(levels(dapc1$grp)))
col.points <- transp(col[as.integer(dapc1$grp)],.7)
pdf(file=paste0(odir, "/genos_modelled/mapping_fam_DAPC.pdf"), height=6, width=6)
# scatter(dapc_pred, col=col, bg="white", scree.da=0, pch="",
#         cstar=0, clab=0, legend=TRUE, posi.leg="bottomright", xlim=c(-7,17), ylim=c(-5,4)) %>% print

scatter(dapc_pred, col=col, bg="white", scree.da=0, pch="",xax=2, yax=3,
        cstar=0, clab=0,  legend=F,xlim=c(-10,10), ylim=c(-1,1)) %>% print
par(xpd=TRUE)
## print COLL
# points(dapc_pred$ind.coord[,2], dapc_pred$ind.coord[,4], pch=15,
#        col=col.points, cex=2)
## print families
col.sup <- col[as.integer(pred.sup$assign)]
points(pred.sup$ind.scores[,2], pred.sup$ind.scores[,3], pch=20,
       col=transp(col.sup,.3), cex=3,xlim=c(-10,10), ylim=c(-1,1))
# points(pred.sup$ind.scores[,2], pred.sup$ind.scores[,4], pch=20,
#        col=transp("black",.3), cex=2)
# add.scatter.eig(dapc_pred$eig,15,1,2, posi="bottomright", inset=.02) %>% print
dev.off()

# export assignements in COL
write.table(clusters, file=paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"),sep="\t", quote=F, row.names=F)
write.table(cluster_fams, file=paste0(odir, "/genos_modelled/assignements_families.txt"), sep="\t", quote=F, row.names=F)

## look at distribution of clusters

clusters<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/genos_modelled/assignments_COLL_DAPC.txt", h=T)
clusters_fam<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/R_output/genos_modelled/assignements_families.txt", h=T)
summary(clusters_fam)
clusters_fam$group=substr(clusters_fam$Name, 1,4)
clusters$group="Collection"
colnames(clusters_fam)<-colnames(clusters)
all_clusters=rbind(clusters, clusters_fam)
all_clusters$Cluster<-as.factor(all_clusters$Cluster)
all_clusters$group<-factor(all_clusters$group, levels=levels(all_clusters$group %>% as.factor) %>% rev)
barplot(table(all_clusters$group,all_clusters$Cluster), col=rainbow(7))

p=ggplot(data=all_clusters, aes(x=Cluster, fill=group) )
pdf(file=paste0(odir, "/genos_modelled/clusters_assignements_histo0.pdf"), height=6, width=6)
pp=p+geom_bar(stat="count",position = "stack",colour="black", size=0.2) +
  # scale_fill_manual(values=c("grey60",'#993404','#d95f0e', '#fe9929','#fec44f','#fee391','#ffffd4'))+
  # scale_fill_manual(values=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494'))+
  # scale_fill_manual(values=c('#f7f7f7','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525'))+
  scale_fill_manual(values=c('#034e7b','#045a8d','#2b8cbe','#74a9cf','#a6bddb','#d0d1e6',"grey60"))+
  labs(x="Cluster", y="Count", fill="Population")+
  theme_minimal()
print(pp)
dev.off()

# Purge obsolete variables
rm(genos_add, Coding.Add, Coding.Dom, Z, genos_round, genos_dom, genos_df)

# Garbage collection
gc()

cat("Genotyped analyzed!\n")
