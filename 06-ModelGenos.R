cat("Analyzing genos...\n")

dir.create(paste0(odir, "/genos_modelled"), showWarnings = FALSE, recursive = TRUE)
# Here we compute different kinship and make some plots

# First we need to round imputed values
genos_round<-apply(genos_ready,2, function(x) as.numeric(x) %>% round(., digits=0))
rownames(genos_round)<-ids
apply(genos_round[,1:20], 2, function(x) as.factor(x) %>%summary) ## now entire values only
cat("sanity check")
genos_round[genos_round==2]<-1
genos_round[genos_round==-2]<- -1
genos_ready[1:6,1:6] ## sanity check

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
cat("How many uninformative SNPs (TRUE)?")
table(InfSnp) 

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
  colnames(Mat) <- row.names(Mat) <- row.names(genos_round)
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
rownames(genos_ready)<-rownames(genos_ready)
A <- A.mat(genos_ready)
colnames(A)<-rownames(A)<-rownames(genos_ready)
png(file=paste0(odir, "/genos_modelled/Ka_Amat.png"), height=1500, width=1500)
heatmap(A, symm=T) %>% print
dev.off()
write.table(A, file=paste0(odir, "/genos_modelled/Ka_Amat.txt"), row.names=T, quote=F, sep="\t")

## Means relatedness (additive) between families and collection
families<-c("FjDe", "FuPi", "FjPL", "GDFj", "GaPi", "GaPL")
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
hetero<-sapply(1:nrow(genos_dom), function(x) mean(as.numeric(genos_dom[x,]), na.rm=T))
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
parents$relatedness<- mapply(function(x,y) K.Add[x,y],as.character(parents$Parent1), as.character(parents$Parent2 ))
parents$het1<- sapply(as.character(parents$Parent1),function(x) hetero[x,])
parents$het2<- sapply(as.character(parents$Parent2),function(x) hetero[x,])
parents$mean_Het<-sapply(1:6, function(x) mean(c(parents[x,"het1"], parents[x,"het2"])))
K.Add["FuMoHo", "Delear"]
write.table(parents, file=paste0(odir, "/genos_modelled/parents_het_relat.txt"), sep="\t", quote=F, row.names=F)

## create IBS for optimisation of training population
# Calculate the number of cores
no_cores <- detectCores() - 5 
# Create structural matrix calc0 to call x and y 
Z<-genos_round
myIBS<- function(x) {
  mysum=summary(Z[calc0[x,1],]==Z[calc0[x,2],]) 
  if (!("FALSE" %in% names(mysum))) {
    ibs<-as.numeric(mysum[["TRUE"]])/as.numeric(mysum[["TRUE"]])
  } else   {
    ibs<-as.numeric(mysum[["TRUE"]])/(as.numeric(mysum[["TRUE"]]) + as.numeric(mysum[["FALSE"]]))
  }
}
calc0<-matrix(c(as.character(sort(rep(rownames(Z), length(rownames(Z))))),   ## create this file to be able to append Ka values between TS and TRS
                as.character(rep(rownames(Z), length(rownames(Z))))), 
              ncol=2, nrow= length(rownames(Z))^2, byrow=F)

start_time <- Sys.time()
cl <- makeCluster(no_cores)
clusterExport(cl, "calc0")
clusterExport(cl, "Z")
IBS_values<-parLapply(cl, 1:nrow(calc0), myIBS)
 # function(x) summary(Z[calc0[x,1],]==Z[calc0[x,2],])[3])

stopCluster(cl)
end_time <- Sys.time()
end_time - start_time
IBS<-matrix(c(unlist(IBS_values)), nrow=length(rownames(Z)), ncol=length(rownames(Z)), byrow=T,
             dimnames=list(rownames(Z),rownames(Z))) 
IBS[1:5,1:5]
write.table(IBS, file=paste0(odir, "/genos_modelled/IBS.txt"), sep="\t", quote=F, row.names=F)
png(file=paste0(odir,"/genos_modelled/IBS_heatmap.png"), height=2000, width=2000)
print(heatmap(IBS))
dev.off()

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

rel_IBS<-mean_rel(IBS,families,"IBS")
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

fviz_pca_ind (res.pca)

geno_dapc<- new("genlight", (genos_round+1)[WhichCOL,]) 
# table(genos_round+1)
grp <- find.clusters(geno_dapc, max.n.clust=15)## keep 300 PCA axes and 5 clusters
dapc1 <- dapc(geno_dapc, grp$grp) ## keep 150 PCA axes and 4 discriminants functions
par(mfrow=c(1,1))
pdf(file=paste0(odir, "/genos_modelled/DAPC_150PC_5clusters.pdf"), height=6, width=8)
scatter(dapc1, ratio.pca=0.3, bg="white", 
        pch=20, cell=0, cstar=0, col=rainbow(5), solid=.4, cex=3, 
        clab=0, mstree=TRUE, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, txt.leg=paste("Cluster",1:5))
dev.off()
# export assignements in COL
saveRDS(dapc1,file=paste0(odir, "/genos_modelled/DAPC_COLL_150PC_4clusters.rds"))
clusters<- cbind(rownames(genos_round)[WhichCOL],dapc1$assign) ## keep, will be used in prediction
colnames(clusters)<-c("Name", "Cluster")
write.table(clusters, file=paste0(odir, "/genos_modelled/assignments_COLL_DAPC.txt"),sep="\t", quote=F, row.names=F)

## add families
sup_id<-c(1:nrow(genos_round))[-WhichCOL]
sup_x<-new("genlight", (genos_round+1)[sup_id,]) 
dapc_pred <- dapc(geno_dapc,grp$grp,n.pca=150,n.da=4)
pred.sup <- predict.dapc(dapc_pred, newdata=sup_x)
names(pred.sup)
cluster_fams<-data.frame(Name=rownames(genos_round)[sup_id], cluster= pred.sup$assign)
head(cluster_fams)
col <- rainbow(length(levels(grp$grp)))
col.points <- transp(col[as.integer(grp$grp)],.7)
pdf(file=paste0(odir, "/genos_modelled/mapping_fam_DAPC.pdf"), height=6, width=8)
scatter(dapc_pred, col=col, bg="white", scree.da=0, pch="",
        cstar=0, clab=0, xlim=c(-10,10), legend=TRUE)

par(xpd=TRUE)
## print COLL
points(dapc_pred$ind.coord[,1], dapc_pred$ind.coord[,2], pch=15,
       col=col.points, cex=2)
## print families
col.sup <- col[as.integer(pred.sup$assign)]
points(pred.sup$ind.scores[,1], pred.sup$ind.scores[,2], pch=20,
       col=transp(col.sup,.3), cex=3)
points(pred.sup$ind.scores[,1], pred.sup$ind.scores[,2], pch=20,
       col=transp("black",.3), cex=2)
add.scatter.eig(dapc_pred$eig,15,1,2, posi="bottomright", inset=.02)
dev.off()

write.table(cluster_fams, file=paste0(odir, "/genos_modelled/assignements_families.txt"), sep="\t", quote=F, row.names=F)
# Purge obsolete variables


# Garbage collection
gc()

cat("Data transformed!\n")
