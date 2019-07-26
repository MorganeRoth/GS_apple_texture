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
# table(InfSnp) 

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

## Means relatedness between families and collection

cat("TO DO: mean relatedness incorporate\n")

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


# Purge obsolete variables
rm(Coding.Add, Coding.Dom)

# Garbage collection
gc()

cat("Data transformed!\n")
