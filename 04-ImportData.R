cat("Importing data...\n")

# Phenos

phenos<-read.table(paste0(idir, "/merging_train_pop_phenos.txt"), h=T, sep="\t")
phenos_raw<-read.table(paste0(idir, "/merging_train_pop_phenos.txt"), h=T, sep="\t")
phenos$Year<-as.factor(phenos$Year)
Names<-sort(levels(phenos$Name))
summary(phenos)
traits=colnames(phenos)[-c(1:7)]
cat("traits:")
traits
par(mfrow=c(3,4))
lapply(traits, function(i) print(hist(phenos[[i]], main="", xlab=i)))
cat("shapiro test on raw data:")
lapply(traits, function(i) print(c(i,shapiro.test(phenos[[i]])$p.value))) %>% unlist ## nothing normally distributed

# Genos

genos_add<-fread(paste0(idir, "/all_genotypes_coll_progenies_additive.txt")) 

new_names<-read.table(paste0(idir, "/new_names_genos_4progenies.txt"), header=T) 
IDs<-merge(genos_add[,1], new_names, by.x="V1", by.y="Geno_name",sort=F )[,2] ## keep order of genos_add
genos_add<-genos_add[,-1]
rownames(genos_add)<-IDs$V1
summary(genos_add[,1:10])
genos_add<-as.matrix(genos_add)
## reshape for snpStats - recode 0,1,2 and NA are 9
genos_fil<-genos_add
genos_fil[genos_fil==NA]<-8
genos_fil<- genos_fil + 1
genos_fil[1:5,1:5]
genos_fil<-new("SnpMatrix", as.matrix(genos_fil) )
rownames(genos_fil)<-IDs$New_geno_name
# if needed: change IDs
# grep("MagGol", IDs)
# rownames(genos_add)[132]<-"Magr"

## filter data

snpsum<-col.summary(genos_fil)
summary(snpsum)
hwe=snpsum$z.HWE
quantile(hwe, probs=seq(0,1,0.25), na.rm=T)
q=quantile(hwe, probs=c(1e-4, (1-1.1e-4)), na.rm=T)
hist(snpsum$Call.rate)
hist(snpsum$z.HWE)
## filter for these values of Z 
snpsel=snpsum[which(snpsum$Call.rate>0.20 & snpsum$MAF>0.05 ),] ## do not filter on HWE
nrow(snpsum)-nrow(snpsel)## filtering out 3020 Mks
idsum<-row.summary(genos_fil)
summary(idsum)

## filters out only 1450 SNPs
# Filter on MAF and call rate
snpsel=rownames(snpsel)
genos_fil<-genos_add[, which(colnames(genos_add) %in% snpsel)]
length( which(colnames(genos_add) %in% snpsel))
cat(ncol(genos_add)-ncol(genos_fil),"SNPs will be removed due to low MAF or low call rate") 
# write.table(genot_sel, file="genotypes_filtered.txt", quote=F, row.names=F, sep="\t")
# genot_sel=read.table("genotypes_filtered.txt", sep="\t", h=T)

## impute
genos_ready<-impute.knn(t(genos_fil)) 
genos_ready<-t(genos_ready$data)
rownames(genos_ready)<-IDs$New_geno_name
saveRDS(genos_ready,file=paste0(idir, "/genos_imputed_for_pred.rds"))

apply(genos_ready[,1:10], 2, function(x) as.factor(x) %>%summary) 
cat("some imputed values visible here, need to round for kinship")

## useful lists

families<-c("FjDe", "GDFj", "FuPi", "FjPL", "GaPi", "GaPL")
cat("Families:")
print(families)
ids<-rownames(genos_ready)
NbInd<- length(ids)

# Purge obsolete variables
rm(genos_fil, genos_add, snpsel)

cat("Data imported!\n")
