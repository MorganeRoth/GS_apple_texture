cat("Importing data...\n")

# Phenos

# phenos<-read.table(paste0(idir, "/merging_train_pop_phenos.txt"), h=T, sep="\t") ## mean values
phenos<-read.table(paste0(idir, "/rep_meas_TEXTURE_4R.txt"), h=T, sep="\t") ## this for analyses
phenos$Year<-as.factor(phenos$Year)
levels(phenos$Name) %>% length
Names<-sort(levels(phenos$Name))
summary(phenos)

traits=colnames(phenos)[-c(1:6)]
cat("traits:")
traits
par(mfrow=c(3,4))
lapply(traits, function(i) print(hist(phenos[[i]], main="", xlab=i)))
## with reps dataset too big for shapiro test
# cat("shapiro test on raw data:")
# lapply(traits, function(i) print(c(i,shapiro.test(phenos[[i]])$p.value))) %>% unlist ## nothing normally distributed

## means per apple
phenos_techrep=phenos ## keep repeated data
phenos<-aggregate(cbind(sapply(traits, function(x) phenos[,x]))~ Year + Name + Apple + Location, data=phenos, mean)
summary(phenos)

# Genos

genos_add<-fread(paste0(idir, "/SNPs_additive_coding_04092019.txt"))[,-1] %>% as.matrix
rownames(genos_add)<-read.table(paste0(idir, "/rownames_SNPs_additive_coding_FuPi_cor.txt"), sep="\t", h=T)[,"name"] %>% as.character

## parents

parents<- read.table(paste0(idir, "/families_parents.txt"), sep="\t", h=T)

## remove outcrossers from both datasets
## all appart from GDFj_024 and FuPi_089 in outcrosser file

outcrossers<-c("FuPi_089", "FjPi_089",
               lapply(c("001","023","048","049","006","007","029","040","050"), function(x) paste0("GaPL_",x ))%>%unlist,
               lapply(c("028","060"), function(x) paste0("FjPL_",x ))%>%unlist,
               lapply(c("024","042","048", "050", "051", "058", "087"), function(x) paste0("GDFj_",x ))%>%unlist,
               lapply(c("065","066", "067", "068", "069", "070","071", "073", "075", "076", "077", "080"), function(x) paste0("FjDe_",x ))%>%unlist)
length(outcrossers)  

lapply(outcrossers, function(x) grep(x,rownames(genos_add))) %>% unlist
to_remove<-lapply(outcrossers, function(x) grep(x,rownames(genos_add))) %>% unlist
genos_add<-genos_add[-to_remove,]
phenos<-phenos[-which(phenos$Name %in% outcrossers),]

## Data already filter MAF >0.05 and call rate > 0.5 for individuals

## reshape for snpStats - recode 0,1,2 and NA are 9
# genos_fil<-genos_add
# genos_fil[which(is.na(genos_fil))]<-8
# genos_fil<- genos_fil + 1
# genos_fil[1:5,1:5]
# table(genos_fil)
# genos_fil<-new("SnpMatrix", as.matrix(genos_fil) )
# rownames(genos_fil)<-IDs$New_geno_name
# if needed: change IDs
# grep("MagGol", IDs)
# rownames(genos_add)[132]<-"Magr"

# ## filter data
# idsum<-row.summary(genos_fil)
# hist(idsum$Call.rate)
# snpsum<-col.summary(genos_fil)
# summary(snpsum)
# hwe=snpsum$z.HWE
# quantile(hwe, probs=seq(0,1,0.25), na.rm=T)
# q=quantile(hwe, probs=c(1e-4, (1-1.1e-4)), na.rm=T)
# hist(snpsum$Call.rate)
# hist(snpsum$z.HWE)
# hist(snpsum$MAF)
# ## filter for these values of Z 
# snpsel=snpsum[which(snpsum$Call.rate>0.20 & snpsum$MAF>0.05 ),] ## do not filter on HWE
# nrow(snpsum)-nrow(snpsel)## filtering out 2263 Mks
# idsum<-row.summary(genos_fil)
# summary(idsum)

# # Filter on MAF and call rate
# snpsel=rownames(snpsel)
# genos_fil<-genos_add[, which(colnames(genos_add) %in% snpsel)]
# length( which(colnames(genos_add) %in% snpsel))
# cat(ncol(genos_add)-ncol(genos_fil),"SNPs will be removed due to low MAF or low call rate") 
# # write.table(genot_sel, file="genotypes_filtered.txt", quote=F, row.names=F, sep="\t")
# # genot_sel=read.table("genotypes_filtered.txt", sep="\t", h=T)

## impute
genos_ready<-impute.knn(t(genos_add)) ## very few to impute due to duplicates
genos_ready<-t(genos_ready$data)
## found some spaces in rownames of genos_ready
grep(" ",rownames(genos_ready))
rownames(genos_ready)[573]<-"Limonc"
rownames(genos_ready)[66]<-"Coop30"

## identify spaces in names
grep(" ", rownames(genos_ready))
phenos$Name[grep(" ", phenos$Name)] ### Gala bukeye will be removed no problem

# ## now identify clones in genos
# identical<-list()
# dim(genos_ready)
# for (i in c(1:(nrow(genos_ready)-1))) {
#   for (j in (i+1):nrow(genos_ready)) {
#     print(i)
#     print(j)
#     if((summary(genos_ready[i,] == genos_ready[j,])["TRUE"]==ncol(genos_ready))){
#       identical<-append(identical, list(c(i,j)))
#     } else {
#       next
#     }
#   }
# }
# 
# lapply(identical, function(x)c(rownames(genos_ready)[x[[1]]], rownames(genos_ready)[x[[2]]]))
# saveRDS(identical, file=paste0(odir, "/genos_modelled/clones_identified_in_geno_ready.rds"))
## need to merge clones at phenotypic levels and to remove them in genotype file
identical=readRDS(file=paste0(odir, "/genos_modelled/clones_identified_in_geno_ready.rds"))
to_remove=lapply(identical, function(x) x[[1]]) %>% unlist
to_replace<-list(lapply(identical, function(x) rownames(genos_ready)[x[[1]]]) %>% unlist, lapply(identical, function(x) rownames(genos_ready)[x[[2]]]) %>% unlist)
genos_ready<-genos_ready[-c(to_remove),]
rownames(genos_ready)
## replace phenos names
phenos_test=phenos
lapply(1:length(to_replace [[1]]), function(x) {
  print(to_replace [[2]][[x]])
  which(phenos_test$Name == to_replace [[2]][[x]])})
## update to  replace: only two changes to make
to_replace[[1]]<-c("BadGol","GaPi_039")
to_replace[[2]]<-c("GolDel","GaPi_040")
phenos_test$Name<-as.character(phenos_test$Name)
mapply(function(x,y) phenos_test[which(phenos_test$Name == x), "Name"] <<- y , to_replace[[1]], to_replace[[2]] ) ## do not forget symbol << recursive to make an actual change!
lapply(to_replace[[1]], function(x) grep(x, phenos$Name) %>% length)
lapply(to_replace[[1]], function(x) grep(x, phenos_test$Name) %>% length)
lapply(to_replace[[2]], function(x) grep(x, phenos$Name) %>% length)
lapply(to_replace[[2]], function(x) grep(x, phenos_test$Name) %>% length)
phenos<-phenos_test
phenos$Name<-as.factor(phenos$Name)
rm(phenos_test)
## keep only what has been phenotyped/genotyped reciprocally
phenos$Name[which(!(phenos$Name %in% rownames(genos_ready)))] %>% as.character %>% unique
length(levels(phenos$Name) )
phenos$Name[which(!(levels(phenos$Name) %in% rownames(genos_ready)))] %>% as.character %>% unique %>% length
phenos$Name[which((phenos$Name %in% rownames(genos_ready)))] %>% as.character %>% unique%>% length
phenos<-phenos[which(phenos$Name %in% rownames(genos_ready)),]
phenos<-droplevels(phenos)
phenos$Trial<-paste0(phenos$Location, phenos$Year)
genos_ready<-genos_ready[which(rownames(genos_ready) %in% levels(phenos$Name)),]
dim(genos_ready)
length(levels(phenos$Name))
saveRDS(phenos, file=paste0(idir, "/phenos_ready_for_pred.rds"))
saveRDS(genos_ready,file=paste0(idir, "/genos_imputed_for_pred.rds"))
saveRDS(phenos_techrep,file=paste0(idir, "/phenos_with_reps_raw_filtered.rds"))
## useful lists

families<-c("FjDe", "GDFj", "FjPi", "FjPL", "GaPi", "GaPL")
lapply(families, function(x) grep(x, rownames(genos_ready)[rownames(genos_ready) %in% phenos$Name]) %>% length)
lapply(families, function(x) grep(x, rownames(genos_ready)[rownames(genos_ready) %in% phenos$Name])) %>% unlist %>% length
WhichFAM<-lapply(families, function(x) grep(x, rownames(genos_ready)[rownames(genos_ready) %in% phenos$Name])) %>% unlist
WhichCOL<-rownames(genos_ready)[-WhichFAM] 
## verify the absence of spaces in rownames
grep(" ", WhichCOL)
grep(" ", rownames(phenos))
length(WhichCOL)
# lapply(families, function(x) grep(x, rownames(genos_ready)[rownames(genos_ready) %in% phenos$Name])) %>% unlist %>% length
cat("Families:")
print(families)
ids<-rownames(genos_ready)
NbInd<- length(ids)

# Purge obsolete variables


cat("Data imported!\n")

