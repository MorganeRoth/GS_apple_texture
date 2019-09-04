
### new geno data 03.09.2019

genew<-fread("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Raw_data/Genotypes/FABRIZIO_110414_FinalReport_design_strand_cut_headers.txt", h=T)

dim(genew)
genew[1:5,1:5]
rownames(genew)<-genew$V1
write.table(genew$V1, file="names.txt")

### markers selected in previous work

genor<-fread("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Raw_data/Genotypes/Training_pop/data_10K_correctMc38_Magr.txt")
corresp<-read.table("/home/evdad.admin.ch/a80848168/projects/GenSel_Italy_no_GIT/Genos_and_kinships/Recoding_SNPs/file_for_recoding_progenies.txt", h=T)
dim(corresp)
dim(genor)
corresp %>% head
corresp[,1:4] %>% head
genor[1:5,1:5]
which(colnames(genor) %in% genew$V1) %>% length
write.table(corresp[,1:4], file="~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Genos_analyzed/corresp_original_new_SNP_names.txt", sep="\t", quote=F, row.names=F)

## names correspond to progeny file
which(corresp$Name_progeny_file %in% genew$V1) %>% length
sel<-which(genew$V1 %in% corresp$Name_progeny_file) 
length(sel) ## select these SNP and work only with those

genew<-genew[sel,] %>% as.matrix %>% t
colnames(genew)<-genew[1,]
genew<-genew[-1,]
genew[1:5,1:5]
IDs<-rownames(genew)
dim(genew)

## we have 994 genotypes here, should include all the ones we need
## find a fast way to recode everything in additive way!

library(snpReady)
genew[which(genew == "NA")] = NA
table(genew)
## remove individuals with more than 50% missing data and impute with knn, recode with -101 (additive)
## keep all individuals, it helps for imputation
x<-raw.data(genew, frame="wide",base=T, maf=0,sweep.sample=0.5,outfile="-101",imput=T, imput.type="knni")
final<-x$M.clean
dim(final)
final[1:5,1:5]

## keep only useful individuals: collection and families
rownames(final)
newnames<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Genos_analyzed/new_names_genos_4progenies.txt", h=T)
head(newnames)
dim(newnames)
rownames(final) %in% newnames$Geno_name %>% summary
rownames(final)[rownames(final) %in% newnames$Geno_name]

# retrieve collection names by hand in excel
write.table(rownames(final), file="names.txt")

corresp2<-read.csv("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Genos_analyzed/corresp_id_names_03092019.csv", h=T, sep=";")
head(corresp2)
corresp2$Name_Fabrizio <- corresp2$Name_Fabrizio  %>% as.character()

final[which(rownames(final) %in% corresp2$Name_Fabrizio),]
which(rownames(final) %in% corresp2$Name_Fabrizio) %>% length

final<-final[which(rownames(final) %in% corresp2$Name_Fabrizio), ]
rownames(corresp2)<-corresp2$Name_Fabrizio %>% as.character
corresp2<-corresp2[rownames(final),]
corresp2 %>% dim
# corresp2<-droplevels(corresp2)
summary(rownames(final) == rownames(corresp2))
rownames(final)<-corresp2$Name_Morgane %>% as.character
final[1:5,1:5]

fwrite(final, file="~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Genos_analyzed/new_genos_03092019.txt")
