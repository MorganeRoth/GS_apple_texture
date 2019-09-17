
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
genew[which(genew == "--")] = NA
# table(genew)
## remove individuals with more than 50% missing data and impute with knn, recode with -101 (additive)
## keep all individuals, it helps for imputation
x<-raw.data(genew, frame="wide",base=T, maf=0.05,sweep.sample=0.5,outfile="-101",imput=T, imput.type = "knni")
final<-x$M.clean
## keep only useful individuals: collection and families
# rownames(final) %>% table %>%  table ## no double name, but with new IDs it should change
# newnames<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Genos_analyzed/new_names_genos_4progenies.txt", h=T)
# rownames(final) %in% newnames$Geno_name %>% summary  ## 417 names included in the data we look for

# retrieve collection names by hand in excel and add families
# write.table(rownames(final), file="names.txt")

corresp2<-read.csv("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Genos_analyzed/corresp_id_names_03092019.csv", h=T, sep=";")
head(corresp2) ## here we have old and new names for each individual
final<-final[which(rownames(final) %in% as.character(corresp2$Name_Fabrizio %>% as.character)),]
dim(final)
## some duplicates
table(corresp2$Name_Morgane) %>% table
double<-names(which(table(corresp2$Name_Morgane) >1))
doubles<-corresp2[which(corresp2$Name_Morgane %in% double),] %>% droplevels
summary(final["2_25_Alkmene",] == final["1_95_ALKMENE",])
## remove duplicates
for (x in levels(doubles$Name_Morgane)) {
  NAMES<-doubles[which (doubles$Name_Morgane==x), "Name_Fabrizio"] %>% as.character()
  comb<-combn(NAMES,2)
  for (y in c(1:choose(length(NAMES),2))){
    print(comb[,y])
    summary(final[comb[1,y],] == final[comb[2,y],]) %>% print
  }  
}

to_update<-matrix(NA, nrow=length(levels(doubles$Name_Morgane)), ncol=ncol(final),
                  dimnames=list(c(levels(doubles$Name_Morgane)), colnames(final)))
## ALKMEN
alkmen<-which(final["2_25_Alkmene",]==final["1_95_ALKMENE",]) ## select markers with exact same genotype
to_update["Alkmen",alkmen]<-final["2_25_Alkmene",alkmen] ## update where it's exact, leave otherwise NA
## COOP 15
Coop15<-which(final["5_59_COOP_15",]==final["1_58_CO-OP_15",]) ## select markers with exact same genotype
to_update["Coop15",Coop15]<-final["5_59_COOP_15",Coop15]
## Pink Lady (CriPin)
PL<-which(final["3_84_Cripps_Pink_(PINK_LADY)",]==final["FxPL_Pink Lady",]) ## select markers with exact same genotype
to_update["CriPin",PL]<-final["FxPL_Pink Lady",PL]
## Delear
Del<-which(final["5_90_DELEARLY",]==final["FxD_Delearly",]) ## select markers with exact same genotype
to_update["Delear",Del]<-final["FxD_Delearly",Del]
## Fuji
which(final["5_68_FUJI_MORI_HOUFU_3a",]!=final["FxD_Fuji",])
which(final["FxG_Fuji",]!=final["5_68_FUJI_MORI_HOUFU_3a",])
which(final["FxG_Fuji",]!=final["FxD_Fuji",])
## I deduce from this that at SNP_FB_0528089 "5_68_FUJI_MORI_HOUFU_3a" and "FxD_Fuji" there is one mismatch 
## this mismatch exists also between "5_68_FUJI_MORI_HOUFU_3a" and "FxG_Fuji"
## which is not found between FxD_Fuji and FxG_Fuji: so they should have the correct genotype
## let's consider the sequence of  "FxD_Fuji" as the good one
to_update["FuMoHo", ]<-final["FxD_Fuji",]

## Gol Del
## "3_68_GOLDEN_DEL_CL_B" and "FxG_Golden"  have the same genotypes we keep one of both
to_update["GolDel", ]<-final["FxG_Golden",]

## Oregon Spur
OS<-which(final["5_18_OREGON_SPUR",]==final["5_19_OREGON_SPUR",]) ## select markers with exact same genotype
to_update["OreSpu",OS]<-final["5_19_OREGON_SPUR",OS]

dim(to_update)

## remove doutful samples and paste updated ones with missing values
final<-final[-which(rownames(final) %in% c(doubles$Name_Fabrizio %>% as.character)),]
keep<-rownames(final)
final<-rbind(final,to_update )

## rename
new<-lapply(rownames(final)[c(1:(nrow(final)-7))], function(x) {
corresp2[which(as.character(corresp2$Name_Fabrizio) == x), "Name_Morgane"] %>% as.character }) %>% unlist
head(new, n=10L)
head(rownames(final), n=10L)

tail(new, n=10L)
tail(rownames(final), n=17L)

rownames(final)[c(1:(nrow(final)-7))]<-new

final<-final[order(rownames(final)),]
dim(final)

fwrite(final, file=paste0(idir, "/SNPs_additive_coding_04092019.txt"), row.names = T)
write.table(rownames(final), file=paste0(idir, "rownames_SNPs_additive_coding_04092019.txt"))

# K=A.mat(final)
# heatmap(K)
# K[1:5,1:5]
# colnames(K[15])

## check where we have high K values 
# k <- arrayInd(which(K>1), dim(K))
# mapply(`[[`, dimnames(K), k) %>% matrix(ncol=2)

