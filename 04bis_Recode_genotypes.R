####################################################
########### Recode genotypes of progenies ##########
########### AA, AB, BB to -1 0 1 as in collection ##
####################################################

# Merge all progeny files
# Update SNP names in progenies file and in collection file
# Create a consensus coding for AA, BB, AB coding for parents
# Use this coding to replace progeny genotypes from AA, AB, BB to additive coding
# In case all parents are heterozygous we cannot determine which allele is AA or BB (-1 or 1) >> missing data

## working place
setwd("~/projects/GenSel_Italy/Genos_and_kinships/Recoding_SNPs/")

## packages

library(readxl)
library(magrittr)
library(reshape2)
library(plyr)
## merge collection family genotypes


## load all raw genos

coll<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Raw_data/Genotypes/Training_pop/data_10K_correctMc38_Magr.txt", h=T)
coll_names_origin<-colnames(coll)[-1]
id_names<-coll[,1] %>% as.character
coll<-gsub(":","",as.matrix(coll))
coll<- gsub("NN", NA, coll)
coll[1:5,1:5]
## need to recode, no filtering
genos_101<-raw.data(coll[,-1], frame = "wide", hapmap=NULL, base=T, 
                    call.rate=0, maf=0.0, imput=F, outfile="-101", plot = T)

genos_101$M.clean[1:5,1:5]
coll<-genos_101$M.clean
######


families<-c("FjDe", "FjGD", "FjPi", "FjPL", "GaPi", "GaPL")
path_BP<-"~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Raw_data/Genotypes/Training_pop/Genotypes/Breeding_pop/"
fams<-lapply(families, function(x) assign(x,read_xlsx(paste0(path_BP,"All_SNP_",x, ".xlsx" ),sheet=1) ))
names(fams)<-families
lapply(families, function(x) fams[[x]][1:5,1:5])### in the 2 first columns we have genotypes of the parents

lapply(fams, dim) ## they have the same length, we can concatenate
lapply(families, function(x) summary(fams[[x]][["Name"]] == fams[["FjDe"]][["Name"]])) ## check that all names are matching
progenies<-fams[["FjDe"]][c("Index","Name")]
head(progenies)
for (fam in families) {
  progenies<-merge(progenies, fams[[fam]][-1], by="Name", no.dups=F, suffixes = c("_x","_y"))
}
dim(progenies)
progenies[1:10,1:10]
colnames(progenies)

## recode SNP names and keep overlap with SNPs from collection - used only one family to do that
## change prefixes

prefixes_coll<-unique(substr(colnames(coll),1,3)) ## we have 16 different types of names in collection
prefixes_FjDe<-unique(substr(fams[["FjDe"]][["Name"]],1,3)) 
prefixes_FjDe[prefixes_FjDe %in% prefixes_coll]
from_FjDe<-c("GDsnp", "SNP_FB")
to_FjDe<-c("GD_", "FB")
replacements<-data.frame(from=from_FjDe, to=to_FjDe)
replacements$from<-as.character(replacements$from)
replacements$to<-as.character(replacements$to)
snp_updated<-gsub("GDsnp","GD_", fams[["FjDe"]][["Name"]]) %>% gsub("SNP_FB","FB", .) %>% gsub("snp","", .)## did not manage to find recursive way

prefixes_updated<-unique(substr(snp_updated,1,3)) 
prefixes_updated
RBlong<-grep("RosBREEDSNP_SNP",snp_updated) ## there replace first 19 characters by RB_
RBshort<-grep("RosBREED_SNP",snp_updated) ## there replace first 19 characters by RB_
length(RBlong)
length(RBshort)
## for RB long replace by RB_ and add the character 20+8
for (x in RBlong){
  snp_updated[x]<-paste0("RB_",substr(snp_updated[x],20,27))
}
## for RB short replace by RB_ and add the character 17+8
for (x in RBshort){
  snp_updated[x]<-paste0("RB_",substr(snp_updated[x],17,24))
}
prefixes_updated<-unique(substr(snp_updated,1,3)) 
prefixes_updated

##### now we need to check all suffixes
to_keep_coll<-matrix(c(10,8,11), nrow=3, ncol=1, dimnames=list(c("FB_", "GD_","RB_"), "numb_char"))
to_rep<-lapply(c("FB_", "GD_", "RB_"), function(x) grep(x,colnames(coll)))
names(to_rep)<-c("FB_", "GD_","RB_")

for (x in 1:length(to_rep)) { 
  nn<-names(to_rep)[[x]]
  colnames(coll)[to_rep[[x]]]<-sapply(to_rep[[x]], function(i) substr(colnames(coll)[i] ,1,to_keep_coll[nn,1] ))
}
colnames(coll)[grep("RB_",colnames(coll))] ### number of characters changes in RB but should be the same as in updated markers
sapply(colnames(coll)[grep("RB_",colnames(coll))], nchar) %>% summary
sapply(snp_updated[grep("RB_",snp_updated)], nchar) %>% summary
snp_updated[grep("RB_",snp_updated)]

summary(snp_updated %in% colnames(coll)) ## 9994 markers

remain<-colnames(coll)[which(!(colnames(coll) %in%snp_updated))]
remain
unique(substr(remain,1,3)) 

## check that the snp names are unique!
length(colnames(coll)) == length(unique(colnames(coll))) ## OK
length(snp_updated) == length(unique(snp_updated)) ## 2 markers to check - there are only two with the same position in family mks, not in collection markers
rrr<-as.data.frame(table(snp_updated))
rrr[which(rrr$Freq==2),1] ## "RB_34439805" "RB_34453757" need to be renamed

snp_updated[which(snp_updated %in% as.character(rrr[which(rrr$Freq==2),1]))] 
fams[["FjDe"]][["Name"]][grep("34439805", fams[["FjDe"]][["Name"]])] ## LG3 and LG1
grep("34439805", fams[["FjDe"]][["Name"]]) == grep("34439805", snp_updated) ## sanity check
snp_updated[grep("34439805", snp_updated)]=c("RB_34439805_LG3","RB_34439805_LG1") ## replacement
fams[["FjDe"]][["Name"]][grep("34453757", fams[["FjDe"]][["Name"]])] ## LG1 and LG5
grep("34453757", fams[["FjDe"]][["Name"]]) == grep("34453757", snp_updated) ## sanity check
snp_updated[grep("34453757", snp_updated)]=c("RB_34453757_LG1","RB_34453757_LG5")

## replace also in colnames(coll)
original<-colnames(read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Raw_data/Genotypes/Training_pop/data_10K_correctMc38_Magr.txt", h=T))
original[grep("34439805", original)] ## LG1
original[grep("34453757", original)] ## LG5
colnames(coll)[grep("34439805", colnames(coll))]="RB_34439805_LG1"
colnames(coll)[grep("34453757", colnames(coll))]="RB_34453757_LG5"

fams[["FjDe"]][["Name_new"]]<-snp_updated
class(fams[["FjDe"]])
fams[["FjDe"]]
summary(fams[["FjDe"]][["Name_new"]] %in% colnames(coll))
length(colnames(coll)) ## so almost all markers where found

###############################################################################################
########## finished! replace new SNP names ####################################################
########## filter on both sides to retain overlap families/collection #########################
###############################################################################################

progenies$Name_new<-snp_updated
write.table(progenies[,c("Index", "Name", "Name_new")], file="new_snps_names27062019_progenies.txt",sep="\t",quote=F, row.names=F )
coll_names_file<-cbind(coll_names_origin, colnames(coll))
colnames(coll_names_file)<-c("Name_origin", "Name_new")
head(coll_names_file)
write.table(colnames(coll), file="new_snps_names26072019_collection.txt",sep="\t",quote=F, row.names=F )
merged<-merge(progenies[,c("Index", "Name", "Name_new")],coll_names_file, by="Name_new"  )
head(merged)
dim(merged)
colnames(merged)<-c("Name_new", "Index", "Name_progeny_file", "Name_collection_file")
write.table(merged, file="overlap_coll_progenies.txt",sep="\t",quote=F, row.names=F )

progenies<-progenies[which(progenies$Name_new %in% colnames(coll) ),] ### now down to 10539  markers

###########################################################################################################
#### Now we need to do recode genotypes - we create a column with the value of AA, AB and BB as -1,0,1 ####
###########################################################################################################

## identifiers parents - load key between SNP progenies files and collection file for parent identifiers
codes<-read.table("~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab/Raw_data/Parents_Codes.txt", h=T, sep="\t")
coll[1:5,1:5]
rownames(coll)=id_names
parents_genos<-coll[which(rownames(coll) %in% codes$Code ), which(colnames(coll) %in% snp_updated )] %>% t 
dim(parents_genos)
## add genos from parents in progeny files
parents_genos<-cbind(rownames(parents_genos), parents_genos)
colnames(parents_genos)[1]<-"Name_new"
selec<-lapply(c("Name_new", as.character(unique(codes$Name.in.All_SNP.file))), function(x) grep(x, colnames(progenies))) %>% unlist 
parents_genos<-merge(parents_genos, progenies[,selec], by="Name_new")
head(parents_genos)
colnames(parents_genos)
### noticed that some genotypes from parents have been sequenced twice so we can impute from each other
consens<-parents_genos[,1:7]
for (i in levels(codes$Code)){
  N<-codes[which(codes$Code == i),"Name.in.All_SNP.file"] %>% as.character %>% unique
  NAM<-lapply(N, function(x)colnames(parents_genos)[grep(x, colnames(parents_genos))]) %>% unlist
  print(NAM)
  consens<-cbind(consens, parents_genos[,NAM[1]])
  colnames(consens)[dim(consens)[2]]=NAM[1]
  WhichNA<-which(consens[,NAM[1]] =="NC")
  print(head(consens))
  for (j in NAM[-1]) {
    if (length(WhichNA)>0 & j!=NAM[length(NAM)]) {
      consens[,NAM[1]][WhichNA]<- parents_genos[,j][WhichNA]
      WhichNA<-which(consens[,NAM[1]] =="NC")
      # print(j)
    } else {
      WhichNA<-WhichNA
    }
  }
  
}

lapply(colnames(consens), function(x) summary(as.factor(consens[,x])))
parents_genos<-as.matrix(parents_genos)
paste(parents_genos[,as.character((codes$Name.in.All_SNP.file)[1])],as.character(parents_genos[,as.character((codes$Code)[1])]), "_" ) %>% as.factor %>% summary

## this confirms that AA can be coded as -1 or 1 and same for BB.. So we need to check marker by marker
### now add a column for AA AB and BB
### need to take into account imputed data

consens$AA <- consens$AB <- consens$BB <- NA
head(consens)
colnames(consens)[c(8,12,13)]=c("FxPL_Pink Lady" , "3_78_PINOVA_M9",  "5_65_ROYAL_GALA")
dim(consens)
counter=0
for (i in 1:dim(consens)[1]){
  geno_AB<- consens[i,8:13] %>% as.matrix %>% as.factor  %>% levels
  for (j in geno_AB){
    if (j !="NC") {
      who<-colnames(consens)[which(consens[i,]==j)] ## select parents with this genotype and find the -1 0 1 coding
      to_rep<-consens[i, as.character(codes[which(codes$Name.in.All_SNP.file  %in% who),"Code" ])]  %>% as.matrix
      to_rep<-to_rep[to_rep %in% c(-1,0,1)] %>% unique
      if (length(to_rep) ==1) { ### need to skip ambiguous results
        consens[i,j]<-to_rep
      } else {
        counter=counter + 1 ## check if there is onyl a single value for a given coding, otherwise we have a problem
      }
    }else{
      j=j
    }
  }
}

counter ## 30 genes with ambiguous coding
write.table(consens, file="consensus_file_parents_progenies.txt", sep="\t", quote=F)
### Now we can just keep the recoding
### quality check
lapply(c("AA", "AB", "BB"), function(x) as.factor(consens[,x]) %>% summary) ## for AB we have only zeros, OK
## some ambigous cases whre AA or BB received 0 > look at the other genotype
consens[which(consens$AA ==0 ), "AA"]<- lapply(consens[which(consens$AA ==0 ), "BB"], function(x) if(!is.na(x)) {-(as.numeric(x))} else NA) %>% unlist
consens[which(consens$BB ==0 ), "BB"]<- lapply(consens[which(consens$BB ==0 ), "AA"], function(x) if(!is.na(x)) {-(as.numeric(x))} else NA) %>% unlist
lapply(c("AA", "AB", "BB"), function(x) as.factor(consens[,x]) %>% summary) ## for AB we have only zeros, OK

### file to use and share
myfile<-merge(merged, consens[,c("Name_new", "AA", "AB", "BB")], by="Name_new")
head(myfile)
dim(myfile)
write.table(myfile, file="file_for_recoding_progenies.txt", )


### recode each marker
progenies_add<-progenies 
progenies_add[,-c(1,2,dim(progenies)[2])]<-NA


for (i in progenies$Name_new){
  ## take each line (marker) and replace AA, BB, AB by corresponding coding (if knows)
  myline<-progenies[which(progenies$Name_new==i),-c(1,2,dim(progenies)[2])] %>% as.character()
  to_rep<-myfile[which(myfile$Name_new==i), c("AA","AB","BB")]
  # print(to_rep)
  myline_new<-mapvalues(myline,from= c("AA", "AB", "BB", "NC"), to= c(to_rep,NA)) %>% unlist
  progenies_add[which(progenies_add$Name_new==i),-c(1,2,dim(progenies)[2])]<-myline_new
}

progenies_add[1:5,1:15]
progenies[1:5,1:15]
#### summary on individuals
prg_names<-colnames(progenies_add)[-c(1,2,dim(progenies))]
res<-lapply(prg_names, function(x)summary(as.factor(progenies_add[,x])))
res<-do.call(rbind,res) %>% as.data.frame()
rownames(res)<-prg_names
res$tot_avail<-do.call(function(x) res[x,"0"]+res[x,"-1"]+res[x,"1"],list(prg_names))
head(res)
write.table(res, file="summary_genotypes_progenies.txt", sep="\t", quote=F)  

#### summary on markers
mk_names<-progenies_add[,c("Name", "Name_new")]
res<-lapply(1:length(mk_names$Name_new), function(x)summary(as.factor(progenies_add[x,-c(1,2,dim(progenies_add)[2])])))
res<-do.call(rbind,res) %>% as.data.frame()
rownames(res)<-mk_names$Name_new
res$tot_avail<-do.call(function(x) res[x,"0"]+res[x,"-1"]+res[x,"1"],list(mk_names$Name_new))
head(res)

write.table(res, file="summary_markers_progenies.txt", sep="\t", quote=F)  

write.table(progenies_add, file="progenies_coding_additive.txt", sep="\t", quote=F)
progenies_add<-read.table("progenies_coding_additive.txt", sep="\t", h=T)

### concatenate with collection
tcoll<-t(coll) %>% as.data.frame()
tcoll$Name_new=rownames(tcoll)
dim(tcoll)
tcoll[1:5,230:234]
progenies_add<-read.table("progenies_coding_additive.txt", sep="\t", h=T)
progenies_add<-progenies_add[,c("Name_new", colnames(progenies_add)[-c(1,2,dim(progenies_add)[2])])]
progenies_add[1:5,1:5]
all<-merge(tcoll, progenies_add, by="Name_new") 

## Last step is to remove outcrossers and additional columns for parents

# remove outcrossers
out<-read_xlsx("~/mnt/agroscope_os/2/2/6/1/1/4/1055/GenSel_Costa_collab/Raw_data/Genotypes/Breeding_pop/outcross.xlsx", sheet=1,col_names=T)
out$ID_modified_Morgane
to_remove=which(colnames(all) %in% out$ID_modified_Morgane ) ## 7 invidividuals
all<-all[,-c(to_remove)]
dim(all)
class(all)
all[1:5,1:5] %>% as.matrix %>% t
rownames(all)<-all[,1]
all<-all[-1,] %>% t

# remove parents in copy
whichPA<-sapply(codes$Name.in.All_SNP.file,function(x) grep(x, rownames(all))) %>% unlist 
rownames(all)[whichPA]
all<-all[-whichPA, ]
rownames(all)
dim(all)

write.table(all, "all_genotypes_coll_progenies_additive.txt", sep="\t", row.names=T, quote=F)

## check that the file is properly saved
dd<-read.table("./all_genotypes_coll_progenies_additive.txt", h=T, sep="\t")
dd[1:5,1:5]
