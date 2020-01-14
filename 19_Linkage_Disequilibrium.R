#################################################################
#### Calculate LD in the population + in subpopulations #########
#################################################################

######### need a map first #########

######## below script used for refpop 

library(snpStats)
library(plink)
library(gridExtra)


#genot.mat2=apply(genot_imp, 2 , as.numeric)
# genot.mat2=cbind(data$fam[,2], genot.mat2)
# colnames(genot.mat2)[1]="MUNK"
data$map$chromosome=as.factor(data$map$chromosome)
bed=data$map
map=data.frame(Name=rownames(bed), Chromosome=bed$chromosome, Position=bed$position)

## clean SNPS without position
map=map[which(is.na(map$Chromosome) ==F &  map$Chromosome != 19  ), ]
summary(map, maxsum=20)
map$Chromosome=droplevels(map$Chromosome)

### calculate LD with preferred data and parameters
## script from tutorial: https://www.bioconductor.org/packages/devel/bioc/vignettes/snpStats/inst/doc/ld-vignette.pdf
set.seed(123)
myLD<- function(nmk,map,filepath.plink,rootfilename.plink, nbValLDchr, nbLDval ) {
  print(paste(nmk, "markers per chromosome", nbLDval, "LD values selected in total", sep=" "))
  p <- list()
  ldl <- list()
  snp <- list()
  spectrum=rainbow(10, start=0, end=1/6)[10:1]
  ## calculate LD for given number of mks per chromosome
  for (i in levels(as.factor(map$Chromosome)) ) { 
    chr_map=map
    snpset=sample(chr_map[chr_map$Chromosome==i,1],nmk) ## enter here number of mks wanted
    data_chr <- read.plink(paste(filepath.plink, rootfilename.plink, sep = "/"), select.snps=as.character(snpset) )
    ldl[[i]]=ld(data_chr$genotypes, stats=c("D.prime", "R.squared"), depth=nmk-1)
    snp[[i]]=chr_map[chr_map$Name %in% snpset,]
    D=ldl[[i]]$D.prime
    name=paste("chr", i, sep="_")
    p[[i]] <- image(D,cuts=9, col.regions=spectrum, colorkey=T, lwd=0)
  }
  print("calculation of LD done")
  # do.call(print(grid.arrange),p)
  D.prime <- vector("list", 17)
  R.squared <-vector("list", 17)
  distance <-vector("list", 17)
  ## ldl is a huge matrix; create new vectors for D and R2
  for (i in (1:17)) {
    pos=snp[[i]]$Position ## take positions of SNPs on this chrosomosome
    # print(head(pos))
    # print(head(snp[[i]]))
    diags <- vector("list",nmk-1) ## all neighbours
    for (j in 1:(nmk-1)) diags[[j]] <- pos[(j+1):nmk]-pos[1:(nmk-j)]### calculate distance between markers 2 by 2
    dist <- bandSparse(nmk, k=1:(nmk-1) , diagonals=diags)
    distance[[i]] <- dist@x
    D.prime[[i]] <- ldl[[i]]$D.prime@x
    R.squared[[i]] <- ldl[[i]]$R.squared@x
  }
  print("LD dataframes ready")
  # plot LD per chromsome 
  for (i in 1:17 ) { 
    print(paste0("number of LD values per chromosome ", nbValLDchr))
    dd=as.data.frame(cbind(distance[[i]], R.squared[[i]]))
    nb=dim(dd)[1]
    ddd= as.data.frame(dd[c(sample(nb, nbValLDchr)),])
    colnames(ddd)=c("Distance", "R2")
    distt=ddd$Distance/1000
    jpeg(filename = paste0("LD_chrom",i,".jpeg"))
    plot(ddd$R2~distt,main=paste("Chromosome", i, sep=" "),
         ylab=expression(paste("R"^"2")),
         xlab="Distance (kb)",
         cex=0.5, pch=20)
    loes <- loess(ddd$R2~distt, span=0.5, degree=2, control = loess.control(surface = "direct"))
    dd=sort(as.numeric(distt))
    pre <- predict(loes, dd)
    frame <- data.frame(dd, pre)
    xval <- frame[which.min(abs(0.2 - frame$pre)),] #find x value where y=0.2
    #plot(file$dist, file$rsq, pch=19, cex=0.9, xlab="distance (Mbp)", ylab="LD (r^2)")
    lines(x=frame$dd, y=frame$pre  , col="gold", lwd=2)
    abline(h=0.2, col="red", lwd=2)
    abline(v=xval[1], col="purple", lwd=2)
    #mtext(round(xval[1,1],2), side=1, line=0.05, at=xval[1,1], cex=0.9, col="purple")
    dev.off()
  }
  print("plotting LD per chromosome done")
  ## gather all value to create a single plot
  ld_tot=data.frame(distance=integer(), R2=numeric())
  for (i in 1:17 ) {
    dd=cbind(distance[[i]], R.squared[[i]])
    colnames(dd)=c("distance", "R2")
    ld_tot=rbind(ld_tot, dd)
  }
  par(mfrow=c(1,1))
  sampled=sample(c(1:dim(ld_tot)[1]), nbLDval) ## there is a huge number of values, select some with nbLDval
  sel=ld_tot[sampled,]
  jpeg(filename = paste0(nmk,"_mks", nbLDval, "_LD_values_LD_decay.jpeg"))
  plot(sel$distance, sel$R2)
  ## plot loess, code found here https://jujumaan.com/2017/07/15/linkage-disequilibrium-decay-plot/
  loes <- loess(sel$R2~sel$distance, span=0.5, degree=2, control = loess.control(surface = "direct"))
  ## important to sort distance before plotting otherwise it will create messy lines
  dd=sort(as.numeric(sel$distance))
  pre <- predict(loes, dd)
  frame <- data.frame(dd, pre)
  xval <- frame[which.min(abs(0.2 - frame$pre)),] #find x value closest to y=0.2
  plot(sel$distance, sel$R2)
  lines(x=frame$dd, y=frame$pre  , col="gold", lwd=2)
  abline(h=0.2, col="red", lwd=2)
  abline(v=xval[1], col="purple", lwd=2)
  mtext(round(xval[1,1],2), side=1, line=0.05, at=xval[1,1], cex=0.9, col="purple")
  dev.off()
  print("LD decay plotted")
  ### density plot of LD based on 1000 randomly selected markers
  ## for reproducible results and publication, use set seed
  jpeg(filename = paste0(nmk,"_mks", nbLDval, "_LD_values_density_LD.jpeg"))
  plot(density(ld_tot$R2, na.rm=T))
  dev.off()
  print("LD density plotted")
}

myLD(1000, map,filepath.plink,rootfilename.plink,10000, 10000  )
