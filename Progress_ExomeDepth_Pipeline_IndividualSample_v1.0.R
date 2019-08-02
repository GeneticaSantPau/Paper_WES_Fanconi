#!Rscript
rm(list=ls())
# Libraries
library(GenomicRanges)
library(IRanges)
library(ExomeDepth)
library(stringr)
#library(XLConnect)
library(png)
library(reshape)
library(bedr)
library(gridExtra)
library(xlsx)
########################### ###########################
# ARGUMENTS from pipeline ###############################
########################### ###########################
args <- commandArgs(TRUE)
ID<- args[1]  # Sample name/ID
data<- args[2]   # Day of thea analysis

outDir<- args[3] # Output folder
captura<-args[4] # Capture ID (Exome-Agilent, Exome-Roche ...)
regionName<-args[5] # Label of geneList to analyse
label<-args[6] # Label for used capture
bedFile<-args[7] # Target bed file of the used capture
candListCheck<-args[8] # File with gene IDs to obtain spcific results from capture panel, one gene per line
######################################################
# DATA needed
######################################################
ExomeCountDir<- paste("/media/gensp/DATA2/BAMs/",captura,sep="")
bamDir<-ExomeCountDir

# Complementary data
fasta<-c("~/NGSpipeline/databases/references/b37/human_g1k_v37.fasta")
exons.hg19<-read.delim("~/databases/UCSC/RefSeqGenes_RefFlat_CodingExonsHg19_31Gener2018_bedtools_joined_RefSeqGenes_RefFlat_RefGeneNMs_AllExonsHg19_31Gener2018_sorted_simplified_FINAL2.merged.bed", head=F)
names(exons.hg19)<-c("chromosome", "start", "end", "name", "strand", "gene","NM","exon")
exons.hg19$chromosome<-gsub("chr", "",exons.hg19$chromosome)
exons.hg19.GRanges <- GRanges(seqnames = exons.hg19$chromosome, IRanges(start=exons.hg19$start, end=exons.hg19$end), names = exons.hg19$name)
# Exon unique intervals:
N.intervals.hg19<-read.delim("~/databases/UCSC/RefSeqGenes_RefFlat_CodingExonsHg19_31Gener2018_bedtools_joined_RefSeqGenes_RefFlat_RefGeneNMs_AllExonsHg19_31Gener2018_sorted_simplified_FINAL2_UniqueIntervalsWithNumber1.bed", head=F)
names(N.intervals.hg19)<-c("chromosome", "start", "end", "name")
N.intervals.hg19$chromosome<-gsub("chr", "",N.intervals.hg19$chromosome)
N.intervals.hg19.GRanges <- GRanges(seqnames = N.intervals.hg19$chromosome, IRanges(start=N.intervals.hg19$start, end=N.intervals.hg19$end), names = N.intervals.hg19$name)

omimData<-read.delim("~/databases/OMIM/31012018/genemap2_hg19_hg20_phenotypes_senseNoConvertits_WithInheritances_Version_31012018.bed")
names(omimData)[1:3]<-c("chr","start", "end")
omimData$chr<-gsub("chr", "", omimData$chr)
omimData2<-omimData[!is.na(omimData$chr),]
omimData.GRanges<-GRanges(seqnames = omimData2$chr, IRanges(start=omimData2$start, end=omimData2$end), names = omimData2$Phenotypes)

# Whole Gene RefSeq # NEW 4 Maig 2017
wholeGenes<-read.delim("~/databases/UCSC/RefSeqGenes_RefFlat_WholeGeneHg19_1Febrer2018.bed", head=F)
names(wholeGenes)<-c("chr","start","end","gene","kk","strand","V7","V8","V9","V10","V11","V12")
wholeGenes$chr<-gsub("chr", "", wholeGenes$chr)
wholeGenes.GRanges<-GRanges(seqnames = wholeGenes$chr, IRanges(start=wholeGenes$start, end=wholeGenes$end), names = wholeGenes$gene)

# Segmental Dups
segdup<-read.delim("~/databases/UCSC/genomicSuperDups_hg19_13Feb2018.bed",head=F)
names(segdup)[1:3]<-c("chr","start", "end")
segdup$chr<-gsub("chr", "", segdup$chr)
segdup.GRanges<-GRanges(seqnames = segdup$chr, IRanges(start=segdup$start, end=segdup$end), names = segdup$V4)

# DGV
dgv<-read.delim("~/databases/UCSC/dgvMerged_hg19_13Feb2018.bed",head=T)
names(dgv)[1:3]<-c("chr","start", "end")
dgv$chr<-gsub("chr", "", dgv$chr)
dgv.GRanges<-GRanges(seqnames = dgv$chr, IRanges(start=dgv$start, end=dgv$end), names = dgv[,c(4:16)])
dgv.GRangesLosses<-GRanges(seqnames = dgv$chr, IRanges(start=dgv$start, end=dgv$end), names = dgv[,14])
dgv.GRangesGains<-GRanges(seqnames = dgv$chr, IRanges(start=dgv$start, end=dgv$end), names = dgv[,13])

# Intervals panel NGS utilitzat.
capture<-read.delim(bedFile,head=F)
names(capture)[1:3]<-c("chromosome", "start", "end")
capture$chromosome<-gsub("chr", "", capture$chromosome)
captureProbes<-capture[,c(1:3)]
capture.GRanges <- GRanges(seqnames = capture$chromosome, IRanges(start=capture$start,end=capture$end))
########################### ###########################
# EN DATA needed
########################### ###########################
setwd(bamDir)
my.bams<-dir(bamDir)[grep(paste(ID,".final.bam$",sep=""), dir(bamDir))] 

####################### BAM COUNTS #######################
####### Recovering data from last existing bamCount
details <- file.info(list.files(pattern="FromExonsExomeCount.dafr*"))
details <- details[with(details, order(as.POSIXct(mtime))), ]
ExomeCount.daframes <- rownames(details)
lastExomeCount<-ExomeCount.daframes[length(ExomeCount.daframes)]
print(paste("LastExomeCount: ", lastExomeCount, sep=""))
ExomeCount.dafr<-read.delim(lastExomeCount)

######## Introducing new data (if any) from bams: ####################### 
#####################################################
for(a in unique(my.bams)){
  a2<-gsub("-","\\.",paste("X",a,sep=""))
  a3<-gsub("-","\\.",a)
  toMatch<-paste(a,a2,a3,sep="|")
  if(length(names(ExomeCount.dafr)[grep(toMatch, names(ExomeCount.dafr))]) < 1 ) {
    print(paste("Adding new ", a, " sample", sep=""))
    my.new.count<-getBamCounts(bed.frame = captureProbes , bam.file = paste(bamDir,  "/",a,  sep="") , include.chr = FALSE, referenceFasta = fasta) 
    test<-as(my.new.count[, colnames(my.new.count)], 'data.frame')
    #names(test)[6]<-a
    ExomeCount.dafr2<-cbind(ExomeCount.dafr, test[names(test)==a | names(test)==a2 ])
    ExomeCount.dafr<-ExomeCount.dafr2
  } else {
    print(paste("Sample ", a, " already added", sep=""))
  }
}

date<-gsub("-","_",Sys.Date())
#write.table(ExomeCount.dafr, paste(bamDir,"/FromExonsExomeCount.dafr.", captura, ".", date, ".txt", sep=""), quote=F, sep="\t", row.names=F)
####################### END BAM COUNTS #######################

##########################################
#### Sample ANALYSIS
##########################################
setwd(outDir)
mostra<-ID
sampleName<-gsub("-",".",mostra)
toMatch<-paste(mostra,sampleName,sep="|")
print(paste("ANALYSING SAMPLE ", mostra, sep=""))

my.test <- ExomeCount.dafr[,grep(toMatch, names(ExomeCount.dafr))]
#########################################################
# Rest of samples as reference
my.ref.samples <- gsub("-", ".", names(ExomeCount.dafr)[grep(toMatch,names(ExomeCount.dafr), invert=TRUE)])[grep("bam$", gsub("-", ".", names(ExomeCount.dafr)[grep(toMatch,names(ExomeCount.dafr), invert=TRUE)]) )]
#my.ref.samples <- gsub("-", ".", my.bams[grep(bamName,my.bams, invert=TRUE)])
my.gender<-c("ContraTotes")
print(paste("Mostres Referència candidates:", sep=""));print(gsub(".bam|.final.bam","",my.ref.samples))

if(length(my.ref.samples)>1){
  my.reference.set <- as.matrix(ExomeCount.dafr[, grep(paste(my.ref.samples, collapse="|"), names(ExomeCount.dafr))])
  my.choice <- select.reference.set (test.counts = my.test,  reference.counts = my.reference.set,  bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000, n.bins.reduced = 10000)
  #print(my.choice[[1]])
  lengthChoice<-length(my.choice[[1]])
    
  # Using the output of this procedure we can construct the reference set:
  if(lengthChoice==1){
  my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE]) #Note that the drop = FALSE option is just used in case the reference set contains a single sample. If this is the case, it makes sure that the subsetted object is a data frame, not a numeric vector.
  print("Reference set només té una mostra")
  }else{
  my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice]) 
  print("Reference set té més d'una mostra")
  }
} else {
  print(paste("Només una altre mostra a la carpeta", sep=""))
  my.ref.sample<-my.ref.samples
  my.reference.set <- as.matrix(ExomeCount.dafr[, grep(paste(my.ref.samples, collapse="|"), names(ExomeCount.dafr))])
    colnames(my.reference.set)<-my.ref.samples
    
  my.choice <- select.reference.set (test.counts = my.test,  reference.counts = my.reference.set,  bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000, n.bins.reduced = 10000)
    #print(paste("Mostres Referència SELECCIONADES: ", my.choice[[1]], sep=""))
  lengthChoice<-length(my.choice[[1]])
  my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE]) #Note that the drop = FALSE option is just used in case the reference set contains a single sample. If this is the case, it mak
  }
  print(paste("Mostres Referència SELECCIONADES:", sep=""));print(my.choice[[1]])
  my.reference.selected <- apply(X = my.matrix,  MAR = 1,  FUN = sum)  

  ################################################################################
  ############# CNV calling ###########################################################
  ################################################################################
  all.exons <- new('ExomeDepth',  test = my.test,  reference = my.reference.selected,  formula = 'cbind(test, reference) ~ 1')
  # We can now call the CNV by running the underlying hidden Markov model:
  #Transition probability of the hidden Markov Chain from the normal copy number state to either a deletion or a duplication. The default (0.0001) expect approximately 20 CNVs genome-wide
  #all.exons <- CallCNVs(x = all.exons,  transition.probability = 10^-4,  chromosome = ExomeCount.dafr$space,  start = ExomeCount.dafr$start,  end = ExomeCount.dafr$end, name = ExomeCount.dafr$names)
  all.exons <- CallCNVs(x = all.exons,  transition.probability = 10^-3,  chromosome = ExomeCount.dafr$space,  start = ExomeCount.dafr$start,  end = ExomeCount.dafr$end, name = ExomeCount.dafr$names)
  #print(all.exons@CNV.calls)
  if(length(all.exons@CNV.calls$start)!=0){
    #all.exons <- AnnotateExtra(x = all.exons, reference.annotation = Conrad.hg19.common.CNVs, min.overlap = 0.5, column.name = 'Conrad.hg19')
    all.exons <- AnnotateExtra(x = all.exons, reference.annotation = wholeGenes.GRanges, min.overlap = 0.0001, column.name = 'Genes')
    all.exons <- AnnotateExtra(x = all.exons, reference.annotation = exons.hg19.GRanges, min.overlap = 0.0001, column.name = 'exons.hg19')
    all.exons <- AnnotateExtra(x = all.exons, reference.annotation = omimData.GRanges, min.overlap = 0.0001, column.name = 'RegionOMIM')
    all.exons <- AnnotateExtra(x = all.exons, reference.annotation = segdup.GRanges, min.overlap = 0.5, column.name = 'SegmDup')
    #all.exons <- AnnotateExtra(x = all.exons, reference.annotation = dgv.GRanges, min.overlap = 0.5, column.name = 'DGV')
    all.exons <- AnnotateExtra(x = all.exons, reference.annotation = dgv.GRangesLosses, min.overlap = 0.5, column.name = 'DGV_Losses')
    all.exons <- AnnotateExtra(x = all.exons, reference.annotation = dgv.GRangesGains, min.overlap = 0.5, column.name = 'DGV_Gains')
    all.exons <- AnnotateExtra(x = all.exons, reference.annotation = N.intervals.hg19.GRanges, min.overlap = 0.0001, column.name = 'N.intervals') 
    #all.exons <- AnnotateExtra(x = all.exons, reference.annotation = dgvGains.GRanges, min.overlap = 0.3, column.name = 'dgvGains.0.3')
    #all.exons <- AnnotateExtra(x = all.exons, reference.annotation = dgvLoss.GRanges, min.overlap = 0.3, column.name = 'dgvLoss.0.3')
    #all.exons <- AnnotateExtra(x = all.exons, reference.annotation = iscaPathogenicGain.GRanges, min.overlap = 0.3, column.name = 'iscaPathGain.0.3')
    #all.exons <- AnnotateExtra(x = all.exons, reference.annotation = iscaPathogenicLoss.GRanges, min.overlap = 0.3, column.name = 'iscaPathLoss.0.3')
    #all.exons <- AnnotateExtra(x = all.exons, reference.annotation = iscaLikelyPathogenic.GRanges, min.overlap = 0.3, column.name = 'iscaLikelyPath.0.3')
    all.exons@CNV.calls$size<-all.exons@CNV.calls$end-all.exons@CNV.calls$start
    outAns<-all.exons@CNV.calls[,c(7,5,6,20,3,9,13,19,4,14:18,10:12)]
    
    names(outAns)[grep("chromosome", names(outAns))]<-c("chr")
    names(outAns)[grep("nexons", names(outAns))]<-c("N.exons")
    outAns$Genes<-sapply(outAns$Genes, function(x) paste(unique(unlist(strsplit(x, ","))), collapse = ","))
    #outAns$exons.hg19<-sapply(outAns$exons.hg19, function(x) paste(unique(unlist(strsplit(x, ","))), collapse = ","))
    outAns$exons.hg19<-sapply(outAns$exons.hg19, function(x) gsub("\\ ","",toString(gsub("\\(|\\)","",unlist(str_split(unique(unlist(strsplit(x,","))),"\\ ",2))[grep("NM_", unlist(str_split(unique(unlist(strsplit(x,","))),"\\ ",2)))]),collapse=",")))
    outAns$N.exons<-sapply(gregexpr(",", outAns$exons.hg19, fixed = TRUE), function(x) sum(x > -1))
    outAns$N.intervals<-sapply(outAns$N.intervals, function(x) sum(as.integer(strsplit(x,",")[[1]])))
    outAns$Genes<-ifelse(outAns$N.exons==0, c("notCoding?"), outAns$Genes)
    outAns$N.intervals<-ifelse(outAns$N.exons==0, c("check"), outAns$N.intervals)
    outAns$exons.hg19<-ifelse(outAns$N.exons==0, c("check"), outAns$exons.hg19)
    
    date<-gsub("-","_",Sys.Date())
    output.file <- paste(outDir,"/",sampleName,".",captura,".CNVs_ExomeDepth.", my.gender,  ".",date, ".txt",sep="")
    write.table(file = output.file,  x = outAns, row.names = FALSE, quote=F, sep="\t")
  } else {
    output<-c("CNVs no trobades per ExomeDepth")
    output.file <- paste(outDir,"/",sampleName,".",captura,".ExomeDepth.SENSE_CNVS.", my.gender, ".",date, ".txt",sep="")
    write.table(file = output.file,  x = all.exons@CNV.calls, row.names = FALSE, quote=F, sep="\t")
}
#######################################################################################
## CNVs-> PLOTS IF there are candidate CNVs
#######################################################################################
if(candListCheck!="SenseCandList" & length(all.exons@CNV.calls$start)!=0){
  candList<-read.delim(candListCheck)
  names(candList)<-c("Gen")
  
  if((regionName==captura) & (captura==label)){  # If ID to analyse is the whole capture, then subset is every gene in the capture
  candCnvs<-outAns
  } else {  # If ID to analyse is NOT the whole capture, then subset genes of interest data
  toMatch<-paste(as.character(candList$Gen),collapse="|")
  candCnvs<-outAns[grep(toMatch, outAns$Genes),]
    
  }
  print(paste("Núm. CNVs en gens candidats ", label, " :", dim(candCnvs)[1], sep=""))
}
  
if(length(candCnvs$chr>0)){
  date<-gsub("-","_",Sys.Date())
  output.file <- paste(outDir,"/",sampleName,".",captura,"_GensCandidats_CNVs_ExomeDepth.", my.gender,  ".",date, ".txt",sep="")
  write.table(file = output.file,  x = candCnvs, row.names = FALSE, quote=F, sep="\t")
  
  outAns<-candCnvs
# PLOT CONSISTENT CNVs, Top 20 based on BF or other conditions:
plots<-outAns[outAns$BF>8,]
if (dim(plots)[1]!=0){
  plots$exonsOK<-c("-")
  for (i in 1:nrow(plots)){
    qq<-as.data.frame(do.call(rbind, strsplit(unlist(strsplit(plots$exons.hg19[i],",")),":")))
    outputExonsOK<-NULL
    if(dim(qq)[1]==0){
      outputExonsOK<-0
    } else {
      for(NM in unique(qq[,1])){
        if(NM=="check"){
        outputExonsOK<-c("check")
        } else {
          #print(NM)
          outLineExons<-paste(qq[qq[,1]==NM,c(2:dim(qq)[2])] , collapse=",")
          outLine<-paste(NM,":",outLineExons,sep="")
          outputExonsOK<-paste(outputExonsOK,outLine,sep="/")
      }
    }
}
      plots$exonsOK[i]<-ifelse(plots$N.exons[i]>6, c("Multiple exons") , gsub("^/","",outputExonsOK))
      nameFile<-gsub("\\?","_check", paste(sampleName,"_CNV_",unlist(strsplit(plots$Genes[i],","))[1],"_chr", plots$chr[i],"_",plots$start[i],".png", sep=""))
      png(paste(outDir,"/",nameFile, sep=""), width = 800, height = 600)
      
      plot (all.exons, sequence = plots$chr[i], xlim = c(plots$start[i] - 10000, plots$end[i] + 10000), count.threshold = 20,  xlab = 'Hg19 coordinates', main = paste(plots$Genes[i], " (", plots$type[i], " ~",plots$end[i]-plots$start[i], "-pb)\nExons:", plots$exonsOK[i], "\n", plots$N.intervals[i], " interval(s) BF=", plots$BF[i], " Ratio=", plots$reads.ratio[i], sep="") , cex.lab = 0.8)
      dev.off()
    }
  } else {
    print(paste(mostra, " : NO HI HA PLOTS PER CNVS GENS CANDIDATS BF>8 per fer PLOT", sep=""))
  }  
  } else {
    print(paste(mostra, " : NO HI HA CNVS GENS CANDIDATS BF>8 per fer PLOT", sep=""))
    date<-gsub("-","_",Sys.Date())
    output.file <- paste(outDir,"/",sampleName,".",captura,"_GensCandidats_ExomeDepth.SENSE_CNVS.", my.gender, ".",date, ".txt",sep="")
    write.table(file = output.file,  x = candCnvs, row.names = FALSE, quote=F, sep="\t")
    
  } 



  

