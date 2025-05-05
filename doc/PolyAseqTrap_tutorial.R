## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6,
  fig.height = 5.5,
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ----message=FALSE, results='hide'--------------------------------------------
library(PolyAseqTrap)
## identify and trim the polyT stretches
file_T <- system.file("extdata", "SRR1168402_T.fastq", package = "PolyAseqTrap")
# The output is 'SRR11837378_A.A.fq', which can be used as input for alignment tools
findTailAT(infile=file_T, odir=NULL,
         poly='T', ml=20, mp=5, mg=10, mm=2, deep=FALSE,
         mtail=6, mper=0.75, mr=3, review=TRUE, debug=TRUE,
          bar=0, reg=1, suf=NULL)


## identify and trim the polyA stretches
file_A <- system.file("extdata", "SRR11837378_A.fastq", package = "PolyAseqTrap")
# The output is 'SRR11837378_A.A.fq', which can be used as input for alignment tools
findTailAT(infile=file_A, odir=NULL, poly='A', 
          ml=20, mp=5, mg=10, mm=2, 
          deep=0, mtail=6, mper=0.75, 
          mr=3, review=TRUE, debug=TRUE, bar=0, reg=1, suf=NULL)
# Example of sequence headers in the FASTQ file
##@SRR1168402.78_TTTCTTTTTTTTTTTT
##@SRR11837378.19_TAACAAATT_AAAAAAAAAAAAA

## ----message=FALSE------------------------------------------------------------
# for Homo sapiens (UCSC genome hg38)
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
bsgenome <-BSgenome.Hsapiens.UCSC.hg38

## ----message=FALSE------------------------------------------------------------
# for Mus musculus (UCSC genome mm39)
library("BSgenome.Mmusculus.UCSC.mm39", quietly = TRUE)
bsgenome <-  BSgenome.Mmusculus.UCSC.mm39

## ----message=FALSE------------------------------------------------------------
# for Arabidopsis thaliana (Ensembl TAIR110)
library("BSgenome.Athaliana.ENSEMBL.TAIR10", quietly = TRUE)
bsgenome <- BSgenome.Athaliana.ENSEMBL.TAIR10

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  install.packages("devtools")
#  require(devtools)
#  install_github("BMILAB/movAPA")
#  library(movAPA)
#  browseVignettes('movAPA')

## ----bsgenome_hg, message=FALSE,eval=FALSE------------------------------------
#  library(movAPA)
#  # for Homo sapiens
#  # download Homo sapiens (hg38) from Ensembl
#  athGFF <- "Homo_sapiens.GRCh38.110.gtf"
#  gff<- parseGff(athGFF)
#  saveRDS(gff,file="Ensembl_Homo_sapiens.GRCh38.110.Rdata")
#  
#  threeUTR.data <-subset(gff$anno.need,type=="three_prime_UTR")
#  threeUTRregion <- makeGRangesFromDataFrame(threeUTR.data,keep.extra.columns = F)
#  saveRDS(threeUTRregion,file="ThreeRegion_Homo_sapiens.Rdata")

## ----bsgenome_mouse, message=FALSE,eval=FALSE---------------------------------
#  library(movAPA)
#  # for Mus musculus
#  # download Mus musculus (mm39) from Ensembl
#  athGFF <- "Mus_musculus.GRCm39.110.gtf"
#  gff<- parseGff(athGFF)
#  saveRDS(gff,file="Ensembl_Mus_musculus.GRCm39.110.Rdata")
#  threeUTR.data <-subset(gff$anno.need,type=="three_prime_UTR")
#  threeUTRregion <- makeGRangesFromDataFrame(threeUTR.data,keep.extra.columns = F)
#  saveRDS(threeUTRregion,file="ThreeRegion_Mus_musculus.Rdata")

## ----bsgenome_tair, message=FALSE,eval=FALSE----------------------------------
#  library(movAPA)
#  # for Arabidopsis thaliana
#  # download Arabidopsis (TAIR10) from Ensembl Plant
#  athGFF <- "Arabidopsis_thaliana.TAIR10.57.gff3"
#  gff<- parseGff(athGFF)
#  saveRDS(gff,file="Ensembl_Arabidopsis_thaliana.TAIR10.57.Rdata")
#  threeUTRregion <- makeGRangesFromDataFrame(threeUTR.data,keep.extra.columns = F)
#  saveRDS(threeUTRregion,file="ThreeRegion_Arabidopsis_thaliana.Rdata")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(PolyAseqTrap,  warn.conflicts = FALSE, quietly=TRUE)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <-  BSgenome.Hsapiens.UCSC.hg38

# load 3'UTR annotation for detecting V8 polyA site
threeUTR_path <- system.file("extdata",
                             "ThreeRegion_Homo_sapiens.Rdata",
                             package = "PolyAseqTrap")
threeUTRregion <- readRDS(threeUTR_path)

# get bam file
bam_T_file <- system.file("extdata",
                          "SRR299116_T_chr22_hg_sorted.bam",
                          package = "PolyAseqTrap")


# identify and quantify PACs, it wouldn't predict V8 polyA site if 
# without providing 3'UTR annotation.
# here "adjust.chr" is set to TRUE to add "chr" prefix
pa.hg.result <- FindPTA(bam=bam_T_file, 
        yieldSize=10^7,
        reverse=F,
        bsgenome=bsgenome,
        d=24,
        poly='T',
        adjust.chr=TRUE,
        threeUTRregion=threeUTRregion,
        cutoffCount = 5,
        ext3UTRlen =   1000 ,
        isDRS = FALSE,
        run.quantify=TRUE,
        detect.repeat=TRUE)
# Display details of alignment and category of aligned reads
rmarkdown::paged_table(head(pa.hg.result$pa.table[,c("readName","cigar","seq",
                                                     "softClipFragment","trimmed_seq",
                                                     "unmapped_seq",
                                                     "reference_seq","is_Arich",
                                                     "chr","strand","coord",
                                                     "level","class","use.as.count")]),
                       options = list(rows.print = 5, cols.print = 5))
#category of aligned reads
t(table(pa.hg.result$pa.table$class))
#subclasses of aligned reads
t(table(pa.hg.result$pa.table$level))

# Display details of PACs
rmarkdown::paged_table(head(pa.hg.result$pa.coord),
                       options = list(rows.print = 5, cols.print = 5)) 
#filter PACs that were supported by at least five reads
pac5.hg <- subset(pa.hg.result$pa.coord,total.count>=5)
# Filter PACs whose representative polyA sites are not located in A-rich repeat regions.For polyA sites located in intronic regions, it is recommended to filter out those in repeat regions.
pac5.hg.clean <- subset(pac5.hg, !repeat_flag)

## ----message=FALSE, eval=FALSE------------------------------------------------
#  library(PolyAseqTrap, warn.conflicts = FALSE, quietly=TRUE)
#  library(BSgenome.Hsapiens.UCSC.hg38)
#  bsgenome <-  BSgenome.Hsapiens.UCSC.hg38
#  
#  # load 3'UTR annotation for detecting V8 polyA site
#  threeUTR_path <- system.file("extdata",
#                               "ThreeRegion_Homo_sapiens.Rdata",
#                               package = "PolyAseqTrap")
#  threeUTRregion <- readRDS(threeUTR_path)
#  
#  # get bam file
#  bam_A_file <- system.file("extdata",
#                            "SRR11837378_A_chr22_hg_sorted.bam",
#                            package = "PolyAseqTrap")
#  
#  # identify and quantify PACs, it wouldn't predict V8 polyA site if
#  #without providing 3'UTR annotation.
#  pa.hg.result <- FindPTA(bam=bam_A_file,
#                          yieldSize=10^7,
#                          reverse=F,
#                          bsgenome=bsgenome,
#                          d=24,
#                          poly='A',
#                          adjust.chr=TRUE,
#                          threeUTRregion=threeUTRregion,
#                          cutoffCount = 5,
#                          ext3UTRlen =   1000 ,
#                          isDRS = FALSE,
#                          run.quantify=TRUE,
#                          detect.repeat=TRUE)
#  # Display details of alignment and category of aligned reads
#  rmarkdown::paged_table(head(pa.hg.result$pa.table[,c("readName","cigar","seq",
#                                                               "softClipFragment","trimmed_seq",
#                                                       "unmapped_seq",
#                                                       "reference_seq","is_Arich",
#                                                       "chr","strand","coord",
#                                                       "level","class","use.as.count")]),
#                         options = list(rows.print = 5, cols.print = 5))
#  #category of aligned reads
#  knitr::kable(t(table(pa.hg.result$pa.table$class)))
#  #subclasses of aligned reads
#  knitr::kable(t(table(pa.hg.result$pa.table$level)))
#  
#  # Display details of PACs
#  rmarkdown::paged_table(head(pa.hg.result$pa.coord),
#                         options = list(rows.print = 5, cols.print = 5))
#  #filter PACs that were supported by at least five reads
#  pac5.hg <- subset(pa.hg.result$pa.coord,total.count>=5)
#  # Filter PACs whose representative polyA sites are not located in A-rich repeat regions.For polyA sites located in intronic regions, it is recommended to filter out those in repeat regions.
#  pac5.hg.clean <- subset(pac5.hg, !repeat_flag)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(PolyAseqTrap,  warn.conflicts = FALSE, quietly=TRUE)
library(BSgenome.Mmusculus.UCSC.mm39)
bsgenome <-  BSgenome.Mmusculus.UCSC.mm39
# load 3'UTR annotation for detecting V8 polyA site
threeUTR_path <- system.file("extdata",
                             "ThreeRegion_Mus_musculus.Rdata", 
                             package = "PolyAseqTrap")
threeUTRregion <- readRDS(threeUTR_path)

# get bam file
bam_T_file <- system.file("extdata", 
                          "SRR766743_T_chr19_mm_sorted.bam", 
                          package = "PolyAseqTrap")
# identify and quantify PACs, it wouldn't predict V8 polyA site if 
#without providing 3'UTR annotation
pa.mm.result <- FindPTA(bam=bam_T_file, 
                        yieldSize=10^7,
                        reverse=F,
                        bsgenome=bsgenome,
                        d=24,
                        poly='T',
                        adjust.chr=TRUE,
                        threeUTRregion=threeUTRregion,
                        cutoffCount = 5,
                        ext3UTRlen =   1000 ,
                        isDRS = FALSE,
                        run.quantify=TRUE,
                        detect.repeat=TRUE)


# Display details of alignment and category of aligned reads
rmarkdown::paged_table(head(pa.mm.result$pa.table[,c("readName","cigar","seq",
                                                     "softClipFragment","trimmed_seq",
                                                     "unmapped_seq",
                                                     "reference_seq","is_Arich",
                                                     "chr","strand","coord",
                                                     "level","class","use.as.count")]),
                       options = list(rows.print = 5, cols.print = 5))
#category of aligned reads
t(table(pa.mm.result$pa.table$class))
#subclasses of aligned reads
t(table(pa.mm.result$pa.table$level))

# Display details of PACs
rmarkdown::paged_table(head(pa.mm.result$pa.coord),
                       options = list(rows.print = 5, cols.print = 5)) 
#filter PACs that were supported by at least five reads
pac5.mm <- subset(pa.mm.result$pa.coord,total.count>=5)
# Filter PACs whose representative polyA sites are not located in A-rich repeat regions.  For polyA sites located in intronic regions, it is recommended to filter out those in repeat regions.
pac5.mm.clean <- subset(pac5.mm, !repeat_flag)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(PolyAseqTrap, warn.conflicts = FALSE, quietly=TRUE)
library(BSgenome.Athaliana.ENSEMBL.TAIR10)
bsgenome <-  BSgenome.Athaliana.ENSEMBL.TAIR10

# load 3'UTR annotation for detecting V8 polyA site
threeUTR_path <- system.file("extdata", 
                             "ThreeRegion_Arabidopsis_thaliana.Rdata", 
                             package = "PolyAseqTrap")
threeUTRregion <- readRDS(threeUTR_path)

# get bam file
bam_T_file <- system.file("extdata", 
                          "SRR5055884_T_ch2_tair_sorted.bam", 
                          package = "PolyAseqTrap")
# identify and quantify PACs, it wouldn't predict V8 polyA site if 
# without providing 3'UTR annotation
pa.tair.result <- FindPTA(bam=bam_T_file, 
                          yieldSize=10^7,
                          reverse=F,
                          bsgenome=bsgenome,
                          d=24,
                          poly='T',
                          adjust.chr=FALSE,
                          threeUTRregion=threeUTRregion,
                          cutoffCount = 5,
                          ext3UTRlen =   1000 ,
                          isDRS = FALSE,
                          run.quantify=TRUE,
                          detect.repeat=TRUE)

# Display details of alignment and category of aligned reads
rmarkdown::paged_table(head(pa.tair.result$pa.table[,c("readName","cigar","seq",
                   "softClipFragment","trimmed_seq",
                                                       "unmapped_seq",
                                                       "reference_seq","is_Arich",
                                                       "chr","strand","coord",
                                                       "level","class","use.as.count")]),
                       options = list(rows.print = 5, cols.print = 5))
#category of aligned reads
t(table(pa.tair.result$pa.table$class))
#subclasses of aligned reads
t(table(pa.tair.result$pa.table$level))

# Display details of PACs
rmarkdown::paged_table(head(pa.tair.result$pa.coord),
                       options = list(rows.print = 5, cols.print = 5)) 
#filter PACs that were supported by at least five reads
pac5.tair<- subset(pa.tair.result$pa.coord,total.count>=5)
# Filter PACs whose representative polyA sites are not located in A-rich repeat regions.For polyA sites located in intronic regions, it is recommended to filter out those in repeat regions.
pac5.tair.clean <- subset(pac5.tair, !repeat_flag)

## ----message=FALSE, eval=FALSE------------------------------------------------
#  library(reticulate)
#  use_condaenv("/path/miniconda3/envs/DeepIP/")
#  # Build your own model
#  py_run_string("trainSeq='train_mini.fa';
#                trainedModel='train_mini.model.hdf5';
#                epoch=10; seqLabel='01'")
#  py_run_file('DeepIP_train.py')
#  

## ----message=FALSE, eval=FALSE------------------------------------------------
#  library(reticulate)
#  use_condaenv("/path/miniconda3/envs/DeepIP/")
#  # to identify whether PACs are internal priming artifacts
#  py_run_string("testSeq='test_mini.fa';
#                trainedModel='train.mini.model.hdf5';
#                outputFile='train.mini.predicted.csv';
#                seqLabel='0/1'")
#  py_run_file('DeepIP_test.py')

## ----message=FALSE, eval=FALSE------------------------------------------------
#  library(PolyAseqTrap, warn.conflicts = FALSE, quietly=TRUE)
#  library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
#  bsgenome <-BSgenome.Hsapiens.UCSC.hg38
#  ## load identified PACs in human genome generated by FindPTA function
#  data(PACs_human)
#  
#  ## extract 200 bp genome sequences surrounding A-rich polyA sites
#  generateFASTA(reads=PACs_human$pa.table,
#                bsgenome=bsgenome,
#                output.name="human_Arich_PACs.fasta")

## ----message=FALSE, eval=FALSE------------------------------------------------
#  ## use DeepIP to classify A-rich polyA sites
#  library(reticulate)
#  use_condaenv("/path/miniconda3/envs/DeepIP/")
#  py_run_string("testSeq='human_Arich_PACs.fasta';
#                trainedModel='human.train.model.hdf5';
#                outputFile='DeepIP_result_hg.csv'")
#  py_run_file("DeepIP_test.py")
#  
#  # Alternatively, run the following Python command outside R:
#  # python DeepIP_test.py -testSeq human_Arich_PACs.fasta  -trainedModel human.train.model.hdf5 -outputFile DeepIP_result_hg.csv

## ----message=FALSE, eval=FALSE------------------------------------------------
#  ### remove internal priming artifacts
#  ip.table <- read.csv("DeepIP_result_hg.csv")
#  head(ip.table)
#  #title                score     true_label predict_label res
#  #>chr22_19354941_-_0 0.4143335          0             0  TN
#  #>chr22_19382772_-_0 0.4143623          0             0  TN
#  #>chr22_19849373_-_0 0.4143335          0             0  TN
#  
#  ## extract and remove internal priming artifacts
#  ip.table <- subset(  ip.table,predict_label==0)
#  ip.table$title <- gsub("^>","",  ip.table$title)
#  PACs_human$pa.table$label <- paste0(
#    PACs_human$pa.table$chr,"_",
#    PACs_human$pa.table$coord,"_",
#    PACs_human$pa.table$strand,"_0")
#  PACs_human$pa.table$level <- as.character(PACs_human$pa.table$level)
#  index <- which(PACs_human$pa.table$label %in% ip.table$title)
#  
#  # remove internal priming artifacts
#  PACs_human$pa.table$level[index] <- "Count"
#  PACs_human$pa.table$level <- factor( PACs_human$pa.table$level,
#          levels=c("V1","V2","V3","V4","V5","V6","V7","V8","Count"))
#  
#  
#  ## regroup nearby cleavage sites
#  PACs_human <- resut.PA(aln.result=PACs_human$pa.table,d=24)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(PolyAseqTrap, warn.conflicts = FALSE, quietly=TRUE)
## load identified PACs in Arabidopsis genome generated by FindPTA function
data("PACs_tair")

#filter PACs that were supported by at least five reads
PACs_tair$pa.coord <- subset(PACs_tair$pa.coord,total.count>=5)

# reduce the impact of microheterogeneity
PACs_tair <-split_pac(pa.data=PACs_tair,d=24,mc.cores=1)

#filter refined PACs that were supported by at least five reads
PACs_tair$split.clusters<- subset(PACs_tair$split.clusters,total.count>=5)


## compare and visualize the difference before and after reducing microheterogeneity
library(dplyr)
library(ggplot2)
width.density<- data.frame(pac.width=c(PACs_tair$pa.coord$width,
                                       PACs_tair$split.clusters$width) ,
                           type=rep(c("before-cluster","after-cluster"),
                                    time=c(length(PACs_tair$pa.coord$width),length(PACs_tair$split.clusters$width)))
)
width.density %>% dplyr::group_by(type) %>% dplyr::summarise(mean=mean(pac.width),median=median(pac.width ))

# density plot
mu <- plyr::ddply(width.density, "type", summarise, grp.mean=mean(pac.width))
ggplot(width.density, aes(x=pac.width,color=type,fill=type)) + 
  geom_density(alpha=0.6)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=type),
             linetype="dashed")+theme_bw()+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+xlim(0,300)+
  labs(x="Width of PACs (nt)",y="Density")+
  guides(fill = guide_legend(title="Category"),color="none")+
  theme(legend.position = c(0.7,0.7))

## ----annotate_PAC, message=FALSE,eval=FALSE-----------------------------------
#  library(PolyAseqTrap, warn.conflicts = FALSE, quietly=TRUE)
#  library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
#  library(movAPA, warn.conflicts = FALSE, quietly=TRUE)
#  bsgenome <-BSgenome.Hsapiens.UCSC.hg38
#  athGFF <- "Homo_sapiens.GRCh38.110.gtf"
#  annotation <- parseGff(athGFF)
#  
#  # load identified PACs in human genome using FindPTA function
#  data(PACs_human)
#  
#  #to check whether the representative polyA site is located in a repeat region
#  #For polyA sites located in intronic regions, it is recommended to filter out those in repeat regions
#  PACs_human$pa.coord <-  check.repeat(pac.result=PACs_human$pa.coord,
#                                       bsgenome =bsgenome )
#  
#  colnames(PACs_human$pa.coord)[1:3] <- c("chr","UPA_start","UPA_end")
#  data.PACds <- readPACds(PACs_human$pa.coord, colDataFile=NULL)
#  # annotate PAC
#  data.PACds <- annotatePAC(data.PACds, aGFF = annotation)
#  # extend 3'UTR region
#  data.PACds <- ext3UTRPACds(data.PACds,ext3UTRlen = 1000)
#  # identify polyA signals
#  data.PACds<- annotateByPAS(data.PACds, bsgenome, grams='AATAAA', from=-50, to=25, label=NULL)
#  data.PACds <- annotateByPAS(data.PACds, bsgenome, grams='V1', from=-50, to=25, label=NULL)
#  data.PACds@anno$pA.signal <- "Others"
#  data.PACds@anno$pA.signal[which(!is.na(data.PACds@anno$V1_dist))] <- "1Variants"
#  data.PACds@anno$pA.signal[which(!is.na(data.PACds@anno$AATAAA_dist))] <- "AATAAA"
#  
#  table(data.PACds@anno$pA.signal)
#  #1Variants    AATAAA    Others
#  #1346       721       126
#  
#  table(data.PACds@anno$ftr)
#  #3UTR       5UTR       exon intergenic     intron
#  #941          1         54        608        589
#  
#  # For intonic PACs, it is recommended to filter out those in repeat regions
#  ipa.result <- subset(data.PACds@anno,ftr=="intron"&repeat_flag==FALSE)
#  save(data.PACds,file="data.PACds.rda")
#  

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  library(ggplot2)
#  library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
#  library(movAPA, warn.conflicts = FALSE, quietly=TRUE)
#  library(ggpubr)
#  bsgenome <-BSgenome.Hsapiens.UCSC.hg38
#  data(data.PACds)
#  faFiles=faFromPACds(data.PACds, bsgenome, what='updn', fapre='updn',
#                      up=-100, dn=100,byGrp='ftr')
#  #608 >>> updn.intergenic.fa
#  #589 >>> updn.intron.fa
#  #941 >>> updn.3UTR.fa
#  #54 >>> updn.exon.fa
#  #1 >>> updn.5UTR.fa
#  
#  ##plot single nucleotide profiles for 3'UTR PACs
#  plotATCGforFAfile("updn.3UTR.fa", ofreq=FALSE, opdf=FALSE,
#                    refPos=101, mergePlots = TRUE)
#  
#  
#  ##calculate a chi-square metric to assess the similarity between the single nucleotide profile of identified polyA sites and the reference profile
#  #load reference PACs profile
#  data(ref.data)
#  
#  #This demonstrates functionality with 10 iterations and 200 random PACs. Larger parameters are recommended (e.g., iteration = 100, use.size = 5000)
#  chisq.result <- cal.chisq(fafile="updn.3UTR.fa",ref.data=ref.data,
#                            iteration=10,use.size=200)
#  statistic <- as.numeric(unlist(chisq.result$statistic))
#  summary(statistic )
#  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  #8.219   8.878   8.957   8.997   9.261   9.589
#  
#  boxplot(statistic ,col="orange",
#          ylab = "Chi-squared metric",
#          border = "brown"
#  )
#  
#  ##plot distributions of PACs across different categories, subclasses, genomic features, and polyA signals
#  data.PACds@anno$coord_class <- factor(data.PACds@anno$coord_level,
#            levels=c("V1","V2","V3","V4","V5","V6","V7","V8","Count"),
#          labels= c("C1","C1","C2","C2","C1","C2","C2","C3","Count"))
#  
#  p<- plot_summary(data=data.PACds@anno)
#  ggarrange(p$p1,p$p2,p$p3,p$p4,    labels = c("A", "B", "C","D"))
#  

## ----echo=FALSE, fig.cap="Distribution of PACs" ,out.width="70%"--------------
knitr::include_graphics("/Users/sevenye/Desktop/Project_Desk/polyAseqTrap/github/refer/features.png")

## -----------------------------------------------------------------------------
sessionInfo()

