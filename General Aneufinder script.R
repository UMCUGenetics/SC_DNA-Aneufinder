# Ewart Kuijk
# 27092021
# The following script is based on "A quick introduction to AneuFinder"

# 3.3 Running Aneufinder
library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg19) #adapt to your experiment
library(Rsubread)

#---------REFERENCE GENOME--------
# Use Input BAM-files, mappability files, variable.width.reference. files, and blacklist files that share the SAME reference genome

#---------MAPPABILITY CORRECTION FILE-------
# if necessary create mappability correction file using the simulateReads() command (see Aneufinder manual/vignette for details)
## Get an example BAM file with single-cell-sequencing reads or a merged BAM file from multiple cells
bamfile <- "filepath/yourfile.bam"
## Simulate 51bp reads for at a distance of every 5000bp
if (require(BSgenome.Hsapiens.UCSC.hg19)) {
  simulateReads(BSgenome.Hsapiens.UCSC.hg19, bamfile=bamfile, readLength=51,
                file="filepath/simulatedreads", every.X.bp=500000)
}

#---------VARIABLE WIDTH REFERENCE AND BLACKLISTING-------
# For good Aneufinder results it is important to have a good variable.width.reference file and blacklist
# Blacklists can be made from a euploid file with Aneufinder using the following script:

## Get a euploid reference (adjust input file to your experiment)
bedfile <- system.file("extdata", "hg19_diploid.bam.bed.gz", package="AneuFinderData")
## Make 100kb fixed-width bins
bins <- binReads(bedfile, assembly='hg19', binsizes=100e3,
                 chromosomes=c(1:22,'X'))[[1]]
## Make a plot for visual inspection and get the blacklist
lcutoff <- quantile(bins$counts, 0.05)
ucutoff <- quantile(bins$counts, 0.999)
p <- plot(bins) + coord_cartesian(ylim=c(0,50))
p <- p + geom_hline(aes(yintercept=lcutoff), color='red')
p <- p + geom_hline(aes(yintercept=ucutoff), color='red')
print(p)

# Alternatively, it is possible to use the ENCODE blacklist: https://github.com/Boyle-Lab/Blacklist/
# Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: 
# Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z
# In my hands, creating your own blacklist gives the best results!

#-------------RUN ANEUFINDER--------
#see Aneufinder manual for setting all the different options
Aneufinder("filepath_to_your_single_cell_seq_BAMS_folder", 
           "/filepath_to_your_output_folder/Aneufinder-output", 
           configfile = NULL, numCPU = 2, reuse.existing.files = TRUE, binsizes = 1e+06, 
           variable.width.reference = "/filepath_to_your_variable.width.reference/variable.width.reference.bam", 
           reads.per.bin = NULL, pairedEndReads = FALSE, chromosomes = c(1:22,'X'), remove.duplicate.reads = TRUE, min.mapq = 10, 
           blacklist = "/filepath_to_your_blacklist/blacklist.file.bed.gz", use.bamsignals = TRUE, 
           reads.store = FALSE, correction.method = 'GC', GC.BSgenome = 'BSgenome.Hsapiens.UCSC.hg19', method = 'edivisive', strandseq = FALSE, R = 10, sig.lvl = 0.1, 
           eps = 0.01, max.time = 60, max.iter = 5000, num.trials = 15, states = c("zero-inflation", paste0(0:10, "-somy")), confint = NULL, 
           refine.breakpoints = FALSE, hotspot.bandwidth = NULL, hotspot.pval = 0.05, cluster.plots = FALSE)



#----------- 3.4 Loading results and plotting single cells---------
files <- list.files('/filepath_to_your_outputfolder/Aneufinder-output/MODELS/method-edivisive', full.names=TRUE) 
models <- loadFromFiles(files)
plot(files[1], type='profile')
plot(files[1], type='histogram')
plot_cell_01 <- plot(files[1], type='karyogram')
plot_cell_02 <- plot(files[2], type='karyogram')

#------------ 3.5 Quality control----------
# We found that simple filtering procedures such as cutoffs on the total number of
# reads etc., are insufficient to distinguish good from bad-quality libraries. 
# Therefore, we have implemented a multivariate clustering approach that works on 
# multiple quality metrics (see ?clusterByQuality for details) to increase 
# robustness of the filtering. Here is an example demonstrating the usage of 
# clusterByQuality().

cl <- clusterByQuality(files, measures=c('spikiness','num.segments','entropy','bhattacharyya','sos'))
plot(cl$Mclust, what='classification')

print(cl$parameters)

# select the best clusters and plot those
selected.files <- unlist(cl$classification[1:10]) 
heatmapGenomewide(selected.files)

#-------------- 3.6 Karyotype measures-------------
# calculating Aneuploidy and heterogeneity scores
# Get karyotype measures: lung and liver examples are from the Aneufinder package

results <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData") 
files.lung <- list.files(results, full.names=TRUE)
results <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData") 
files.liver <- list.files(results, full.names=TRUE)
## Get karyotype measures
k.lung <- karyotypeMeasures(files.lung) 
k.liver <- karyotypeMeasures(files.liver)                    


df <- rbind(lung = k.lung$genomewide, liver = k.liver$genomewide) 
print(df)

plotHeterogeneity(hmms.list = list(lung=files.lung, liver=files.liver))

## Get karyotype measures using your own data
files.yourcells <- list.files('/filepath_to_your_outputfolder/Aneufinder-output/MODELS/method-edivisive', full.names=TRUE)
k.c2NI <- karyotypeMeasures(files.yourcells)
df <- rbind(lung = k.lung$genomewide, yourcells = k.c2NI$genomewide) 
print(df)

plotHeterogeneity(hmms.list = list(lung=files.lung, yourcells=files.yourcells))


# 3.7 Principal component analysis
lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData") 
lung.files <- list.files(lung.folder, full.names=TRUE)

## Get results from the liver metastasis of the same patient
liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData") 
liver.files <- list.files(liver.folder, full.names=TRUE)

## Plot a clustered heatmap
classes <- c(rep('lung', length(lung.files)), rep('liver', length(liver.files))) 
labels <- c(paste('lung',1:length(lung.files)), paste('liver',1:length(liver.files))) 
plot_pca(c(lung.files, liver.files), colorBy=classes, PC1=2, PC2=4)

## PCA from your data 
yourcells.files <- list.files('/filepath_to_your_outputfolder/Aneufinder-output/MODELS/method-edivisive', full.names=TRUE)

yourcells.filesa <- yourcells.files[1:5]
yourcells.filesb <- yourcells.files[6:10]

## Plot a clustered heatmap
classes <- c(rep('a', length(yourcells.filesa)), rep('b', length(yourcells.filesb))) 
labels <- c(paste('a',1:length(yourcells.filesa)), paste('b',1:length(yourcells.filesb))) 
plot_pca(c(yourcells.filesa, yourcells.filesb), colorBy=classes, PC1=2, PC2=4)

