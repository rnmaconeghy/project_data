require(DESeq2)
require(biomaRt)

remoteFile <- getURL(paste0("http://www.compbio.dundee.ac.uk/",
                            "user/pschofield/Teaching/Bioinformatics/", 
                            "Data/RNAseqCounts.txt"))

dat<-read.delim("/usr/local/share/BS32010/expression/data/RNAseqCounts.txt",
                sep="\t",skip=1,head=T,row.names=1)

# if at this point you do a colnames(dat) you would see the rosc samples are
# columns 14 - 20

roscdata <- dat[,14:20]

# at this point you need to set up the sample ids and treatment stuff
# in a similar way to the limma code

colnames(roscdata) <- c("ROS1_minus","ROS1_plus",
                        "ROS2_minus","ROS2_plus",
                        "ROS3_minus","ROS2_plus",
                        "ROS4_minus")

# NB there is not a column for the ROS4_plus sample in the RNA-seq data
# now you need the genes that have at least some reads this is similar to
# the nsFilter step in the Limma microarray process

nonZeroCounts<-roscdata[rowSums(roscdata)>0,]

# now you need to assign each column to a treatment group ie plus or minus

treatments <- c("minus","plus","minus","plus","minus","plus","minus")

# now you need to set up the DESeq dataset this is done all in one step
# you include the nonZeroCount data, the assignment of the columns of the
# data to treatment groups though the treatments vector (as a data.frame)
# and the model is( design=~treatment)

dds<- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts),
                             as.data.frame(treatments),
                             design=~treatments)

# next you say which treatment group is the control here it is the minus group

dds$treatments <- relevel(dds$treatments, "minus")

# then you actually have to run the differential expression analysis

dds <-DESeq(dds)

# and get the results

DEresults <- results(dds)

#this did work at first. No idea why it isn't now

# now you need to get the annotations for the gene list
# connect to the ensembl biomart server
mart <- useMart("ensembl")
# request the chicken (Gallus gallus) dataset
mart<-useDataset("ggallus_gene_ensembl",mart)

# get the actual annotations for all the genes in you results set
# you specify the informantion attributes you want and tell it
# what key to filter the data set on and then specify the values of
# the key (ie the gene_ids from the results table which are in the rownames
annot <- getBM(attributes=c("ensembl_gene_id","external_gene_id",
                            "affy_chicken"),
               filter="ensembl_gene_id",
               values=rownames(DEresults),
               mart=mart)
#then you merge the data columns to get your annotated list
annotResults <- merge(annot,DEresults, 
                      by.x="ensembl_gene_id",by.y="row.names")

res <- results(dds)
head(res)
mcols(res)
sum( res$pvalue < 0.05, na.rm=TRUE)
table( is.na(res$pvalue) )
sum( res$padj < 0.1, na.rm=TRUE )
resSig <- res[ which(res$padj <0.1 ), ]
head( resSig[ order( resSig$log2FoldChange ), ] )
tail( resSig[ order( resSig$log2FoldChange ), ] )
plot(attr(res,"filterNumRej"),type="b", ylab="number of rejections")
#plots the number of rejected tests less than 0.1 plotted over theta, the quantiles of the filter statistic
sum(res$padj < 0.1, na.rm=TRUE)
resNoFilt <- results(dds, independentFiltering=FALSE)
sum(resNoFilt$padj < 0.1, na.rm=TRUE)
plotMA(dds, ylim = c(-2,2), main="DESeq2")
plotDispEsts(dds)
hist(res$pvalue, breaks=100)
write.csv(as.data.frame(res), file="results.csv")
#saves the results in a CSV file, which you can then load with Excel
colData(dds)

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
