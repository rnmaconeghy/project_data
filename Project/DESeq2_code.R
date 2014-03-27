require(DESeq2)
require(RCurl)
require(biomaRt)

remoteFile <- getURL(paste0("http://www.compbio.dundee.ac.uk/",
                            "user/pschofield/Teaching/Bioinformatics/", 
                            "Data/RNAseqCounts.txt"))

dat<-read.delim("/usr/local/share/BS32010/expression/data/RNAseqCounts.txt", 
                sep="\t",skip=1,head=T,row.names=1)

geneCounts<-read.delim(textConnection(remoteFile), 
                       head=T, sep="\t",skip=1)


nonZeroCounts<-geneCounts[rowSums(geneCounts[,7:12])>0,7:12]

treatments <- as.factor(substr(colnames(nonZeroCounts),1,1))

dds<- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts), 
                             as.data.frame(treatments),
                             design=~treatments)

dds$treatments <- relevel(dds$treatments, "U")

dds <-DESeq(dds)
DEresults <- results(dds)


if(interactive()){
  mart <- useMart("ensembl")
  datasets <- listDatasets(mart)
  mart<-useDataset("ggallus_gene_ensembl",mart)
  getBM(attributes=c("ensembl_gene_id","external_gene_id",
                     "affy_chicken"),
        filter="ensembl_gene_id",
        values=c("ROS1-_9", "ROS1+_10", "ROS2-_11", "ROS2+_12", "ROS3-_13", "ROS3+_14", "ROS4-_15", "ROS4+_16"), 
        mart=mart)
}

#mart <- useMart (biomart="ENSEMBL_MART_ENSEMBL",
 #                host="www.ensembl.org",
  #               path="/biomart/martservice")
#returns null results

gg4 <- useDataset("ggallus_gene_ensembl", mart=mart)

annot <- getBM(attributes=c("ensembl_gene_id",
                            "external_gene_id",
                            "affy_chicken"), 
               filter="ensembl_gene_id", 
               values=rownames(DEresults), 
               mart=gg4)
#something is wrong as this gives 0 obs. of 3 variables

annotResults <- merge(annot,DEresults, 
                      by.x="ensembl_gene_id",by.y="row.names")
#all reults are null



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

#design<-model.matrix(~0 + sample)
 #fit<-lmFit(eset,design),
  #design(dds)<-formula(~type+condition),
  #dds<-DESeq(dds),
