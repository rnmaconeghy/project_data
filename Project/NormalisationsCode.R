if(! require(affy)){
  source("http://bioconductor.org/biocLite.R")
  biocLite("affy")
  require(affy)
}

library(GEOquery)
library(arrayQualityMetrics)
library(chicken.db)
library(chickenprobe)
library(annotate)
library(genefilter)
library(cluster)
library(limma) #Two-colour pre-processing; differential expression
library(affyPLM)


workingDir="/homes/rnmaconeghy/BS32010/Project/GCOS_Cel"

filenames <- c("ROS1-_9.CEL", "ROS1+_10.CEL", "ROS2-_11.CEL", "ROS2+_12.CEL", "ROS3-_13.CEL", "ROS3+_14.CEL", "ROS4-_15.CEL", "ROS4+_16.CEL")
samplenames <- c("ROS1-_9.CEL", "ROS1+_10.CEL", "ROS2-_11.CEL", "ROS2+_12.CEL", "ROS3-_13.CEL", "ROS3+_14.CEL", "ROS4-_15.CEL", "ROS4+_16.CEL")
targets <- c("minus", "plus", "minus", "plus", "minus", "plus", "minus", "plus")



#import data, describing the experimental design
phenodata<-as.data.frame(cbind(filenames,samplenames,targets))
#write.table(phenodata,paste(workingDir,))
celRAW <- ReadAffy(celfile.path=workingDir, phenoData=phenodata)

.plotDensity <- function(exps, filename)
{
  pdf(filename)
  d <- apply(exps,2,
             function(x){
               density(x)
             })
  xmax <- max(sapply(d,function(x)max(x$x)))
  xmin <- min(sapply(d,function(x)min(x$x)))
  ymax <- max(sapply(d,function(x)max(x$y)))
  plot(0,pch='',ylab='',xlab='',
       xlim=c(xmin,round(xmax+1)),ylim=c(0,ymax))
  lapply(1:length(d),function(x) lines(d[[x]],col=x))
  dev.off()
}

.doPLM <- function(celRAW)
{
  pdf("celRAWqc.pdf")
  celRAWqc <- fitPLM(celRAW)
  image(celRAWqc, which=1, add.legend=TRUE)
  image(celRAWqc, which=4, add.legend=TRUE)
  RLE(celRAWqc, main="RLE")
  NUSE(celRAWqc, main="NUSE")
  dev.off()
}

.doCluster <- function(celRMA)
{
  eset <- exprs(celRMA)
  distance <- dist(t(eset), method="maximum")
  clusters <- hclust(distance)
  plot(clusters)
}

.doFilter <- function(celRMA)
{
  celfiles.filtered <- nsFilter(celRMA, 
                                require.entrez=FALSE,
                                remove.dupEntrez=FALSE)
}

.doDE <- function(eset)
{
 samples <- eset$targets
 samples <- as.factor(samples)
 design <- model.matrix(~0 + samples)
 colnames(design) <- c("plus", "minus")
 fit <- lmFit(exprs(eset), design)
 ebFit <- eBayes(fit)
 ttab <- topTable(ebFit, number=10000)
 annotation <- as.data.frame(select(chicken.db,
                                    rownames(ttab),
                                    c("ENSEMBL", "SYMBOL")))
 colnames(annotation) <- c("probeId", "ensemblId", "geneSymbol")
 results <- merge(annotation, ttab,by.x="probeId",by.y="row.names")
 
 head(results)
 write.table(results, "results.txt", sep="\t", quote=FALSE)
 return(results)
}

if(!exists("celResults"))
{
  eset <- exprs(celRAW)
  celGCRMA <- gcrma(celRAW)
  celRMA <- rma(celRAW)
  .plotDensity(log2(exprs(celRAW)),"densityRAW.pdf")
  .plotDensity(log2(exprs(celRMA)),"densityRMA.pdf")
  celFilt <- .doFilter(celRMA)
  celResults <- .doDE(celFilt$eset)
}

