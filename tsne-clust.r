#!/usr/bin/env Rscript
library(Rtsne)
library(ggplot2)
library(dbscan)
library(getopt)
library(RColorBrewer)

spec = matrix(c(
	'help','h', 0, "logical","Print short help message",
	'log','L', 0, "logical","Perform logarythmization (not really recommended)",
	'clustering','c', 2, "character","Clustering mode: h - HDBSCAN (default), k- k-means",
	'inputfile','i', 1, "character", "Data file with TNF, coverage and contig length",
	'contignumber','n', 1, "integer","Total number of analyzed contigs <10000>",
	'minimumlength','m',1,"integer","Minimum length of contig <1000>"
), byrow=TRUE, ncol=5)
opt = getopt(spec)


if ( !is.null(opt$help) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

if ( is.null(opt$inputfile) ) {
	cat("No input file specified!\n")
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# Clustering mode
if ( is.null( opt$clustering    ) ) { opt$clustering    = "h" }
# Max contig number
if ( is.null( opt$contignumber  ) ) { opt$contignumber  = 10000 }
# Minimum contig length
if ( is.null( opt$minimumlength ) ) { opt$minimumlength = 1000  }
# Logarythm
if ( is.null( opt$log           ) ) { opt$log           = FALSE }


rdata <- read.delim(opt$inputfile,header=F)
lastcol = ncol(rdata)
rdata <- rdata[as.matrix(rdata[,lastcol])>=opt$minimumlength,]

n <- min(opt$contignumber, nrow(rdata))

mdata <- as.matrix(rdata[1:n,2:(lastcol-1)])
meanByCol <- apply(mdata,2,mean)
normdata <- sweep(mdata,2,meanByCol,"/")

# Do logarythmization of data
if(opt$log){
	normdata <- apply(normdata + unique(sort(as.vector(normdata)))[2]/4,2,log)
}

cat(paste0("Loaded data for ", n, " sequences.\n"))


nn <- round(n / 100)/10

#Perplexity for tSNE:
perplexities = ceiling(c(20,40,60) * nn**0.3)
#K-means cluster number
clustersvector = ceiling(c(5,10,15) * nn**0.4)


mytsne <- function(i){
	tsne <- Rtsne(normdata, perplexity=i)
	return(tsne$Y)
}

# For testing purposes only
mykmean <- function(tsne,nc,i){
	kx <- kmeans(tsne,centers=nc,nstart=100,iter.max=10000)
	cldata <- data.frame(rdata[1:n,1],tsne,kx$cluster,rdata[1:n,(lastcol-1):lastcol])
	colnames(cldata)=c("Contig ID","X","Y","cluster","Ncoverage","length")
	cat(paste0("Clustering results for ",nc," clusters:\n"))
	cat("Bin#\tbin size, nt.\n")
	for (x in 1:nc){
		cat(paste0(x,"\t",sum(cldata[cldata$cluster==x,6]),"\n"))
	}
	datafile <- paste0("yamb-pp-",i,"-cl-",nc,".csv")
	write.table(cldata,file=datafile,sep="\t",row.names = FALSE, quote=FALSE)
	p <- ggplot(cldata) + 
		geom_point(aes(x=X,y=Y,size=length,alpha=Ncoverage,colour=as.factor(cluster))) + 
		scale_size(range=c(0.2,10)) +
		geom_label(data = as.data.frame(kx$centers),aes(x=V1,y=V2,label=rownames(kx$centers))) + 
		theme_bw()
	return(p)
}

myhdbscan <- function(tsne){
	hx <- hdbscan(tsne,minPts=20)
	nc <- max(hx$cluster)
	cldata <- data.frame(rdata[1:n,1],tsne,hx$cluster,rdata[1:n,(lastcol-1):lastcol])
	colnames(cldata)=c("Seq ID","X","Y","cluster","Ncoverage","length")
	cat(paste0("Clustering results for ",nc," clusters:\n"))
	cat("Bin#\tbin size, nt.\n")
	for (x in 0:nc){
		cat(paste0(x,"\t",sum(cldata[cldata$cluster==x,6]),"\n"))
	}
	datafile <- paste0("yamb-pp-",i,"-hdbscan.csv")
	write.table(cldata,file=datafile,sep="\t",row.names = FALSE, quote=FALSE)
	p <- ggplot(cldata) + 
		geom_point(aes(x=X,y=Y,size=length,alpha=Ncoverage,colour=as.factor(cluster))) + 
		scale_size(range=c(0.2,10)) + theme_bw() + 
		scale_color_manual(values=c("black",colorRampPalette(brewer.pal(12,"Paired"))(nc)))
	return(p)
}

for(i in perplexities){
	tsne <- mytsne(i)
	if(opt$clustering == "k"){
		for(clusters in clustersvector){
			pngfile <- paste0("yamb-",nn,"K-pp-",i,"-cl-",clusters,".png")
			mykmean(tsne,clusters,i)
			ggsave(file=pngfile)
		}
	}else if(opt$clustering == "h"){
		pngfile <- paste0("yamb-pp-",i,"-hdbscan.png")
		myhdbscan(tsne)
		ggsave(file=pngfile)
	}
}

cat("R analysis completed!\n")
q()
