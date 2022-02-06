library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'rds',	'r',	1,	"character",
	'gene',	'g',	1,	"character",
	'outdir',	'o',	1,	"character"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	prefix:the output prefix of files,such as pictures and excel.
	rds: the rds file use for plot vln
	gene: the gene file name ,one gene per line
	outdir:outdir  of outputs,we will setwd(opt$outdir)
	Usage example:
	Rscript this.r -p pref_out -r sample.rds -g gene.name.xls  -o outdir 
	Options:
	--help		h	NULL		get this help
	--rds	r	character	input rds file[forced]
	--gene	g	character	gene file [forced]
	--outdir	o	character	The	result of outdir for analysis [forced]
	--prefix	p	character	the prefix for outputfiles [forced]
	\n")
	q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$rds) )	{ cat("Please input the data rds file ...\n\n") ; print_usage(para)}
if ( is.null(opt$gene) )	{ cat("Please gene file ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )	{ cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para) }

require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)

mkdirs <- function(outdir,fp) {
	if(!file.exists(file.path(outdir,fp))) {
		dir.create(file.path(outdir,fp))
	}else{
			print(paste(fp,"Dir already exists!",sep="     "))
			unlink(file.path(outdir,fp), recursive=TRUE)
			dir.create(file.path(outdir,fp))
		}
}

genes= read.table(opt$gene,header=F,stringsAsFactors=F)$V1
immune.combined = readRDS(opt$rds)
DefaultAssay(immune.combined) <- "RNA"
out_dir=opt$outdir

pre_fix=opt$prefix
pdf(paste(paste(out_dir,pre_fix,sep="/"),".pdf",sep=""),w=6,h=2)
plots=VlnPlot(immune.combined,features=genes ,pt.size = 0, combine = FALSE)
print(plots)
dev.off()

tiff(paste(paste(out_dir,pre_fix,sep="/"),".tiff",sep=""),w=6,h=2*length(genes),unit='in', compression="lzw", res=300)
VlnPlot(immune.combined,features=genes ,pt.size = 0, combine = TRUE, ncol=1)
dev.off()
