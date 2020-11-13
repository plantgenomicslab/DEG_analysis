###==================================
#	Modular Tool Set For Differential 
#		Analysis For Genes
#	@author Froilan Luna-Lopez
###==================================
library("argparse")
suppressMessages({library("edgeR")
					library("DESeq2")
				})

getTPM <- function(featureCounts){
	featureCounts[,2:ncol(featureCounts)] <- sweep(featureCounts[,2:ncol(featureCounts)], 1, featureCounts[,1], FUN="/")
	size_ratio <- colSums(featureCounts[,2:ncol(featureCounts)]) / 1000000
	featureCounts[,2:ncol(featureCounts)] <- sweep(featureCounts[,2:ncol(featureCounts)], 2, size_ratio, FUN="/")
	return(featureCounts)

}

getRPKM <- function(featureCounts){
	size_ratio <- colSums(featureCounts[,6:ncol(featureCounts)]) / 1000000
	featureCounts[,6:ncol(featureCounts)] <- sweep(featureCounts[,6:ncol(featureCounts)], 2, size_ratio, FUN="/")
	featureCounts[,6:ncol(featureCounts)] <- sweep(featureCounts[,6:ncol(featureCounts)], 1, featureCounts[,5], FUN="/")
	return(featureCounts)
}

getPairs <- function(samplefile){
	return(levels(samplefile[,1][1]))
}

getRepetitions <- function(samplefile){
	# Create empty data frame, two rows two cols
	reps <- data.frame("Reps"=c(0,0), "test"=c(5,5), row.names=getPairs(samplefile))
	reps$Reps <- as.numeric(as.character(reps$Reps))
	# Loop through rows in samplefile
	for (i in 1:nrow(samplefile)){
		# If row's condition is already in data frame
		reps$Reps[samplefile[i,1]] <- reps$Reps[samplefile[i,1]] + 1
	}
	# Return second column of data frame (should be a vector)
	return(reps$Reps)
}

deseqAnalysis <- function(ftc_matrix, samplefile, args){
	print("Conducting deseq analysis...")
	ftc_matrix <- ftc_matrix[,2:ncol(ftc_matrix)]
	samples <- read.table(samplefile, sep='\t')
	colnames(samples) <- c("condition", "batch")
	# PERFORM DESEQ
	fpkm_dds <- DESeqDataSetFromMatrix(countData = ftc_matrix, colData = samples, design= ~ condition) # generate DESeq matrix with count data against conditions
	fpkm_dds <- DESeq(fpkm_dds, test = args$type) #analyze count data gainst different conditions
	keep <- rowSums(counts(fpkm_dds)) >= 1
	fpkm_dds <- fpkm_dds[keep,]
	cat('\n')
	resultsNames(fpkm_dds)
	res <- results(fpkm_dds, name="condition_WT_vs_CEB")
	return(res)
}

edgeRAnalysis <- function(ftc_matrix, samplefile, methd, args){
	ftc_matrix <- ftc_matrix[,2:ncol(ftc_matrix)]
	sample_conditions <- read.delim(samplefile, header=F, sep='\t')
	conditions <- factor(as.vector(sample_conditions[,1]))
	reps <- getRepetitions(sample_conditions)

	dge_ftcnts_list <- DGEList(counts=ftc_matrix, group=conditions)
	keep <- filterByExpr(dge_ftcnts_list)
	dge_ftcnts_list <- dge_ftcnts_list[keep, keep.lib.sizes=F]
	dge_ftcnts_list <- calcNormFactors(dge_ftcnts_list, method=methd)
	design <- model.matrix(~conditions)
	if(reps[1] > 1 && reps[2] > 1){
		dge_ftcnts_list <- estimateDisp(dge_ftcnts_list, design)
	} else if(is.na(args$dispersion) == 0){
		dge_ftcnts_list <- exactTest(dge_ftcnts_list, pair = getPairs(sample_conditions), dispersion = args$dispersion)
	} else{
		warning("Dispersion value needed ( --dispersion )")
	}
	fit <- glmFit(dge_ftcnts_list, design)
	lrt <- glmLRT(fit, coef=ncol(fit$design))
	toptags <- topTags(lrt, n=nrow(lrt))$table
	return(toptags)
}

# Command Line Argument Parser
#============
# TODO:
# - Separate groups for DESeq, EdgeR Options
# - Add dispersion, min_reps_min_cpm, type for deseq, htseq options
parser <- ArgumentParser(description="Command Line Arguments:")
req_parser <- parser$add_argument_group("required arguments")
req_parser$add_argument("-i", "--input", help = "Input file containing read count matrix.")
req_parser$add_argument("-m", "--method", help = "Method to perform DE analysis [deseq, edger].")
req_parser$add_argument("-sf", "--samplefile", help = "Input file containing sample and treatment.")

edgerparser <- parser$add_argument_group("edgeR optional arguments")
edgerparser$add_argument("-d", "--dispersion", help = "Dispersion value for edgeR to utilize.",
					default = NULL)
edgerparser$add_argument("-t", "--tpm", help = "Optional output of tpm normalize read counts.",
					default = 0)
edgerparser$add_argument("-f", "--fpkm", help = "Optional output of fpkm normalized read counts.",
					default = 0)
edgerparser$add_argument("-r", "--rpkm", help = "Optional output of rpkm normalized read counts.",
					default = 0)
edgerparser$add_argument("-tm", "--tmm", help = "Optional output of TMM normalized read counts.",
					default = 1)

deseqparser <- parser$add_argument_group("deseq2 optional arguments")
deseqparser$add_argument("-T", "--type", help = "Type of hypothesis test to utilize [LRT, Wald].",
					default = "LRT")

parser$add_argument("-mrc", "--min_reps_min_cpm", help = "Minimum repetitions and counts per million allowed. format: min_reps,min_cpm",
					default = "1,1")
parser$add_argument("-s", "--sep", help = "Separation delimiter to use.",
					default = "\t")
parser$add_argument("-o", "--output", help = "File prefix to output resultant data to.",
					default = "DEG_results")

args <- parser$parse_args()
# End parser

# Test if matrix file exists
if(file.exists(args$input)==0){
	stop(paste("Input file ", args$input, " cannot be found."))
}
# End test

matrix_file <- read.table(args$input, header=T, row.names=1, sep = args$sep)
samplefile <- args$samplefile

#=============
# TODO:
# -Modify output file names to match some options
#=============
# DESeq2 Analysis Pathway
if(args$method == "deseq"){
	deseq_results <- deseqAnalysis(matrix_file, samplefile, args)
	write.table(deseq_results, file=args$output, sep="\t", quote=F)
#end DESeq2 pathway
# EdgeR Analysis Option Pathway
} else if(args$method == "edger"){
	if (args$tpm == 1){
		deseq_results <- edgeRAnalysis(getTPM(matrix_file), samplefile, "none", args)
		write.table(deseq_results, file=paste(args$output, ".TPM", sep=""), sep="\t", quote=F)
	}
	if (args$rpkm == 1){
		deseq_results <- edgeRAnalysis(getRPKM(matrix_file), samplefile, "none", args)
		write.table(deseq_results, file=paste(args$output, ".RPKM", sep=""), sep="\t", quote=F)
	}
	if (args$tmm == 1){
		deseq_results <- edgeRAnalysis(matrix_file, samplefile, "TMM", args)
		write.table(deseq_results, file=paste(args$output, ".TMM", sep=""), sep="\t", quote=F)
	}
}
# end edgeR pathway
