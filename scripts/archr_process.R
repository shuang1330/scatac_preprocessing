library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(Matrix)

num_threads <- as.numeric(Sys.getenv("OMP_NUM_THREADS", unset = 1))
addArchRThreads(threads = num_threads)

args <- commandArgs(trailingOnly = TRUE)
fragment_filepath <- args[1]
samplename <- args[2]


# setwd(paste("arrowfiles/", samplename, sep=""))

save_genematrix <- function(fragment_filepath, samplename){
    ArrowFiles <- createArrowFiles(
        inputFiles = fragment_filepath, 
        sampleNames = samplename,
        minTSS = 4, 
        minFrags = 1000, 
        addTileMat = FALSE,
        addGeneScoreMat = TRUE
    )
    proj <- ArchRProject(
      ArrowFiles = ArrowFiles, 
      outputDirectory = paste(samplename, "_proj/", sep = ""),
      copyArrows = FALSE
    )
    genematrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
    mat <- assay(genematrix)  # assuming this is already a dgCMatrix
    Matrix::writeMM(mat, file = paste(samplename, "_genescore.mtx", sep = ""))
    write.table(as.data.frame(rowData(genematrix)), file = paste(samplename, "_features.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
    write.table(as.data.frame(colData(genematrix)), file = paste(samplename, "_barcodes.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}

create_arrowfile <- function(fragment_filepath, samplename){
    ArrowFiles <- createArrowFiles(
        inputFiles = fragment_filepath, 
        sampleNames = samplename,
        minTSS = 4, 
        minFrags = 1000, 
        addTileMat = FALSE,
        addGeneScoreMat = TRUE
    )
}

# print(paste(fragment_filepath, samplename))
save_genematrix(fragment_filepath, samplename)
# create_arrowfile(fragment_filepath, samplename)