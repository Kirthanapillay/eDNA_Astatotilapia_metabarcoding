##analysis of reads in R using DADA2 

library(dada2)

miseq_path <- "/Users/gkpillay/Desktop/Write_up/Florida_gut/Data_analysis/COI/"
list.files(miseq_path)

setwd("/Users/gkpillay/Desktop/Write_up/Florida_gut/Data_analysis/")

fnFs <- sort(list.files(miseq_path,pattern="_COI_1.fastq.gz",full.names = TRUE))
fnRs <- sort(list.files(miseq_path,pattern="_COI_2.fastq.gz",full.names = TRUE))

length(fnFs)
length(fnRs)

sampleNames <- sapply(strsplit(basename(fnFs), "_COI_cut.1.fastq.gz"), `[`, 1)
head(sampleNames)

plotQualityProfile(fnRs[102:104])

filtFs <- file.path(miseq_path, "filtered", paste0(sampleNames, "_F_filt.fastq"))
filtRs <- file.path(miseq_path, "filtered", paste0(sampleNames, "_R_filt.fastq"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,truncLen = c(200,200),
                         maxN=0, rm.phix=TRUE, truncQ=2,maxEE = c(2,2),
                         compress=TRUE, multithread=TRUE, minLen = 100)


#only keep reads that have passed the previous step
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

head(out)
plotQualityProfile(filtFs[102:104])

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

head(dadaFs)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))
# only keep AVS 305-316 bp


seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(305,316)]
dim(seqtab2)

table(nchar(getSequences(seqtab2)))

seqtab2.nonchim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 7890 bimeras out of 12408 input sequences.
dim(seqtab2.nonchim)

sum(seqtab2.nonchim)/sum(seqtab2)

table(nchar(getSequences(seqtab2.nonchim)))


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2.nonchim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
head(track)

write.table(seqtab2.nonchim, "seqtab2.nonchim.txt")
write.csv(seqtab2.nonchim, "seqtab2.nonchim_COI.csv")

#make into fasta format
#Get ASV sequences
asv_seqs <- colnames(seqtab2.nonchim)
asv_headers <- vector(dim(seqtab2.nonchim)[2], mode="character")

for (i in 1:dim(seqtab2.nonchim)[2]) {
  asv_headers[i] <- paste(">ASV", i , sep="")
  
}

#Write to FASTA file
asv_fasta<-c(rbind(seqtab2.nonchim, asv_seqs))
write(asv_fasta, "ASVs_COI.fa")

#import ASVs_COI.fa into linux to assign taxonomy using blast
#repeat above for 18S sequences



