#Load packages
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

#pheno_data_file <- "mouse.csv"
pheno_data = read.csv("mouse.csv")
bg_mouse = ballgown(dataDir = "ballgown", samplePattern = "GRCm38.p6_gtf_", pData = pheno_data)
bg_mouse

#pheno_data <- read.table(pheno_data_file, header=TRUE, colClasses = rep("character", 2))
pheno_data

#Make sure all samples are included in the analysis.
list.files("ballgown")

# this table needs to be alphabetically ordered...if you didn't create it this way use the following line to organize it.
#pheno_data = pheno_data[order(pheno_data$ids),]

# this test should return TRUE
all(pheno_data$ids == list.files("ballgown"))

# filter -> only take rows with ...
bg_mouse_filt = subset(bg_mouse,"rowVars(texpr(bg_mouse)) > 1.5",genomesubset=TRUE)
bg_mouse_filt

# find de genes
# results_genes = stattest(bg_tick_filt, feature="gene", covariate="X..inf_status",  getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_mouse_filt, feature="gene", covariate="inf_status",  getFC=TRUE, meas="FPKM")

results_genes = arrange(results_genes,pval)
head(results_genes, n=30)
write.csv(results_genes, "mouse_de_gene_pval_results.csv", row.names = FALSE)

results_genes = arrange(results_genes,desc(fc))
head(results_genes, n=30)
write.csv(results_genes, "mouse_de_gene_fc_results.csv", row.names = FALSE)

results_genes = arrange(results_genes,qval)
head(results_genes, n=30)
write.csv(results_genes, "mouse_de_gene_qval_results.csv", row.names = FALSE)

#find de transcripts
results_transcripts = stattest(bg_mouse_filt, feature="transcript",covariate="inf_status", getFC=TRUE, meas="FPKM")
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_mouse_filt), geneIDs=ballgown::geneIDs(bg_mouse_filt), results_transcripts)

results_transcripts = arrange(results_transcripts,pval)
head(results_transcripts, n=30)
write.csv(results_transcripts, "mouse_de_transcript_pval_results.csv", row.names = FALSE)

results_transcripts = arrange(results_transcripts,desc(fc))
head(results_transcripts, n=30)
write.csv(results_transcripts, "mouse_de_transcript_fc_results.csv", row.names = FALSE)

results_transcripts = arrange(results_transcripts,qval)
head(results_transcripts, n=30)
write.csv(results_transcripts, "mouse_de_transcript_qval_results.csv", row.names = FALSE)

table(results_transcripts$qval < 0.05)

#plot results
results_transcripts_qvals <- results_transcripts$qval
hist(results_transcripts_qvals)
head(results_transcripts_qvals, n=30)

library(ggplot2)
ggplot(results_transcripts, aes(x = qval, color=fc)) + geom_histogram(bins=100) + scale_y_log10()
ggplot(results_transcripts, aes(x = qval, color=fc)) + geom_histogram(bins=100)

library(cowplot)
results_transcripts$mean <- rowMeans(texpr(bg_mouse_filt))

ggplot(results_transcripts, aes(log2(mean), log2(fc), colour = qval<0.05)) +
  scale_color_manual(values=c("#999999", "#FF0000")) +
  geom_point() +
  geom_hline(yintercept=0)

#Below are some things I tried to get to work with variable success and much of this has been improved and incorporated above already.
#transcriptNames(bg_mouse_filt)

#options(max.print=10000)

#get the full table of filtered 
#full_table <- texpr(bg_mouse_filt, 'all')
#full_table


#genename_ids <- dplyr::filter(results_transcripts, geneNames!=".") %>% dplyr::distinct(geneNames, geneIDs)

#Convert archaic MSTG ids to gene names
#results_genes1 <- dplyr::left_join(results_genes, genename_ids, by=c("id"="geneIDs"))
#results_genes1

###tried to get transcript names but not working yet
#results_transcripts1 <- dplyr::left_join(results_transcripts_qvals, transcriptNames, by=c("id"=transcriptIDs))

#export gene results with gene names
#write.table(results_genes1, "/Users/jlee337/izzotemp/results_genes2.txt", sep="\t")
