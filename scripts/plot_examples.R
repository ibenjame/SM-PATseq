path<-getwd()
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(stringr)

#Read In Table
sample_table<-read_delim("samples.txt", delim="\t")
samples<-sample_table$ID
genotypes<-sample_table$Genotype

#Read in datasets
for(x in samples) {
  assign(paste0("sample",x),read_delim(paste0(x,"/",x,".full.txt.gz"), delim="\t", col_names = FALSE))
}

#Column description: ReadID, Length, Direction, C, G, T, TerminalUs, Gene, AlignmentScore, GencodeLongString
#simplify tables; remove artifacts
for(x in 1:length(samples)) {
  sam<-samples[x]
  assign(paste0("lengths",sam),data.frame("Gene"=get(paste0("sample",sam))$X8,"Length"=get(paste0("sample",sam))$X2, "Cs"=get(paste0("sample",sam))$X4, "Gs"=get(paste0("sample",sam))$X5, "Ts"=get(paste0("sample",sam))$X6, "Us"=get(paste0("sample",sam))$X7, "Sample"=samples[x], "Genotype"=genotypes[x]))
  assign(paste0("lengths",sam),filter(get(paste0("lengths",sam)),Length>10,Length<500,!is.na(Gene)))
}

#Create combined results table
#Set levels to factor in desired plot order
LengthsTable<-data.frame("Gene"=c(),"Length"=c(),"Cs"=c(),"Gs"=c(),"Ts"=c(), "Us"=c(),"Sample"=c(), "Genotype"=c())
for (x in samples) {
  LengthsTable<-rbind(LengthsTable,get(paste0("lengths",x)))
}
LengthsTable$Sample<-factor(LengthsTable$Sample, levels = samples)
LengthsTable$Genotype<-factor(LengthsTable$Genotype)

#cleanup of _mRNA from IDs (for Cerevisiae)
LengthsTable$Gene<- gsub('_mRNA', '', LengthsTable$Gene)

#create plots subdirectory
dir.create(file.path(path,"plots"))

#Cumulative Distribution plots
pdf(file="plots/Cumulative_by_sample.pdf")
ggplot(LengthsTable, aes(Length, color=Sample)) + stat_ecdf(geom = "smooth") +
  labs(title = "All Gene Cumulative Plot", y = "Cumulative Fraction", x="PAT Length") +
  xlim(0,300)
dev.off()

pdf(file="plots/Cumulative_by_genotype.pdf")
ggplot(LengthsTable, aes(Length, color=Genotype)) + stat_ecdf(geom = "smooth") +
  labs(title = "All Gene Cumulative Plot", y = "Cumulative Fraction", x="PAT Length") +
  xlim(0,300)
dev.off()

#Violin Plots
pdf(file="plots/Violin_by_sample.pdf")
ggplot(LengthsTable, aes(x=Sample, y=Length, color=Sample)) + geom_violin() +
  geom_boxplot(width=0.1) +
  labs(title = "All Genes Violin Plot", y = "PAT Length (nt)", x="Sample") +
  ylim(0,300)
dev.off()

pdf(file="plots/Violin_by_genotype.pdf")
ggplot(LengthsTable, aes(x=Genotype, y=Length, color=Genotype)) + geom_violin() +
  geom_boxplot(width=0.1) +
  labs(title = "All Genes Violin Plot", y = "PAT Length (nt)", x="Genotype") +
  ylim(0,300)
dev.off()

#Raw Histogram (No density smoothing)
pdf(file="plots/Histogram_plot_by_sample.pdf")
ggplot(LengthsTable, aes(Length, color=Sample, fill=Sample)) + stat_bin(aes(y=..count../sum(..count..)), binwidth = 1) +
  facet_grid(Sample ~ .) +
  labs(title = "All Gene Tail Lengths", y = "Fraction", x="PAT Length") +
  xlim(0,300)
dev.off()

pdf(file="plots/Histogram_plot_by_genotype.pdf")
ggplot(LengthsTable, aes(Length, color=Genotype, fill=Genotype)) + stat_bin(aes(y=..count../sum(..count..)), binwidth = 1) +
  facet_grid(Genotype ~ .) +
  labs(title = "All Gene Tail Lengths", y = "Fraction", x="PAT Length") +
  xlim(0,300)
dev.off()
