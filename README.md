# Pop-Seq
Identification of QTL for BSA (Bulk Degregant Analysis) and for diverse populations. 

Pheotype based bulks are alternative to population based approaches like Bi-Parental QTL mapping and Association mapping.
SNP frequency differences of these bulks helps identify candidate regions/genes. FIXED sites offer the best resolution for identification of QTL.

There are 4 snp types based on the frequency in each of the bulk.
Fixed site - Allele frequency is 100% for each of the pop, but Alternate alleles in each of the population.
Shared site - Both alleles present in both the populations.
Unique_pop (Uniq_pop1 and Uniq_pop2) - Allele frequency of one allele is 100% in only one population, but other population has both alleles. 


Usage:
The first awk code reads a vcf file (SNP from bulks) and determines the site type. The output is then read into R to create a distribution using sliding windows either for all sites or specific set or based on delta 
(difference in frequencies of the two bulks)

### Step 1: Determine SNP Type
Uncompressed file : cat my.vcf | awk -f identifySnpTypePools.awk > my.vcf.sitetype
bz2 compressed file : bzip2 -dck my.vcf.bz2 | awk -f identifySnpTypePools.awk > my.vcf.sitetype

### Step 2:read the awk input from previous awk step
setwd("/Users/myname/myAnalysis")
If Fixed sites (Defaults 10k/2k): sliding_windows_frequency("my.vcf.sitetype",type="fixed")
If all Sites and 100k/10k slides: sliding_windows_frequency("my.vcf.sitetype",window=100000, slide=10000, type="all")

### Step 3: determine cutoffs for the site type using boostraps (Output from Above)
mydata <- read.table("my.vcf.sitetype_10000_2000_freq_dist.txt", header=TRUE)
source("/Users/myname/myAnalysis/R_functions")

nsamples <- sapply(1:1000,function(i) sample(mydata$S, replace =T))
cutoffs <- quantile_cutoffs(nsamples, probs = c(0.95,0.99,0.999),myname="Sites_of_Imp")
rm(nsamples)

#### Extract values form bootstrap table
sign95 <- cutoffs[1,2]
sign99 <- cutoffs[2,2]

### Step 4: Manhattan plot using qqman package
if("qqman" %in% rownames(installed.packages()) == FALSE)
{install.packages("qqman")} else
{library(qqman)}

par(mar=c(1,1,1,1))
mydata$chr = as.numeric(gsub("Chr","",mydata$Chrom))

png("Sites_distribution.png")
manhattan(mydata, chr = "chr", bp = "Start_Mbp", p = "S", snp = (paste0("Chrom","Start_Mbp", sep="_")),
col = c("dodgerblue", "orange"), chrlabs = c("1","2","3","4", "5", "6","7","8","9","10","11"),
suggestiveline = sign95, genomewideline = sign99,
highlight = NULL, logp = FALSE,main="Sites Distribution",
xlab ="Chromosome", ylab ="Frequency", cex=0.2)
dev.off()



