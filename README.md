# metagenomics_halophytes
Miguel Camacho Sanchez, José M. Barcia-Piedras, Susana Redondo-Gómez, Maria Camacho (2020) Mediterranean seasonality and the halophyte Arthrocnemum macrostachyum determine the bacterialcommunity in salt marsh soils in Southwest Spain. Applied Soil Ecology
---
Supplementary methods
author: Miguel Camacho-Sánchez
date: 10/4/2018
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE,eval=F)
```

## 1. Trimming of sequences

Cutadapt 1.8.3 was used to trim F341 (CCTACGGGNGGCWGCAG) and R785 (GACTACHVGGGTATCTAATCC) primers (Klindworth et al. 2013) from R1 and R2 reads respectively.
List all files in directory `ls`:
```{bash, eval = T, include=T,echo=F}
ls | grep R[1-2].fastq
```

Run loop to (1) trim forward primer anchored to 5' end in R1 reads and reverse primer anchored to 5' end in R2 reads, and (2) write to ouput with only paired sequences for which both primers were found:
```{bash eval=FALSE,include=T}
fastq1=($(find *1.fastq))
fastq2=($(find *2.fastq))
long=($(echo ${#fastq1[@]}))
for i in `seq 0 $((long-1))`
do
cutadapt -g ^CCTACGGGNGGCWGCAG -G ^GACTACHVGGGTATCTAATCC --trimmed-only -o ${fastq1[i]%%.*}_cutadapt.fastq -p ${fastq2[i]%%.*}_cutadapt.fastq ${fastq1[i]} ${fastq2[i]}
done >cutadatp_filtering.txt
```

## 2. Determination of Amplicon Sequence Variants (ASVs)

We followed the strategy described in [Callahan et al. 2016](http://dx.doi.org/10.12688/f1000research.8986.2) in its latest [update](https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html).

Load libraries into R:
```{r echo=T, eval=T, include=T,message=FALSE}
library("knitr")
#Start session
.cran_packages <- c("ggplot2", "gridExtra")
  #load libraries
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
  # Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
```

Set `path` to path directory with trimmed sequences from cutadapt:
`path<-"path_to_sequences"`

Store paths for sequences to be read and for filtered sequences:
```{r eval=F,include=T}
fnFs <- sort(list.files(path, pattern="_R1_cutadapt.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cutadapt.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#create path for sequences that will be filtered downstream
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
```

Check quality plots for R1 and R2:
```{r eval=F}
qualplt<-plotQualityProfile(c(fnFs,fnRs),n=5000,aggregate=F)#for 5000 random reads
```

```{r eval=T,fig.cap="Figure S1.1. Quality plot"}
qualplt
```

Based on the quality plots we decided to truncate R1 reads to 275 nt and R2 reads to 180 nt. Considering the amplicon length is ~425 nt (without primers), we would still have have ~30 nt of overlap for merging R1 and R2, downstream. Due to the lower quality of R2, we set less stringent parameters for expected maximum errors (maxEE=6).
```{r eval=F}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,maxN=0,maxEE=c(3,6), truncQ=2,compress=TRUE, multithread=TRUE,truncLen=c(275,180))
```

Around 75% of the paired reads passed the filter:
```{r, eval=T}
cbind(out,"Retained"=out[,2]/out[,1])
```

In dada2 the amount of computing power needed growths quadratically with sequencing depth, [big data](https://benjjneb.github.io/dada2/bigdata.html). For that reason, we decided to calculate the error matrices from a random subsample of 10,000 reads from each sample:
```{r eval=F}
library(ShortRead)
n=10000#seq to subsample
filtFs10000<-NA#empty vector to paste filepath for R1 reads
for(i in 1:length(filtFs)){#loop to subsample n seq
  sampler<-FastqSampler(filtFs[i],n)
  temp<- yield(sampler)
  pathx<-file.path(filt_path,paste0(sample.names[i],"_F_filt_10000_fastq.gz"))#name for new file
  writeFastq(temp,pathx)
  filtFs10000[i]<-pathx#filepath for newly created files
}
#the same loop for R2 reads
filtRs10000<-NA
for(i in 1:length(filtRs)){
  sampler<-FastqSampler(filtRs[i],n)
  temp<- yield(sampler)
  pathx<-file.path(filt_path,paste0(sample.names[i],"_R_filt_10000_fastq.gz"))
  writeFastq(temp,pathx)
  filtRs10000[i]<-pathx
}
```

Learn the error rates for R1 and R2 reads independently:
```{r, eval=F}
errF <- learnErrors(filtFs10000, multithread=TRUE)
errR <- learnErrors(filtRs10000, multithread=TRUE)
```

Errors increase with low-quality:
```{r, eval=T, warning=F,fig.cap="Figure S1.2. Error plots"}
plotErrors(errF, nominalQ=TRUE)
```

Derreplication removes sequences which are the same thus lowering the computing burden for downstream analysis:
```{r, eval=F}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

Infer ASVs for each sample:
```{r eval=F}
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

Merge paired reads. Spurious sequence variants are further reduced by merging overlapping reads:
```{r eval=F}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,minOverlap = 20,maxMismatch=0)
```

Construct sequence table (analogous to an OTU table) from the provided list of samples:
```{r eval=F}
seqtab <- makeSequenceTable(mergers)
```
Number of samples and variants:
```{r eval=T}
dim(seqtab)
```
Distribution of sequence lengths:
```{r eval=T}
table(nchar(getSequences(seqtab)))
```

Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing:
```{r eval=F}
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(387,431)]
```

Remove bimeras:
```{r eval=F}
seqtabNoC <- removeBimeraDenovo(seqtab2)
```
Proportion of sequences kept:
```{r eval=T}
seqkept<-sum(seqtabNoC)/sum(seqtab2)
```
Proportion of variants kept:
```{r eval=T}
varkept<-ncol(seqtabNoC)/ncol(seqtab2)
```
After removing bimeras `r seqkept` of the sequences were kept, which implied `r varkept` of the variants.

###2.1 Track reads through the pipeline
```{r eval=T}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtabNoC))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
```
```{r eval=T}
track
```
Variants per sample:
```{r eval=T}
apply(seqtabNoC,1,function(x){sum(x>0)})
```

###2.2 Assign taxonomy
Assign taxonomy with formatted [Silva v128 database](https://zenodo.org/record/801832#.WoVUHZOdW34):
```{r eval=F}
fastaRef <- file.path("path_to/silva_nr_v128_train_set.fa.gz")
taxTabSilva <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)#assigns taxonomy down to Genus
taxTabSilva_esp <- addSpecies(taxTabSilva, "path_to/silva_species_assignment_v128.fa.gz")#adds a Species column
```

###2.3 Taxonomy filtering

Phyloseq package allows integrating taxonomy table, ASVs table, metadata, phylogenetic trees into an unique phyloseq object for further manipulation.
Import metadata:

```{r, eval=F}
meta<-read.csv("path_to/metadata.txt",colClasses=c("ID"="character"))
```
Do all samples in the ASV table (seqtabNoC) are present in the metadata?
```{r eval=T}
all(rownames(seqtabNoC) %in% meta$ID)
rownames(meta) <- meta$ID
```
Create phyloseq object:
```{r eval=T}
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxTabSilva_esp))
ps
```
Show taxonomic ranks in the dataset:
```{r eval=T}
rank_names(ps)
```

ASVs per Kingdom and Phylum:
```{r eval=T}
table(tax_table(ps)[, "Kingdom"], exclude = NULL)
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```
Remove ASVs with NO Phylum assigned and only keep Kingdom Bacteria:
```{r eval=T}
ps <- subset_taxa(ps, !is.na(Phylum))
ps <- subset_taxa(ps, Kingdom=='Bacteria')
```
###2.4 Negative control

After having a look at the ASVs in the negative control they seem to have the most abundant variants in their composition (compare to Figxxx) likely as a result of  some cross contamination. Furthermore, the its low number of sequences (see "track" object) and variants seem to point to the fact that contamination has not been an issue. Therefore the risk of removing important biological diversity in the metacommunity by removing the ASVs in the negative control does not compensate for the possible removal of contaminant variants (https://github.com/benjjneb/dada2/issues/114). For that reason the negative control was directly pruned from the pyloseq object without any further processing:
```{r eval=T, include=T,warning=F,fig.cap="Figure S1.3. per-Phylum ASV prevalence for the negative control"}
#Function to have a look at the abundance plots per phylum(also crates prevalence table [[1]],and plots[[3]])
plot_abundances<-function(phylo){
  prev = apply(X = otu_table(phylo),
               MARGIN = ifelse(taxa_are_rows(phylo), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prev = data.frame(Prevalence = prev,
                    TotalAbundance = taxa_sums(phylo),
                    tax_table(phylo))
  
  ##Get prevalence per phylum
  
  prevalence_per_phylum<-plyr::ddply(prev, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
  
  ##Plot abundances vs prevalence
  
  plot1<-ggplot(prev, aes(TotalAbundance, Prevalence / nsamples(phylo),color=Phylum)) +
    # Include a guess for parameter
    geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")
  return(list(prev,prevalence_per_phylum,plot1))
}

ps_neg<-prune_samples(sample.names=="negativo",ps)#create phyloseq object with the negative control
neg_summary<-plot_abundances(ps_neg)
neg_summary[[3]]

ps<-prune_samples(sample_names(ps) != "negativo",ps)#remove negative sample from dataset
```

###2.5 Abundance/prevalence filtering
Secondly, we set filters based on abundances and prevalence.
We evaluated replicas to estimate abundance/prevalence thresholds for ASVs.
Count in how many samples all given ASVs are present (prevalence) and plot in an histogram the abundance of ASVs with Prevalence = 1 vs 2.
```{r eval=T,include=T, message=F, fig.cap="Figure S1.4. Shared vs non-shared ASV for replicas"}
library(data.table)
ps_rep<-subset_samples(ps,ID %like% "L2-3")#keeps only samples L2-3 and L2-3rep.

prevdf_rep = apply(X = otu_table(ps_rep),#creates vector with prevalence
                MARGIN = ifelse(taxa_are_rows(ps_rep), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf_rep = data.frame(Prevalence = prevdf_rep,#adds abundance and prevalence vector to data frame
                     TotalAbundance = taxa_sums(ps_rep),
                     tax_table(ps_rep))

# Histogram Grey Color
hist(log10(prevdf_rep$TotalAbundance[prevdf_rep$Prevalence=="2"]), col=rgb(1,0,0,0.5),main="L2-3/L2-3-rep; Prevalence 1, blue; Prev. 2, red",xlab="Log10 Abundance")#plot the abundance for varians with a prevalence of 2

hist(log10(prevdf_rep$TotalAbundance[prevdf_rep$Prevalence=="1"]), col=rgb(0,0,1,0.5), add=T)#plot the abundance for variants with a prevalence of 1
rm(ps_rep)
rm(prevdf_rep)
```
The histograms seem to indicate that for ASVs with an abundance >100 (log10(100)=2), the chances of having variants only present in one of the replicas is very low.

We remove replicas L2-3-rep and L1-2a.
```{r eval=T, include=T}
ps<-prune_samples(!(sample_names(ps) %in% c("L1-2","L2-3-rep")),ps)
ps
```

Based on the histogram above we decided to remove all ASVs for which the product of the abundances for the samples in which they were present was below 100:
```{r eval=T, include=T}
prevdf = apply(X = otu_table(ps),#vector with prevalence
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
mult = apply(X = otu_table(ps),#vector with product of all non-0 abundances
              MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
              FUN = function(x){prod(x[which(x > 0)])})
prevdf = data.frame(Prevalence = prevdf,#data frame with above information
                    TotalAbundance = taxa_sums(ps),
                    Product = mult,
                    tax_table(ps))

#Set abundance/prevalence threshold for filter:
prev_abundance_Threshold = 100#abundance threshold

#Execute filter:
keepTaxa = rownames(prevdf)[(prevdf$Product > prev_abundance_Threshold)]
ps_filt = prune_taxa(keepTaxa, ps)
```
### 2.6 Remove mitochondrial and chloroplastic variants:
```{r eval=T,include=T}
ps_filt <- subset_taxa(ps_filt, !Class %in% "Chloroplast" & !Family %in% "Mitochondria")
ps_filt
```
Plot ASVs abundances per Phylum:
```{r eval=T,include=T,fig.cap="Figure S1.5.per-Phylum ASV prevalence"}
plot_abundances(ps_filt)[[3]]
```

Final filtered sequences and variants per sample. There seems to be no correlation between number of sequences and number of variants:
```{r eval=T,include=T}
final_Sequences=sample_sums(ps_filt)
final_Variants=apply(ps_filt@otu_table,1,function(x){sum(x>0)})
data.frame(Final_seq=final_Sequences,Final_Variants=final_Variants)
```
```{r eval=T,include=T,fig.cap="Figure S1.6.Correlation between number of variants and number of sequences"}
plot(final_Sequences,final_Variants)
```

(Why not rarefy?: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)

## 3. Analysis of diversity
###3.1 Merge phylogenetic tree to phyloseq object

Export fasta from variants:
```{r eval=F,include=T}
seqs <- getSequences(otu_table(ps_filt)@.Data)
names(seqs) <- seqs # names variant  as its sequence

library(ShortRead)
writeFasta(DNAStringSet(seqs),"path_to/seq.fa")
```

Align-reconstruct tree iteratively with PASTA:

```{bash eval=F,include=T}
python run_pasta.py \
-i seq.fa \
-d dna \
--aligner=mafft \
--merger=OPAL \
--tree-estimator fasttree \
--num-cpus=2 \
--iter-limit=3 \
-j metarizo \
-o pasta_folder
# Note that the final alignment outputted by PASTA is NOT masked, but masked versions of the output are also saved as temporary files.
```
Open "tree.tre" in the text editor and remove lines 2 and 3.
Import resulting tree from PASTA into phyloseq object:
```{r eval=F, include=T}
tree_metarizo<-phyloseq::read_tree("path_to_tree.tre")
tree_metarizo<-midpoint(tree_metarizo)#root tree
```

Merge phylogenetic tree to phyloseq object:
```{r eval=T,include=T}
ps_filt_t<-merge_phyloseq(ps_filt,tree_metarizo)
```

Add sample names:
```{r eval=T}
sample_names(ps_filt_t)<-c("jul1","jul2","jul3","jul4","oct1","oct2","oct3","oct4")
```

###3.2 Estimate alfa diversity

```{r,eval=T,include=T, fig.cap="Figure S1.7. Alpha diversity"}
alpha_div<-estimate_richness(ps_filt_t,measures=c("Observed","Shannon"))
alpha_div<-data.frame(alpha_div,alpha_div$Shannon/log(alpha_div$Observed))
names(alpha_div)<-c("Observed","Shannon","Pielou")
alpha_div
```
```{r eval=T,include=T, }
plot_richness(ps_filt_t,x="sample", measures=c("Shannon", "Observed"), color="rhizosphere",shape="site") + theme_bw()
```

####3.2.1 Change of alpha diversity with agglomeration
We agglomerated variants using cophenetic distances from the tree to check how diversity changes with agglomeration and to detect potential overstimation of ASVs, which would be evident after shallow agglomeration.
```{r eval=T,include=T,warning=F}
h=0
shannon_div<-data.frame(c(rep(NA,nsamples(ps_filt_t))))
observed_div<-data.frame(c(rep(NA,nsamples(ps_filt_t))))

glom0<-ps_filt_t
for (i in 1:8){
  assign(paste0("glom",i),tip_glom(get(paste0("glom",i-1)),h=h[i]))
  h[i+1]=0.05+h[i]
  shannon_div[,i]<-estimate_richness(get(paste0("glom",i)),measures=c("Shannon"))$Shannon
  observed_div[,i]<-estimate_richness(get(paste0("glom",i)),measures=c("Observed"))$Observed
  names(shannon_div)[i]<-paste0("glom",i,"_",h[i-1])
  names(observed_div)[i]<-paste0("glom",i,"_",h[i-1])
}
```
```{r eval=T,include=T, fig.cap="Figure S1.8. Change of alpha diversity with agglomeration (h1, tree distance)"}
h1<-head(h,-1)
par(mfrow=c(1,2))
plot(h1,as.numeric(shannon_div[1,]),type="l",ylim = c(3,max(shannon_div)),ylab="Shannon diversity")
for (i in 2:8){
  lines(h1,shannon_div[i,],col=i)
}
legend(0.2,6,legend=sample_names(ps_filt_t),col=c(1:length(h1)),lty=1,bty="n",cex=0.6)

plot(h1,as.numeric(observed_div[1,]),type="l",ylim = c(3,max(observed_div)),ylab="Observed diversity")
for (i in 2:8){
  lines(h1,observed_div[i,],col=i)
}
legend(0.2,800,legend=sample_names(ps_filt_t),col=c(1:length(h1)),lty=1,bty="n",cex=0.6)
par(mfrow=c(1,1))
```

####3.2.2 Create Qiime-like OTU tables

Transform abundance to relative values:
```{r eval=T,include=T}
ps_filt_tra = transform_sample_counts(ps_filt_t, function(x){x / sum(x)})
```
Create Qiime-like OTU tables for each taxonomic level:
```{r eval=F,include=T}
taxa_levels<-colnames(tax_table(ps_filt_tra)@.Data)[2:6]#names of taxonomic levels
Qotu_tables<-list()#empty list
library(xlsx)
write.xlsx("Qiime-like OTU tables", file="path_to/otu_tables.xlsx",row.names=F,col.names = F)#empty Excel file

for (i in 1:length(taxa_levels)){#loop
  ps3 = tax_glom(ps_filt_tra, taxa_levels[i], NArm = F)
OTUs<-t(otu_table(ps3)@.Data)
taxa<-tax_table(ps3)@.Data[,c(1:(1+i))]
level_tax<-assign(paste0("qiimeOTU_",taxa_levels[i]),cbind(taxa,OTUs))
Qotu_tables[[i]]<-level_tax
write.xlsx(level_tax,file="path_to/otu_tables.xlsx",sheetName=paste0(taxa_levels[i]),row.names=FALSE, append=TRUE)
}
```

Proportion of ASVs assigned to different taxonomic levels:
```{r eval=T}
total_taxa<-ntaxa(ps_filt_t)
nClass<-ntaxa(subset_taxa(ps_filt_t, !is.na(Class)))
nOrder<-ntaxa(subset_taxa(ps_filt_t, !is.na(Order)))
nFamily<-ntaxa(subset_taxa(ps_filt_t, !is.na(Family)))
nGenus<-ntaxa(subset_taxa(ps_filt_t, !is.na(Genus)))
data.frame("nClass"=nClass/total_taxa, "nOrder"=nOrder/total_taxa, "nFamily"=nFamily/total_taxa,"nGenus"= nGenus/total_taxa)
```

####3.2.3 Most abundant phyla

Plot cumulative abundance per Phylum:
```{r eval=T, include=T}
#calculate culumative abundance per phylum
phylum_ps<-tax_glom(ps_filt_t, "Phylum")#agglomerate abundances per phylum
phylum_psm<-psmelt(phylum_ps)#melt ps
sum_phylum<-tapply(phylum_psm$Abundance,phylum_psm$Phylum,sum)/sum(phylum_psm$Abundance)#sum abundances per phylum
sum_phylum<-sort(sum_phylum)#sort
cum_sum_phylum<-cumsum(sum_phylum)#calculate cummulative sum

#calculate diversity per phylum
asv_phyla1<-table(tax_table(ps_filt_t)[, "Phylum"], exclude = NULL)#ASV per phylum
asv_phyla<-as.numeric(asv_phyla1)/sum(asv_phyla1)#sort and make relative
names(asv_phyla)<-names(asv_phyla1);rm(asv_phyla1)

cum_div_abun<-merge(data.frame("rel_ASV"=asv_phyla),data.frame(cum_sum_phylum),by=0,all=T)#combine diversity and abundance
cum_div_abun<-data.frame(cum_div_abun,row.names = 1)
cum_div_abun<-merge(cum_div_abun,data.frame("rel_abund"=sum_phylum),by=0,all=T)
cum_div_abun<-cum_div_abun[order(cum_div_abun[,3]),]#order by cumulative abundance
cum_div_abun<-cbind(cum_div_abun,"ASV"=cumsum(cum_div_abun$rel_ASV))
cum_div_abun<-data.frame(cum_div_abun,row.names = 1)
cum_div_abun
```
```{r eval=T, include=T, fig.cap="Figure S1.9. Most dominant/diverse phyla"}
#plot cumulative abundance/diversity
plot(cum_div_abun$cum_sum_phylum,type="n",ylab="Cumulative abundance/diversity per phylum",xlab="Phyla",yaxt="n")#Plot
axis(2, las=2)
segments(c(1:length(cum_div_abun$cum_sum_phylum)),c(0,head(cum_div_abun$cum_sum_phylum,-1)),c(1:length(cum_div_abun$cum_sum_phylum)),cum_div_abun$cum_sum_phylum,lwd=10,col = "grey" ,lend=1)#add phyla
abline(h=0.1,lty=2)#0.1 line
text(c(20:26),cum_div_abun$ASV[20:26],rownames(cum_div_abun)[20:26],pos=2)#phyla names

segments(c(1:length(cum_div_abun$ASV)),c(0,head(cum_div_abun$ASV,-1)),c(1:length(cum_div_abun$ASV)),cum_div_abun$ASV,lwd=5,col = 1,lend=1)#overplot
```

###3.2 Estimate beta diversity
Transformation of the abundances to natural logarithmic scale:
```{r eval=T,include=T}
ps_filt_tlog <- transform_sample_counts(ps_filt_t, function(x) log(1 + x))
```

MDS with weighted UNIFRAC distances:
```{r eval=T, include=T,fig.cap="Figure S1.10. Beta diversity: MDS from weighted UNIFRAC distances"}
unif_meta <- ordinate(ps_filt_tlog, method = "MDS", distance = "wunifrac")
evals <- unif_meta$values$Eigenvalues
library("ggrepel")#has function to "repel" overlaping labels in plot
plot_ordination(ps_filt_tlog, unif_meta,axes = c(1,2),shape="site", color = "rhizosphere") +
  labs("site","rhizosphere") +
  geom_point(size=5)+
  coord_fixed(sqrt(evals[2] / evals[1]))+
  geom_text_repel(aes(label=sample_names(ps_filt_tlog)),size=4,min.segment.length=1,point.padding=0.5)

unif_meta$values
```

Exploration of other dimensions. The scaling of y axis is proportional to the variance explained:
```{r eval=T, include=T,fig.cap="Figure S1.11. Beta diversity: MDS from weighted UNIFRAC distances for dimensions 3 to 5"}
for (i in 3:5){
  assign(paste0("plot",i),
  plot_ordination(ps_filt_tlog, unif_meta,axes = c(1,i),shape="site", color = "rhizosphere") +
  labs("site","rhizosphere") +
  geom_point(size=5)+
  coord_fixed(sqrt(evals[i] / evals[1]))+
  geom_text_repel(aes(label=sample_names(ps_filt_tlog)),size=4,min.segment.length=1,point.padding=0.5))
}
library(grid)
library(gridExtra)
grid.arrange(plot3,plot4,plot5)
```

##4. Differential abundances of ASVs between conditions

Determine differential abundant ASVs per condition:
```{r eval=T, message=F,include=T,fig.cap="Figure S1.12. Dispersion estimates for month, rhizosphere and site"}
library(dplyr)
library(DESeq2)

variable<-c("month","rhizosphere","site")#vector with variables to determine differential abundances

gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}#geometric means prior to estimate size factors
par(mfrow=c(1,3))
tax <- tax_table(ps_filt_t) %>% data.frame();tax$seq <- rownames(tax);tax$uniqueID<-as.character(1:nrow(tax))#dataframe with taxonomic data
tax$compact_names<-paste(tax$Phylum,paste0("O.",tax$Order),sep="-",paste0(tax$Genus,tax$Species),tax$uniqueID)%>%gsub(pattern = "NA",replacement = "")
options(digits=3)#for digits of p-values
otut<-t(otu_table(ps_filt_t))%>%data.frame();otut$seq<-rownames(otut)#dataframe with abundance data

for (i in 1:length(variable)){
meta_deseq<-phyloseq_to_deseq2(ps_filt_t, as.formula(paste("~ ",variable[i],sep = "")))###phyloseq to deseq2
geoMeans = apply(counts(meta_deseq), 1, gm_mean)
meta_deseq = estimateSizeFactors(meta_deseq, geoMeans = geoMeans)
meta_deseq = DESeq(meta_deseq)

plotDispEsts(meta_deseq)

sigtab<-as(results(meta_deseq), "data.frame"); sigtab<-cbind("seq"=rownames(sigtab),sigtab,stringsAsFactors=F)

#ASV with most contribute to differences among the given condition
taxa_diff<-sigtab[which(sigtab$padj<0.05),]%>%arrange(padj)
if (all(taxa_diff$seq %in% phy_tree(ps_filt_t)$tip.label,T)){
  
 ##join taxonomic assignations + frequencites
merged_results<-left_join(taxa_diff,tax,by=c("seq","seq"))%>%left_join(otut,by=c("seq","seq"))
assign(paste("merged_results",variable[i],sep="_"),merged_results) 
}
}
par(mfrow=c(1,1))
```


Write result tables from differential abundances to Excel file:

```{r eval=F,include=T}
pathxls<-"path_to/differential_abundances.xlsx"
```
```{r eval=T, include=T}
write.xlsx("Differential_abundances", file=pathxls,row.names=F,col.names = F)#empty Excel file
for (i in 1:length(variable)){
  write.xlsx(get(ls(pattern = "merged_results_")[i]),file=pathxls,sheetName=paste0(variable[i]),row.names=FALSE, append=TRUE)
  }
```

```{r eval=T, include=T}
x=c();for (i in 1:length(variable)){x<-c(x,get(ls(pattern = "merged_results_")[i])$compact_names)}
ids<-unique(x)
```

Plot tree with non-significant ASVs pruned and corresponding fold changes (see Figure 6 from main text):
```{r eval=F, include=T}
library(phytools)
tree<-phy_tree(ps_filt_t)
x=c();for (i in 1:length(variable)){x<-c(x,get(ls(pattern = "merged_results_")[i])$seq)};ids<-unique(x)#determine ASVs to retain

to_remove<-tree$tip.label[-match(ids, tree$tip.label)]
pruned.tree<-drop.tip(tree,to_remove)#remove non-significant taxa from tree
pruned.tree$tip.label<-tax$compact_names[match(pruned.tree$tip.label,tax$seq)]

nt <- normTransform(meta_deseq) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)
matching<-match(tax$seq[which(tax$compact_names %in% pruned.tree$tip.label)],rownames(log2.norm.counts))

log2.norm.counts<-log2.norm.counts[matching,]
rownames(log2.norm.counts)<-tax$compact_names[matching]

#fix(phylo.heatmap)#add a small sepparation between labels and image (START+.1, END+.1), change "value" label of bar to "log2 abundance",  return START value and replace "phylogram" by a modified version "phylogram2",wich returns to the plotting environment a dataframe with the tiplabels arranged by y, so it can later be used to map other data to the tree.

phylo.heatmap(pruned.tree,log2.norm.counts,tip.labels=FALSE,fsize=c(0.3,1,1),pts=F)
lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
ordered_tips<-lastPP$label.tips
ordered_tips<-arrange(desc(ordered_tips))

for (i in 1:length(variable)){#add squares with transparency correlated to their p-values
  if(i==1){xx=0}
  merged_results<-get(ls(pattern = "merged_results_")[i])
  matching_tips<-merged_results$compact_names%>%match(ordered_tips$name)
  color_pvalues<-merged_results$padj[order(matching_tips,decreasing = T)]*-20+1
  points(x =rep(START-xx,Ntip(pruned.tree))[matching_tips],y=sort(ordered_tips$y[matching_tips]),pch=15,cex=0.6,col=hsv((i/3)-1/3,alpha=color_pvalues))
  xx<-xx+0.025
}
```

Abundances of the significant AVSs:
```{r eval=T, include=T}
melted_ps<-psmelt(ps_filt_t)
melted_ps[melted_ps$OTU==taxa_diff$seq[1],]
```
