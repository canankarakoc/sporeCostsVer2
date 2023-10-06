# Spore Costs v2
# 6 October 2023 - last update
# Author: C. Karakoc
# Global expression data from SporeWeb: https://sporeweb.molgenrug.nl/
# Alternative expression data for validation: https://doi.org/10.3390/ijms22179345
# Germination data: https://doi.org/10.3390/ijms232113614
# Germination expression data: https://doi.org/110.1128/mSphere.00463-20   
# Single gene deletion library from: https://doi.org/10.1016/j.cels.2016.12.013
# Global protein abundance data: https://pax-db.org/
# List of gene categories & annotation: https://subtiwiki.uni-goettingen.de/
# Protein sequence: Uniprot https://www.uniprot.org/taxonomy/224308
# Amino acid costs: https://doi.org/10.1073/pnas.1701670114

######################
# Packages & Plotting 
######################
##########################################################################
library(tidyverse)       
library(patchwork) # for nls
library(ggpmisc) # for the formulas of fits
library(stringr) 
#devtools::install_github("dgrtwo/fuzzyjoin")
library(fuzzyjoin) # for merging protein sequences
# ggplot theme
library(ggsci) # nature publishing colors 

mytheme <- theme_bw()+
  theme(axis.text = element_text(size=16, color = "black"),
        axis.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=14),
        legend.background = element_blank(),
        legend.title = element_text(size=14,face="bold"),
        plot.title = element_text(size=16, face="bold", hjust=0.5),
        strip.text = element_text(size=16, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.length = unit(-2, "mm"))

# Color palette inspired by plots in Nature Reviews Cancer

paletteNature <- c("Cinnabar" = "#E64B35", "Shakespeare" = "#4DBBD5",
  "PersianGreen" = "#00A087", "Chambray" = "#3C5488",
  "Apricot" = "#F39B7F", "WildBlueYonder" = "#8491B4",
  "MonteCarlo" = "#91D1C2", "Monza" = "#DC0000",
  "RomanCoffee" = "#7E6148", "Sandrift" = "#B09C85")

setwd("~/GitHub/sporeCostsVer2")

# Data 
######################
# Data
######################
##########################################################################
# Gene&protein length from SubtiWiki
annotationData <- read.table("./otherData/subtiwiki.gene.export.2022-05-11.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Regulator list - sigma factors #Daniel
sigma          <- read.table("./otherData/regulations.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Data with synonyms # different sources
nameMap        <- read.table("./otherData/nameMap.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Protein abundances from PaxDB
protAbun       <- read.table("./otherData/protAbunData.csv", sep = ',', dec = ".", header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Protein sequences from Uniprot
protSeq        <- read.delim("./otherData/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.02.08-14.23.07.59.txt")
#protSeq       <- read.table("./otherData/uniprotkb_proteome_UP000001570_2023_09_01.csv", sep = ",", header = T)

# Nucleotide & Amino acid costs 
aaCosts        <- read.table("./otherData/aaCosts_pnas.1701670114.sd01.csv", sep = ",", dec = "." , header = T)
nucCosts       <- read.table("./otherData/nucleotideCosts_pnas.1701670114.sd01.csv", sep = ",", dec = "." , header = T)

# Single gene deletion library from Koo et al., 2018
deletionLib    <- read.table("./mutantLibraryData/mmc7_Koo_etal_2017_TableS6A.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))
mutantLib      <- read.table("./mutantLibraryData/mmc7_Koo_etal_2017_TableS6B.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))
conservedLib   <- read.table("./mutantLibraryData/mmc7_Koo_etal_2017_TableS6E.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Other traits from SubtiWiki lists 
otherTraits    <- read.table("./otherData/traits.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Germination time course - Gao et. al. 2022
germination    <- read.table("./germinationData/proteinTimeCourseGaoEtAl2022.csv", sep = ",", header = T) 

# Germination time course - Swarge et al. 2020 
germination1   <- read.table("./germinationData/Swarge_DifferentiallyExpressedGenes.csv", sep = ",", header = T)
germination2   <- read.table("./germinationData/Swarge_DifferentialProteinExp.csv", sep = ",", header = T)
germination3   <- read.table("./germinationData/Swarge_NonDifferentiallyExpressedProt.csv", sep = ",", header = T)
germination4   <- read.table("./germinationData/Swarge_ProteinScore.csv", sep = ",", header = T) 
germination5   <- read.table("./germinationData/Swarge_SporeNewproteinRatio_geoMean.csv", sep = ",", header = T) 
germination6   <- read.table("./germinationData/SwargeEtAl_onlyNewProteins.csv", sep = ",", header = T) 
  
# Spore time course alternative test
sporulation1   <- read.table("./otherData/TuEtAll_Sporulation_ProteinExpression.csv", sep = ",", header = T) 

# Sporulation efficiency
efficiency     <- read.table("./otherData/efficiency.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

                           ###########
###########################SPORULATION#####################################
                           ###########
#####################
# Time series prep
#####################
##########################################################################

# Expression data from SporeWeb
files  <- list.files(path = "./SporeWebHeatmaps/" , pattern = "*.csv", full.names = T)
exp_files <-  list()
for (i in 1:length(files)){
  exp_files[[i]] <- read.table(files[i], header = T, sep = ",", dec = ".")
}
merged_exp                   <- bind_rows(exp_files, .id = 'sourceID')
merged_exp$sourceID          <- as.factor(merged_exp$sourceID)
levels(merged_exp$sourceID ) <- c("1.vegetative", "2.starvation", "3.onset", "4.commitment", "5.engulfment")

# Wrangle 
expressionLong <- merged_exp[,c(1, 3:5, 7:14)] %>%
  pivot_longer(cols = c('t1','t2','t3','t4','t5','t6','t7','t8'),
               names_to =  "time", values_to = "expression") %>%
  group_by(time, sourceID, regulators, gene, locus_tag) %>%
  summarise(meanexp = mean(expression)) %>%
  ungroup() %>%
  filter(meanexp > 0) %>%
mutate(time_h = gsub('t', '', time)) 


# Abundance data prep
#####################
# Protein abundance
#####################
##########################################################################

# Protein abundance data does not include locus tags
# It is often easier to merge data sets with locus tags, because genes have a lot of synonyms

# A gene name assigned two different locus tag, corrected based on PaxDB
which(nameMap$gene1 =="nrgB")
nameMap[246,4] <- ""

# Logical map with matching gene names
MAP <- as.data.frame(t(apply(nameMap, 1, function(x) x %in% protAbun$gene)))
colnames(MAP) <- colnames(nameMap)

# Convert logical map into a column mapped with the matching column names
MAP$index = apply(MAP, 1, function(x) paste(names(x)[x], collapse=", "))
nameMap$index = MAP$index

# Match the column names by the rows
id <- lapply(seq(nrow(nameMap)), function(i) nameMap[i,nameMap$index[i]])

# Convert the list output to vector column, but prevent loosing NULLs, convert them
# First to NAs, so that we know which genes do not have abundance info
nameMap$gene <- do.call(rbind, lapply(id, function(x) if(is.null(x)) NA else x))
nameMap$gene <- ifelse(is.na(nameMap$gene), nameMap$geneP, nameMap$gene)

# Now use this collapsed matching column to merge abundances, genes and locus tag
# Merge data
mergedAbunData <- protAbun  %>%
  right_join(nameMap, by = "gene", multiple = "all") %>%
  left_join(annotationData[-2], by = "locus_tag", multiple = "all") %>%
  left_join(sigma[-4], by = "locus_tag", multiple = "all") %>%
  distinct(locus_tag, .keep_all =TRUE) %>%  #left join duplicates 
  dplyr::select(locus_tag, gene, regulator, abundance, protein_length, gene_length)

# Fill NAs with median values

# Average protein abundance length
protMed <- as.numeric(as.vector(protAbun$abundance))
median(protMed, na.rm = T) #17.6 #round

mergedAbunData$protein_length <- as.numeric(mergedAbunData$protein_length)
mergedAbunData$gene_length    <- as.numeric(mergedAbunData$gene_length)

# Median protein and gene length 
median(mergedAbunData$protein_length, na.rm = T) #254
median(mergedAbunData$gene_length, na.rm = T) #765

# This function takes a while to run 
protSeqTidy <- protSeq %>%
  regex_inner_join(mergedAbunData, by = "gene") %>%
  distinct(sequence, .keep_all = TRUE) #remove the duplicated rows bases on unique sequences
 
# Amino acid alphabet 
alphabet = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

seqCount <- protSeqTidy %>%
  rowwise() %>%
  reframe(aac = str_count(sequence, pattern = alphabet)) %>%
  mutate(symbol = rep(alphabet, times = 4243)) %>%
  mutate(gene = rep(1:4243, each = 20)) %>%
  left_join(aaCosts, by = "symbol") %>%
  mutate(aa_opportunity = aac*opportunity_costs, aa_direct = aac*direct_costs) %>%
  group_by(gene) %>%
  summarize(aa_opportunitySum = sum(aa_opportunity), aa_directSum = sum(aa_direct)) %>%
  mutate(locus_tag = protSeqTidy$locus_tag) %>%
  distinct(locus_tag, .keep_all = TRUE)
  
protSeqTidyAbun <-  protSeqTidy %>%
  left_join(seqCount[,-1], by = "locus_tag") %>%
  select("gene.y", "protID", "locus_tag", "abundance", "protein_length", "gene_length", "aa_opportunitySum", "aa_directSum") %>%
  distinct(locus_tag, .keep_all = TRUE)


# Accounting data prep
#####################
# Expression data
#####################
##########################################################################

# Merge with expression data
# Here I'll create two different data sets. One for calculating transcription
# costs, another for translation. Since proteins are degraded much slower, I
# I will account for only repolimerization costs of transcripts 
# Here genes are accounted once base on first appearance

mergedExpData_time_distinct <-  expressionLong %>%
  distinct(locus_tag, .keep_all = TRUE) %>% # genes are accounted only once
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

mergedExpData_time <-  expressionLong %>%
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))


# Replication costs (Whole genome) 
#####################
# Genome costs
#####################
##########################################################################
# Genome size = (https://www.nature.com/articles/36786)
# Opportunity 
# 2 * Genome size * (34+1)
genome_opp <- 2 * 4214810 * 35 #295036700
# Direct 
# 2 * Genome size * (11+2)
genome_dir <- 2 * 4214810 * 14 #118014680
# Total
genome_tot <- 295036700 + 118014680 #413051380

# Membrane costs 
#####################
# Membrane costs
#####################
##########################################################################

# Number of lipid molecules = Cellular membrane areas/head-group areas of membrane lipid molecules
# Head group a1 = 0.65 nm2 (Nagle and Tristram-Nagle 2000; Petrache et al. 2000; Kucerka et al. 2011).

# Thickness of the bilayer (h):
# The thickness of a bilayer is approximately twice the radius of the head-group area, which 0.5 nm in all cases,
# plus the total length of the internal hydrophobic tail domains (Lewis and Engelman 1983; Mitra et al. 2004), 
# generally 3.0 nm

# Summed bilayer area 4π[r2 + (r − h)2] 
# Average Bacillus diameter =  	0.87 µm (BNID100211)

# Outer #4πr^2
4*3.141593*(0.87*1000)^2 #9511487

# Inner membrane 4π(r − h)^2
4*3.141593*((0.87*1000)-3.5)^2 #9435112

# Sum 
9511487+9435112 #18946599 molecules 

# Cost of lipid head & tail 
# Opportunity = 212, Direct = 18

# Opportunity 
18946599*212 # 4016678988
# Direct 
18946599*18 # 341038782

4016678988+341038782 #4357717770

# 50% discount, because of proteins
4357717770/2 #2178858885

# Septum is 1/6 of the total membrane
2178858885/6 #363143148

lipid_opp <- (4016678988/2)/6
lipid_dir <- (341038782/2)/6
lipid_tot <- lipid_opp+lipid_dir

# Replication costs of expressed genes
#####################
# Spore Replication 
#####################
##########################################################################
# 878 genes
# 4429 genes whole genome
sum(mergedExpData_time_distinct$gene_length.filled)
# total length = 722878
# genome length = 4214810
# 722878/4214810 %17
# 878/4429 %20

# Opportunity costs 
sporeRep  <- mergedExpData_time_distinct %>%
  mutate(opportunity = 2*gene_length.filled*35) %>% #2 = doublestring  35 = nucleotide costs 
  mutate(direct = 2*gene_length.filled*14)

# Sum 
sporeRepSum <- sporeRep %>%
  summarise(sumOpp = sum(opportunity, na.rm =T), sumDir = sum(direct, na.rm =T))
sporeRepTotal <- sporeRepSum$sumOpp+sporeRepSum$sumDir
# opportunity 50601460
# direct 20240584
# total 70842044

# percentage compared to total genome
70842044/413051380  #17%

# Transcription costs
#####################
# Spore Transcription
#####################
##########################################################################

# 1 mRNA can yield to 100-1000 proteins (Cell biology by the numbers)
# I will count opportunity costs separately, so I can consider repolimerization costs

# Opportunity costs 
sporeTranscriptOpp <- mergedExpData_time_distinct %>%
  mutate(estimation = (((abundance.filled/1e2)*1774445)/1e6)*gene_length.filled) %>% #protein abundance/1000 X 1.8 X gene length
  # the reason I multiply with 1.8 is that the protein abundance is reported as parts per million. An average size bacteria has about 3 million protein molecules
  # this was reported experimentally as average 1774445 in Bacillus (Maass et.al. 2011)
  mutate(opportunity = estimation*31) 

# Sum 
sporeTranscriptOppSum <- sporeTranscriptOpp %>%
  group_by(time) %>%
  summarise(sumOpp = sum(opportunity, na.rm =T))

# Sum stages
sporeTranscriptOppSum2 <- sporeTranscriptOpp %>%
  group_by(sourceID) %>%
  summarise(sumOpp = sum(opportunity, na.rm =T))
  
# Direct costs 
sporeTranscriptDir <- mergedExpData_time %>%
  mutate(estimation = (((abundance.filled/1e2)*1774445)/1e6)*gene_length.filled) %>% 
  mutate(direct = estimation*(10+(2*12*1)))#hours
  # average sporulation time is 8hours, median mRNA degradation rate of Bacillus is 12 per hour 
  # (DOI: 10.1007/s00438-003-0883-6)
  # = 12 re-polymerization events per hour
  # assuming nucleotides are well recycled and it only affects polymerization costs

# Sum 
sporeTranscriptDirSum <-sporeTranscriptDir %>% 
  group_by(time_h) %>%
  summarise(sumDir = sum(direct, na.rm =T))

# Stages 
sporeTranscriptDirSum2 <-sporeTranscriptDir %>% 
  group_by(sourceID) %>%
  summarise(sumDir = sum(direct, na.rm =T))

# Cumulative costs
transcriptCosts <- cbind.data.frame(sporeTranscriptDirSum$time_h, opportunity = sporeTranscriptOppSum$sumOpp, 
                                    direct = sporeTranscriptDirSum$sumDir, 
                                    total = sporeTranscriptOppSum$sumOpp + sporeTranscriptDirSum$sumDir)

transcriptSum <- colSums(transcriptCosts[,-1])

#Cost of genes 
sporeTranscriptDirDist <- sporeTranscriptDir %>%
              group_by(locus_tag) %>%
              summarise(sumDist = sum(direct, na.rm =T))

# Translation costs
#####################
# Spore Translation
#####################
##########################################################################

# Opportunity and direct costs 
sporeTranslationOppDir <- mergedExpData_time_distinct %>%
  mutate(estimation = (abundance.filled*1774445)/1e6) %>% 
  mutate(direct = estimation*aa_directSum) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum) %>% 
  mutate(total = direct + opportunity)

# Sum
sporeTranslationOppDirSum <- sporeTranslationOppDir %>%
  group_by(time_h) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T), 
            total = sum(total, na.rm = T))

translationSum <- colSums(sporeTranslationOppDirSum[,-1])

# Opportunity 
opp <- translationSum[1]+genome_opp+sporeRepSum$sumOpp+lipid_opp
# Direct
dir <- translationSum[2]+genome_dir+sporeRepSum$sumDir+lipid_dir

(opp/(opp+dir))*100 #80
(dir/(opp+dir))*100 #20

# Total costs, plots, pie, bars, model - TODO: Remove membrane
#####################
# Total, model, plot
#####################
##########################################################################

# Total costs 
cost_rep_all    <- genome_opp+genome_dir
cost_rep_part   <- sporeRepSum$sumOpp+sporeRepSum$sumDir
cost_rep_rest   <- cost_rep_all-cost_rep_part
cost_transcript <- sum(transcriptCosts$total)
cost_translation<- sum(sporeTranslationOppDirSum$total)
cost_membrane   <- lipid_opp+lipid_dir
all_pie_costs   <- cost_rep_all+cost_transcript+cost_translation+cost_membrane 
total_spore     <- cost_rep_part+cost_transcript+cost_translation+cost_membrane
total_spore_nm  <- cost_rep_part+cost_transcript+cost_translation   

# % proportions 
rep_partial     <- cost_rep_part/all_pie_costs*100 #1.7%
rep_rest        <- cost_rep_rest/all_pie_costs*100 #16.1%
transcript      <- cost_transcript/all_pie_costs*100 #%14.9 
translation     <- cost_translation/all_pie_costs*100 #51.7
membrane        <- cost_membrane/all_pie_costs*100 #15.6%
  
proportion <- c(rep_partial, rep_rest, transcript, translation, membrane)
pieCost    <- c("replication_partial", "replication_rest", "transcription", "translation", "membrane")
pieData    <- cbind.data.frame(pieCost, proportion)

ggplot(pieData, aes(x = "", y = proportion, fill = pieCost))+
  geom_bar(width = 1, stat = "identity")+ 
  coord_polar("y", start = 0)+
  mytheme+ 
  scale_fill_npg()

# spore percentage
(cost_rep_part/cost_rep_all)*100

### Total costs of sporulation ###
### Figure  ###
time        <- rep(c(1:8), times = 2)
opportunity <- transcriptCosts$opportunity + sporeTranslationOppDirSum$opportunity
direct      <- transcriptCosts$direct + sporeTranslationOppDirSum$direct
costs       <- c(opportunity, direct)
type        <- rep(c("opportunity", "direct"), each = 8) 

sporulationCosts <- cbind.data.frame(time, costs, type)
sum(sporulationCosts$costs) #1547593909

sporulationTotal <- cbind.data.frame(time = c(1:8), 
                                     costs = direct + opportunity)

# Figure 1 Model 
fit <- nls(costs ~ SSasymp(time, yf, y0, log_alpha), data = sporulationTotal)
coef(fit) 

label    <- coef(fit)[3]

tt       <- seq(1,8, by = 0.1)
pred     <- predict(fit, list(time = tt))
preddata <- cbind.data.frame(pred, tt)

RSS.p <- sum(residuals(fit)^2)
y     <- as.numeric(as.character(sporulationTotal$costs))
TSS   <- sum((log10(y) - mean(y))^2)
R2    <- 1 - (RSS.p/TSS)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Plot 
ggplot(NULL, aes(x = x, y = y))+
  #geom_hline(yintercept = 295036700., color='grey50', linetype = "dashed")+
  #geom_hline(yintercept = 118014680, color='grey50', linetype = "dashed")+
  geom_bar(data = sporulationCosts, 
           aes(x = time, y = costs, fill = type), stat = "identity", color = 'grey25')+
  ylab("Costs in units of ATP")+
  xlab("Time (hrs)")+
  geom_line(data = preddata, aes(x = tt, y = pred))+
  mytheme+
  scale_y_continuous(labels = scientific_10)+
  scale_fill_npg()+
  annotate(geom = "text", x = 8.5, y = 2e8, label = paste0("-\U03BB", "=" , round(label, 3)), hjust = "right", size = 4, fontface = 'italic')+
  annotate(geom = "text", x = 8.5, y = 1e8, label = paste0("R2", "=" , round(R2, 3)), hjust = "right", size = 4, fontface = 'italic')+
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())

# END Figure 1 

# Costs until full commitment 
line_integral <- function(x, y) {
  dx <- diff(x)
  end <- length(y)
  my <- (y[1:(end - 1)] + y[2:end]) / 2
  sum(dx *my)
} 

x <- preddata$tt[1:11]
y <- preddata$pred[1:11]
plot(x,y,"l")
line_integral(x,y)
(442154587/sum(sporulationTotal))*100 #28.6%. %only transcription and translation 

# New figure for replication 
type_rep  <- c("opportunity", "direct", "opportunity", "direct")
cost_rep  <- c(genome_opp, genome_dir, sporeRepSum$sumOpp, sporeRepSum$sumDir)
genes     <- c("whole", "whole", "partial", "partial")

repCosts <- cbind.data.frame(type_rep, cost_rep, genes)

ggplot(NULL, aes(x = x, y = y))+
  geom_bar(data = repCosts, 
           aes(x = genes, y = cost_rep, color = type_rep, fill = type_rep), stat = "identity", color = "black")+
  ylab("ATP molecules")+
  xlab("Genes")+
  mytheme+
  scale_y_continuous(labels = scientific_10)+
  scale_fill_npg()+
  theme(legend.position = "none")


# Regulons
#selected <- c("SigE", "SigF", "SigG", "SigH", "SigK")
#SigE: early mother cell-specific sporulation sigma factor
#SigF: early forespore-specific sporulation sigma factor
#SigG: late forespore-specific sporulation sigma factor
#SigH: sigma factor that controls genes of the transition phase
#SigK: late mother cell-specific sporulation sigma factor
### Regulons
                          
                          ###########
##########################GERMINATION#####################################
                          ###########
#####################
# Time course of germination 
#####################
##########################################################################
protSeqTidyAbun$gene <- protSeqTidyAbun$gene.y
germinationLong <- germination %>%
  pivot_longer(cols = c("t30", "t60", "t90", "t120", "t150"), names_to = "time", values_to = "presence") %>%
  left_join(protSeqTidyAbun, by = "gene") %>%
  filter(!presence == 0) %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

germinationLongDIST <- germination %>%
  pivot_longer(cols = c("t30", "t60", "t90", "t120", "t150"), names_to = "time", values_to = "presence") %>%
  filter(!presence == 0) %>%
  distinct(gene, .keep_all = TRUE) %>% # genes are accounted only once
  left_join(protSeqTidyAbun, by = "gene") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

#####################
# Replication costs
#####################
##########################################################################

germRep  <- germinationLongDIST %>%
  mutate(opportunity = 2*gene_length.filled*35) %>% #2 = doublestring  35 = nucleotide costs 
  mutate(direct = 2*gene_length.filled*14)

# Sum 
germRepSum <- germRep %>%
  group_by(time) %>%
  summarise(sumOpp = sum(opportunity, na.rm =T), sumDir = sum(direct, na.rm =T))
germRepTotal <- germRepSum$sumOpp+germRepSum$sumDir


# Transcription costs
######################
# Transcription costs
######################
##########################################################################

# Opportunity costs 
germTranscriptOpp   <- germinationLongDIST %>%
  mutate(estimation  = (abundance.filled/1e2)*(1774445)/1e6) %>% 
  mutate(opportunity = estimation*31) 

# Sum
germTranscriptOppSum <- germTranscriptOpp %>%
  group_by(time) %>% 
  summarise(value  = sum(opportunity, na.rm = T)) %>%
  mutate(source = rep("transcription")) %>%
  mutate(name = rep("opportunity")) 

# Direct costs 
germTranscriptDir    <- germinationLong %>%
  mutate(estimation  = (abundance.filled/1e2)*(1774445/1e6)*gene_length.filled) %>% 
  mutate(direct      = estimation*(10+(2*6))) # 6 repolimerization events every 30 min.

# Sum
germTranscriptDirSum <- germTranscriptDir %>% 
  group_by(time) %>% 
  summarise(value = sum(direct, na.rm = T)) %>%
  mutate(source = rep("transcription")) %>%
  mutate(name = rep("direct")) 
  
######################
# Translation costs
######################
##########################################################################

# Opportunity and direct costs 
germTranslationOppDir <- germinationLongDIST %>%
  mutate(estimation = ((abundance.filled)*(1774445))/1e6) %>% 
  mutate(direct = estimation*aa_directSum) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum)

# Sum
germTranslationOppDirSum <- germTranslationOppDir %>%
  group_by(time) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T)) %>%
  pivot_longer(cols = 2:3) %>%
  mutate(source = rep("translation")) %>%
  mutate(name = rep(c("direct", "opportunity"), times = 5))

sumAll_Germ <- rbind.data.frame(germTranscriptOppSum, germTranscriptDirSum, germTranslationOppDirSum)

######################
# Total costs, plots 
######################
##########################################################################

all_Germ_sum <- sumAll_Germ %>%
  group_by(time, name) %>%
  summarise(costs = sum(value)) %>%
  mutate(hours = case_when(
    time == "t30" ~ 0.5,
    time == "t60" ~ 1,
    time == "t90" ~ 1.5,
    time == "t120" ~ 2,
    time == "t150" ~ 2.5))

# Figure Model

fit2 <- nls(costs ~ SSasymp(hours, yf, y0, log_alpha), data = all_Germ_sum)
coef(fit2) 

label2    <- coef(fit2)[3]

tt2       <- seq(0.25, 2.5, by = 0.01)
pred2     <- predict(fit2, list(hours = tt2))
preddata2 <- cbind.data.frame(pred2, tt2)

RSS.p2 <- sum(residuals(fit2)^2)
y2     <- as.numeric(as.character(all_Germ_sum$costs))
TSS2   <- sum((log10(y2) - mean(y2))^2)
R22    <- 1 - (RSS.p2/TSS2)

# Figure 
ggplot(NULL, aes(x = x, y = y))+
  geom_bar(data = all_Germ_sum, 
           aes(x = hours, y = costs, fill = name), stat = "identity", color = 'grey25')+
  ylab("Costs in units of ATP")+
  xlab("Time (hrs)")+
  geom_line(data = preddata2, aes(x = tt2, y = pred2))+
  mytheme+
  scale_y_continuous(labels = scientific_10)+
  scale_fill_npg()+
  annotate(geom = "text", x = 2.5, y = 2e9, label = paste0("-\U03BB", "=" , round(label, 3)), hjust = "right", size = 4, fontface = 'italic')+
  annotate(geom = "text", x = 2.5, y = 1.5e9, label = paste0("R2", "=" , round(R2, 3)), hjust = "right", size = 4, fontface = 'italic')+
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())+
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5), limits = c(0.25,2.8))
# END Figure  

# Pie 
all_Germ_sum_pie<- sumAll_Germ %>%
  group_by(source) %>%
  summarise(costs = sum(value))

replicationGerm   = sum(germRepTotal)/sum(all_Germ_sum$costs+germRepTotal)*100 #germination genes 0.6%
transcriptionGerm = 681377848/sum(all_Germ_sum$costs+germRepTotal)*100 #8.1%
translationGerm   = 6317061229/sum(all_Germ_sum$costs+germRepTotal)*100 #74.8%

# Removing the membrane, couln't find a reference for the size of the membrane synthesized. 
# Previously I calculated the whole mebrane

proportion2 <- c(round(replicationGerm, 2), round(transcriptionGerm, 2), round(translationGerm, 2))
pieCost2    <- c("replication", "transcription", "translation")
pieData2    <- cbind.data.frame(pieCost2, proportion2)

ggplot(pieData2, aes(x = "", y = proportion2, fill = pieCost2))+
  geom_bar(width = 1, stat = "identity", color = "black")+ 
  coord_polar("y", start=0)+
  mytheme+ 
  scale_fill_manual(values = c("#4DBBD5","#3C5488","#F39B7F"))


# OTHER TRAITS
######################
# Other traits
######################
##########################################################################

# Merge trait data
mergedTraitData <- otherTraits %>%
  filter(!category == "sporulation") %>%
  filter(!category == "germination") %>%
  left_join(annotationData, by = "gene") %>%
  left_join(protSeqTidyAbun[,c("locus_tag","abundance", "aa_opportunitySum", "aa_directSum")], by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  group_by(category) %>%
  mutate(no_genes = length(category))
  
# All costs
totalCosts_traits <- mergedTraitData %>%
  mutate(translationAll = abundance.filled*(1774445/1e6)) %>%
  mutate(translationDirect = translationAll*aa_directSum) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*aa_opportunitySum) %>%
  mutate(translationTotal = translationDirect + translationOpportunity) %>%
  mutate(transcriptionAll = (abundance.filled/1e3)*(1774445/1e6)*as.numeric(gene_length)) %>%
  mutate(transcriptionDirect = transcriptionAll*(10+(2*12*1))) %>% #assuming that mRNAs transcribed at least 1 hour
  mutate(transcriptionOpportunity = transcriptionAll*31) %>%
  mutate(transcriptionTotal = transcriptionDirect + transcriptionOpportunity) %>%
  mutate(costs = transcriptionTotal + translationTotal)
  
# Cumulative costs
costs_sum_traits <- totalCosts_traits %>%
  group_by(category) %>%
  summarize(sumCosts = sum(costs, na.rm = T),
            number_genes = mean(no_genes)) %>%
  select(category, sumCosts) %>%
  add_row(category = "total_budget", sumCosts = 2.6e10) %>%
  add_row(category = "basal_metabolism", sumCosts = 3.5e8) %>%
  add_row(category = "sporulation", sumCosts = 1483262323) 

# Plot cumulative costs
ggplot(costs_sum_traits, aes(x = log10(sumCosts), y = reorder(category, sumCosts)))+
  geom_col(fill="grey", color="black")+
  xlab("Log(Total translation costs in units of ATP)")+
  ylab("Complex traits")+
  coord_cartesian(xlim = c(7.8,10.5))+
  mytheme

# Add a second axis 
costs_sum_traits_rel <- costs_sum_traits %>%
  mutate(relative = (sumCosts/2.6e10)*100)

### Figure 2 ###
ggplot(costs_sum_traits_rel, aes(x = log10(sumCosts), y = reorder(category, sumCosts)))+
  geom_col(fill="grey", color="black")+
  #geom_vline(xintercept = log10(1737135887), linetype = 'dashed')+
  
  # Custom the Y scales:
  scale_x_continuous(
    
    # Features of the first axis
    name = "Total costs in units of log(ATP)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans = ~log10((10^./2.6e10)*100), name = "% Costs relative to sporulation"))+
  
  coord_cartesian(xlim = c(7.8,10.5))+
  mytheme
### Figure 2 ###

#Basal metabolism an hour 350000000 
#350000000*8 During spore formation
(350000000*8 /2178858885)*100 


# ESSENTIAL GENES #
######################
# Essential genes
######################
##########################################################################
# Distribution of essential genes and others 
# Costs of genes genome wide 

deletionLibMerge <- deletionLib %>%
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

# Average sporulation score under 95%
mean(deletionLibMerge$rSS_average)*0.95 #0.88

deletionLibMergeRelevant <- deletionLibMerge %>%
  mutate(success = case_when(rSS_average < 0.88 & rSS_average > 0.31 ~ "auxilary",
                             rSS_average < 0.31  ~ "essential",
                             rSS_average > 0.88 ~ "nonessential"))

# All costs 
delAllCosts <- deletionLibMergeRelevant %>%
  mutate(transcript = (((abundance.filled/1e3)*1774445)/1e6)*gene_length.filled) %>%
  mutate(protein = (abundance.filled*1774445)/1e6) %>% 
  mutate(replicationOpp = gene_length.filled*35*2) %>% 
  mutate(replicationDir = gene_length.filled*14*2) %>%
  mutate(transcriptOpp = transcript*31) %>%
  mutate(transcriptDir = transcript*(10+(2*12*1))) %>%
  mutate(translationOpp = protein*aa_opportunitySum) %>% # ignoring protein degradation
  mutate(translationDir = protein*aa_directSum) %>% 
  group_by(success) %>% 
  summarise_at(vars(replicationOpp:translationDir),
               sum, na.rm = TRUE) %>% 
  pivot_longer(!success, names_to = 'source', values_to = 'costs') %>% 
  mutate(process = rep(c("replication", "transcription", "translation"), each = 2, times = 3)) %>%
  mutate(source = rep(c("opportunity", "direct"), times = 9)) %>% 
  filter(!success == "nonessential")

delAllCosts$ratio <- delAllCosts$costs*100/sum(delAllCosts$costs)
sum(delAllCosts$ratio[7:12])

### Figure 4 ###
ggplot(delAllCosts, aes(y = costs, x = source,  fill = process))+
  mytheme+
  facet_grid(~success)+
  geom_bar(stat="identity", width=1, color="white")+
  scale_fill_npg()+
  coord_cartesian(ylim = c(1e8, 1.7e9))

### Figure 4 ###

#totals 

totals <- delAllCosts %>% 
  group_by(success, source) %>% 
  summarise(sum = sum(costs))

######################
# Quality 
#####################
###########################################################################
delAllCosts2 <- deletionLibMergeRelevant %>%
  mutate(transcript = (((abundance.filled/1e2)*1774445)/1e6)*gene_length.filled) %>%
  mutate(protein = (abundance.filled*1774445)/1e6) %>% 
  mutate(replicationOpp = gene_length.filled*35*2) %>% 
  mutate(replicationDir = gene_length.filled*14*2) %>%
  mutate(transcriptOpp = transcript*31) %>%
  mutate(transcriptDir = transcript*(10+(2*12*1))) %>%
  mutate(translationOpp = protein*aa_opportunitySum) %>% # ignoring protein degradation
  mutate(translationDir = protein*aa_directSum) %>% 
  mutate(sum = translationOpp+translationDir+replicationOpp+replicationDir+transcriptOpp+transcriptDir) %>% 
  right_join(mergedExpData_time_distinct, by = "locus_tag") 
  
  ggplot(delAllCosts2, aes(x = log10(sum), y = rSS_average))+
  geom_point(alpha = .1)+
  ylab("spore quality")+
  xlab("gene costs")

# SPORULATION EFFICIENCY #
#####################
# Efficiency
#####################
##########################################################################

ggplot(efficiency, aes(efficiency))+
  geom_density()+
  xlab("Sporulation efficiency")+
  ylab("Frequency")+
  geom_vline(xintercept = 34.7, color = "#A42820", linetype = "dashed" )+
  coord_cartesian(ylim = c(0.0025, 0.0125))+
  mytheme

library(truncnorm)
fit   <- density(efficiency$efficiency, from = 0, to = 100)
N     <- 1e6
x.new <- rtruncnorm(a = 0, b = 100, N, sample(efficiency$efficiency, size = N, replace = TRUE))
plot(density(x.new, bw = fit$bw))

densityPlot <- as.data.frame(x.new)

means = densityPlot %>%
  summarise(M = median(x.new), SD = sd(x.new), N = n())

# Sample size of 50. T Distribution intervals
means$error = qt(0.975, df = 50-1)*means$SD/sqrt(means$N)
means$upper = means$M+means$error
means$lower = means$M-means$error

### Figure Efficiency ###
ggplot(densityPlot, aes(x = x.new))+
  geom_density(alpha = 0.1, size = 0.8, bw = 11.35555)+
  geom_vline(xintercept = c(means$M), linetype = 'dashed', color = "#A42820")+
  annotate("text", x = 36, y = 0.013, label= "Median = 30.6%", color = "#A42820")+
  labs(x="Sporulation efficiency", y = "Density")+
  mytheme+
  scale_x_continuous(limits = c(-30,130))
### Figure Efficiency ###


####################
# EVOLUTION
####################
########################################################################

# Merge with mutant library
# Essential genes (66 of them)
mergedMutDataTranslation <- sporeTranslationOppDir %>%
  left_join(mutantLib, by = "locus_tag") %>%
  filter(mutant == "known sporulation mutant") %>%
  select(locus_tag)

# Evolution of the genes
singleGeneTranscript  <- sporeTranscriptOpp$opportunity + sporeTranscriptDirDist$sumDist
singleGeneTranslation <- sporeTranslationOppDir$direct + sporeTranslationOppDir$opportunity
singleGeneRepOpp      <- sporeTranslationOppDir$gene_length*35*2
singleGeneRepDir      <- sporeTranslationOppDir$gene_length*13*2

scReplicationOpp <- (singleGeneRepOpp/2.6e10)*0.69
scReplicationDir <- (singleGeneRepDir/2.6e10)*0.69
scTranscriptOpp  <- (sporeTranscriptOpp$opportunity/2.6e10)*0.69
scTranscriptDir  <- (sporeTranscriptDirDist$sumDist/2.6e10)*0.69
scTranslationOpp <- (sporeTranslationOppDir$direct/2.6e10)*0.69
scTranslationDir <- (sporeTranslationOppDir$opportunity/2.6e10)*0.69

sc     <- c(scReplicationOpp, scReplicationDir, scTranscriptOpp, scTranscriptDir, 
                           scTranslationOpp, scTranslationDir)
source <- rep(c("replication", "transcription", "translation"), each = 878*2)
type   <- rep(c("opportunity", "direct"), each = 878, times = 3)

scDataAll <- cbind.data.frame(gene.lenght = sporeTranscriptOpp$gene_length.filled, 
                           protein.lenght = sporeTranscriptOpp$protein_length.filled,
                           protein.abundance = sporeTranscriptOpp$abundance.filled,
                           sc, source, type, 
                           locus_tag = sporeTranscriptOpp$locus_tag, 
                           gene = sporeTranscriptOpp$gene)

#scData <- scData[complete.cases(scData), ]

scDataMut <- scDataAll %>%
  mutate(Ne = log10(sc)*-1) %>%
  full_join(mutantLib, by = "locus_tag") 
#(Ne from Masel and Maugan 2007)

mutant <- scDataMut %>%
  filter(mutant == "known sporulation mutant")
  
  
# All sporulation genes 
ggplot(scDataMut, aes(x = log10(sc)))+
  geom_vline(xintercept = -log10(1e8), color = "#E64B35", linetype = "dashed", size=.7)+
  geom_vline(xintercept = -log10(1e5), color = "#E64B35", linetype = "dashed", size=.7)+
  #geom_histogram(aes(y=..density.., fill = source), alpha = 0.1) +
  geom_density(aes(color=source), size = 1)+
  mytheme+
  scale_color_manual(values = c( "#4DBBD5", "#00A087", "#3C5488"))+
  scale_fill_manual(values = c( "#4DBBD5", "#00A087", "#3C5488"))+
  xlab("Fractional costs (Sc)")+
  ylab("Frequency")+
  #geom_histogram(data = scDataTotal, aes(x = log10(mean), y=..density..), alpha = 0.1)+
  geom_density(data = scDataTotal, aes(x = log10(mean)), size = 1, linetype = "dashed") 


ggplot(scDataMut, aes(x = Ne, y = sc))+
  geom_vline(xintercept = 8, color = "#E64B35", linetype = "dashed", size=.7)+
  geom_vline(xintercept = 5, color = "#E64B35", linetype = "dashed", size=.7)+
  geom_point()+
  mytheme+
  scale_color_manual(values = c( "#4DBBD5", "#00A087", "#3C5488"))+
  scale_fill_manual(values = c( "#4DBBD5", "#00A087", "#3C5488"))+
  xlab("Ne")+
  ylab("Sc")+
  geom_point(data = mutant, aes(x = Ne, y = sc), color = "red")
  

##############################
# CONSERVED SPORULATION GENES
##############################
########################################################################

# Conserved sporulation genes
# After Acidamicoccus fermentans DSM 20731 non-sporulating species, add an id column

col <- rep(c('sporulating', 'non-sporulating'), times = c(21, 19))

mergedOrtData <- conservedLib %>%
  pivot_longer(cols = Alicyclobacillus_acidocaldarius.subsp._acidocaldariusDSM446:
                 Thermodesulfobium_narugenseDSM14796, names_to = "species", values_to = "presence")%>%
  mutate(col = rep(col, times = 149)) %>%
  inner_join(delAllCostsEvo, by = "locus_tag", multiple = "all") %>%
  group_by(locus_tag, process, source, presence) %>% 
  summarise_at(c("costs", "selection"), mean, na.rm = T)
  

                          ##################
##########################RETAKE GERMINATION##############################
                          ##################
####################################
# 1.Differential gene expression 
####################################
#########################################################################
germLong1 <- germination1 %>%
  pivot_longer(cols = c("Tminus30_Tminus15","Tminus15_T0","T0_T15","T15_T30","T30_T60","T60_T150","T150_T330"),
               names_to =  "time_interval", values_to = "expression") %>%
  filter(!expression == 0) %>%
  filter(!expression == -1) %>%
  mutate(time_h = case_when(
    time_interval == "Tminus30_Tminus15" ~ 0.25,
    time_interval == "Tminus15_T0" ~ 0.25,
    time_interval == "T0_T15" ~ 0.25,
    time_interval == "T15_T30" ~ 0.25,
    time_interval == "T30_T60" ~ 0.5,
    time_interval == "T60_T150" ~ 1.5,
    time_interval == "T150_T330" ~ 1)) %>%
  mutate(ordered_interval = case_when(
    time_interval == "Tminus30_Tminus15" ~ 'H00',
    time_interval == "Tminus15_T0" ~ 'H00',
    time_interval == "T0_T15" ~ 'H0.25',
    time_interval == "T15_T30" ~ 'H0.5',
    time_interval == "T30_T60" ~ 'H1',
    time_interval == "T60_T150" ~ 'H2.5',
    time_interval == "T150_T330" ~ 'H5.5'))

germLong1_merged <- germLong1 %>%
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

germLong1_merged_DIST <-  germLong1 %>%
  distinct(locus_tag, .keep_all = TRUE) %>% # genes are accounted only once
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

### Transcription costs ###
# Opportunity costs 
germination1opp   <- germLong1_merged_DIST %>%
  mutate(estimation  = (abundance.filled/1e2)*(1774445)/1e6) %>% 
  mutate(opportunity = estimation*31) %>%
  group_by(ordered_interval) %>% 
  summarise(value  = sum(opportunity, na.rm = T))%>%
  mutate(source = rep("translation"))%>%
  mutate(name = rep("direct")) 

# Direct costs 
germination1dir    <- germLong1_merged_DIST  %>%
  mutate(estimation  = (abundance.filled/1e2)*(1774445/1e6)*gene_length.filled) %>% 
  mutate(direct      = estimation*(10+(2*12*time_h)))%>%
  group_by(ordered_interval) %>% 
  summarise(value  = sum(direct, na.rm = T)) %>%
  mutate(source = rep("translation")) %>%
  mutate(name = rep("direct")) 

### Translation costs ###

# Opportunity and direct costs 
germination1translation <- germLong1_merged_DIST  %>%
  mutate(estimation = ((abundance.filled)*(1774445))/1e6) %>% 
  mutate(direct = estimation*aa_directSum) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum) %>% 
  mutate(total = direct + opportunity)%>% 
  group_by(ordered_interval) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T)) %>%
  pivot_longer(cols = 2:3) %>%
  mutate(source = rep("translation"))

germination1_all <- rbind.data.frame(germination1dir, germination1opp, germination1translation)

germination1sum <- germination1_all %>%
  group_by(ordered_interval, name) %>%
  summarise(costs = sum(value))
ggplot(data = germination1sum, aes(x = ordered_interval, y = costs, fill = name))+
  geom_bar(stat = "identity")+
  ggtitle("differential gene exp. + prot database ")

#####################################
# 2.Differential protein expression
#####################################
#########################################################################

germLong2 <- germination2 %>%
  pivot_longer(cols = c("Tminus30_T0","T0_T15","T15_T30","T30_T45","T45_T60","T60_T90","T90_T150","T150_T210", "T210_T330"),
               names_to =  "time_interval", values_to = "expression") %>%
  filter(!expression == 0) %>%
  filter(!expression == -1) %>%
  mutate(time_h = case_when(
    time_interval == "Tminus30_T0" ~ 0.5,
    time_interval == "T0_T15" ~ 0.25,
    time_interval == "T15_T30" ~ 0.25,
    time_interval == "T30_T45" ~ 0.25,
    time_interval == "T45_T60" ~ 0.25,
    time_interval == "T60_T90" ~ 0.5,
    time_interval == "T90_T150" ~ 1,
    time_interval == "T150_T210" ~ 1,
    time_interval == "T210_T330" ~ 2)) %>%
  mutate(ordered_interval = case_when(
    time_interval == "Tminus30_T0" ~ 'H00',
    time_interval == "T0_T15" ~ 'H0.25',
    time_interval == "T15_T30" ~ 'H0.5',
    time_interval == "T30_T45" ~ 'H0.75',
    time_interval == "T45_T60" ~ 'H1',
    time_interval == "T60_T90" ~ 'H1.5',
    time_interval == "T90_T150" ~ 'H2.5',
    time_interval == "T150_T210" ~ 'H3.5',
    time_interval == "T210_T330" ~ 'H5.5')) 

protSeqTidyAbun$gene <- protSeqTidyAbun$gene.y
germLong2_merged <- germLong2 %>%
  left_join(protSeqTidyAbun, by = "gene") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

germLong2_merged_DIST <-  germLong2 %>%
  distinct(gene, .keep_all = TRUE) %>% # genes are accounted only once
  left_join(protSeqTidyAbun, by = "gene") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

### Transcription costs ###
# Opportunity costs 
germination2opp   <- germLong2_merged_DIST %>%
  mutate(estimation  = (abundance.filled/1e2)*(1774445)/1e6) %>% 
  mutate(opportunity = estimation*31) %>%
  group_by(ordered_interval) %>% 
  summarise(value  = sum(opportunity, na.rm = T))%>%
  mutate(source = rep("translation"))%>%
  mutate(name = rep("direct")) 

# Direct costs 
germination2dir    <- germLong2_merged  %>%
  mutate(estimation  = (abundance.filled/1e2)*(1774445/1e6)*gene_length.filled) %>% 
  mutate(direct      = estimation*(10+(2*12*time_h)))%>%
  group_by(ordered_interval) %>% 
  summarise(value  = sum(direct, na.rm = T)) %>%
  mutate(source = rep("translation")) %>%
  mutate(name = rep("direct")) 

### Translation costs ###

# Opportunity and direct costs 
germination2translation <- germLong2_merged_DIST  %>%
  mutate(estimation = ((abundance.filled)*(1774445))/1e6) %>% 
  mutate(direct = estimation*aa_directSum) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum) %>% 
  mutate(total = direct + opportunity)%>% 
  group_by(ordered_interval) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T)) %>%
  pivot_longer(cols = 2:3) %>%
  mutate(source = rep("translation"))

germination2_all <- rbind.data.frame(germination2dir, germination2opp, germination2translation)

germination2sum <- germination2_all %>%
  group_by(ordered_interval, name) %>%
  summarise(costs = sum(value))
ggplot(data = germination2sum, aes(x = ordered_interval, y = costs, fill = name))+
  geom_bar(stat = "identity")+
  ggtitle("differential protein exp. + prot. database")


# 3. Newly synthesized proteins - Goes to manuscript
#####################################
# 3.New Proteins-Goes to the manuscript
#####################################
#########################################################################
germination6_interval <- cbind.data.frame(
  protID = germination6$ProtID, H0.25 = germination6$T15-germination6$T0,
  H0.5 = germination6$T30-germination6$T15, H1 = germination6$T60-germination6$T30, 
  H1.5 = germination6$T90-germination6$T60, H2 = germination6$T90-germination6$T60,
  H3 = germination6$T150-germination6$T90, H4 = germination6$T210-germination6$T150, 
  H5.5 = germination6$T330-germination6$T210)


germLong6_interval <- germination6_interval %>%
  pivot_longer(cols = c(H0.25:H5.5),
               names_to =  "time_interval", values_to = "score") %>%
  filter(!score <= 0)

protSeqTidyAbun$gene <- protSeqTidyAbun$gene.y

germLong6_merged_interval <- germLong6_interval %>%
  left_join(protSeqTidyAbun, by = "protID") %>%
  mutate(score.filled = replace_na(score, 18)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(time_h = case_when(
    time_interval == "H0.25" ~ 0.25,
    time_interval == "H0.5" ~ 0.25,
    time_interval == "H1" ~ 0.5,
    time_interval == "H1.5" ~ 0.5,
    time_interval == "H2" ~ 0.5,
    time_interval == "H3" ~ 1,
    time_interval == "H4" ~ 1,
    time_interval == "H5.5" ~ 1.5))

### Replication costs ###

germRep6  <- germLong6_merged_interval %>%
  mutate(opportunity = 2*gene_length.filled*35) %>% #2 = doublestring  35 = nucleotide costs 
  mutate(direct = 2*gene_length.filled*14)

# Sum 
germRepSum6 <- germRep6 %>%
  group_by(time_interval) %>%
  summarise(sumOpp = sum(opportunity, na.rm =T), sumDir = sum(direct, na.rm =T))
  germRepTotal <- germRepSum$sumOpp+germRepSum$sumDir
  
  germinationRepCost <- germRepSum6$sumOpp[1] + germRepSum6$sumDir[1]
  outgrowthRepCost   <- sum(germRepSum6$sumOpp[2:7]) + sum(germRepSum6$sumDir[2:7])  
  
### Transcription costs ###
# Opportunity costs 
germination6opp_interval   <- germLong6_merged_interval %>%
  mutate(estimation  = (score.filled/1e2)*(1774445)/1e6) %>% 
  mutate(opportunity = estimation*31) %>%
  group_by(time_interval) %>% 
  summarise(value  = sum(opportunity, na.rm = T))%>%
  mutate(source = rep("translation"))%>%
  mutate(name = rep("direct")) 

# Direct costs 
germination6dir_interval    <- germLong6_merged_interval  %>%
  mutate(estimation  = (score.filled/1e2)*(1774445/1e6)*gene_length.filled) %>% 
  mutate(direct      = estimation*(10+(2*12*time_h)))%>%
  group_by(time_interval) %>% 
  summarise(value  = sum(direct, na.rm = T)) %>%
  mutate(source = rep("translation")) %>%
  mutate(name = rep("direct")) 

### Translation costs ###

# Opportunity and direct costs 
germination6translation_interval <- germLong6_merged_interval  %>%
  mutate(estimation = ((score.filled)*(1774445))/1e6) %>% 
  mutate(direct = estimation*aa_directSum) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum) %>% 
  mutate(total = direct + opportunity)%>% 
  group_by(time_interval) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T)) %>%
  pivot_longer(cols = 2:3) %>%
  mutate(source = rep("translation"))

germination6_all_interval <- rbind.data.frame(germination6dir_interval, germination6opp_interval, germination6translation_interval)

germination6sum_interval <- germination6_all_interval %>%
  group_by(time_interval, name) %>%
  summarise(costs = sum(value)) %>%
  mutate(hour = gsub("[A-Z]", "", time_interval))

ggplot(data = germination6sum_interval, aes(x = time_interval, y = costs))+
  geom_bar(stat = "identity")+
  ggtitle("newly synthesized and counted proteins")
  #coord_cartesian(ylim = c(7.5,9.5))

germinationTT <- sum(germination6sum_interval$costs[1:2])
outgrowthTT   <- sum(germination6sum_interval$costs[2:14])

#####################################
# Model, plot, pie -Remove membrane  
#####################################
#########################################################################
germination6sum_interval$hours <- as.numeric(germination6sum_interval$hour)

germination6sum_int_sum <- germination6sum_interval %>%
  group_by(hours) %>%
  summarize(costs = sum(costs)) %>%
  filter(!hours == 0.25) %>%
  filter(!hours == 5.5)

fit3 <- nls(costs ~ SSasymp(hours, yf, y0, log_alpha), data = germination6sum_int_sum)
coef(fit3) 

label3    <- coef(fit3)[3]

tt3       <- seq(0.5, 4, by = 0.01)
pred3     <- predict(fit3, list(hours = tt3))
preddata3 <- cbind.data.frame(pred3, tt3)

RSS.p3 <- sum(residuals(fit3)^2)
y3     <- as.numeric(as.character(germination6sum_int_sum$costs))
TSS3   <- sum((log10(y3) - mean(y3))^2)
R23    <- 1 - (RSS.p3/TSS3)

# Plot 

germination6sum_interval_plotSum <- germination6sum_interval %>%
  #filter(!hours == 0.25) %>%
  #filter(!hours == 5.5) %>%
  group_by(hours) %>%
  summarize(sum = sum(costs))

ggplot(NULL, aes(x = x, y = y))+
  geom_bar(data = germination6sum_interval_plotSum, 
           aes(x = min(hours)+ seq(0.5, by = mean(diff(hours)), length = length(hours)), 
                                   y = sum), stat = "identity", color = 'grey25')+
  ylab("ATP molecules")+
  xlab("Time (h)")+
  geom_line(data = preddata3, aes(x = tt3, y = pred3))+
  mytheme+
  scale_y_continuous(labels = scientific_10)+
  scale_fill_npg()+
  #annotate(geom = "text", x = 4.1, y = 1.2e9, label = paste0("\U03BB", "=" , round(label, 3)), hjust = "right", size = 4, fontface = 'italic')+
  #annotate(geom = "text", x = 4.1, y = 1.5e9, label = paste0("r2", "=" , round(R2, 3)), hjust = "right", size = 4, fontface = 'italic')+
  #annotate(geom = "text", x = 4.1, y = 1.8e9, label = "outgrowth", hjust = "right", size = 4, fontface = 'italic')+
  #annotate(geom = "text", x = 4.1, y = 1e9, label = "germination", hjust = "right", size = 4, fontface = 'italic' )+
  theme(legend.position = "none")
  
# Pie
all_Germ_sum_pie<- sumAll_Germ %>%
  group_by(source) %>%
  summarise(costs = sum(value))

replicationGerm   = germRepTotal/sum(all_Germ_sum$costs+germRepTotal+29367246728)*100 #germination genes 0.6%
transcriptionGerm = 681377848/sum(all_Germ_sum$costs+germRepTotal+29367246728)*100 #8.1%
translationGerm   = 6317061229/sum(all_Germ_sum$costs+germRepTotal+29367246728)*100 #74.8%
membraneGerm      = 29367246728/sum(all_Germ_sum$costs+germRepTotal+29367246728)*100 #%1

proportion2 <- c(round(replicationGerm, 2), round(transcriptionGerm, 2), round(translationGerm, 2), round(membraneGerm, 2) )
pieCost2    <- c("replication", "transcription", "translation", "membrane")
pieData2    <- cbind.data.frame(pieCost2, proportion2)

ggplot(pieData2, aes(x = "", y = proportion2, fill = pieCost2))+
  geom_bar(width = 1, stat = "identity", color = "black")+ 
  coord_polar("y", start=0)+
  mytheme+ 
  scale_fill_manual(values = c("#E64B35","#4DBBD5","#3C5488","#F39B7F"))


# New take on Figure 2 

#########################
# TRAITS NEW FIGURE
#########################
#########################################################################

# Merge trait data
mergedTraitData <- otherTraits %>%
  filter(!category == "sporulation") %>%
  filter(!category == "germination") %>%
  left_join(annotationData, by = "gene") %>%
  left_join(protSeqTidyAbun[,c("locus_tag","abundance", "aa_opportunitySum", "aa_directSum")], by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  group_by(category) %>%
  mutate(no_genes = length(category))

# All costs
totalCosts_traits <- mergedTraitData %>%
  mutate(translationAll = abundance.filled*(1774445/1e6)) %>%
  mutate(translationDirect = translationAll*aa_directSum) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*aa_opportunitySum) %>%
  mutate(translationTotal = translationDirect + translationOpportunity) %>%
  mutate(transcriptionAll = (abundance.filled/1e3)*(1774445/1e6)*as.numeric(gene_length)) %>%
  mutate(transcriptionDirect = transcriptionAll*(10+(2*12*1))) %>% #assuming that mRNAs transcribed at least 1 hour
  mutate(transcriptionOpportunity = transcriptionAll*31) %>%
  mutate(transcriptionTotal = transcriptionDirect + transcriptionOpportunity) %>%
  mutate(RepOpportunity = 2*as.numeric(gene_length)*35) %>% #2 = doublestring  35 = nucleotide costs 
  mutate(RepDirect = 2*as.numeric(gene_length)*14) %>%
  mutate(RepTotal = RepOpportunity+RepDirect) %>%
  mutate(costs = transcriptionTotal + translationTotal + RepTotal)


# Cumulative costs
costs_sum_traits <- totalCosts_traits %>%
  group_by(category) %>%
  summarize(sumCosts = sum(costs, na.rm = T),
            number_genes = mean(no_genes)) %>%
  select(category, sumCosts) %>%
  add_row(category = "growth requirements", sumCosts = 2.6e10) %>%
  add_row(category = "basal metabolism", sumCosts = 3.5e8) %>%
  #add_row(category = "maintenance oligotrophic", sumCosts = 3.5e8/96) %>%
  add_row(category = "membrane", sumCosts = cost_membrane) %>%
  add_row(category = "developmental programm", sumCosts = total_spore_nm+
            outgrowthRepCost+outgrowthTT+germinationRepCost+germinationTT)%>%
  add_row(category = "genome replication", sumCosts = genome_tot) 
  
 
category = factor(c("sporulation", "germination", "outgrowth"))
sumCosts = c(total_spore_nm,germinationRepCost+germinationTT, outgrowthRepCost+outgrowthTT)
costs_sum_dev <- cbind.data.frame(category, sumCosts)

# Add a second axis 
costs_sum_traits_rel <- costs_sum_traits %>%
  mutate(relative = (sumCosts/2.6e10)*100)

### Figure 2 ###

ggplot(costs_sum_traits_rel, aes(x = log10(sumCosts), y = reorder(category, sumCosts)))+
  geom_col(fill="#F39B7F", color="black")+
  geom_vline(xintercept = log10(costs_sum_dev[1,2]))+
  geom_vline(xintercept = log10(costs_sum_dev[1,2]+costs_sum_dev[2,2]))+

  # Custom the Y scales:
  scale_x_continuous(
    
    # Features of the first axis
    name = "Costs (ATP molecules)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans = ~log10((10^./2.6e10)*100), name = "% Costs relative to sporulation"))+
  
  coord_cartesian(xlim = c(7.8,10.5))+
  mytheme+
  theme(axis.ticks.length = unit(2, "mm"))+
  theme(axis.ticks.y = element_blank())+
  theme(axis.title.y = element_blank())
  

ggplot(costs_sum_dev, aes(y = log10(sumCosts), x = reorder(category, sumCosts)))+
  geom_col(fill="grey", color="black")+
  scale_y_continuous(name = "Costs (ATP molecules)")+
  coord_cartesian(ylim = c(7.8,10.5))+
  mytheme

### Figure 2 ###

#########################
# ALTERNATIVE SPORULATION 
#########################
#########################################################################
# Alternative expression data

# Data = sporulation1
# Add new columns for differences

#for (col in 3:ncol(sporulation1)) {
#  new_col_name <- paste0("diff_", col)
#  sporulation1[[new_col_name]] <- 0  
  
#  for (i in 1:nrow(sporulation1)) {
#    sporulation1[[new_col_name]][i] <- sporulation1[[col+1]][i] - sporulation1[[col]][i]

#    }
#}

# Since the other data set is normalized with 0 

sporulation1_long <- sporulation1 %>%
  mutate(across(X0min:X480min, ~ .x - X0min)) %>%
  mutate(across(everything(), function(x){replace(x, which(x<0), 0)})) %>%
  mutate(h1 = X15min+X30min+X45min+X60min, 
         h2 = X75min+X90min+X105min+X120min,
         h3 = X135min+X150min+X120min+X180min,
         h4 = X210min+X240min,
         h5 = X270min+X300min,
         h6 = X330min+X360min,
         h7 = X390min+X420min,
         h8 = X450min+X480min) %>% 
  select(c(1, 2, starts_with("h"))) %>%
  pivot_longer(cols = starts_with("h"), names_to = 'expression', values_to = 'value') %>%
  filter(!value == 0)
  
sporulation1_distinct <-  sporulation1_long %>%
  distinct(locus_tag, .keep_all = TRUE) %>% # genes are accounted only once
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

sporulation1_time <-  sporulation1_long %>%
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

# Replication

# Opportunity costs 
sporulation1Rep  <- sporulation1_distinct%>%
  mutate(opportunity = 2*gene_length.filled*35) %>% #2 = doublestring  35 = nucleotide costs 
  mutate(direct = 2*gene_length.filled*14)

# Sum 
sporulation1RepSum <- sporulation1Rep %>%
  summarise(sumOpp = sum(opportunity, na.rm =T), sumDir = sum(direct, na.rm =T))
sporulation1ReTotal <- sporulation1RepSum$sumOpp+sporulation1RepSum$sumDir

# Transcription costs

# Opportunity costs 
sporulation1TranscriptOpp <- sporulation1_distinct %>%
  mutate(estimation = (((abundance.filled/1e2)*1774445)/1e6)*gene_length.filled) %>% #protein abundance/1000 X 1.8 X gene length
  # the reason I multiply with 1.8 is that the protein abundance is reported as parts per million. An average size bacteria has about 3 million protein molecules
  # this was reported experimentally as average 1774445 in Bacillus (Maass et.al. 2011)
  mutate(opportunity = estimation*31) 

# Sum 
sporulation1TranscriptOppSum <- sporulation1TranscriptOpp %>%
  group_by(expression) %>%
  summarise(sumOpp = sum(opportunity, na.rm =T))

# Direct costs 
sporulation1TranscriptDir <- sporulation1_time %>%
  mutate(estimation = (((abundance.filled/1e2)*1774445)/1e6)*gene_length.filled) %>% 
  mutate(direct = estimation*(10+(2*12*1))) #hours

# Sum 
sporulation1TranscriptDirSum <-sporulation1TranscriptDir %>% 
  group_by(expression) %>%
  summarise(sumDir = sum(direct, na.rm =T))

#cbind.fill<-function(...){
#  nm <- list(...) 
#  nm<-lapply(nm, as.matrix)
#  n <- max(sapply(nm, nrow)) 
#  do.call(cbind, lapply(nm, function (x) 
#    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
#}

# Cumulative costs
transcriptCostsspor1 <- as.data.frame(cbind.fill(sporulation1TranscriptDirSum$expression, opportunity = sporulation1TranscriptOppSum$sumOpp, 
                                    direct = sporulation1TranscriptDirSum$sumDir))
transcriptCostsspor1[is.na(transcriptCostsspor1)] <- 0
names(transcriptCostsspor1) <- c("hours", "opportunity", "direct")

# Translation costs

# Opportunity and direct costs 
sporeTranslationOppDirspor1 <- sporulation1_distinct %>%
  mutate(estimation = (abundance.filled*1774445)/1e6) %>% 
  mutate(direct = estimation*aa_directSum) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum) %>% 
  mutate(total = direct + opportunity)

# Sum
sporeTranslationOppDirspor1Sum <- sporeTranslationOppDirspor1 %>%
  group_by(expression) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T), 
            total = sum(total, na.rm = T))


# Total costs, plots, pie, bars, model - TODO: Remove membrane

# Total costs 
cost_rep_all    <- genome_opp+genome_dir
cost_rep_part2   <- sporulation1RepSum$sumOpp+sporulation1RepSum$sumDir
cost_rep_rest2   <- cost_rep_all-cost_rep_part
cost_transcript2  <- sum(sporulation1TranscriptOppSum$sumOpp, sporulation1TranscriptDirSum$sumDir)
cost_translation2 <- sum(sporeTranslationOppDirspor1Sum$total)
cost_membrane    <- lipid_opp+lipid_dir
all_pie_costs2    <- cost_rep_all2+cost_transcript2+cost_translation2+cost_membrane 
total_spore2      <- cost_rep_part2+cost_transcript2+cost_translation2+cost_membrane
total_spore_nm2   <- cost_rep_part2+cost_transcript2+cost_translation2   

# % proportions 
rep_partial2     <- cost_rep_part2/all_pie_costs2*100 
rep_rest2        <- cost_rep_rest2/all_pie_costs2*100 
transcript2      <- cost_transcript2/all_pie_costs2*100 
translation2     <- cost_translation2/all_pie_costs2*100 
membrane2        <- cost_membrane2/all_pie_costs2*100 

proportion2 <- c(rep_partial2, rep_rest2, transcript2, translation2, membrane)
pieCost2    <- c("replication_partial", "replication_rest", "transcription", "translation", "membrane")
pieData2    <- cbind.data.frame(pieCost2, proportion2)

ggplot(pieData2, aes(x = "", y = proportion, fill = pieCost))+
  geom_bar(width = 1, stat = "identity")+ 
  coord_polar("y", start = 0)+
  mytheme+ 
  scale_fill_npg()

### Total costs of sporulation ###
### Figure  ###
time2        <- rep(1:8, times = 2)
copp2        <- as.data.frame(cbind.fill(transcriptCostsspor1$opportunity, sporeTranslationOppDirspor1Sum$opportunity))
copp2[is.na(copp2)] <- 0
opportunity2 <- as.numeric(copp2$V1)+as.numeric(copp2$V2)
cdir2        <- as.data.frame(cbind.fill(transcriptCostsspor1$direct, sporeTranslationOppDirspor1Sum$direct))
cdir2[is.na(cdir2)] <- 0
direct2      <- as.numeric(cdir2$V1)+ as.numeric(cdir2$V2)
costs2       <- c(opportunity2, direct2)
type2        <- rep(c("opportunity", "direct"), each = 8, times=2) 

sporulationCosts2 <- cbind.data.frame(time2, costs2, type2)
sum(sporulationCosts2$costs) #8572993254

ggplot(sporulationCosts2, aes(x = time2, y = costs2))+
  geom_bar(stat = "identity")
