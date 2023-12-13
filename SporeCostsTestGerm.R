# Spore Costs v2
# 26 October 2023 - last update
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
# COGs: https://doi.org/10.1128/jb.00079

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
library(scales)
# for COGs
library(vegan)
library(ggfortify)
library(ggrepel)

mytheme <- theme_bw()+
  theme(axis.ticks.length = unit(.25, "cm"))+
  theme(legend.text = element_text(size=16))+
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 18))+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.border = element_rect(fill=NA, colour = "black", 
size=1))+
theme(strip.text.x = element_text(size = 18))+
  theme(legend.title=element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", 
                                    linewidth=1)) +
theme(axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
        axis.text.y.right = element_blank(), axis.title.y.right = element_blank())+
   theme(axis.title.x = element_text(margin=margin(10,0,0)),
   axis.title.y = element_text(margin=margin(0,10,0,0)),
   axis.text.x = element_text(margin=margin(10,0,0,0)),
   axis.text.y = element_text(margin=margin(0,10,0,0)))

mytheme2 <- theme_bw()+
  theme(axis.ticks.length = unit(.25, "cm"))+
  theme(legend.text = element_text(size=14))+
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", 
                                    size=1))+
  theme(strip.text.x = element_text(size = 18))+
  theme(legend.title=element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", 
                                    linewidth=1)) +
  #theme(axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
  #     axis.text.y.right = element_blank(), axis.title.y.right = element_blank())+
  theme(axis.title.x = element_text(margin=margin(10,0,0)),
        axis.title.y = element_text(margin=margin(0,10,0,0)),
        axis.text.x = element_text(margin=margin(10,0,0,0)),
        axis.text.y = element_text(margin=margin(0,10,0,0)))+
  theme(axis.title.x.top = element_text(margin=margin(10,10,0)),
        axis.title.y.right = element_text(margin=margin(0,10,0,0)),
        axis.text.x.top = element_text(margin=margin(10,0,10,0)),
        axis.text.y.right = element_text(margin=margin(0,10,0,0)))

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

# Regulator list - sigma factors - annotations
sigma          <- read.table("./otherData/regulations.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))
functions      <- read.table("./otherData/subtiwiki.gene.export.2023-11-17.csv", sep = ",", header = T, fill = TRUE,  stringsAsFactors = F, na.strings=c(" ","NA"))
descriptions   <- read.table("./otherData/SWxWWxRS_sporulation_genes.csv", sep = ",", header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

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
efficiency     <- read.table("./otherData/efficiency.csv", sep = ",", dec = "," , header = T, stringsAsFactors = F)

# COGs - Galperin
cogs_dat_gene  <- read.table("./COGs_Galperin/geneMatrix.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)
cogs_dat_sp    <- read.table("./COGs_Galperin/speciesMatrix.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)
#replaced empty cells with 0


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

mergedExpData_time_distinct <- expressionLong %>%
  distinct(locus_tag, .keep_all = TRUE) %>% # genes are accounted only once
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))

mergedExpData_time <- expressionLong %>%
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

# DNA unwinding 1ATP per base pair #https://doi.org/10.1016/j.cell.2006.10.049
unwind <- 4214810

# Opportunity 
# 2 * Genome size * (34+1)
genome_opp <- 2 * 4214810 * 35 #295036700
# Direct 
# 2 * Genome size * (11+2)
genome_dir <-( 2 * 4214810 * 14 ) + unwind #122229490
# Total
genome_tot <- genome_opp + genome_dir #417266190

#####################
# Membrane costs
#####################
##########################################################################

# Number of lipid molecules = Cellular membrane areas/head-group areas of membrane lipid molecules
# Head group area is a1 = 0.65 nm2 (Nagle and Tristram-Nagle 2000; Petrache et al. 2000; Kucerka et al. 2011).

# Thickness of the bilayer (h):
# The thickness of a bilayer is approximately twice the radius of the head-group area, which 0.5 nm in all cases,
# plus the total length of the internal hydrophobic tail domains (Lewis and Engelman 1983; Mitra et al. 2004), 
# generally 3.0 nm, so total is 4nm

# Bacillus average length (a) and width (b) Barak et al. 2018 
a1 <- 0.65
h  <- 4
a  <- 2.5*1000 #convert to nm
b  <- 1*1000 #also septum
# height/width=2.5
c  <- 0.4 

# Outer #4πa*b
outer <- 4*pi*a*b 
moleculesOut <- outer/a1 
# Inner membrane 4π(a − h)(b - h)
inner <- 4*pi*((a-h)*(b-h))
moleculesInn <- inner/a1 
# 50% discount, because of protein
totalMol <- (moleculesOut+moleculesInn)/2 #total lipid molecules

# Cost of lipid head & tail 
# Opportunity = 212, Direct = 18
# costs
membraneOpp <- totalMol*212 
membraneDir <- totalMol*18 
membraneTot <- membraneOpp+membraneDir

# Septum should be 1µm 1000 nm (as the width of the cell) 
septumOut <- ((4*b)/a1)/2
septumInn <- ((4*(b-h))/a1)/2 
septumTot <- septumOut + septumInn

# costs
septumOpp     <- septumTot*212 
septumDir     <- septumTot*18 

# Germination assuming that they recycle membrane of the endospore
# Whole membrane - (endospore sphere + septum)
# Septum stretches 
# Endospore size is 1/6 of the total cell
membraneGerm  <- membraneTot/6

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
sporeRepTotal/413051380  #17%

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
# Fill missing protein sequence cost estimations
median(mergedExpData_time_distinct$aa_opportunitySum, na.rm = T) #5723
median(mergedExpData_time_distinct$aa_directSum, na.rm = T) #1351

# Opportunity and direct costs 
sporeTranslationOppDir <- mergedExpData_time_distinct %>%
  mutate(estimation = (abundance.filled*1774445)/1e6) %>% 
  mutate(aa_opportunitySum.filled = replace_na(aa_opportunitySum, 5723)) %>%
  mutate(aa_directSum.filled = replace_na(aa_directSum, 1351)) %>%
  mutate(direct = estimation*aa_directSum.filled) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum.filled) %>% 
  mutate(total = direct + opportunity)

# Sum
sporeTranslationOppDirSum <- sporeTranslationOppDir %>%
  group_by(time_h) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T), 
            total = sum(total, na.rm = T))

translationSum <- colSums(sporeTranslationOppDirSum[,-1])

# Opportunity 
opp <- translationSum[1]+genome_opp+sporeRepSum$sumOpp+septumOpp
# Direct
dir <- translationSum[2]+genome_dir+sporeRepSum$sumDir+septumDir

(opp/(opp+dir))*100 #78
(dir/(opp+dir))*100 #22

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
cost_membrane   <- septumOpp+septumDir
all_pie_costs   <- cost_rep_all+cost_transcript+cost_translation+cost_membrane

# % proportions 
rep_partial     <- cost_rep_part/all_pie_costs*100 
rep_rest        <- cost_rep_rest/all_pie_costs*100 
((cost_rep_part+cost_rep_rest)/all_pie_costs)*100
transcript      <- cost_transcript/all_pie_costs*100
translation     <- cost_translation/all_pie_costs*100
membrane        <- cost_membrane/all_pie_costs*100
  
proportion <- c(rep_partial, rep_rest, transcript, translation, membrane)
pieCost    <- c("replication_partial", "replication_rest", "transcription", "translation", "membrane")
pieData    <- cbind.data.frame(pieCost, proportion)

pieData$labels <- paste(pieData$pieCost, round(pieData$proportion, 1), "%")

pieSpore <- ggplot(pieData, aes(x = "", y = proportion, fill = pieCost))+
  geom_bar(width = 1, stat = "identity")+ 
  coord_polar("y", start = 0)+
  mytheme+ 
  scale_fill_npg()+
  geom_text_repel(aes(label = labels), size = 4.5, show.legend = FALSE)+
  theme_void()+
  theme(legend.position = "none")

ggsave("~/GitHub/sporeCostsVer2/figures/pieSpore.pdf", pieSpore, height = 5, width = 6)

# spore percentage
(cost_rep_part/cost_rep_all)*100 #17%

### Total costs of sporulation ###
### Figure  ###
time        <- rep(c(1:8), times = 2)
opportunity <- transcriptCosts$opportunity + sporeTranslationOppDirSum$opportunity
direct      <- transcriptCosts$direct + sporeTranslationOppDirSum$direct
costs       <- c(opportunity, direct)
type        <- rep(c("opportunity", "direct"), each = 8) 

sporulationCosts <- cbind.data.frame(time, costs, type)
sum(sporulationCosts$costs) #1554146274 (transcription and translation)

cost_rep_part+cost_rep_rest+sum(sporulationCosts$costs) #1971412464

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

#scientific_10 <- function(x) {
#  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
#}

sporulationCosts_lay1 <- sporulationCosts %>%
  group_by(time) %>%
  summarize(sum = sum(costs))

sporulationCosts$type <- factor(sporulationCosts$type, c("opportunity","direct"))
# Plot 
# Note that some adjustments are done in Adobe Illustrator

my_lab <- c(expression(P['D']),
            expression(P['O']), 
            expression(P['T']))


f1 <- ggplot(NULL, aes(x = x, y = y))+
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_bar(data = sporulationCosts_lay1, 
           aes(x = time, y = sum), stat = "identity", color = "grey90", fill = "grey75", alpha = 0.5)+
  geom_bar(data = sporulationCosts, 
           aes(x = time, fill = type, y = costs), stat = "identity", color = "grey25", position = position_dodge(width=1))+
  ylab("ATP molecules")+
  xlab("Time (h)")+
  geom_line(data = preddata, aes(x = tt, y = pred), linewidth = 1)+
  mytheme+
  scale_y_continuous(breaks = c(2*10^8, 4*10^8, 6*10^8, 8*10^8, 10*10^8), 
                     labels = c(2,4,6,8,10), sec.axis=dup_axis())+
  scale_x_continuous(sec.axis=dup_axis())+
  annotate("text",x=-0.7,y=1.2*10^9,label=paste("(x10^8)"), parse =T, size = 18/.pt)+
  coord_cartesian(xlim = c(0.5, 8.5), clip="off")+
  
  annotate(geom = "text", x = 8, y = 2e8, label = paste("-\U03BB==", round(label, 3)), hjust = "right", size = 6, fontface = 'italic', parse = T)+
  annotate(geom = "text", x = 8, y = 1.2e8, label = paste("R^2==", round(R2, 3)), hjust = "right", size = 6, fontface = 'italic', 
           parse=TRUE)+
  theme(legend.position = c(0.32, 0.83), legend.title = element_blank())+
  scale_fill_npg(labels=c(my_lab[1], 
                               my_lab[2],
                               my_lab[3]))

ggsave("~/GitHub/sporeCostsVer2/figures/sporeCostsTime.pdf", f1, height = 5, width = 6)
# END Figure 1 

# Costs until full commitment 
line_integral <- function(x, y) {
  dx <- diff(x)
  end <- length(y)
  my <- (y[1:(end - 1)] + y[2:end]) / 2
  sum(dx *my)
} 

x <- preddata$tt[1:11] #2 hours
y <- preddata$pred[1:11]
plot(x,y,"l")
commitment <- line_integral(x,y)
(commitment/sum(sporulationTotal))*100 #34.9%. %only transcription and translation 

# New figure for replication 
type_rep  <- c("opportunity", "direct", "opportunity", "direct")
cost_rep  <- c(genome_opp, genome_dir, sporeRepSum$sumOpp, sporeRepSum$sumDir)
genes     <- c("whole", "whole", "partial", "partial")

repCosts <- cbind.data.frame(type_rep, cost_rep, genes)

# Regulons
#selected <- c("SigE", "SigF", "SigG", "SigH", "SigK")
#SigE: early mother cell-specific sporulation sigma factor
#SigF: early forespore-specific sporulation sigma factor
#SigG: late forespore-specific sporulation sigma factor
#SigH: sigma factor that controls genes of the transition phase
#SigK: late mother cell-specific sporulation sigma factor
### Regulons
                          
                ############################
################GERMINATION - UPDATED BELOW##################################
                ############################
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
# Fill empty cost estimations
median(germinationLongDIST$aa_opportunitySum, na.rm = T) #6337
median(germinationLongDIST$aa_directSum, na.rm = T) #1527

germTranslationOppDir <- germinationLongDIST %>%
  mutate(estimation = ((abundance.filled)*(1774445))/1e6) %>% 
  mutate(aa_opportunitySum.filled = replace_na(aa_opportunitySum, 6337)) %>%
  mutate(aa_directSum.filled = replace_na(aa_directSum, 1527)) %>%
  mutate(direct = estimation*aa_directSum.filled) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum.filled) 

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

# Figure 1 - germination 
spore <- ggplot(NULL, aes(x = x, y = y))+
  geom_bar(data = all_Germ_sum, 
           aes(x = hours, y = costs, fill = name), stat = "identity", color = 'grey25')+
  ylab("ATP molecules")+
  xlab("Time (h)")+
  geom_line(data = preddata2, aes(x = tt2, y = pred2))+
  mytheme+
  scale_y_continuous(labels = scientific_10)+
  scale_fill_npg()+
  annotate(geom = "text", x = 2.5, y = 2e9, label = paste0("-\U03BB", "=" , round(label, 3)), hjust = "right", size = 4, fontface = 'italic')+
  annotate(geom = "text", x = 2.5, y = 1.5e9, label = paste0("R2", "=" , round(R2, 3)), hjust = "right", size = 4, fontface = 'italic')+
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())+
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5), limits = c(0.25,2.8))
# END Figure 1 - germination  

# Pie 
all_Germ_sum_pie<- sumAll_Germ %>%
  group_by(source) %>%
  summarise(costs = sum(value))

replicationGerm   = sum(germRepTotal)/sum(all_Germ_sum$costs+germRepTotal)*100 #germination genes 0.6%
transcriptionGerm = sum(germTranscriptDirSum$value)/sum(all_Germ_sum$costs+germRepTotal)*100 #9.5%
translationGerm   = sum(germTranslationOppDirSum$value)/sum(all_Germ_sum$costs+germRepTotal)*100 #74.8%

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
# Other traits- UPDATED BELOW
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
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) %>%
  mutate(translationDirect = translationAll*aa_directSum) %>% #ignoring protein degradation
  mutate(translationOpportunity = translationAll*aa_opportunitySum.filled) %>%
  mutate(translationTotal = translationDirect + aa_opportunitySum.filled) %>%
  mutate(transcriptionAll = (abundance.filled/1e2)*(1774445/1e6)*as.numeric(gene_length)) %>%
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
    name = "ATP molecules",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans = ~log10((10^./2.6e10)*100), name = "Costs relative to the total budget (%)"))+
  
  coord_cartesian(xlim = c(7.8,10.5))+
  mytheme
### Figure 2 ###

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

# Average sporulation score under 50% - %90
mean(deletionLibMerge$rSS_average)*0.5 
mean(deletionLibMerge$rSS_average)*0.9

deletionLibMergeRelevant <- deletionLibMerge %>%
  mutate(success = case_when(rSS_average < 0.46 & rSS_average > 0.31 ~ "50% failure",
                             rSS_average < 0.83 & rSS_average > 0.46 ~ "10% failure",
                             rSS_average < 0.31 ~ "essential",
                             rSS_average > 0.83 ~ "90% success"))

# All costs 
delAllCosts <- deletionLibMergeRelevant %>%
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) %>%
  mutate(transcript = (((abundance.filled/1e2)*1774445)/1e6)*gene_length.filled) %>%
  mutate(protein = (abundance.filled*1774445)/1e6) %>% 
  mutate(replicationOpp = gene_length.filled*35*2) %>% 
  mutate(replicationDir = gene_length.filled*14*2) %>%
  mutate(transcriptOpp = transcript*31) %>%
  mutate(transcriptDir = transcript*(10+(2*12*1))) %>%
  mutate(translationOpp = protein*aa_opportunitySum.filled) %>% # ignoring protein degradation
  mutate(translationDir = protein*aa_directSum.filled) %>% 
  group_by(success) %>% 
  summarise_at(vars(replicationOpp:translationDir),
               sum, na.rm = TRUE) %>% 
  pivot_longer(!success, names_to = 'source', values_to = 'costs') %>% 
  mutate(process = rep(c("replication", "transcription", "translation"), each = 2, times = 4)) %>%
  mutate(source = rep(c("opportunity", "direct"), times = 12)) %>% 
  filter(!success == "90% success")

delAllCosts$ratio <- delAllCosts$costs*100/sum(delAllCosts$costs)
sum(delAllCosts$ratio[7:12])

#totals 
totals <- delAllCosts %>% 
  group_by(success, source) %>% 
  summarise(sum = sum(costs))

### Figure S ###
my_lab <- c(expression(P['D']),
            expression(P['O']), 
            expression(P['T']))

essential <- ggplot(totals, aes(y = sum, x = success,  fill = source))+
  mytheme+
  #facet_wrap(~source, scales = "free")+
  geom_bar(stat="identity", width = 1, color="black")+
  scale_fill_npg(labels = my_lab)+
  theme(legend.position = c(0.85, 0.85))+
  theme(axis.title.x = element_blank())+
  ylab("ATP molecules")+
  scale_y_continuous(breaks = c(0, 0.5*10^9, 1*10^9, 1.5*10^9), 
                     labels = c(0,0.5,1,1.5), sec.axis=dup_axis())
  
### Figure S ###
ggsave("~/GitHub/sporeCostsVer2/figures/essentialGenes.pdf", essential, height = 4, width = 5)

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
  mutate(aa_opportunitySum.filled = replace_na(median(aa_opportunitySum, na.rm = T))) %>%
  mutate(aa_directSum.filled = replace_na(median(aa_directSum, na.rm = T))) %>%
  mutate(translationOpp = protein*aa_opportunitySum.filled) %>% 
  mutate(translationDir = protein*aa_directSum.filled) %>% 
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
efficiency$efficiency <- as.numeric(efficiency$efficiency)
  
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
efficiencyPlot <- ggplot(densityPlot, aes(x = x.new))+
  geom_density(alpha = 0.1, linewidth = 0.8, bw = 11.35555)+
  geom_vline(xintercept = c(means$M), linetype = 'dashed', color = "#A42820")+
  annotate(geom = "label", x = 36, y = 0.013, color = "#A42820", fill = "white", label = "Median = 30.6%")+
  labs(x="Sporulation efficiency", y = "Density")+
  mytheme+
  scale_x_continuous(limits = c(-30,130), sec.axis = dup_axis())+
  scale_y_continuous(sec.axis = dup_axis())
### Figure Efficiency ###

ggsave("~/GitHub/sporeCostsVer2/figures/efficiency.pdf", efficiencyPlot, height = 5, width = 6)
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
  ylab("Frequency")

ggplot(scDataMut, aes(x = Ne, y = sc))+
  #geom_vline(xintercept = 8, color = "#E64B35", linetype = "dashed", size=.7)+
  #geom_vline(xintercept = 5, color = "#E64B35", linetype = "dashed", size=.7)+
  geom_point()+
  mytheme+
  scale_color_manual(values = c( "#4DBBD5", "#00A087", "#3C5488"))+
  scale_fill_manual(values = c( "#4DBBD5", "#00A087", "#3C5488"))+
  xlab("Ne")+
  ylab("Sc")+
  geom_point(data = mutant, aes(x = Ne, y = sc), color = "red")
  

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
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) %>%
  mutate(direct = estimation*aa_directSum.filled) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum.filled) %>% 
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
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) %>%
  mutate(direct = estimation*aa_directSum.filled) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum.filled) %>% 
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
# 3.New Proteins - Goes to the manuscript
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
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) %>%
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
  mutate(direct = estimation*aa_directSum.filled) %>% # 
  mutate(opportunity = estimation*aa_opportunitySum.filled) %>% 
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

#Add germination to the first data point 
germination6sum_interval[3,3] <- germination6sum_interval[3,3] + germination6sum_interval[1,3]
germination6sum_interval[4,3] <- germination6sum_interval[4,3] + germination6sum_interval[2,3]

germinationTT <- sum(germination6sum_interval$costs[1:2])
outgrowthTT   <- sum(germination6sum_interval$costs[2:14])
outTT <- (outgrowthTT+membraneGerm)-germinationTT
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

tt3       <- seq(0.5, 4.9, by = 0.01)
pred3     <- predict(fit3, list(hours = tt3))
preddata3 <- cbind.data.frame(pred3, tt3)

RSS.p3 <- sum(residuals(fit3)^2)
y3     <- as.numeric(as.character(germination6sum_int_sum$costs))
TSS3   <- sum((log10(y3) - mean(y3))^2)
R23    <- 1 - (RSS.p3/TSS3)

# Plot 

germination6sum_interval_plotSum <- germination6sum_interval %>%
  filter(!hours == 0.25) %>%
  filter(!hours == 5.5) %>%
  group_by(hours) %>%
  summarize(sum = sum(costs))

germination6sum_interval$name <- factor(germination6sum_interval$name, c("opportunity", "direct"))
  
germination6sum_interval_plotSum_n <- germination6sum_interval %>%
  filter(!hours == 0.25) %>%
  filter(!hours == 5.5)

germination6sum_interval_plotSum_out <- germination6sum_interval %>%
  filter(!hours == 0.25) %>%
  group_by(hours) %>%
  summarize(sum = sum(costs))

germination6sum_interval_plotSum_out_all <- germination6sum_interval %>%
  filter(!hours == 0.25) %>%
  group_by(hours, name) %>%
  summarize(sum = sum(costs))

min15_sep <- germination6sum_interval %>%
  filter(hours == 0.25)

# Figure 1 - germination - new 
germs <- ggplot(NULL, aes(x = x, y = y))+
  geom_vline(xintercept = 4.92, linetype = "dashed")+
  geom_bar(data = germination6sum_interval_plotSum_out, 
           aes(x = seq(0.5, by = mean(diff(hours)), length = length(hours)), 
                                   y = sum), stat = "identity", color = "grey90", fill = "grey75", alpha = .5)+
  
  annotate("rect", xmin = 0.11, xmax = 0.85, ymin = 0, ymax = 84970721+323997094,
           alpha = .5,fill = "grey25")+
  geom_segment(aes(x = 0.11, xend = 0.85, y = 84970721+323997094, yend = 84970721+323997094), color = "grey25")+
  geom_bar(data = germination6sum_interval_plotSum_out_all, 
           aes(x = seq(0.5, by = mean(diff(hours)), length = length(hours)), fill = name, 
               y = sum), position = position_dodge(0.9), stat = "identity", color = "grey25")+

  ylab("ATP molecules")+
  xlab("Time (h)")+
  geom_line(data = preddata3, aes(x = tt3, y = pred3), linewidth = 1)+
  annotate(geom = "text", x = 4.8, y = 5.1e8, label = paste("-\U03BB==", round(label3, 3)), hjust = "right", size = 5, fontface = 'italic', parse = T)+
  annotate(geom = "text", x = 4.8, y = 3.8e8, label = paste("R^2==", round(R23, 3)), hjust = "right", size = 5, fontface = 'italic', 
           parse=TRUE)+
  mytheme+
  scale_y_continuous(breaks = c(4*10^8, 8*10^8, 12*10^8, 16*10^8, 20*10^8), 
                     labels = c(4,8,12,16,20), sec.axis=dup_axis())+
  scale_x_continuous(sec.axis=dup_axis())+
  annotate("text",x=-0.55,y=2.3*10^9,label=paste("(x10^8)"), parse =T, size = 16/.pt)+
  coord_cartesian(xlim = c(0.25, 6), clip="off")+
  scale_fill_npg(labels=c(my_lab[1], 
                          my_lab[2],
                          my_lab[3]))+
  annotate("rect", xmin = 0.33, xmax = 0.67 , ymin = 0, ymax = 323997094,
           alpha = .6,fill = "grey25")+
  geom_segment(aes(x = 0.33, xend = 0.67, y = 323997094, yend = 323997094), color = "grey25")+
  annotate("rect", xmin = 0.71, xmax = 1.05, ymin = 0, ymax = 84970721,
                                              alpha = .6,fill = "grey25")+
  geom_segment(aes(x = 0.71, xend = 1.05, y = 84970721, yend = 84970721), color = "grey25")+
  theme(legend.position = c(0.22, 0.88), legend.title = element_blank())

ggsave("~/GitHub/sporeCostsVer2/figures/germCostsTime.pdf", germs, height = 5, width = 6)
# End

# Pie
all_Germ_sum_pie<- sumAll_Germ %>%
  group_by(source) %>%
  summarise(costs = sum(value))

all_Germ_cost <- (sum(all_Germ_sum$costs))+ (sum(germRepTotal))+membraneGerm

replicationGerm   = (sum(germRepTotal))/all_Germ_cost*100 #germination genes 0.6%
transcriptionGerm = all_Germ_sum_pie$costs[1]/all_Germ_cost*100 #7.8%
translationGerm   = all_Germ_sum_pie$costs[2]/all_Germ_cost*100 #70.9%
membrGerm         = membraneGerm/all_Germ_cost*100 #%20.7

proportion2 <- c(round(replicationGerm, 2), round(transcriptionGerm, 2), round(translationGerm, 2), round(membrGerm, 2))
pieCost2    <- c("replication", "transcription", "translation", "membrane")
pieData2    <- cbind.data.frame(pieCost2, proportion2)

pieData2$labels <- paste(pieData2$pieCost, round(pieData2$proportion, 1), "%")

pieGerm <- ggplot(pieData2, aes(x = "", y = proportion2, fill = pieCost2))+
  geom_bar(width = 1, stat = "identity")+ 
  coord_polar("y", start=0)+
  mytheme+ 
  scale_fill_manual(values = c("#E64B35","#4DBBD5","#3C5488","#F39B7F"))+
geom_text_repel(aes(label = labels), size = 4.5, show.legend = FALSE)+
  theme_void()+
  theme(legend.position = "none")

ggsave("~/GitHub/sporeCostsVer2/figures/pieGerm.pdf", pieGerm, height = 5, width = 5)

germinationTT <- sum(germination6sum_interval$costs[1:2])
outgrowthTT   <- sum(germination6sum_interval$costs[3:14])
outTT <- (outgrowthTT+membraneGerm)-germinationTT

# New take on Figure 2 

#########################
# TRAITS UPDATED FIGURE
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
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) %>%
  mutate(translationDirect = translationAll*aa_directSum.filled) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*aa_opportunitySum.filled) %>%
  mutate(translationTotal = translationDirect + translationOpportunity) %>%
  mutate(transcriptionAll = (abundance.filled/1e2)*(1774445/1e6)*as.numeric(gene_length)) %>%
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
  add_row(category = "membrane", sumCosts = membraneOpp+membraneDir) %>%
  add_row(category = "developmental programm", sumCosts = all_pie_costs+
            germinationTT+outTT)%>%
  add_row(category = "genome replication", sumCosts = genome_tot) 
  
 
category = factor(c("sporulation", "germination", "outgrowth"))
sumCosts = c(all_pie_costs,germinationTT,outTT)
costs_sum_dev <- cbind.data.frame(category, sumCosts)

# Add a second axis 
costs_sum_traits_rel <- costs_sum_traits %>%
  mutate(relative = (sumCosts/2.6e10)*100)

### Figure 2 ###
sumAxis = all_pie_costs+outTT+germinationTT

all_pie_costs/sumAxis*100
outgrowthTT/sumAxis*100
germinationTT/sumAxis*100

log10(sumAxis)
log10(outgrowthTT+germinationTT)

costs_sum_traits_rel$category2 <- c("biofilm structure", "chemotaxis", "competence", 
                                    "essential genes", "flagella", "heat-shock proteins", 
                                    "homeostasis", "swarming", "total cell budget", "basal metabolism h-1", 
                                    "membrane lipid synthesis", "developmental programm", "genome replication")
  
f2 <- ggplot(costs_sum_traits_rel, aes(x = log10(sumCosts), y = reorder(category2, sumCosts)))+
  geom_col(fill="grey75", color="grey25")+
  geom_vline(xintercept = log10(outgrowthTT))+
  geom_vline(xintercept = log10(all_pie_costs))+
  geom_vline(xintercept = log10(germinationTT))+
  mytheme2+
  # Custom the Y scales:
  scale_x_continuous(
                  
    # Features of the first axis
    name = "ATP molecules",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans = ~log10((10^./2.6e10)*100), name = "Costs relative to the cell budget (%) "))+
  
  coord_cartesian(xlim = c(7.8,10.5))+
  
  theme(axis.ticks.y = element_blank())+
  theme(axis.title.y = element_blank())
ggsave("~/GitHub/sporeCostsVer2/figures/otherTraits.pdf", f2, height = 5, width = 6)
### Figure 2 ###  

# germination, sporulation, outgrowth 
ggplot(costs_sum_dev, aes(y = log10(sumCosts), x = reorder(category, sumCosts)))+
  geom_col(fill="grey", color="black")+
  scale_y_continuous(name = "Costs (ATP molecules)")+
  coord_cartesian(ylim = c(7.8,10.5))+
  mytheme

#########################
# ALTERNATIVE SPORULATION 
#########################
#########################################################################
# Alternative expression data for validation

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
  mutate(protein_length.filled = replace_na(protein_length, 254)) %>%
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) 

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

cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# Cumulative costs
transcriptCostsspor1 <- as.data.frame(cbind.fill(sporulation1TranscriptDirSum$expression, opportunity = sporulation1TranscriptOppSum$sumOpp, 
                                    direct = sporulation1TranscriptDirSum$sumDir))
transcriptCostsspor1[is.na(transcriptCostsspor1)] <- 0
names(transcriptCostsspor1) <- c("hours", "opportunity", "direct")

# Translation costs

# Opportunity and direct costs 
sporeTranslationOppDirspor1 <- sporulation1_distinct %>%
  mutate(estimation = (abundance.filled*1774445)/1e6) %>% 
  mutate(direct = estimation*aa_directSum.filled) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum.filled) %>% 
  mutate(total = direct + opportunity)

# Sum
sporeTranslationOppDirspor1Sum <- sporeTranslationOppDirspor1 %>%
  group_by(expression) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T), 
            total = sum(total, na.rm = T))

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


#########################
# COGS
#########################
#########################################################################
cogs_dat_gene_cop <- cogs_dat_gene
cogs_dat_sp_cop   <- cogs_dat_sp

gene_mat <- cogs_dat_gene_cop[,10:246]
sp_mat   <- cogs_dat_sp[,19:189]

distance_matrix1 <- vegdist(gene_mat, method = "jaccard", binary = TRUE, na.rm = TRUE)
pco1 <- cmdscale(distance_matrix1, eig = T)

coordinates1 <- as.data.frame(pco1$points)
combined_dist1 <- cbind.data.frame(cogs_dat_gene_cop[,1:9], PCoA1 = coordinates1$V1, 
                                  PCoA2 = coordinates1$V2)

ggplot(na.omit(combined_dist1), aes(x = PCoA1, y = PCoA2, color = sporulating, shape = Spo0A))+
  geom_point(alpha = .5, size = 3)+
  mytheme

# Orientation with genes
distance_matrix2 <- vegdist(sp_mat, method = "jaccard", binary = TRUE, na.rm = TRUE)
pco2 <- cmdscale(distance_matrix2, eig = T)

coordinates2 <- as.data.frame(pco2$points)
combined_dist2 <- cbind.data.frame(cogs_dat_sp_cop[,1:18], PCoA1 = coordinates2$V1, 
                                   PCoA2 = coordinates2$V2)

ggplot(na.omit(combined_dist2), aes(x = PCoA1, y = PCoA2, color = role))+
  geom_point(alpha = .5, size = 3)+
  mytheme

mean(protSeqTidyAbun$aa_directSum, na.rm = T)
mean(protSeqTidyAbun$aa_opportunitySum, na.rm = T)

merge_species <- combined_dist2 %>%
  left_join(protSeqTidyAbun, by = "protID") %>%
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254))%>%
  mutate(translationAll = abundance.filled*(1774445/1e6)) %>%
  mutate(translationDirect = translationAll*aa_directSum.filled) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*aa_opportunitySum.filled) %>%
  mutate(translationTotal = translationDirect + translationOpportunity) %>%
  mutate(transcriptionAll = (abundance.filled/1e2)*(1774445/1e6)*as.numeric(gene_length.filled)) %>%
  mutate(transcriptionDirect = transcriptionAll*(10+(2*12*1))) %>% #assuming that mRNAs transcribed at least 1 hour
  mutate(transcriptionOpportunity = transcriptionAll*31) %>%
  mutate(transcriptionTotal = transcriptionDirect + transcriptionOpportunity) %>%
  mutate(RepOpportunity = 2*as.numeric(gene_length.filled)*35) %>% #2 = doublestring  35 = nucleotide costs 
  mutate(RepDirect = 2*as.numeric(gene_length.filled)*14) %>%
  mutate(RepTotal = RepOpportunity+RepDirect) 
  
merge_species$costs = as.numeric(merge_species$transcriptionTotal)+ 
    as.numeric(merge_species$translationTotal)+as.numeric(merge_species$RepTotal)
  
ggplot(na.omit(merge_species), aes(x = PCoA1, y = PCoA2, color = role, size = log10(costs)))+
  geom_point(alpha = .2)+
  mytheme 

ggplot(na.omit(merge_species), aes(x = PCoA1, y = PCoA2, color = role, size = Per_in_spore.formers))+
  geom_point(alpha = .2)+
  mytheme 

# Corelation costs vs. PCoA1

ggplot(na.omit(merge_species), aes(x = log10(costs), y = PCoA1))+
  geom_point(alpha = .2)+
  mytheme

# Subset spore formers 
cogs_dat_gene_copS <- cogs_dat_gene_cop %>%
  filter(sporulating == "N")
species <- c(cogs_dat_gene_copS$regulator)

names.useN <- sp_mat[sp_mat %in% species]

sp_mat.subsetN <- sp_mat[sp_mat %in% names.useN]


distance_matrixN <- vegdist(sp_mat.subsetN, method = "jaccard", binary = TRUE, na.rm = TRUE)

pco3 <- cmdscale(distance_matrixN)



coordinates2 <- as.data.frame(pco2$points)
combined_dist2 <- cbind.data.frame(cogs_dat_sp_cop[,1:18], PCoA1 = coordinates2$V1, 
                                   PCoA2 = coordinates2$V2)

ggplot(na.omit(combined_dist2), aes(x = PCoA1, y = PCoA2, color = role))+
  geom_point(alpha = .5, size = 3)+
  mytheme