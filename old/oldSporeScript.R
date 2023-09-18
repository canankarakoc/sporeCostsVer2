# Spore Costs v3
# 24 January 2023 - last update
# Author: C. Karakoc
# Global expression data from SporeWeb: https://sporeweb.molgenrug.nl/
# Single gene deletion library from: https://doi.org/10.1016/j.cels.2016.12.013
# Global protein abundance data: https://pax-db.org/
# List of gene categories & annotation: subtiwiki.uni-goettingen.de/

# packages
library(tidyverse)       
library(patchwork)

setwd("~/GitHub/sporeCostsVer2")

# ggplot theme
library(ggsci)

mytheme<- theme_bw()+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=18,face="bold"),
        legend.text = element_text(size=14),
        legend.background = element_blank(),
        legend.title = element_text(size=14,face="bold"),
        plot.title = element_text(size=18, face="bold", hjust=0.5),
        strip.text = element_text(size=14, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Data & proteins

# Gene&protein length from SubtiWiki
annotationData <- read.table("./Other_data/subtiwiki.gene.export.2022-05-11.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Regulator list - sigma factors #Daniel
sigma          <- read.table("./Other_data/regulations.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Data with synonyms
nameMap        <- read.table("./Other_data/nameMap.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Protein abundances from PaxDB
protAbun       <- read.table("./Other_data/protAbunData.csv", sep = ',', dec = ".", header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Single gene deletion library from Zoo et al. 2018
deletionLib    <- read.table("./Other_data/mmc7_Koo_etal_2017_TableS6A.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))
mutantLib      <- read.table("./Other_data/mmc7_Koo_etal_2017_TableS6B.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))
conservedLib   <- read.table("./Other_data/mmc7_Koo_etal_2017_TableS6E.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Other traits and subsetted sporulation genes #SubtiWiki lists 
otherTraits           <- read.table("./Other_data/traits.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Germination time course - Swarge et al. 2020 - Protein presence
germinationTimecourse <- read.table("./Other_data/mmc2ProteomeSwargeEtAl2020.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Expression data from SporeWeb
files  <- list.files(path = "./SporeWebHeatmaps/" , pattern = "*.csv", full.names = T)
exp_files <-  list()
for (i in 1:length(files)){
  exp_files[[i]] <- read.table(files[i], header = T, sep = ",", dec = ".")
}
merged_exp <- bind_rows(exp_files, .id = 'sourceID')

merged_exp$sourceID <- as.factor(merged_exp$sourceID)
levels(merged_exp$sourceID ) <- c("1.vegetative", "2.starvation", "3.onset", "4.commitment", "5.engulfment")

# Expression time series
expressionLong <- merged_exp[,c(1, 3:5, 7:14)] %>%
  pivot_longer(cols = c('t1','t2','t3','t4','t5','t6','t7','t8'),
               names_to =  "time", values_to = "expression") %>%
  group_by(time, sourceID, regulators, gene, locus_tag) %>%
  summarise(meanexp = mean(expression)) %>%
  ungroup() %>%
  filter(meanexp>=0.1) %>%
  distinct(locus_tag, .keep_all = TRUE) %>%
  mutate(hours = gsub('t', '', time))

#*******************************************************************************************************************
# Merge abundances
# Protein abundance data does not include locus tags
# It is often easier to merge data sets with locus tags, because genes have
# a lot of synonyms

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

# Now use this collapsed matching column to merge abundances, genes and locus tag
# Merge data
mergedAbunData <- protAbun  %>%
  right_join(nameMap, by = "gene") %>%
  left_join(annotationData[-2], by = "locus_tag") %>%
  left_join(sigma[-4], by = "locus_tag") %>%
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

# Merge with expression data
mergedExpData_time <-  expressionLong %>%
  left_join(mergedAbunData[,-2], by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(as.numeric(gene_length), 765)) %>%
  mutate(protein_length.filled = replace_na(as.numeric(protein_length), 254))

# Alternative germination timecourse data
mergedGermData_time <-  germinationTimecourse %>%
  left_join(mergedAbunData, by = "gene") %>%
  pivot_longer(cols = c('t30','t60','t90','t120','t150'),
               names_to =  "time", values_to = "presence") %>%
  filter(presence > 0) %>% #present genes
  distinct(locus_tag, .keep_all = TRUE) %>% #based on first appearance 
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(as.numeric(gene_length, 765))) %>%
  mutate(protein_length.filled = replace_na(as.numeric(protein_length, 254)))%>%
  mutate(minutes = gsub('t', '', time))

#********************************************************************************************
# Transcription costs
# 1 mRNA can yield to 100-1000 proteins (Cell biology by the numbers)

transcriptCosts <- mergedExpData_time %>%
  mutate(estimation = (abundance.filled/1000)*(1774445/1000000)*gene_length.filled) %>% #protein abundance/100 X 1.8 X gene length
  # the reason I multiply with 1.8 is that the protein abundance is reported as parts per million. An average size bacteria has about 3 million protein molecules
  # this was reported experimentally as average 1774445 in Bacillus (Maass et.al. 2011)
  group_by(time) %>%
  mutate(direct = estimation*(10+(2*12))) %>% # average sporulation time is 8hours, we would expect 60 minutes / 5 minutes (degradation rate) 
  # = 12 re-polymerization events per hour
  # assuming nucleotides are well recycled and it only affects polymerization costs
  # however the cells do not divide, neither grow so the proteins dilute. I am not sure if we have to
  # include degradation rate.
  mutate(opportunity = estimation*31) %>%
  mutate(total = direct + opportunity)%>%
  mutate(numGenes = length(time))

# Cumulative costs
transcriptCosts_sum <- transcriptCosts %>%
  group_by(time) %>%
  summarise(sumTotal = sum(total, na.rm =T)) %>%
  mutate(type = rep('transcription'))

# Cumulative costs - Regulons
transcriptCosts_sum_reg<- transcriptCosts %>%
  group_by(time, regulator) %>%
  summarise(sumTotal = sum(total, na.rm =T)) %>%
  mutate(type = rep('transcription'))

sum(transcriptCosts_sum$sumTotal[1:8]) #8874479


# Translation costs
translationCosts <- mergedExpData_time %>%
  mutate(estimation = (abundance.filled)*(1774445/1000000)*protein_length.filled) %>% #protein abundance/100 X 1.8 X gene length
  # the reason I multiply with 1.8 is that the protein abundance is reported as parts per million. An average size bacteria has about 3 million protein molecules
  # this was reported experimentally as average 1774445  in Maass et.al. 2011 Bacillus
  mutate(direct = estimation*(4+2)) %>% # ignoring protein degradation
  mutate(opportunity = estimation*24) %>%
  mutate(total = direct + opportunity) %>%
  group_by(sourceID) %>%
  mutate(num_genes = length(time)) 

# Cumulative costs
translationCosts_sum <- translationCosts %>%
  group_by(time) %>%
  summarise(sumTotal = sum(total, na.rm = T)) %>%
  mutate(type = rep('translation'))

# Cumulative costs - Regulons
translationCosts_sum_reg<- translationCosts %>%
  group_by(time, regulator) %>%
  summarise(sumTotal = sum(total, na.rm = T)) %>%
  mutate(type = rep('transcription'))

sum(translationCosts_sum$sumTotal[1:8]) #1359813753


# Total costs 
combined <- cbind.data.frame(time = transcriptCosts_sum$time, 
                             sumTotal = transcriptCosts_sum$sumTotal+translationCosts_sum$sumTotal)

#SigE: early mother cell-specific sporulation sigma factor
#SigF: early forespore-specific sporulation sigma factor
#SigG: late forespore-specific sporulation sigma factor
#SigH: sigma factor that controls genes of the transition phase
#SigK: late mother cell-specific sporulation sigma factor

combined_reg <- cbind.data.frame(time  = transcriptCosts_sum_reg$time, 
                                 regulator = transcriptCosts_sum_reg$regulator,
                                 sumTotal  = transcriptCosts_sum_reg$sumTotal+
                                   translationCosts_sum_reg$sumTotal)

selected <- c("SigE", "SigF", "SigG", "SigH", "SigK")
sigmaFactors <- combined_reg[combined_reg$regulator %in% selected,]

#SigE: early mother cell-specific sporulation sigma factor
#SigF: early forespore-specific sporulation sigma factor
#SigG: late forespore-specific sporulation sigma factor
#SigH: sigma factor that controls genes of the transition phase
#SigK: late mother cell-specific sporulation sigma factor

# add the replication costs to the time t1
combined$sumTotal[1] <- combined$sumTotal[1]+388671200
sum(combined$sumTotal[1:8]) #1757359431

#Figure 1B Model addition
time        <- c(1:8)
sporulation <- as.numeric(combined$sumTotal)
model.data  <- cbind.data.frame(time, sporulation)
fit         <- nls(sporulation ~ SSasymp(time, yf, y0, log_alpha), data = model.data)

coef(fit) 

tt <- seq(1,8, by = 0.1)
pred <- predict(fit, list(time = tt))
preddata <- cbind.data.frame(pred, tt)

RSS.p <- sum(residuals(fit)^2)
y = as.numeric(as.character(sumTotal))
TSS   <- sum((log10(y) - mean(y))^2)
1 - (RSS.p/TSS)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Plot without germination
ggplot(NULL, aes(x = x, y = y))+
  geom_bar(data = combined, 
           aes(x = time, y = sumTotal), position="dodge", stat="identity", color = 'grey25')+
  ylab("Costs in units of ATP")+
  xlab("Time (hrs)")+
  geom_line(data = preddata, aes(x = tt, y = pred))+
  mytheme+
  scale_y_continuous(labels = scientific_10)
### END Figure 1B ### 

### Regulons ###

# Cumulative costs - Regulons
sigmaFactors_spo_sum  <- sigmaFactors  %>%
  group_by(regulator) %>%
  summarise(sumTotal = sum(sumTotal))

### Supplementary Figure ###
ggplot(sigmaFactors_spo_sum, aes(x = reorder(regulator, -sumTotal), y = as.numeric(as.character(sumTotal))),)+
  geom_bar(position="dodge", stat="identity")+
  ylab("Costs in units of log(ATP)")+
  xlab("Time(h)")+
  mytheme+
  scale_y_continuous(labels = scientific_10)
### END Supplementary Figure ###

#******************************************************************************************
# Timecourse of germination 
# Protein presence data 
# Transcription costs
# 1 mRNA can yield to 100-1000 proteins (Cell biology by the numbers)

transcriptCosts_ger <- mergedGermData_time %>%
  mutate(estimation = (abundance.filled/1000)*(1774445/1000000)*gene_length.filled) %>% #protein abundance/100 X 1.8 X gene length
  # the reason I multiply with 1.8 is that the protein abundance is reported as parts per million. An average size bacteria has about 3 million protein molecules
  # this was reported experimentally as average 1774445 in Maass et.al. 2011 Bacillus
  group_by(minutes) %>%
  mutate(direct = estimation*(10+(2*12*0.5))) %>% # we would expect 60 minutes / 5 minutes (degradation rate). Germination has been measured every 30 min.
  # = 12 re-polymerization events per hour.
  mutate(opportunity = estimation*31) %>%
  mutate(total = direct + opportunity)%>%
  mutate(numGenes = length(minutes))

# Cumulative costs
transcriptCosts_sum_ger <- transcriptCosts_ger %>%
  group_by(minutes) %>%
  summarise(sumTotal = sum(total, na.rm = T)) %>%
  mutate(type = rep('transcription'))

# Translation costs
translationCosts_ger <- mergedGermData_time %>%
  mutate(estimation = (abundance.filled)*(1774445/1000000)*protein_length.filled) %>% #protein abundance/100 X 1.8 X gene length
  # the reason I multiply with 1.8 is that the protein abundance is reported as parts per million. An average size bacteria has about 3 million protein molecules
  # this was reported experimentally as average 1774445  in Maass et.al. 2011 Bacillus
  group_by(minutes) %>%
  mutate(direct = estimation*(4+2)) %>% # ignoring protein degradation
  mutate(opportunity = estimation*24) %>%
  mutate(total = direct + opportunity) 

# Cumulative costs
translationCosts_sum_ger <- translationCosts_ger %>%
  group_by(minutes) %>%
  summarise(sumTotal = sum(total, na.rm = T)) %>%
  mutate(type = rep('translation'))

# Total costs 
combined_ger <- cbind.data.frame(time = transcriptCosts_sum_ger$minutes, 
                                 sumTotal = transcriptCosts_sum_ger$sumTotal+translationCosts_sum_ger$sumTotal)

### Figure 1C ###
ggplot(combined_ger, aes(x = as.numeric(time), y = sumTotal))+
  geom_bar(position="dodge", stat="identity", color = "grey25")+
  ylab("Costs in units of ATP")+
  xlab("Time(min)")+  
  mytheme+
  scale_y_continuous(labels = scientific_10)+
  scale_fill_npg()
### END Figure 1C ###

#***************************************************************************************************
# Now merge this data with deletion library

mergedAbunData$gene_length    <- as.numeric(mergedAbunData$gene_length)
mergedAbunData$protein_length <- as.numeric(mergedAbunData$protein_length)

mergedDelData <- deletionLib[-2] %>%
  left_join(mergedAbunData, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254)) %>%
  mutate(success = case_when(
    rSS_average < 0.31 ~ "essential", # assume these genes are needed for successful sporulation
    rSS_average > 0.31 & rSS_average < 1.32 ~ "non-essential")) %>% # these genes are affecting sporulation (75%), but not esential
  filter(!is.na(success))

# Transcription costs
transcriptCosts_del <- mergedDelData %>%
  mutate(estimation = (abundance.filled/1000)*(1774445/1000000)*gene_length.filled) %>%
  mutate(direct = estimation*(10+(2*12*1))) %>% #assuming that mRNAs transcribed at least 1 hour
  mutate(opportunity = estimation*31) %>%
  mutate(total = direct + opportunity) %>%
  mutate(type = rep('transcription')) %>% 
  group_by(success) %>% 
  mutate(no_genes = length(success)) 

# Translation costs
translationCosts_del <- mergedDelData %>%
  mutate(estimation = (abundance.filled)*(1774445/1000000)*protein_length.filled) %>% # protein abundance X protein length
  # again the reason I multiply with 3 is that the protein abundance is reported as parts per million, and an average size bacteria has about 3 million protein molecules
  mutate(direct = estimation*(4+2)) %>% # ignoring protein degradation
  mutate(opportunity = estimation*24) %>%
  mutate(total = direct + opportunity) %>%
  mutate(type = rep('translation')) %>% 
  group_by(success) %>% 
  mutate(no_genes = length(success)) 

combined_del_sum <- transcriptCosts_del %>%
  bind_rows(translationCosts_del) %>%
  group_by(success, type) %>%
  
  
  summarize(direct = sum(direct, na.rm = T),
            opportunity = sum(opportunity, na.rm = T),
            total = sum(total, na.rm = T),
            no_genes = mean(no_genes)) %>%
  mutate(norm = total/no_genes)

### Figure 2 ###
ggplot(combined_del_sum, aes(x = success, y = total, color = type, fill = type))+
  geom_point(color = "grey25")+
  ylab("Costs in units of log(ATP)")+
  xlab("Genes")+  
  mytheme+
  scale_fill_npg()+
  coord_cartesian(ylim = c(6,10))
### Figure 2 ###

ggplot(translationCosts_del, aes(y = rSS_average, x = log10(total), color = success))+
  geom_point(alpha = .5)+
  geom_violin(alpha =.5)+
  xlab("Costs in units of log(ATP)")+
  ylab("Relative sporulation success")+  
  mytheme+
  scale_color_npg()

ggplot(translationCosts_del, aes(x = log10(total), color = success))+
  geom_density()+
  ylab("Frequency")+
  xlab("Relative sporulation success")+  
  mytheme+
  scale_color_npg()
#*******************************************************************************************
# Sporulation efficiency
efficiency   <- read.table("./Other_data/efficiency.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

ggplot(efficiency, aes(efficiency))+
  geom_density()+
  xlab("Sporulation efficiency")+
  ylab("Frequency")+
  geom_vline(xintercept = 34.7, color = "#A42820", linetype = "dashed" )+
  coord_cartesian(ylim=c(0.0025, 0.0125))+
  mytheme

library(truncnorm)
fit <- density(efficiency$efficiency, from = 0, to = 100)
N <- 1e6
x.new <- rtruncnorm(a = 0, b = 100, N, sample(efficiency$efficiency, size = N, replace = TRUE))
plot(density(x.new, bw = fit$bw))

densityPlot <- as.data.frame(x.new)

means = densityPlot %>%
  summarise(M = mean(x.new), SD = sd(x.new), N = n())

#Sample size of 50. T Distribution intervals
means$error=qt(0.975, df=50-1)*means$SD/sqrt(means$N)
means$upper=means$M+means$error
means$lower=means$M-means$error

### Figure 3 ###
ggplot(densityPlot, aes(x = x.new))+
  geom_density(alpha=0.1, size=0.8, bw = 11.35555)+
  geom_vline(xintercept = c(means$M), linetype = 'dashed', color = "#A42820")+
  annotate("text", x = 36, y=0.013, label= "Mean = 34.8%", color = "#A42820")+
  labs(x="Sporulation efficiency", y = "Density")+
  mytheme+
  scale_x_continuous(limits=c(-30,130))
### Figure 3 ###

#********************************************************************************************
# Merge trair data
mergedTraitData <- otherTraits %>%
  left_join(annotationData, by = "gene") %>%
  left_join(mergedAbunData[,c("locus_tag", "abundance")], by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 558)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 251)) %>%
  group_by(category) %>%
  mutate(no_genes = length(category))

# All costs
mergedTraitData$protein_length.filled <- as.numeric(mergedTraitData$protein_length.filled)
mergedTraitData$gene_length.filled <- as.numeric(mergedTraitData$gene_length.filled)

totalCosts_traits <- mergedTraitData %>%
  mutate(translationAll = abundance.filled*(1774445/1000000)*protein_length.filled) %>%
  mutate(translationDirect = translationAll*(4+2)) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*24) %>%
  mutate(translationTotal = translationDirect + translationOpportunity) %>%
  mutate(transcriptionAll = (abundance.filled/1000)*(1774445/1000000)*gene_length.filled) %>%
  mutate(transcriptionDirect = transcriptionAll*(10+(2*12*1))) %>% #assuming that mRNAs transcribed at least 1 hour
  mutate(transcriptionOpportunity = transcriptionAll*31) %>%
  mutate(transcriptionTotal = transcriptionDirect + transcriptionOpportunity) %>%
  mutate(costs = transcriptionTotal + translationTotal)

# Cumulative costs
costs_sum_traits <- totalCosts_traits %>%
  group_by(category) %>%
  summarize(sumCosts = sum(costs, na.rm = T),
            number_genes = mean(no_genes)) %>%
  select(category, sumCosts, number_genes) %>%
  add_row(category = "growth", sumCosts = 26000000000, number_genes = 0) %>%
  add_row(category = "basal_metabolism", sumCosts = 350000000, number_genes = 0)

# Add replication costs to sporulation 
costs_sum_traits$sumCosts[8] <- costs_sum_traits$sumCosts[8]+388671200

# Plot cumulative costs
ggplot(costs_sum_traits, aes(x = log10(sumCosts), y = reorder(category, sumCosts)))+
  geom_col(fill="grey", color="black")+
  xlab("Log(Total translation costs in units of ATP)")+
  ylab("Complex traits")+
  coord_cartesian(xlim = c(7.8,10.5))+
  mytheme

# Add a second axis 

costs_sum_traits_rel <- costs_sum_traits %>%
  mutate(relative = (932926184/sumCosts)*100)

### Figure4 ###
ggplot(costs_sum_traits_rel, aes(x = log10(sumCosts), y = reorder(category, sumCosts)))+
  geom_col(fill="grey", color="black")+
  geom_vline(xintercept = log10(932926184), linetype = 'dashed')+
  # Custom the Y scales:
  scale_x_continuous(
    
    # Features of the first axis
    name = "Total costs in units of log(ATP)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~log10((10^./932926184)*100), name="% Costs realative to sporulation")
  )+
  coord_cartesian(xlim = c(7.8,10.5))+
  mytheme
### Figure4 ###

#**********************************************************************************************
# Conserved sporulation genes
# After Acidamicoccus fermentans DSM 20731 non-sporulating species, add an id column
col = rep(c('sporulating', 'non-sporulating'), times = c(21, 19))

mergedOrtData <- conservedLib %>%
  left_join(mergedAbunData, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 558)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 251)) %>%
  mutate(presence_occurrence = str_count(text, present)) 

names(mergedOrtData[,4:43]) <- col

# Absence/presence of essential genes

#mergedOrtData$presence <- apply(mergedOrtData[,4:43], 1, table)
mergedOrtData$present <- apply(mergedOrtData[,4:43], 1, function(x) length(which(x=='present')))
mergedOrtData$absent  <- apply(mergedOrtData[,4:43], 1, function(x) length(which(x=='absent')))
mergedOrtData$freq    <- (mergedOrtData$presence/40)*100

### Supplementary Figures ###
ggplot(mergedOrtData, aes(y = freq, x = gene_length.filled))+
  geom_jitter(alpha=0.5)+
  geom_smooth(method = 'lm')+
  ylab("Percentage of presence")+
  xlab("Gene length")+
  mytheme+
  scale_fill_npg()

ggplot(mergedOrtData, aes(y = freq, x = abundance.filled))+
  geom_jitter(alpha=0.5)+
  geom_smooth(method = 'lm')+
  ylab("Percentage of presence")+
  xlab("Protein abundance")+
  mytheme+
  scale_fill_npg()+
  scale_y_continuous(limits=c(0,40))


mergedOrtData_costs <- mergedOrtData %>%
  mutate(translationAll = abundance.filled*(1774445/1000000)*protein_length.filled) %>%
  mutate(translationDirect = translationAll*(4+2)) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*24) %>%
  mutate(translationTotal = translationDirect + translationOpportunity) %>%
  mutate(transcriptionAll = (abundance.filled/1000)*(1774445/1000000)*gene_length.filled) %>%
  mutate(transcriptionDirect = transcriptionAll*(10+(2*12*1))) %>% #assuming that mRNAs transcribed at least 1 hour
  mutate(transcriptionOpportunity = transcriptionAll*31) %>%
  mutate(transcriptionTotal = transcriptionDirect + transcriptionOpportunity) %>%
  mutate(costs = transcriptionTotal + translationTotal)

ggplot(mergedOrtData_costs, aes(y = freq, x = log10(costs)))+
  geom_jitter(alpha= 0.5)+
  geom_smooth(method="lm")+
  ylab("Percentage of presence among taxa")+
  xlab("Costs of sporulation genes log(ATP)")+
  mytheme+
  scale_fill_npg()

fit.lm <- lm(data = mergedOrtData_costs, freq~log10(costs))

summary(fit.lm)

#***********************************************************************************
# Cost of being a spore 
# Too high, estimated for new spores
# Supplementary if necessary 

mergedSporeData <-  sporeTranscript %>%
  left_join(mergedAbunData, by = "gene") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 558)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 251))

# Transcription costs
transcriptCosts_spore <- mergedSporeData %>%
  mutate(estimation = (abundance.filled/1000)*(73411/1000000)*gene_length.filled) %>%
  mutate(direct = estimation*(10+(2))) %>% #assuming that mRNAs transcribed at least 1day long 
  mutate(opportunity = estimation*31) %>%
  mutate(total = direct + opportunity) %>%
  mutate(type = rep('transcription'))

# Translation costs
translationCosts_spore <- mergedSporeData %>%
  mutate(estimation = (abundance.filled)*(73413/1000000)*protein_length.filled) %>% # protein abundance X protein length
  # again the reason I multiply with 3 is that the protein abundance is reported as parts per million, and an average size bacteria has about 3 million protein molecules
  mutate(direct = estimation*(4+2)) %>% # ignoring protein degradation
  mutate(opportunity = estimation*24) %>%
  mutate(total = direct + opportunity) %>%
  mutate(type = rep('translation'))

combined_spore_sum <- transcriptCosts_spore %>%
  bind_rows(translationCosts_spore) %>%
  summarize(direct = sum(direct, na.rm = T),
            opportunity = sum(opportunity, na.rm = T),
            total = sum(total, na.rm = T)) 

#***********************************************************************************
