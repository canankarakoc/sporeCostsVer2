# Spore Costs v3
# 20 January 2022 - last update
# Author: C. Karakoc
# Global expression data from SporeWeb https://sporeweb.molgenrug.nl/
# Single gene deletion library from https://doi.org/10.1016/j.cels.2016.12.013
# Global protein abundance data https://pax-db.org/

library(tidyverse)
library(RColorBrewer)

setwd("~/OneDrive - Indiana University/SporeCost/Test_calculations")

# ggplot theme

mycolors    <- c("#74828f", "#c25b56", "#d2b274", "#74867b", "grey")
mycolors2   <- c("#E1BD6D", "#EABE94", "#0B775E", "#35274A" ,"#F2300F")
mycolors3   <- c("#798E87", "#C27D38", "#CCC591")
mycolors3.1 <- c("#C27D38", "#CCC591")
mycolors4   <- c("grey10", "grey55")
mycolors4.1 <- c("grey25", "grey75")

n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

mytheme<- theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12,face="bold"),
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.title = element_text(size=12,face="bold"),
        plot.title = element_text(size=14, face="bold", hjust=0.5),
        strip.text = element_text(size=12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Data & proteins

# Gene&protein length from SubtiWiki
annotationData <- read.table("./Other_data/subtiwiki.gene.export.2021-06-25.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Data with synonyms
nameMap        <- read.table("./Other_data/nameMap.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Protein abundances from PaxDB
protAbun       <- read.table("./Other_data/protAbunData.csv", sep = ',', dec = ".", header = T, stringsAsFactors = F)

# Single gene deletion library from Zoo et al. 2018
deletionLib    <- read.table("./Other_data/mmc7_Koo_etal_2017_TableS6A.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)
mutantLib      <- read.table("./Other_data/mmc7_Koo_etal_2017_TableS6B.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)
conservedLib   <- read.table("./Other_data/mmc7_Koo_etal_2017_TableS6E.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Found synonyms to match the databases
# synonyms01    <- read.table("./Found_synonyms/synonyms.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Expression data from SporeWeb

files  <- list.files(path = "./SporeWebHeatmaps/" , pattern = "*.csv", full.names = T)
exp_files        <-  list()
for (i in 1:length(files)){
  exp_files[[i]] <- read.table(files[i], header = T, sep = ",", dec = ".")
}
merged_exp <- bind_rows(exp_files, .id = 'sourceID')

merged_exp$sourceID <- as.factor(merged_exp$sourceID)
levels(merged_exp$sourceID ) <- c("1.vegetative", "2.starvation", "3.onset", "4.commitment", "5.engulfment")

# Regulators time series

merged_exp2 <- merged_exp
names(merged_exp2) <- c("sourceID", "Class","regulators","gene","locus_tag","0","1",
                        "2","3","4","5","6","7","8","135","150","180")

expressionLong <- merged_exp2[,c(1, 3, 6:14)] %>%
  pivot_longer(cols = c('0','1','2','3','4','5','6','7','8'),
             names_to =  "time", values_to = "expression") %>%
  group_by(time, sourceID, regulators) %>%
  summarise(meanexp = mean(expression))

# Group genes based on regulons

expressionLong_vegetative   <- expressionLong[expressionLong$sourceID == "1.vegetative",]
expressionLong_starvation   <- expressionLong[expressionLong$sourceID == "2.starvation",]
expressionLong_onset        <- expressionLong[expressionLong$sourceID == "3.onset",]
expressionLong_commitment   <- expressionLong[expressionLong$sourceID == "4.commitment",]
expressionLong_engulfment   <- expressionLong[expressionLong$sourceID == "5.engulfment",]

# Plot time series óf expression
#vegetation
ggplot(expressionLong_vegetative, aes(x = as.numeric(time) , y = meanexp, group = regulators, color = regulators))+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_line()+
  ylab("normalized expression")+
  xlab("time (hours after initiation)")+
  ggtitle('vegetative')+
  mytheme+
  scale_color_manual(values = col_vector)

#starvation
ggplot(expressionLong_starvation, aes(x = as.numeric(time) , y = meanexp, group = regulators, color = regulators))+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_line()+
  ylab("normalized expression")+
  xlab("time (hours after initiation)")+
  ggtitle('starvation')+
  mytheme+
  scale_color_manual(values = col_vector)


#onset
ggplot(expressionLong_onset, aes(x = as.numeric(time) , y = meanexp, group = regulators, color = regulators))+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_line()+
  ylab("normalized expression")+
  xlab("time (hours after initiation)")+
  ggtitle('onset')+
  mytheme

#commitment
ggplot(expressionLong_commitment, aes(x = as.numeric(time) , y = meanexp, group = regulators, color = regulators))+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_line()+
  ylab("normalized expression")+
  xlab("time (hours after initiation)")+
  ggtitle('commitment')+
  mytheme+
  scale_color_manual(values = col_vector)

#engulfment
ggplot(expressionLong_engulfment, aes(x = as.numeric(time) , y = meanexp, group = regulators, color = regulators))+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_line()+
  ylab("normalized expression")+
  xlab("time (hours after initiation)")+
  ggtitle('engulfment')+
  mytheme+
  scale_color_manual(values = col_vector)

# Subset the upregulated genes, remove the duplicate genes, keep the first one

merged_exp_positive <- merged_exp[,c(1, 3, 5, 7:14)] %>%
  pivot_longer(names_to = "time", values_to = "expression_relative",
               cols = t1:t8) %>%
  filter_at(vars(contains("t")), all_vars(.>0)) %>%
  group_by(sourceID, regulators, locus_tag) %>%
  summarise(mean_pos_exp = mean(expression_relative))%>%
  distinct(locus_tag, .keep_all = TRUE)

# Both upregulated & downregulated genes

merged_exp_all <- merged_exp[,c(1, 3, 5, 7:14)] %>%
  pivot_longer(names_to = "time", values_to = "expression_relative",
               cols = t1:t8) %>%
  group_by(sourceID, regulators, locus_tag) %>%
  summarise(mean_pos_exp = mean(expression_relative)) %>%
  distinct(locus_tag, .keep_all = TRUE)

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

# number of NAs

# Now use this collapsed matching column to merge abundances, genes and locus tag

# Merge data
mergedAbunData <- nameMap %>%
  left_join(protAbun, by = "gene") %>%
  left_join(annotationData, by = "locus_tag") %>%
  select(locus_tag, abundance, protein_length, gene_length)

# Fill NAs with average values
# Find better estimates?

# average protein length:
# median 256aa, BNID 106444, Milo et al. 2010

# average gene length
# mean 924±9, BNID 111922, Milo et al. 2010

# average protein abundance length
mean(protAbun$abundance) #252
median(protAbun$abundance) #17

# Merge with expression data
mergedExpData <-  merged_exp_positive %>%
  left_join(mergedAbunData, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 252)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 924)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 256))

# 99 missing abundance data out of 1494
# 18 missing length data out of 1494

#********************************************************************************************
# Transcription costs

# 1 mRNA can yield to 100-1000 proteins (Cell biology by the numbers)

transcriptCosts <- mergedExpData %>%
  mutate(transcriptAll = (abundance.filled/1000)*(1774445/1000000)*gene_length.filled) %>% #protein abundance/100 X 1.8 X gene length
  # the reason I multiply with 1.8 is that the protein abundance is reported as parts per million. An average size bacteria has about 3 million protein molecules
  # this was reported experimentally as average 1774445 in Maass et.al. 2011 Bacillus
  group_by(sourceID) %>%
  mutate(transcriptDirect = transcriptAll*(10+(2*12*8))) %>% # average sporulation time is 8hours, we would expect 8*60 minutes / 5 minutes (degradation rate) = 96 re-polymerization events
  # assuming nucleotides are well recycled and it only affects polymerization costs
  # however the cells do not divide, neither grow so the proteins dilute. I am not sure if we have to
  # include degradation rate.
  mutate(transcriptOpportunity = transcriptAll*31) %>%
  mutate(transcriptTotal = transcriptDirect + transcriptOpportunity)

# Divide into regulons

transcriptCosts_vegetative  <- transcriptCosts[transcriptCosts$sourceID == "1.vegetative",]
transcriptCosts_starvation  <- transcriptCosts[transcriptCosts$sourceID == "2.starvation",]
transcriptCosts_onset       <- transcriptCosts[transcriptCosts$sourceID == "3.onset",]
transcriptCosts_commitment  <- transcriptCosts[transcriptCosts$sourceID == "4.commitment",]
transcriptCosts_engulfment  <- transcriptCosts[transcriptCosts$sourceID == "5.engulfment",]

# Cumulative costs

typesCost <- c('transcriptDirect', 'transcriptOpportunity', 'transcriptTotal')
tv <- apply(transcriptCosts_vegetative[,typesCost], 2, sum)
ts <- apply(transcriptCosts_starvation[,typesCost], 2, sum)
to <- apply(transcriptCosts_onset[,typesCost], 2, sum)
tc <- apply(transcriptCosts_commitment[,typesCost], 2, sum)
te <- apply(transcriptCosts_engulfment[,typesCost], 2, sum)

# Dataset to plot cumulative costs

costs     <- c(tv, ts, to, tc, te)
type      <- rep(c('direct', 'opportunity', 'total'), times= 5)
stages    <- rep(c('1.vegetative', '2.starvation', '3.onset', '4.commitment', '5.engulfment'), each =3)
plotCosts <- cbind.data.frame(costs, type, stages)

# Plot cumulative costs

ggplot(plotCosts, aes(x = stages, y = log10(costs), color = type))+
  geom_point(size = 3, alpha = 0.8)+
  ylab("log(Costs) in units of ATP")+
  xlab("stages")+
  ggtitle('cumulative costs transcription')+
  mytheme+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_color_manual(values = mycolors)

# Plot cumulative costs

ggplot(plotCosts, aes(x = stages, y = costs, fill = type))+
  geom_bar(position="dodge", stat="identity")+
  #ylab("log(Costs) in units of ATP")+
  #xlab("stages")+
  ggtitle('cumulative costs transcription')+
  mytheme+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_manual(values = mycolors)

#************************************************************************************************
# Translation costs

translationCosts <- mergedExpData %>%
  mutate(translationAll = (abundance.filled)*(1774445/1000000)*protein_length.filled) %>% #protein abundance/100 X 1.8 X gene length
  # the reason I multiply with 1.8 is that the protein abundance is reported as parts per million. An average size bacteria has about 3 million protein molecules
  # this was reported experimentally as average 1774445  in Maass et.al. 2011 Bacillus
  mutate(translationDirect = translationAll*(4+2)) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*24) %>%
  mutate(translationTotal = translationDirect + translationOpportunity)


translationCosts_vegetative   <- translationCosts[translationCosts$sourceID == "1.vegetative",]
translationCosts_starvation   <- translationCosts[translationCosts$sourceID == "2.starvation",]
translationCosts_onset        <- translationCosts[translationCosts$sourceID == "3.onset",]
translationCosts_commitment   <- translationCosts[translationCosts$sourceID == "4.commitment",]
translationCosts_engulfment   <- translationCosts[translationCosts$sourceID == "5.engulfment",]

# Cumulative costs

typesCost_p <- c('translationDirect', 'translationOpportunity', 'translationTotal')
pv <- apply(translationCosts_vegetative[,typesCost_p], 2, sum)
ps <- apply(translationCosts_starvation[,typesCost_p], 2, sum)
po <- apply(translationCosts_onset[,typesCost_p], 2, sum)
pc <- apply(translationCosts_commitment[,typesCost_p], 2, sum)
pe <- apply(translationCosts_engulfment[,typesCost_p], 2, sum)

# Dataset to plot cumulative costs

costs_p     <- c(pv, ps, po, pc, pe)

type_p      <- rep(c('direct', 'opportunity','total'), times= 5)
stages_p    <- rep(c('1.vegetative', '2.starvation','3.onset','4.commitment', '5.engulfment'), each =3)
plotCosts_p <- cbind.data.frame(costs_p, type_p, stages_p)

# Cumulative costs

ggplot(plotCosts_p, aes(x = stages, y = log10(costs_p), color = type))+
  geom_point(size = 3, alpha = 0.8)+
  ylab("log(Costs) in units of ATP")+
  xlab("stages")+
  ggtitle('cumulative costs translation')+
  mytheme+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_color_manual(name = "", values = mycolors)

ggplot(plotCosts_p, aes(x = stages, y = costs, fill = type))+
  geom_bar(position="dodge", stat="identity")+
  #ylab("log(Costs) in units of ATP")+
  #xlab("stages")+
  ggtitle('cumulative costs translation')+
  mytheme+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_manual(values = mycolors)

# Individual donut plots for the figure
# Total costs transcription+translation+replication (replication only during starvation)

total_costs <- c(419318396+9844847.6, 5107896.0+214241478+388671200, 9931002.6+424155650,
                 27030889.6+1151896186, 21317677.2+907671830)

stages <- c("vegetative", "starvation", "onset", "commitment", "engulfment")

donut_plot <- cbind.data.frame(costs=total_costs, stages=stages)

# Compute percentages
donut_plot$fraction = (donut_plot$costs / sum(donut_plot$costs))*100

# Compute the cumulative percentages (top of each rectangle)
donut_plot$ymax = cumsum(donut_plot$fraction)

# Compute the bottom of each rectangle
donut_plot$ymin = c(0, head(donut_plot$ymax, n=-1))

# Compute label position
donut_plot$labelPosition  <- (donut_plot$ymax + donut_plot$ymin) / 2

# Compute a good label
donut_plot$label  <- paste0(donut_plot$stages, " : ", round(donut_plot$fraction,1))


ggplot(donut_plot, aes(ymax = ymax , ymin = ymin, xmax=4, xmin=3, fill = stages))+
  geom_rect()+
  coord_polar(theta="y")+
  xlim(c(2, 4))+
  geom_label( x = 3.5, aes(y = labelPosition, label = label), size = 4) +
  #ylab("log(Costs) in units of ATP")+
  #xlab("stages")+
  ggtitle('Total costs of sporulation stages')+
  mytheme+
  scale_fill_manual(values = mycolors2)+
  theme_void()+
  theme(legend.position = "none")

# Breakdown of central dogma costs

plotCosts_sum   <- sum((subset(plotCosts,  plotCosts$type=="total")$costs))
plotCosts_p_sum <- sum((subset(plotCosts_p,  plotCosts_p$type=="total")$costs))

cent_costs <- c(388671200, plotCosts_sum, plotCosts_p_sum)
stages <- c("replication", "transcription", "translation")

donut_plot2 <- cbind.data.frame(costs=cent_costs, stages=stages)

# Compute percentages
donut_plot2$fraction = (donut_plot2$costs / sum(donut_plot2$costs))*100

# Compute the cumulative percentages (top of each rectangle)
donut_plot2$ymax = cumsum(donut_plot2$fraction)

# Compute the bottom of each rectangle
donut_plot2$ymin = c(0, head(donut_plot2$ymax, n=-1))

# Compute label position
donut_plot2$labelPosition  <- (donut_plot2$ymax + donut_plot2$ymin) / 2

# Compute a good label
donut_plot2$label  <- paste0(donut_plot2$stages, " : ", round(donut_plot2$fraction,1))

ggplot(donut_plot2, aes(ymax = ymax , ymin = ymin, xmax=4, xmin=3, fill = stages))+
  geom_rect()+
  coord_polar(theta="y")+
  xlim(c(2, 4))+
  geom_label( x = 3.5, aes(y = labelPosition, label = label), size = 4)+
  ggtitle('Total costs of central dogma')+
  mytheme+
  scale_fill_manual(values = mycolors3)+
  theme_void()+
  theme(legend.position = "none")

# Breakdown direct opportunity costs

transcription_sum <- plotCosts %>%
  group_by(type)%>%
  summarize(costs = sum(costs))

translation_sum <- plotCosts_p %>%
  group_by(type_p)%>%
  summarize(costs = sum(costs_p))

# replication costs:
# direct = 107504800
# opportunity = 281166400
# total = 388671200

type_costs <- c(634889580+623456708+107504800, 97433549+2493826832+281166400)
type <- c("direct", "opportunity")

donut_plot3 <- cbind.data.frame(costs = type_costs, type)

# Compute percentages
donut_plot3$fraction = (donut_plot3$costs / sum(donut_plot3$costs))*100

# Compute the cumulative percentages (top of each rectangle)
donut_plot3$ymax = cumsum(donut_plot3$fraction)

# Compute the bottom of each rectangle
donut_plot3$ymin = c(0, head(donut_plot3$ymax, n=-1))

# Compute label position
donut_plot3$labelPosition  <- (donut_plot3$ymax + donut_plot3$ymin) / 2

# Compute a good label
donut_plot3$label  <- paste0(donut_plot3$type, " : ", round(donut_plot3$fraction,1))

ggplot(donut_plot3, aes(ymax = ymax , ymin = ymin, xmax=4, xmin=3, fill = type))+
  geom_rect()+
  coord_polar(theta="y")+
  xlim(c(2, 4))+
  geom_label( x = 3.5, aes(y = labelPosition, label = label), size = 4) +
  ggtitle('Type of costs')+
  mytheme+
  scale_fill_manual(values = mycolors4)+
  theme_void()+
  theme(legend.position = "none")

#***************************************************************************************************
# Now merge this data with deletion library

mergedDelData <- deletionLib[-2] %>%
left_join(mergedAbunData, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 252)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 924)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 256))

# Average sporulation score of all genes tested

SS0   <- min(deletionLib$rSS_average)
SS100 <- max(deletionLib$rSS_average)
SS50  <- mean(deletionLib$rSS_average) # assume these genes are needed for successful sporulation
SS75  <- SS100*0.75 # assume these genes are not essential
SS25  <- SS100*0.25 # assume these genes are  essential
SS1   <- SS100*0.01 # assume these genes are  essential

# Transcription costs

transcriptCosts <- mergedDelData %>%
  mutate(transcriptAll = (abundance.filled/1000)*(1774445/1000000)*gene_length.filled) %>%
  mutate(success = case_when(
    rSS_average > SS75  ~ "75% above",
    rSS_average < SS75 &  rSS_average > SS50 ~ "50% - 75%",
    rSS_average < SS50 &  rSS_average > SS25 ~ "25% - 50%",
    rSS_average < SS25 &  rSS_average > SS1  ~ "1%  - 25%",
    rSS_average < SS1  ~ "0% - 1% ")) %>%
  mutate(transcriptDirect = transcriptAll*(10+(2*12*8))) %>%
  mutate(transcriptOpportunity = transcriptAll*31) %>%
  mutate(transcriptTotal = transcriptDirect + transcriptOpportunity)

# Plot costs of single genes

ggplot(transcriptCosts, aes(x = success, y = log10(transcriptTotal)))+
  geom_boxplot(position = position_dodge(width = 0.5))+
  geom_jitter(alpha = 0.2)+
  ylab("log(Costs) in units of ATP")+
  xlab("average sporulation score (rSS)")+
  ggtitle('transcription costs - single gene deletion')+
  mytheme

# Cumulative costs

transcriptCosts_sum <- transcriptCosts %>%
  group_by(success) %>%
  summarize(transcriptDirect = sum(transcriptDirect, na.rm = T),
            transcriptOpportunity = sum(transcriptOpportunity, na.rm = T),
            transcriptTotal = sum(transcriptTotal, na.rm = T)) %>%
  pivot_longer(cols = c('transcriptDirect', 'transcriptOpportunity', 'transcriptTotal'), values_to = 'costs',
               names_to = 'type')

# Plot cumulative costs

ggplot(transcriptCosts_sum, aes(x = success, y = log10(costs), color = type))+
  geom_point(size = 3, alpha = 0.8)+
  ylab("log(Costs) in units of ATP")+
  xlab("average sporulation score (rSS)")+
  ggtitle('transcription costs - single gene deletion')+
  scale_color_manual("type", values = mycolors)+
  mytheme

# Translation costs

translationCosts <- mergedDelData %>%
  mutate(translationAll = (abundance.filled)*(1774445/1000000)*protein_length.filled) %>% # protein abundance X protein length
  # again the reason I multiply with 3 is that the protein abundance is reported as parts per million, and an average size bacteria has about 3 million protein molecules
  mutate(success = case_when(
    rSS_average > SS75  ~ "75% above",
    rSS_average < SS75 &  rSS_average > SS50 ~ "50% - 75%",
    rSS_average < SS50 &  rSS_average > SS25 ~ "25% - 50%",
    rSS_average < SS25 &  rSS_average > SS1  ~ "1%  - 25%",
    rSS_average < SS1  ~ "0% - 1%")) %>%
  mutate(translationDirect = translationAll*(4+2)) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*24) %>%
  mutate(translationTotal = translationDirect + translationOpportunity)

# Plot costs of single genes

ggplot(translationCosts, aes(x = success, y = log10(translationTotal)))+
  geom_boxplot(position = position_dodge(width = 0.5))+
  geom_jitter(alpha = 0.2)+
  ylab("log(Costs) in units of ATP")+
  xlab("average sporulation score (rSS)")+
  ggtitle('translation costs - single gene deletion')+
  mytheme

# Cumulative costs

translationCosts_sum <- translationCosts %>%
  group_by(success) %>%
  summarize(translationDirect = sum(translationDirect, na.rm = T),
            translationOpportunity = sum(translationOpportunity, na.rm = T),
            translationTotal = sum(translationTotal, na.rm = T)) %>%
  pivot_longer(cols = c('translationDirect', 'translationOpportunity', 'translationTotal'), values_to = 'costs',
               names_to = 'type')

# Plot cumulative costs

ggplot(translationCosts_sum, aes(x = success, y = log10(costs), color = type))+
  geom_point(size = 3, alpha = 0.8)+
  ylab("log(Costs) in units of ATP")+
  xlab("average sporulation score (rSS)")+
  ggtitle('translation costs - single gene deletion')+
  scale_color_manual("", values = mycolors)+
  mytheme

#***************************************************************************************************

# Known mutants

mergedMutData <- mutantLib[-2] %>%
  left_join(mergedAbunData, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 252)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 924)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 256))

# Transcription costs

transcriptCosts_mut <- mergedMutData %>%
  mutate(transcriptAll = (abundance.filled/1000)*(1774445/1000000)*gene_length.filled) %>%
  mutate(success = case_when(
    rSS_average > SS75  ~ "75% above",
    rSS_average < SS75 &  rSS_average > SS50 ~ "50% - 75%",
    rSS_average < SS50 &  rSS_average > SS25 ~ "25% - 50%",
    rSS_average < SS25 &  rSS_average > SS1  ~ "1%  - 25%",
    rSS_average < SS1  ~ "0% - 1% ")) %>%
  mutate(transcriptDirect = transcriptAll*(10+(2*12*8))) %>% # average sporulation time is 8hours, we would expect 8*60 minutes / 5 minutes (degradation rate) = 96 re-polymerization events, assuming nucleotides are well recycled and it only affects polymerization costs
  mutate(transcriptOpportunity = transcriptAll*31) %>%
  mutate(transcriptTotal = transcriptDirect + transcriptOpportunity)

# Plot costs of single genes

# Sum of all
sum_t <- sum(transcriptCosts_mut$transcriptTotal)

ggplot(transcriptCosts_mut, aes(x = success, y = log10(transcriptTotal)))+
  geom_boxplot(position = position_dodge(width = 0.5))+
  geom_jitter(alpha = 0.2)+
  ylab("log(Costs) in units of ATP")+
  xlab("average sporulation score (rSS)")+
  ggtitle('transcription costs - known mutants')+
  annotate("text", x=1, y=6, label = paste('total costs =', format(sum_t, digits = 2)))+
  mytheme

# Cumulative costs

transcriptCosts_sum_mut <- transcriptCosts_mut %>%
  group_by(success) %>%
  summarize(transcriptDirect = sum(transcriptDirect, na.rm = T),
            transcriptOpportunity = sum(transcriptOpportunity, na.rm = T),
            transcriptTotal = sum(transcriptTotal, na.rm = T)) %>%
  pivot_longer(cols = c('transcriptDirect', 'transcriptOpportunity', 'transcriptTotal'), values_to = 'costs',
               names_to = 'type')

# Plot cumulative costs

ggplot(transcriptCosts_sum_mut, aes(x = success, y = log10(costs), color = type))+
  geom_point(size = 3, alpha = 0.8)+
  ylab("log(Costs) in units of ATP")+
  xlab("average sporulation score (rSS)")+
  ggtitle('transcription costs - known mutants')+
  scale_color_manual("type", values = mycolors)+
  mytheme

# Translation costs

translationCosts_mut <- mergedMutData %>%
  mutate(translationAll = abundance.filled*(1774445/1000000)*protein_length.filled) %>%
  mutate(success = case_when(
    rSS_average > SS75  ~ "75% above",
    rSS_average < SS75 &  rSS_average > SS50 ~ "50% - 75%",
    rSS_average < SS50 &  rSS_average > SS25 ~ "25% - 50%",
    rSS_average < SS25 &  rSS_average > SS1  ~ "1%  - 25%",
    rSS_average < SS1  ~ "0% - 1%")) %>%
  mutate(translationDirect = translationAll*(4+2)) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*24) %>%
  mutate(translationTotal = translationDirect + translationOpportunity)

# Plot costs of single genes

# Sum of all

sum_p <- sum(translationCosts_mut$translationTotal)

ggplot(translationCosts_mut, aes(x = success, y = log10(translationTotal)))+
  geom_boxplot(position = position_dodge(width = 0.5))+
  geom_jitter(alpha = 0.2)+
  ylab("log(Costs) in units of ATP")+
  xlab("average sporulation score (rSS)")+
  ggtitle('translation costs - known mutants')+
  scale_color_manual("type", values = mycolors)+
  annotate("text", x=1, y=8, label = paste('total costs =', format(sum_p, digits = 2)))+
  mytheme

# Cumulative costs

translationCosts_sum_mut <- translationCosts_mut %>%
  group_by(success) %>%
  summarize(translationDirect = sum(translationDirect, na.rm = T),
            translationOpportunity = sum(translationOpportunity, na.rm = T),
            translationTotal = sum(translationTotal, na.rm = T)) %>%
  pivot_longer(cols = c('translationDirect', 'translationOpportunity', 'translationTotal'), values_to = 'costs',
               names_to = 'type')

# Plot cumulative costs

ggplot(translationCosts_sum_mut, aes(x = success, y = log10(costs), color = type))+
  geom_point(size = 3, alpha = 0.8)+
  ylab("log(Costs) in units of ATP")+
  xlab("average sporulation score (rSS)")+
  ggtitle('translation costs - known mutants')+
  scale_color_manual("type", values = mycolors)+
  mytheme


# New plot
transcriptCosts_mut
translationCosts_mut

plotData <- cbind.data.frame(transcriptDirect       = transcriptCosts_mut$transcriptDirect,
                             transcriptOpportunity  = transcriptCosts_mut$transcriptOpportunity,
                             transcriptTotal        = transcriptCosts_mut$transcriptTotal,
                             translationDirect      = translationCosts_mut$translationDirect,
                             translationOpportunity = translationCosts_mut$translationOpportunity,
                             translationTotal       = translationCosts_mut$translationTotal)

plotDataLong <- pivot_longer(plotData, cols = everything() , values_to = 'costs', names_to = 'type')
plotDataLong$source <- rep(c("transcription","translation"), each = 3, times = 174)
plotDataLong$type   <- rep(c("direct","opportunity", "total"), times = 348)

plotDataLong_sum <- plotDataLong %>%
  group_by(type,source) %>%
  summarize(sum(costs, na.rm = T))

plotDataLong_dir_oppt <- subset(plotDataLong, type!="total")

ggplot(plotDataLong_dir_oppt, aes(x = source, y = log10(costs), color = type))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(), alpha = .2)+
  ylab("log(Costs) in units of ATP")+
  scale_color_manual("type", values = mycolors4)+
  mytheme


# Breakdown central dogma costs

plotDataLong_central_sum <- plotDataLong_dir_oppt %>%
  group_by(source) %>%
  summarize(sum(costs, na.rm = T))

cent_costs_mut  <- c(14983974, 638793533)
stages_mut      <- c("transcription", "translation")

donut_plot2_mut <- cbind.data.frame(costs=cent_costs_mut, stages=stages_mut)

# Compute percentages
donut_plot2_mut$fraction = (donut_plot2_mut$costs / sum(donut_plot2_mut$costs))*100

# Compute the cumulative percentages (top of each rectangle)
donut_plot2_mut$ymax = cumsum(donut_plot2_mut$fraction)

# Compute the bottom of each rectangle
donut_plot2_mut$ymin = c(0, head(donut_plot2_mut$ymax, n=-1))

# Compute label position
donut_plot2_mut$labelPosition  <- (donut_plot2_mut$ymax + donut_plot2_mut$ymin) / 2

# Compute a good label
donut_plot2_mut$label  <- paste0(donut_plot2_mut$stages, " : ", round(donut_plot2_mut$fraction,1))

ggplot(donut_plot2_mut, aes(ymax = ymax , ymin = ymin, xmax=4, xmin=3, fill = stages))+
  geom_rect()+
  coord_polar(theta="y")+
  xlim(c(2, 4))+
  geom_label( x = 3.5, aes(y = labelPosition, label = label), size = 4) +
  #ylab("log(Costs) in units of ATP")+
  #xlab("stages")+
  ggtitle('Total costs of central dogma')+
  mytheme+
  scale_fill_manual(values = mycolors3.1)+
  theme_void()+
  theme(legend.position = "none")

# Breakdown of cost types

plotDataLong_type_sum <- plotDataLong_dir_oppt %>%
  group_by(type) %>%
  summarize(sum(costs, na.rm = T))

type_type_mut   <- c(140749105, 513028402)
stages_type_mut <- c("direct", "opportunity")

donut_plot3_mut <- cbind.data.frame(costs=type_type_mut, stages=stages_type_mut)

# Compute percentages
donut_plot3_mut$fraction = (donut_plot3_mut$costs / sum(donut_plot3_mut$costs))*100

# Compute the cumulative percentages (top of each rectangle)
donut_plot3_mut$ymax = cumsum(donut_plot3_mut$fraction)

# Compute the bottom of each rectangle
donut_plot3_mut$ymin = c(0, head(donut_plot3_mut$ymax, n=-1))

# Compute label position
donut_plot3_mut$labelPosition  <- (donut_plot3_mut$ymax + donut_plot3_mut$ymin) / 2

# Compute a good label
donut_plot3_mut$label  <- paste0(donut_plot3_mut$stages, " : ", round(donut_plot3_mut$fraction,1))

ggplot(donut_plot3_mut, aes(ymax = ymax , ymin = ymin, xmax=4, xmin=3, fill = stages))+
  geom_rect()+
  coord_polar(theta="y")+
  xlim(c(2, 4))+
  geom_label( x = 3.5, aes(y = labelPosition, label = label), size = 4) +
  ggtitle('Direct and opportunity costs')+
  mytheme+
  scale_fill_manual(values = mycolors4.1)+
  theme_void()+
  theme(legend.position = "none")

#********************************************************************************************

# Costs of conserved sporulation genes
# After Acidamicoccus fermentans DSM 20731 non-sporulating species, add an id column

col = rep(c('sporulating', 'non-sporulating'), times = c(21, 19))

mergedOrtData <- conservedLib %>%
  pivot_longer(cols = names(conservedLib[,4:43]), values_to = 'presence',
               names_to = 'species') %>%
  mutate(sporulation = rep(col, times = 149)) %>%
  left_join(mergedAbunData, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 252)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 924)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 256))

# Transcription costs

transcriptCosts_ort <- mergedOrtData %>%
  mutate(transcriptAll = (abundance.filled/1000)*(1774445/1000000)*gene_length.filled) %>%
  mutate(transcriptDirect = transcriptAll*(10+(2*12*8))) %>%
  mutate(transcriptOpportunity = transcriptAll*31) %>%
  mutate(transcriptTotal = transcriptDirect + transcriptOpportunity)

# Plot costs of single genes

ggplot(transcriptCosts_ort, aes(x = sporulation, y = log10(transcriptTotal), color = presence))+
  geom_jitter(alpha = 0.2)+
  geom_boxplot(position = position_dodge(width = 0.5))+
  ylab("log(Costs) in units of ATP")+
  xlab("species")+
  ggtitle('transcription costs - across species')+
  scale_color_manual(name = 'gene presence', values = mycolors)+
  mytheme


# Translation costs

translationCosts_ort <- mergedOrtData %>%
  mutate(translationAll = abundance.filled*(1774445/1000000)*protein_length.filled) %>%
  mutate(translationDirect = translationAll*(4+2)) %>%
  mutate(translationOpportunity = translationAll*24) %>%
  mutate(translationTotal = translationDirect + translationOpportunity)

# Plot costs of single genes

ggplot(translationCosts_ort, aes(x = sporulation, y = log10(translationTotal), color = presence))+
  geom_jitter(alpha = 0.2)+
  geom_boxplot(position = position_dodge(width = 0.5))+
  ylab("log(Costs) in units of ATP")+
  xlab("species")+
  ggtitle('translation costs - across species')+
  scale_color_manual(name = 'gene presence', values = mycolors)+
  mytheme



