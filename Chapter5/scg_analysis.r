  ##load libraries####################################################################################################
  library(plyr)
  library(data.table)
  library(dplyr)
  library(vegan)
  library(reshape2)
  library(ggplot2)
  library(ggvegan)
  library(metacoder)
##setup functions###################################################################################################
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



##multimerge function###
multbind = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){((read.csv(file=x,header=T, sep = ",")))})
  Reduce(function(x,y) {full_join(x,y, by= "Gene")}, datalist)}


###replace NA with zeroes###
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}


##merge files######################################################################################################

###point file path to folder containing final annotations files######
all_data_RNA =  multbind("/home/jonny/Documents/data/2015_data/bg/mapping/tot_reads/rna/")
all_data_DNA =  multbind("/home/jonny/Documents/data/2015_data/bg/mapping/tot_reads/dna/")

all_data_no_zero <- na.zero(all_data)
##save file###

write.csv(all_data_RNA, file="all_bg_rna.csv", quote=FALSE, row.names = FALSE)
write.csv(all_data_DNA, file="all_bg_dna.csv", quote=FALSE, row.names = FALSE)

rna_tots <- c(11627682,23106560,15635254,19673146,22284836,20884710,20389216,16063534,24953610,31093464,26014696,27473168,26445188,25569736,21347046,20568324,27537876,24282974,26893654,24495798,25263310)
dna_tots <- c(178728628,287225804,233779812,198068956,197243492,24816830,26340810,18915710,28260112,24928364,27488284,23905258,21247590,22449786,27860216,284810538,259674992,221777750,238261726,206798312,23912724,25756704,21757908,25824124,22854534,29122676,30062108,32796018,22145254,25676230,195609520,183941758,202565964,159621730,202690830,25658524,26820300,23849806,27962464,25071848,28521492,23497974,24222736,20536298,20453460)


#normalise RNA

a <- sweep(as.matrix(all_data_RNA[-1]), 2, rna_tots, '/')
b <- a*1000000
d <- as.data.frame(b)
rownames(d) <- all_data_RNA$Gene
d
write.csv(d, file="all_bg_rna_norm.csv", quote = FALSE)


#normalsie DNA

a <- sweep(as.matrix(all_data_DNA[-1]), 2, dna_tots, '/')
b <- a*1000000
d <- as.data.frame(b)
rownames(d) <- all_data_DNA$Gene
d
write.csv(d, file="all_bg_dna_norm.csv", quote = FALSE)



##############################################################################################
# scg counts
#############################################################################################

setwd("/home/jonny/Documents/data/2015_data/scg/")
aln <- read.csv("alneberg_counts.txt", header=FALSE, sep = " ")
aln$set <- "Alneberg_et_al"
cam <- read.csv("campbell_counts.txt", header=FALSE, sep = " ")
cam$set <- "Campbell_et_al"
dup <- read.csv("dupont_counts.txt", header=FALSE, sep = " ")
dup$set <- "Dupont_et_al"
rin <- read.csv("rinke_counts.txt", header=FALSE, sep = " ")
rin$set <- "rinke_et_al"

all_scg <- rbind(cam,aln,dup,rin)

tail(as.data.frame(unique(cam$V2)))
tail(as.data.frame(unique(aln$V2)))
tail(as.data.frame(unique(dup$V2)))
tail(as.data.frame(unique(rin$V2)))
library(ggplot2)

ggplot(all_scg, aes(x= set, y= V1, colour = set)) +
  scale_color_brewer(palette = "Set1") +
  geom_violin() + geom_jitter(alpha=0.6)


ggplot(all_scg, aes(x=V1)) + 
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_density(aes(group=set, colour=set, fill=set), alpha=0.3) +
  facet_grid( ~ set, scales = "free_y")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

mean(cam$V1)
median(cam$V1)
getmode(cam$V1)

mean(aln$V1)
median(aln$V1)
getmode(aln$V1)

mean(dup$V1)
median(dup$V1)
getmode(dup$V1)

mean(rin$V1)
median(rin$V1)
getmode(rin$V1)




#all_scg_DNA =  multbind("/home/jonny/Documents/data/2015_data/scg/coverage/dna/")
#all_scg_RNA =  multbind("/home/jonny/Documents/data/2015_data/scg/coverage/rna/")
#write.csv(all_scg_DNA, file="all_scg_dna", quote = FALSE)
#write.csv(all_scg_RNA, file="all_scg_rna", quote = FALSE)

all_scg_DNA <- read.csv("all_scg_dna", header = TRUE, row.names = 1)
all_scg_RNA <- read.csv("all_scg_rna", header = TRUE, row.names = 1)

rna_tots <- c(11627682,23106560,15635254,19673146,22284836,20884710,20389216,16063534,24953610,31093464,26014696,27473168,26445188,25569736,21347046,20568324,27537876,24282974,26893654,24495798,25263310)
dna_tots <- c(178728628,287225804,233779812,198068956,197243492,24816830,26340810,18915710,28260112,24928364,27488284,23905258,21247590,22449786,27860216,284810538,259674992,221777750,238261726,206798312,23912724,25756704,21757908,25824124,22854534,29122676,30062108,32796018,22145254,25676230,195609520,183941758,202565964,159621730,202690830,25658524,26820300,23849806,27962464,25071848,28521492,23497974,24222736,20536298,20453460)


#normalise RNA

norm_rna = function(x){
a <- sweep(as.matrix(x[-1]), 2, rna_tots, '/')
b <- a*1000000
d <- as.data.frame(b)
rownames(d) <- x$Gene
d
}

norm_dna = function(x){
  a <- sweep(as.matrix(x[-1]), 2, dna_tots, '/')
  b <- a*1000000
  d <- as.data.frame(b)
  rownames(d) <- x$Gene
  d
}

scg_rna_norm <- norm_rna(all_scg_RNA)
scg_dna_norm <- norm_dna(all_scg_DNA)

#write.csv(scg_rna_norm, file="all_scg_rna_norm.csv", quote = FALSE)
#write.csv(scg_dna_norm, file="all_scg_dna_norm.csv", quote = FALSE)

############################################################
#ordinations
############################################################

library(plyr)
library(data.table)
library(dplyr)
library(vegan)
library(reshape2)
library(ggplot2)
library(ggvegan)

setwd("/home/jonny/Documents/data/2015_data/scg/")
scg_rna_norm <- read.csv("all_scg_rna_norm_norep.csv", header = TRUE)
scg_dna_norm <- read.csv("all_scg_dna_norm.csv", header = TRUE)
scg_annotations <- read.csv("scg_kaiju.csv", header = TRUE, sep = "\t")
rna_factors <-read.table("RNAfactors.csv", sep=',', header=TRUE)
dna_factors <-read.table("factors.csv", sep=',', header=TRUE)

taxbactprop$mergetax<-paste(taxbact$Kingdom,"|",taxbact$Phylum,"|",taxbact$Class,"|",taxbact$Order,"|",taxbact$Family,"|",taxbact$Genus,"|",taxbact$Species)


scg_annotations$merge <- paste(scg_annotations$phylum,"|",scg_annotations$class,"|",scg_annotations$order,"|",scg_annotations$family,"|",scg_annotations$genus,"|",scg_annotations$species,"|")
scg_annotations2 <- scg_annotations[c(1,9)]
scg_dna_full <- full_join(scg_dna_norm, scg_annotations2, by = "Gene")
scg_rna_full <- full_join(scg_rna_norm, scg_annotations2, by = "Gene")
head(scg_dna_full)
head(scg_rna_full)
scg_rna_sum <- aggregate(scg_dna_full[2:46], by=list(scg_dna_full$merge), FUN=sum)
scg_rna_sum <- aggregate(scg_rna_full[2:16], by=list(scg_dna_full$merge), FUN=sum)
write.csv(scg_dna_sum,"test.csv",quote = FALSE)
write.csv(scg_rna_sum,"testrna.csv",quote = FALSE)


### DNA ####

bacttaxmap<-parse_taxonomy_table("test4.csv", taxon_col = c("class" = 1), 
                                 other_col_type="obs_info", sep=",",
                                 class_sep = "\\|")

data <- mutate_taxa(bacttaxmap, 
            mean_percentage = sapply(obs(bacttaxmap), function(i) sum(bacttaxmap$obs_data$av_abund[i])))

data_filt <- filter_taxa(data,mean_percentage > 0.05, supertaxa = TRUE)

heat_tree(data_filt, node_size=mean_percentage,
          edge_color = mean_percentage,
          node_color = mean_percentage,
          node_color_range = c("cyan", "magenta", "green"),
          edge_color_range = c("cyan", "magenta", "green"),
          initial_layout = "reingold-tilford", layout = "davidson-harel",
          node_label = ifelse(mean_percentage > 0.001, name, NA))


### RNA ###

bacttaxmap<-parse_taxonomy_table("test4rna.csv", taxon_col = c("class" = 1), 
                                 other_col_type="obs_info", sep=",",
                                 class_sep = "\\|")

data <- mutate_taxa(bacttaxmap, 
                    mean_percentage = sapply(obs(bacttaxmap), function(i) sum(bacttaxmap$obs_data$av_abund[i])))

data_filt <- filter_taxa(data,mean_percentage > 0.05, supertaxa = TRUE)

heat_tree(data_filt, node_size=mean_percentage,
          edge_color = mean_percentage,
          node_color = mean_percentage,
          node_color_range = c("cyan", "magenta", "green"),
          edge_color_range = c("cyan", "magenta", "green"),
          initial_layout = "reingold-tilford", layout = "davidson-harel",
          node_label = ifelse(mean_percentage > 0.001, name, NA))


heat_tree(filter_taxa(data_filt, name == "Archaea", subtaxa = TRUE),
          node_size=n_obs,
          edge_color = mean_percentage,
          node_color = mean_percentage,
          node_color_range = c("cyan", "magenta", "green"),
          edge_color_range   = c("#555555", "#EEEEEE"),
          initial_layout = "reingold-tilford", layout = "davidson-harel",
          node_label = ifelse(mean_percentage > 50, name, NA))

phy <- scg_annotations[c(1,2)]
class <- scg_annotations[c(1,3)]
gen <- scg_annotations[c(1,6)]
spec <- scg_annotations[c(1,7)]

head(scg_annotations)
scg_dna_perc <- decostand(scg_dna_norm[2:46], method = "total", MARGIN = 2)
scg_dna_perc <- scg_dna_perc * 100
scg_dna_perc <- cbind(scg_dna_norm[1], scg_dna_perc)
scg_rna_perc <- decostand(scg_rna_norm[2:16], method = "total", MARGIN = 2)
scg_rna_perc <- scg_rna_perc * 100
scg_rna_perc <- cbind(scg_rna_norm[1], scg_rna_perc)

phy_dna <- full_join(scg_dna_perc, phy, by = "Gene")
phy_rna <- full_join(scg_rna_perc, phy, by = "Gene")
class_dna <- full_join(scg_dna_perc, class, by = "Gene")
class_rna <- full_join(scg_rna_perc, class, by = "Gene")
gen_dna <- full_join(scg_dna_perc, gen, by = "Gene")
gen_rna <- full_join(scg_rna_perc, gen, by = "Gene")
spec_dna <- full_join(scg_dna_perc, spec, by = "Gene")
spec_rna <- full_join(scg_rna_perc, spec, by = "Gene")

phy_dna_f <- aggregate(phy_dna[2:46], by=list(phy_dna$phylum), FUN=sum, na.rm=TRUE)
phy_rna_f <- aggregate(phy_rna[2:16], by=list(phy_rna$phylum), FUN=sum, na.rm=TRUE)
class_dna_f <- aggregate(class_dna[2:46], by=list(class_dna$class), FUN=sum, na.rm=TRUE)
class_rna_f <- aggregate(class_rna[2:16], by=list(class_rna$class), FUN=sum, na.rm=TRUE)
gen_dna_f <- aggregate(gen_dna[2:46], by=list(gen_dna$genus), FUN=sum, na.rm=TRUE)
gen_rna_f <- aggregate(gen_rna[2:16], by=list(gen_rna$genus), FUN=sum, na.rm=TRUE)
spec_dna_f <- aggregate(spec_dna[2:46], by=list(spec_dna$species), FUN=sum, na.rm=TRUE)
spec_rna_f <- aggregate(spec_rna[2:16], by=list(spec_rna$species), FUN=sum, na.rm=TRUE)

phy_dna_f1 <- phy_dna_f[2:46]
rownames(phy_dna_f1) <- phy_dna_f$Group.1

phy_rna_f1 <- phy_rna_f[2:16]
rownames(phy_rna_f1) <- phy_rna_f$Group.1
phy_rna_f2 = phy_rna_f1[ rowSums(phy_rna_f1)!=0, ] 

class_dna_f1 <- class_dna_f[2:46]
rownames(class_dna_f1) <- class_dna_f$Group.1

class_rna_f1 <- class_rna_f[2:16]
rownames(class_rna_f1) <- class_rna_f$Group.1
class_rna_f2 = class_rna_f1[ rowSums(class_rna_f1)!=0, ] 

gen_dna_f1 <- gen_dna_f[2:46]
rownames(gen_dna_f1) <- gen_dna_f$Group.1
gen_dna_f2 = gen_dna_f1[ rowSums(gen_dna_f1)!=0, ]

gen_rna_f1 <- gen_rna_f[2:16]
rownames(gen_rna_f1) <- gen_rna_f$Group.1
gen_rna_f2 = gen_rna_f1[ rowSums(gen_rna_f1)!=0, ] 

spec_dna_f1 <- spec_dna_f[2:46]
rownames(spec_dna_f1) <- spec_dna_f$Group.1
spec_dna_f2 = spec_dna_f1[ rowSums(spec_dna_f1)!=0, ]

spec_rna_f1 <- spec_rna_f[2:16]
rownames(spec_rna_f1) <- spec_rna_f$Group.1
spec_rna_f2 = spec_rna_f1[ rowSums(spec_rna_f1)!=0, ] 

eco_dist <- vegdist((as.data.frame(t(class_dna_f1))), dist = "euc")
class_dna_dis <- capscale(eco_dist~1)
eco_dist2 <- vegdist((as.data.frame(t(class_rna_f2))), dist = "euc")
class_rna_dis <- capscale(eco_dist2~1)
eco_dist3 <- vegdist((as.data.frame(t(gen_dna_f2))), dist = "euc")
gen_dna_dis <- capscale(eco_dist3~1)
eco_dist4 <- vegdist((as.data.frame(t(gen_rna_f2))), dist = "euc")
gen_rna_dis <- capscale(eco_dist4~1)

f1 <- fortify(class_dna_dis)
f2 <- fortify(gen_dna_dis)
f3 <- fortify(class_rna_dis)
f4 <- fortify(gen_rna_dis)




##subset and add taxrank descriptor##

fps_1 <- subset(f1, Score =="sites")
fps_2 <- subset(f2, Score =="sites")
fps_3 <- subset(f3, Score =="sites")
fps_4 <- subset(f4, Score =="sites")


fps_1$Rank <- "DNA_class" 
fps_1$Dim2 <- (fps_1$Dim2*-1)
fps_2$Rank <- "DNA_gen"
fps_3$Rank <- "RNA_class"
fps_4$Rank <- "RNA_gen"


allm <- rbind(fps_1, fps_2, fps_3, fps_4)
all_factors <- rbind(dna_factors, rna_factors)
pca_final <- full_join(allm, all_factors, by = "Label")
pca_final <- subset(pca_final, Rank != "NA")
#pca_final$Rank_f = factor(pca_final$Rank, levels=c('Phylum','Class','Order','Family', 'Genus', 'Species', 'SEED level 1', 'SEED level 2', 'SEED level 3'))



# plot by horizon "using hulls ################################################################
ggplot(pca_final, aes(x = Dim1, y = Dim2, colour = Depth, label = Site, shape = Site)) +
  geom_point(size = 2, stroke = 1.2) +
  #geom_text(hjust = 1.1, vjust = -1) +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(name = "Site", values = 1:nlevels(pca_final$Site)) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 11),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5)) +
  facet_wrap(~ Rank, scales = "free") +
  xlab("PCoA 1") + ylab("PCoA 2")


extract_importance = function(x){
  a <- (summary(eigenvals(x)))
  b <- t(a$importance)
  b[1:2,]
}

extract_importance(class_dna_dis)
extract_importance(gen_dna_dis)
extract_importance(class_rna_dis)
extract_importance(gen_rna_dis)



################################################
# boxplots
################################################

class_dna_f$avg <- rowMeans(class_dna_f[2:46])
sub_class_dna <- head(arrange(class_dna_f,desc(avg)), n = 20)
m_class <- melt(sub_class_dna[1:46])
colnames(m_class) <- c("Taxa", "Label", "value")
final_dna_m <- full_join(m_class, dna_factors, by= "Label")


class_rna_f$avg <- rowMeans(class_rna_f[2:16])
sub_class_rna <- head(arrange(class_rna_f,desc(avg)), n = 20)
m_class_rna <- melt(sub_class_rna[1:16])
colnames(m_class_rna) <- c("Taxa", "Label", "value")
final_rna_m <- full_join(m_class_rna, rna_factors, by= "Label", na.rm = T)
final_rna_m <- subset(final_rna_m, Taxa != "NA")

# RNA # 
ggplot(final_rna_m, aes(x=reorder(Taxa, -value, FUN=median), y= value)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color = final_rna_m$Depth, alpha = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size=12),
        strip.text = element_text(size = 11),
        panel.border = element_rect(fill = NA, size=0.5))
        

# DNA # 
ggplot(final_dna_m, aes(x=reorder(Taxa, -value, FUN=median), y= value)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color = final_dna_m$Depth, alpha = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size=12),
        panel.border = element_rect(fill = NA, size=0.5))
        

reorder(Taxa, value, FUN=median)

###############################################
# alpha div
##############################################

#shannon
adiv_dna_g <- diversity(as.data.frame(t(gen_dna_f2)), index = "shannon")
dna_df <- (as.data.frame(adiv_dna_g))
setDT(dna_df, keep.rownames = TRUE)[]
names(dna_df)[2]<-"adiv"
names(dna_df)[1]<-"Label"

adiv_rna_g <- diversity(as.data.frame(t(gen_rna_f2)), index = "shannon")
rna_df <- (as.data.frame(adiv_rna_g))
setDT(rna_df, keep.rownames = TRUE)[]
names(rna_df)[2]<-"adiv"
names(rna_df)[1]<-"Label"

adiv_all_sh <- rbind(dna_df,rna_df)
adiv_plots <- full_join(adiv_all_sh, all_factors, by = "Label")
plot_final <- subset(adiv_plots, adiv != "NA")

adiv_plot_f_m <- melt(plot_final)

adiv_plot_se <- summarySE(adiv_plot_f_m, measurevar="value", groupvars=c("type","Site","Depth"))
adiv_dna_plot <- subset(adiv_plot_se, type == "DNA")
adiv_rna_plot <- subset(adiv_plot_se, type == "RNA")

ggplot(adiv_dna_plot, aes(x=Depth, y = value, colour = Site)) +
  geom_point(size=3, shape=21, fill="white") +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black") +
  geom_line(aes(group = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap( ~ Site) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5))

ggplot(adiv_rna_plot, aes(x=Depth, y = value, colour = Site)) +
  geom_point(size=3, shape=21, fill="white") +
  #ylim(5, 6.25) +
  geom_line(aes(group = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap( ~ Site) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5))


##############################################
# beta div
##############################################

bdis_prep = function(x){ 
  b <- vegdist(x, method = "euclidean")
  with(dna_factors, betadisper(b, Depth, type = "centroid"))
}

rna_bdis_fact <- subset(rna_factors, Rep =="A")

bdis_rna_prep = function(x){ 
  b <- vegdist(x, method = "euclidean")
  with(rna_bdis_fact, betadisper(b, Depth, type = "centroid"))
}

bdis_boxplot <- function(x) {
  a <- cbind((as.data.frame(x$distances)), (as.data.frame((x$group))))
  colnames(a) <- c("Dist", "group")
  a
}

bdis_dna <- bdis_prep(as.data.frame(t(gen_dna_f2)))
bdis_rna <- bdis_rna_prep(as.data.frame(t(gen_rna_f2)))


anova(bdis_dna)
anova(bdis_rna)



bdis_plot1 <- bdis_boxplot(bdis_dna)
bdis_plot1$rank = "DNA"
bdis_plot2 <- bdis_boxplot(bdis_rna)
bdis_plot2$rank = "RNA"


bdis_all <- rbind(bdis_plot1, bdis_plot2)
setDT(bdis_all, keep.rownames = TRUE)[]
names(bdis_all)[2]<-"bdiv"
names(bdis_all)[1]<-"Label"
all_factors <- rbind(dna_factors, rna_factors)
bdis_plots <- full_join(bdis_all, all_factors, by = "Label")
plot_final <- subset(bdis_plots, bdiv != "NA")
bdiv_plot_se <- summarySE(bdis_plots, measurevar="bdiv", groupvars=c("type","Site","Depth"))
bdiv_dna_plot <- subset(bdiv_plot_se, type == "DNA")
bdiv_rna_plot <- subset(plot_final, type == "RNA")

ggplot(bdiv_dna_plot, aes(x=Depth, y = bdiv, colour = Site)) +
  geom_point(size=3, shape=21, fill="white") +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=bdiv-se, ymax=bdiv+se), colour = "black") +
  geom_line(aes(group = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap( ~ Site) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5))

ggplot(bdiv_rna_plot, aes(x=Depth, y = bdiv, colour = Site)) +
  geom_point(size=3, shape=21, fill="white") +
  #ylim(5, 6.25) +
  geom_line(aes(group = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap( ~ Site) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5))


ggplot(bdis_all, aes(x = group, y = bdiv, colour = group)) +
  geom_boxplot() +
  xlab("") +
  ylab("Distance to centroid") +
  scale_color_brewer(name = "", palette = "Set1") +
  stat_compare_means(method = "anova",label= "p.format", label.y.npc = 0.9) +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none",
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5)) +
  facet_wrap( ~ rank, scales = 'free')



###############################################
# anosim / adonis
#############################################

make_dist = function(x){
  helling_Taxa<-decostand(x,method="hellinger")
  vegdist(helling_Taxa,dist="euclidean")
}

do_anosim = function(x){
  dist <- make_dist(x)
  anosim(dist, dna_factors$Depth, permutations = 9999)
}

anosim_plots <- function(x){
  a <-cbind((as.data.frame(x$dis.rank)),(as.data.frame(x$class.vec)))
  colnames(a) <- c("rank", "cat")
  a
}


dna_dist <- make_dist(as.data.frame(t(gen_dna_f2)))
a <- do_anosim(dna_dist)
b <- anosim_plots(a)
anosim_plot <- ggplot(b, aes(x = cat, y = rank, colour = cat)) +
  geom_boxplot()

rna_dist <- make_dist(as.data.frame(t(gen_rna_f2)))



phy_ano <- do_anosim(a_phy_f)
class_ano <- do_anosim(a_class_f)
ord_ano <- do_anosim(a_ord_f)
fam_ano <- do_anosim(a_fam_f)
gen_ano <- do_anosim(a_gen_f)
spec_ano <- do_anosim(a_spec_f)
s1_ano <- do_anosim(a_seed1_f)
s2_ano <- do_anosim(a_seed2_f)
s3_ano <- do_anosim(a_seed3_f)

summary(phy_ano)
summary(class_ano)
summary(ord_ano)
summary(fam_ano)
summary(gen_ano)
summary(spec_ano)
summary(s1_ano)
summary(s2_ano)
summary(s3_ano)


class_ano_p <- anosim_plots(class_ano)
gen_ano_p <- anosim_plots(gen_ano)
seed1_ano_p <- anosim_plots(s1_ano)
seed3_ano_p <- anosim_plots(s3_ano)

class_ano_p$class <- "Class"
gen_ano_p$class <- "Genus"
seed1_ano_p$class <- "SEED level 1"
seed3_ano_p$class <- "SEED level 3"

ano_merged <- rbind(class_ano_p, gen_ano_p, seed1_ano_p, seed3_ano_p)


anosim_ploty <- ggplot(ano_merged, aes(x = cat, y = rank, colour = cat)) +
  geom_boxplot() +
  xlab("") +
  ylab("Ranked dissimilarity") +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none",
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 11),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5)) +
  facet_wrap( ~ class)

#######################################################
# adonis
#######################################################

adon_hor = function(x){
  dist <- make_dist(x)
  adonis2(dist ~ Horizon , data=alps_fact_anosim, permutations = 99999)
}

adon_site = function(x){
  dist <- make_dist(x)
  adonis2(dist ~ Site , data=alps_fact_anosim, permutations = 99999)
}


# depth no strata # 
adonis2(dna_dist ~ Depth , data=dna_factors, permutations = 9999)
# site no strata # 
adonis2(dna_dist ~ Site , data=dna_factors, permutations = 9999)
# depth stratified by site #
adonis2(dna_dist ~ Depth , data=dna_factors, strata= Site, permutations = 9999)
#site stratefied by depth
adonis2(dna_dist ~ Site , data= dna_factors, strata= Depth, permutations = 9999)
# depth * rep stratified by site # 
adonis2(dna_dist ~ Depth * Rep , data=dna_factors, strata= Site, permutations = 9999)
# site * rep stratified by depth #
adonis2(dna_dist ~ Site * Rep , data=dna_factors, strata= Depth, permutations = 9999)
# site * depth stratefied by rep # 
adonis2(dna_dist ~ Site * Depth , data=dna_factors, strata= Rep, permutations = 9999)


rna_factors_sub <- subset(rna_factors, Rep =="A")

# depth no strata # 
adonis2(rna_dist ~ Depth , data=rna_factors_sub, permutations = 9999)
# site no strata # 
adonis2(rna_dist ~ Site , data=rna_factors_sub, permutations = 9999)
# depth stratified by site #
adonis2(rna_dist ~ Depth , data=rna_factors_sub, strata= Site, permutations = 9999)
#site stratefied by depth
adonis2(rna_dist ~ Site , data= rna_factors_sub, strata= Depth, permutations = 9999)
# depth * rep stratified by site # 
adonis2(rna_dist ~ Depth * Rep , data=rna_factors_sub, strata= Site, permutations = 9999)
# site * rep stratified by depth #
adonis2(rna_dist ~ Site * Rep , data=rna_factors_sub, strata= Depth, permutations = 9999)
# site * depth stratefied by rep # 
adonis2(rna_dist ~ Site * Depth , data=rna_factors_sub, strata= Rep, permutations = 9999)



ad_phy_hor <- adon_hor(a_phy_f)
ad_phy_site <- adon_site(a_phy_f)
ad_class_hor <- adon_hor(a_class_f)
ad_class_site <- adon_site(a_class_f)
ad_ord_hor <- adon_hor(a_ord_f)
ad_ord_site <- adon_site(a_ord_f)
ad_fam_hor <- adon_hor(a_fam_f)
ad_fam_site <- adon_site(a_fam_f)
ad_gen_hor <- adon_hor(a_gen_f)
ad_gen_site <- adon_site(a_gen_f)
ad_spec_hor <- adon_hor(a_spec_f)
ad_spec_site <- adon_site(a_spec_f)
ad_seed1_hor <- adon_hor(a_seed1_f)
ad_seed1_site <- adon_site(a_seed1_f)
ad_seed2_hor <- adon_hor(a_seed2_f)
ad_seed2_site <- adon_site(a_seed2_f)
ad_seed3_hor <- adon_hor(a_seed3_f)
ad_seed3_site <- adon_site(a_seed3_f)



##############################################
# simper
##############################################

do_simp = function(x){
  simp <- with(dna_factors, simper(x, Depth, permutations = 9999))
  summary(simp, ordered = TRUE)
}

a <- with(dna_factors, simper(gen_dna_f2, Depth, permutations = 999))

simp_class <- do_simp(as.data.frame(t(class_dna_f1)))
ab <- simp_class$Ah_Bw
ac <- simp_class$Ah_Cox
bc <- simp_class$Bw_Cox
write.csv(ab, file = "phy_a_b_simper", quote = FALSE)
write.csv(ac, file = "phy_a_c_simper", quote = FALSE)
write.csv(bc, file = "phy_b_c_simper", quote = FALSE)


##############################################
# heatmaps
###############################################

perc_filter = function(x){
  maxab <- apply(x, 2, max)
  n1 <- names(which(maxab < 100))        ##### control % filter
  x[, -which(names(x) %in% n1)]
}

t_dna_c <- (as.data.frame(t(class_dna_f1)))
t_dna_g <- (as.data.frame(t(gen_dna_f2)))
t_rna_c <- (as.data.frame(t(class_rna_f2)))
t_rna_g <- (as.data.frame(t(gen_rna_f2)))


h1 <- decostand(t_dna_c, method = "total")
h2 <- decostand(t_dna_g, method = "total")
h3 <- decostand(t_rna_c, method = "total")
h4 <- decostand(t_rna_g, method = "total")

library(scales)

h1 <- perc_filter(t_dna_c)
h1_s <- decostand(h1, method = "standardize", MARGIN = 1)
hm1 <- melt(t(h1))
h2 <- perc_filter(t_rna_c)
h2_s <- decostand(h2, method = "standardize", MARGIN = 1)
hm2 <- melt(t(h2))
h3 <- perc_filter(t_dna_g)
h3_s <- decostand(h3, method = "standardize", MARGIN = 1)
hm3 <- melt(t(h3))
h4 <- perc_filter(t_rna_g)
h4_s <- decostand(h4, method = "standardize", MARGIN = 1)
hm4 <- melt(t(h4))

library(RColorBrewer)
scaleyellowred <- colorRampPalette()(1000)

ggplot(hm1, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient2(low = "darkgreen",mid = "lightyellow", high = "darkred") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust =1))

ggplot(hm2, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient2(low = "darkgreen",mid = "lightyellow", high = "darkred") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust =1))

ggplot(hm3, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient2(low = "darkgreen",mid = "lightyellow", high = "darkred") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust =1))

ggplot(hm4, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient2(low = "darkgreen",mid = "lightyellow", high = "darkred") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust =1))



###############################################
# sig diffs
###############################################

library(reshape)
library(ggplot2)
library(ggforce)



head(class_dna_f)

#Load the abundance table 

a <- decostand(class_rna_f1, method = "total", MARGIN = 2)

round(class_dna_f1, 3)

class_dna_p <- percent_convert(class_dna_f1)
gen_dna_p <- percent_convert(gen_dna_f2)
spec_dna_p <- percent_convert(spec_dna_f2)

class_rna_p <- percent_convert(class_rna_f2)
gen_rna_p <- percent_convert(gen_rna_f2)
spec_rna_p <- percent_convert(spec_rna_f2)

perc_filter = function(x){
  a <- as.data.frame(t(x))
  maxab <- apply(a, 2, max)
  n1 <- names(which(maxab < 0.1))        ##### control % filter
  a[, -which(names(a) %in% n1)]
}

colSums(phy_dna_f)
phy_dna_p_f <- perc_filter(phy_dna_f1)
class_dna_p_f <- perc_filter(class_dna_f1)
gen_dna_p_f <- perc_filter(gen_dna_f2)
spec_dna_p_f <- perc_filter(spec_dna_f2)

phy_rna_p_f <- perc_filter(phy_rna_f1)
class_rna_p_f <- perc_filter(class_rna_f1)
gen_rna_p_f <- perc_filter(gen_rna_f2)
spec_rna_p_f <- perc_filter(spec_rna_f2)

class_rna_p_f_log <- decostand(class_rna_p_f, method = "log")

abund_table <- as.data.frame(phy_dna_p_f)



#abund_table<-t(abund_table)
#env <- as.data.frame(t(env))
#Use countries as grouping information
groups <- factor(dna_factors$Depth) # levels=c('Ah','Bw','Cox')
#groups <- as.factor(env[,2], levels=c('5-10','10-15','15-20','20-25', '25-30'))
#Apply normalisation (either use relative or log-relative transformation)
#data<-abund_table/rowSums(abund_table)
#data<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
data<-as.data.frame(abund_table)

kruskal.wallis.alpha=0.1
kruskal.wallis.table <- data.frame()
for (i in 1:dim(data)[2]) {
  ks.test <- kruskal.test(data[,i], g=groups)
  # Store the result in the data frame
  kruskal.wallis.table <- rbind(kruskal.wallis.table,
                                data.frame(id=names(data)[i],
                                           p.value=ks.test$p.value
                                ))
  # Report number of values tested
  cat(paste("Kruskal-Wallis test for ",names(data)[i]," ", i, "/", 
            dim(data)[2], "; p-value=", ks.test$p.value,"\n", sep=""))
}


kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]

kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, 
                                    size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)

kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value,
                                                   decreasing=FALSE), ]
kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
kruskal.wallis.table$q.value <- p.adjust(kruskal.wallis.table$p.value, method = "fdr")

last.significant.element <- max(which(kruskal.wallis.table$q.value <= kruskal.wallis.alpha))
selected <- 1:last.significant.element
diff.cat.factor <- kruskal.wallis.table$id[selected]
diff.cat <- as.vector(diff.cat.factor)

print(kruskal.wallis.table[selected,])
dim(kruskal.wallis.table[selected,])
write.csv(kruskal.wallis.table, "kw_rna_phy_depth_tab", quote =FALSE, row.names = FALSE)


#Now we plot taxa significantly different between the categories
df<-NULL
for(i in diff.cat){
  tmp<-data.frame(data[,i],groups,rep(paste(i," q = ",round(kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"],5),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Type","Taxa")

n_pages <- ceiling(length(levels(df$Taxa))/12)
n_pages

for (i in seq_len(n_pages)){ 
  print(ggplot(df,aes(Type,Value,colour=Type)) +
          ylab("Percent abundnace") +
          xlab("Depth (cm)")+
          geom_boxplot(outlier.shape=NA) + geom_jitter(alpha = 0.7)+
          theme_bw() +
          scale_color_brewer(name = "", palette = "Set1") +
          theme(axis.text.x=element_text(angle=45,hjust=1, size = 12),
                legend.position = "none",
                axis.text.y = element_text(size = 12),
                strip.text = element_text(size = 12),
                axis.title = element_text(size = 14, face = "bold")) +
          facet_wrap_paginate( ~ Taxa , scales="free_y", ncol=3, nrow=4 ,page = i))
}

library(ggforce)



library(data.table)
setDT(class_rna_p_f, keep.rownames = TRUE)[]
colnames(class_rna_p_f)[1] <- "Label"
m_class <- melt(class_rna_p_f)
df <- m_class
df <- full_join(df, rna_factors_sub, by = "Label")

df$variable <- factor(df$variable, levels = c("Actinobacteria", "Bacilli", "Nitrospira", "Verrucomicrobiae", "Clostridia", "Acidobacteriia", "Ktedonobacteria", "Solibacteres", "Alphaproteobacteria", "Gammaproteobacteria", "Spartobacteria", "Deltaproteobacteria", "Blastocatellia", "Eurotiomycetes", "Gemmatimonadetes", "Leotiomycetes", "Nitrososphaeria", "Planctomycetia", "Sphingobacteriia", "Thermoleophilia"))

n_pages <- ceiling(length(levels(df$variable))/12)
n_pages


for (i in seq_len(n_pages)){ 
  print(ggplot(df, aes(x=Site, y = value, group = Depth , fill = Site)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap_paginate( ~ variable , scales="free_y", ncol=3, nrow=4 ,page = i) +
  theme_bw() +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none"))
}

ggplot(df, aes(x=Site, y = value, group = Depth , fill = Site)) +
        geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
        #ylim(5, 6.25) +
        scale_fill_brewer(palette = "Set1") +
        facet_wrap_paginate( ~ variable , scales="free_y", ncol=3, nrow=4 ,page = 1) +
        theme_bw() +
        theme(legend.title = element_text(size=14, face="bold"),
              legend.text = element_text(size=14),
              axis.title = element_text(size = 14),
              axis.text.x = element_text(size = 12, angle = 45, hjust =1),
              axis.text.y = element_text(size = 12),
              strip.text = element_text(size = 14),
              panel.border = element_rect(fill = NA, size=0.5),
              strip.background = element_rect(color = "black", size = 0.5),
              legend.position="none")
