##load libraries####################################################################################################
library(plyr)
library(data.table)
library(dplyr)
library(vegan)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggvegan)
library(metacoder)
##setup functions###################################################################################################


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

setwd("/home/jonny/Documents/data/2015_data/bg/foam_set/final_counts/")



e10names <- read.csv("all_e10_names_final_sorted", header = TRUE)
foam_hierarchy <- read.csv("foam_hierarchy_for_r.csv", sep = "\t", header = TRUE)
foam_hierarchy_u <- unique(foam_hierarchy)

annotations <- full_join(e10names, foam_hierarchy_u, by = "KO")
annotations_filt <- subset(annotations, Gene != "NA")
head(annotations_filt)
annotations_top_lvl <- annotations[1:4]
annotations_top_lvl <- unique(annotations_top_lvl)
annotations_top_lvl[!duplicated(annotations_top_lvl$Gene),]
annotations_final <- annotations_top_lvl[!duplicated(annotations_top_lvl$Gene),]
head(annotations_final)
tail(e10names)

annotations_no_na <- annotations_final[complete.cases(annotations_final), ]

missing <- subset(annotations_final, is.na(annotations_final$L1))
missing_kos <- as.data.frame(missing)
missing_kos_u<- unique(missing_kos)
write.csv(missing_kos_u, file="missing_kos", quote = FALSE, row.names = FALSE)



rna_counts <- read.csv("rna_counts.final", header = TRUE)
head(rna_counts)
dna_counts <- read.csv("dna_counts_final_cut", header = TRUE)
head(dna_counts)


rna_count_table <- full_join(annotations_final, rna_counts, by = "Gene")
dna_count_table <- full_join(annotations_final, dna_counts, by = "Gene")


######restart from here ##########
write.csv(rna_count_table, "rna_count_table", quote = FALSE, row.names = FALSE)
write.csv(dna_count_table, "dna_count_table", quote = FALSE, row.names = FALSE)
rna_count_table <- read.csv("rna_count_table", header = TRUE)
dna_count_table <- read.csv("dna_count_table", header = TRUE)

rna_count_table_no_na <- na.omit(rna_count_table)
dna_count_table_no_na <- na.omit(dna_count_table)

write.csv(rna_count_table_no_na, "rna_count_table_no_na", quote = FALSE, row.names = FALSE)
write.csv(dna_count_table_no_na, "dna_count_table_no_na", quote = FALSE, row.names = FALSE)
rna_c_nona_2 <- read.csv("rna_count_table_no_na", header =TRUE)
dna_c_nona_2 <- read.csv("dna_count_table_no_na", header =TRUE)





### do rpkm here #####


c_fix_dna <- read.csv("c_fix_dna", header = TRUE)
c_fix_rna <- read.csv("c_fix_rna", header = TRUE)
c_gene_lengths <- read.csv("c_fix_gene_lengths", header = TRUE)

gene_lengths <- as.vector(c_gene_lengths$length)



gene_lengths <- read.csv("gene_lengths", header = TRUE)
tail(c_fix_dna)
colnames(gene_lengths)[1] <- "Gene"

a <- sweep(as.matrix(c_fix_dna[2:46]), 1, gene_lengths, '/')
b <- a*1000
d <- as.data.frame(b)
head(d)
final_c_dna <- cbind(c_fix_dna[47:50], d)

dna_cbb_fix <- aggregate(c_fix_dna[2:46], by=list(c_fix_dna$Cat), 
          FUN=sum, na.rm=TRUE)

a <- sweep(as.matrix(c_fix_rna[2:16]), 1, gene_lengths, '/')
b <- a*1000
d <- as.data.frame(b)
head(d)
final_c_dna <- cbind(c_fix_rna[17:20], d)

rna_cbb_fix <- aggregate(c_fix_rna[2:16], by=list(c_fix_rna$Cat), 
                         FUN=sum, na.rm=TRUE)


test <- full_join(gene_lengths, as.data.frame(rna_c_nona_2[1:2]), by ="Gene")
test_no_na <- test[complete.cases(test), ]
gene_lengths <- as.vector(test_no_na$Length)


a <- sweep(as.matrix(rna_c_nona_2[5:25]), 1, gene_lengths, '/')
b <- a*1000
d <- as.data.frame(b)
head(d)
final_rna_norm <- cbind(rna_c_nona_2[1:4], d)
colSums(final_rna_norm[5:25])

write.csv(final_dna_norm, file="final_dna_norm_l", quote = FALSE, row.names = FALSE)
write.csv(final_rna_norm, file="final_rna_norm_l", quote = FALSE, row.names = FALSE)

final_rna_norm_l <- read.csv(file="final_rna_norm_l", header = T)
final_rna_norm_l_reduced <- subset(final_rna_norm_l, select = -c(G11B.1.rna,G11C.1.rna,G1B.1.rna,G1C.1.rna,G1C.1.rna,V10B.1.rna,V10C.1.rna))
final_dna_norm_l <- read.csv(file="final_dna_norm_l", header = T)

agg_dna_kos_rpkb <- aggregate(final_dna_norm_l[5:49], by=list(final_dna_norm_l$KO), 
                                            FUN=sum, na.rm=TRUE)
colnames(agg_dna_kos_rpkb)[1] <- "KO"
agg_rna_kos_rpkb <- aggregate(final_rna_norm_l_reduced[5:19], by=list(final_rna_norm_l_reduced$KO), 
                              FUN=sum, na.rm=TRUE)
colnames(agg_rna_kos_rpkb)[1] <- "KO"

write.csv(agg_dna_kos_rpkb, "agg_dna_kos_rpkb", quote = FALSE, row.names = FALSE)
write.csv(agg_rna_kos_rpkb, "agg_rna_kos_rpkb", quote = FALSE, row.names = FALSE)


dna_h <- final_dna_norm_l[2:4]
dna_h_u <- unique(dna_h)

dna_kos_w_h <- full_join(dna_h_u, agg_dna_kos_rpkb, by ="KO")
rna_kos_w_h <- full_join(dna_h_u, agg_rna_kos_rpkb, by ="KO")

write.csv(dna_kos_w_h, "dna_kos_w_h", quote = FALSE, row.names = FALSE)
write.csv(rna_kos_w_h, "rna_kos_w_h", quote = FALSE, row.names = FALSE)

dna_kos_w_h <- read.csv("dna_kos_w_h", header= TRUE)
rna_kos_w_h <- read.csv("rna_kos_w_h", header= TRUE)

(colSums(as.data.frame(colSums(dna_kos_w_h[4:48])))/45)/1000
(colSums(as.data.frame(colSums(rna_kos_w_h[4:18])))/15)/1000

dna_kos_filt <- dna_kos_w_h[rowSums(dna_kos_w_h[4:48] > 300) >= 3, ]
rna_kos_filt <- rna_kos_w_h[rowSums(rna_kos_w_h[4:18] > 300) >= 3, ]

rna_tots <- c(11627682,23106560,15635254,19673146,22284836,16063534,24953610,31093464,26014696,27473168,21347046,20568324,27537876,24282974,26893654)
dna_tots <- c(178728628,287225804,233779812,198068956,197243492,24816830,26340810,18915710,28260112,24928364,27488284,23905258,21247590,22449786,27860216,284810538,259674992,221777750,238261726,206798312,23912724,25756704,21757908,25824124,22854534,29122676,30062108,32796018,22145254,25676230,195609520,183941758,202565964,159621730,202690830,25658524,26820300,23849806,27962464,25071848,28521492,23497974,24222736,20536298,20453460)

rna_scaling <- as.vector(colSums(all_kos_rna[4:18]))
rna_scaling <- rna_scaling/1000000
dna_scaling <- as.vector(colSums(all_kos_dna[4:48]))
dna_scaling <- dna_scaling/1000000

rna_cbb_fix <- write.csv(rna_cbb_fix, "rna_cbb_fix", quote = FALSE, row.names = FALSE)
dna_cbb_fix <- write.csv(dna_cbb_fix, "dna_cbb_fix", quote = FALSE, row.names = FALSE)

rna_cbb_fix <- read.csv("rna_cbb_fix", header= TRUE)
dna_cbb_fix <- read.csv("dna_cbb_fix", header =TRUE)


all_kos_rna <- rbind(rna_kos_filt, rna_cbb_fix)
all_kos_dna <- rbind(dna_kos_filt, dna_cbb_fix)

a <- sweep(as.matrix(all_kos_dna[4:48]), 2, dna_scaling, '/')
d <- as.data.frame(a)
final_dna_kos_tpm <- cbind(all_kos_dna[1:3], d)


a <- sweep(as.matrix(all_kos_rna[4:18]), 2, rna_scaling, '/')
d <- as.data.frame(a)
final_rna_kos_tpm <- cbind(all_kos_rna[1:3], d)

write.csv(final_dna_kos_tpm, file = "final_dna_kos_tpm_filt_0.1", quote = FALSE, row.names=FALSE)
write.csv(final_rna_kos_tpm, file = "final_rna_kos_tpm_filt_0.1", quote = FALSE, row.names=FALSE)


dna_tmp <- read.csv("final_dna_tpm", header = TRUE)
rna_tmp <- read.csv("final_rna_tpm", header = TRUE)

test <- aggregate(final_dna_kos_tpm[4:48], by=list(final_dna_kos_tpm$L1), 
                                            FUN=median, na.rm=TRUE)
test2 <- aggregate(final_rna_kos_tpm[4:18], by=list(final_rna_kos_tpm$L1), 
                  FUN=median, na.rm=TRUE)

head(dna_tmp)
colSums(dna_tmp[5:49], na.rm = TRUE)

#Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#Divide the RPK values by the “per million” scaling factor. This gives you TPM.

final_rna_tpm[final_rna_tpm == 0] <- NA
final_dna_tpm[final_dna_tpm == 0] <- NA

agg_dna_counts_median_l1_test1 <- aggregate(final_dna_tpm[5:49], by=list(final_dna_tpm$L1), 
                                      FUN=median, na.rm=TRUE)


agg_rna_counts_sum_l1 <- aggregate(final_rna_tpm[5:25], by=list(final_rna_tpm$L1), 
                         FUN=sum, na.rm=TRUE)

agg_dna_counts_sum_l1 <- aggregate(final_dna_tpm[5:49], by=list(final_dna_tpm$L1), 
                            FUN=sum, na.rm=TRUE)

agg_rna_counts_mean_l1 <- aggregate(final_rna_tpm[5:25], by=list(final_rna_tpm$L1), 
                                FUN=mean, na.rm=TRUE)

agg_dna_counts_mean_l1 <- aggregate(final_dna_tpm[5:49], by=list(final_rna_tpm$L1), 
                                FUN=mean, na.rm=TRUE)

agg_rna_counts_median_l1 <- aggregate(final_rna_tpm[5:25], by=list(final_rna_tpm$L1), 
                                FUN=median, na.rm=TRUE)

agg_dna_counts_media_l1 <- aggregate(final_dna_tpm[5:49], by=list(final_rna_tpm$L1), 
                                FUN=median, na.rm=TRUE)

agg_rna_counts_sum_l2 <- aggregate(final_rna_tpm[5:25], by=list(final_rna_tpm$L2), 
                                   FUN=sum, na.rm=TRUE)

agg_dna_counts_sum_l2 <- aggregate(final_dna_tpm[5:49], by=list(final_dna_tpm$L2), 
                                   FUN=sum, na.rm=TRUE)

agg_rna_counts_mean_l2 <- aggregate(final_rna_tpm[5:25], by=list(final_rna_tpm$L2), 
                                    FUN=mean, na.rm=TRUE)

agg_dna_counts_mean_l2 <- aggregate(final_dna_tpm[5:49], by=list(final_rna_tpm$L2), 
                                    FUN=mean, na.rm=TRUE)

agg_rna_counts_median_l2 <- aggregate(final_rna_tpm[5:25], by=list(final_rna_tpm$L2), 
                                      FUN=median, na.rm=TRUE)

agg_dna_counts_media_l2 <- aggregate(final_dna_tpm[5:49], by=list(final_rna_tpm$L2), 
                                     FUN=median, na.rm=TRUE)

write.csv(agg_dna_counts_median_l1_test1,"dna_tot_norm_test_l1", quote = FALSE, row.names=FALSE)

write.csv(agg_rna_counts_sum_l1,"rna_sum_l1", quote = FALSE, row.names=FALSE)
write.csv(agg_rna_counts_sum_l2,"rna_sum_l2", quote = FALSE, row.names=FALSE)
write.csv(agg_rna_counts_mean_l1,"rna_mean_l1", quote = FALSE, row.names=FALSE)
write.csv(agg_rna_counts_mean_l2,"rna_mean_l2", quote = FALSE, row.names=FALSE)
write.csv(agg_rna_counts_median_l1,"rna_median_l1", quote = FALSE, row.names=FALSE)
write.csv(agg_rna_counts_median_l2,"rna_median_l2", quote = FALSE, row.names=FALSE)

write.csv(agg_dna_counts_sum_l1,"dna_sum_l1", quote = FALSE, row.names=FALSE)
write.csv(agg_dna_counts_sum_l2,"dna_sum_l2", quote = FALSE, row.names=FALSE)
write.csv(agg_dna_counts_mean_l1,"dna_mean_l1", quote = FALSE, row.names=FALSE)
write.csv(agg_dna_counts_mean_l2,"dna_mean_l2", quote = FALSE, row.names=FALSE)
write.csv(agg_dna_counts_media_l1,"dna_median_l1", quote = FALSE, row.names=FALSE)
write.csv(agg_dna_counts_media_l2,"dna_median_l2", quote = FALSE, row.names=FALSE)


write.csv(agg_dna_counts,"length_normalised_dna", quote = FALSE, row.names=FALSE)

read.csv("length_normalised_rna", header = FALSE)


colnames(agg_dna_counts_sum_l2)[1] <- "L2"
foam_h <- foam_hierarchy_u[1:2]
foam_h_2 <- unique(foam_h)
write.csv(foam_h_2, file = "foam_h_2", quote = FALSE, row.names = FALSE)
foam_h_2 <- read.csv("foam_h_2", header = TRUE)

final_dna <- full_join(foam_h_2, agg_dna_counts_sum_l2, by = "L2")
final_dna[is.na(final_dna)] <- 0
final_rna <- full_join(foam_h_2, agg_rna_counts, by = "L2")
final_rna[is.na(final_rna)] <- 0

write.csv(final_dna, file = "dna_sums_l2_w_names", quote = FALSE, row.names = FALSE)
write.csv(agg_rna_counts, file = "agg_rna_counts", quote = FALSE, row.names = FALSE)


write.csv(final_rna, file = "final_rna", quote = FALSE, row.names = FALSE)
write.csv(final_dna, file = "final_dna", quote = FALSE, row.names = FALSE)

colnames(agg_rna_counts)[1] <- "L2"
unique(final_rna_counts)

colnames(agg_dna_counts)[1] <- "L2"



########################################
# actual graphs
########################################

# load functions


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



# load data

bg_dna <- read.csv("final_dna", header = TRUE)
bg_rna <- read.csv("final_rna2.csv", header = TRUE)
dna_meta <- read.csv("dna_meta.csv", header=TRUE)
rna_meta <- read.csv("rna_meta.csv", header=TRUE)

rna_tots <- c(11627682,23106560,15635254,19673146,22284836,16063534,24953610,31093464,26014696,27473168,21347046,20568324,27537876,24282974,26893654)
dna_tots <- c(178728628,287225804,233779812,198068956,197243492,24816830,26340810,18915710,28260112,24928364,27488284,23905258,21247590,22449786,27860216,284810538,259674992,221777750,238261726,206798312,23912724,25756704,21757908,25824124,22854534,29122676,30062108,32796018,22145254,25676230,195609520,183941758,202565964,159621730,202690830,25658524,26820300,23849806,27962464,25071848,28521492,23497974,24222736,20536298,20453460)


#normalise RNA
a <- sweep(as.matrix(bg_rna[4:18]), 2, rna_tots, '/')
b <- a*1000000
d <- as.data.frame(b)
final_rna_norm <- cbind(bg_rna[1:3], d)
write.csv(final_rna_norm, file="final_rna_norm", quote = FALSE, row.names = FALSE)

#normalsie DNA

a <- sweep(as.matrix(bg_dna[4:48]), 2, dna_tots, '/')
b <- a*1000000
d <- as.data.frame(b)
final_dna_norm <- cbind(bg_dna[1:3], d)
write.csv(final_dna_norm, file="final_dna_norm", quote = FALSE, row.names = FALSE)


#################################
# start here # nope
#################################
setwd("/home/jonny/Documents/data/2015_data/bg/foam_set/final_counts/")

bg_dna <- read.csv("final_dna_kos_tpm_filt_0.1", header = TRUE)
bg_rna <- read.csv("final_rna_kos_tpm_filt_0.1", header = TRUE)


foam_h_2 <- read.csv("foam_h_2", header= TRUE)

test <- read.csv("dna_tot_norm_test_l1", header = T)
test[is.na(test)] <- 0

dna_sum_l1 <- read.csv("dna_sum_l1", header = TRUE)
dna_sum_l1[is.na(dna_sum_l1)] <- 0

dna_sum_l2 <- read.csv("dna_sum_l2", header = TRUE)
dna_sum_l2 <- full_join(foam_h_2, dna_sum_l2, by = "Group.1")
dna_sum_l2[is.na(dna_sum_l2)] <- 0

dna_mean_l1 <- read.csv("dna_mean_l1", header = TRUE)
dna_mean_l1[is.na(dna_mean_l1)] <- 0

dna_mean_l2 <- read.csv("dna_mean_l2", header = TRUE)
dna_mean_l2 <- full_join(foam_h_2, dna_mean_l2, by = "Group.1")
dna_mean_l2[is.na(dna_mean_l2)] <- 0

dna_median_l1 <- read.csv("dna_median_l1", header = TRUE)
dna_median_l1[is.na(dna_median_l1)] <- 0


dna_median_l2 <- read.csv("dna_median_l2", header = TRUE)
dna_median_l2 <- full_join(foam_h_2, dna_median_l2, by = "Group.1")
dna_median_l2[is.na(dna_median_l2)] <- 0


rna_sum_l1 <- read.csv("rna_sum_l1", header = TRUE)
rna_sum_l1[is.na(rna_sum_l2)] <- 0

rna_sum_l2 <- read.csv("rna_sum_l2", header = TRUE)
rna_sum_l2 <- full_join(foam_h_2, rna_sum_l2, by = "Group.1")
rna_sum_l2[is.na(rna_sum_l2)] <- 0

rna_mean_l1 <- read.csv("rna_mean_l1", header = TRUE)
rna_mean_l1[is.na(rna_mean_l1)] <- 0

rna_mean_l2 <- read.csv("rna_mean_l2", header = TRUE)
rna_mean_l2 <- full_join(foam_h_2, rna_mean_l2, by = "Group.1")
rna_mean_l1[is.na(rna_mean_l1)] <- 0

rna_median_l1 <- read.csv("rna_median_l1", header = TRUE)
rna_median_l1[is.na(rna_median_l1)] <- 0

rna_median_l2 <- read.csv("rna_median_l2", header = TRUE)
rna_median_l2 <- full_join(foam_h_2, rna_median_l2, by = "Group.1")
rna_median_l2[is.na(rna_median_l2)] <- 0

dna_meta <- read.csv("dna_meta.csv", header=TRUE)
rna_meta <- read.csv("rna_meta.csv", header=TRUE)

#####################
# start here
#####################


setwd("/home/jonny/Documents/data/2015_data/bg/foam_set/")

dna_tmp <- read.csv("final_dna_tpm_e50", header = TRUE)
rna_tmp <- read.csv("final_rna_tpm_e50", header = TRUE)



dna_h <- dna_tmp[2:4]
dna_h <- unique(dna_h)
test <- aggregate(dna_tmp[5:49], by=list(dna_tmp$KO), 
                  FUN=sum, na.rm=TRUE)
colnames(test)[1] <- "KO"
head(test)
test_w_h <- full_join(dna_h, test, by = "KO")
write.table(test_w_h, "test3", quote = FALSE, row.names = FALSE, sep= "\t")
dna_h2 <- dna_tmp[3:4]
dna_h2 <- unique(dna_h2)

test_3 <- full_join(dna_h2, test2, by = "L2")
write.table(test_3, "dna_tmp_e50", quote = FALSE, row.names = FALSE, sep= "\t")


rna_h <- rna_tmp[2:4]
rna_h <- unique(rna_h)
test_rna <- aggregate(rna_tmp[5:25], by=list(rna_tmp$KO), 
                  FUN=sum, na.rm=TRUE)
colnames(test_rna)[1] <- "KO"
head(test_rna)
test_rna_w_h <- full_join(rna_h, test_rna, by = "KO")
write.table(test_rna_w_h, "rna_tmp_e50", quote = FALSE, row.names = FALSE, sep= "\t")

colnames(test)[1] <- "KO"
head(test)
test_w_h <- full_join(dna_h, test, by = "KO")
write.table(test_w_h, "test3", quote = FALSE, row.names = FALSE, sep= "\t")
dna_h2 <- dna_tmp[3:4]
dna_h2 <- unique(dna_h2)




setwd("/home/jonny/Documents/data/2015_data/bg/foam_set/final_counts/")
bg_dna <- read.csv("final_dna_kos_tpm_filt_fixed.csv", header = TRUE)
bg_rna <- read.csv("final_rna_kos_tpm_filt_fixed.csv", header = TRUE)
dna_meta <- read.csv("dna_meta.csv", header=TRUE)
rna_meta <- read.csv("rna_meta.csv", header=TRUE)


head(bg_dna)
##################
## carbon cycling
##################
#bg_rna <- read.csv("test3", sep= ",", header = T)

c_dna <- subset(bg_dna, Cycle =="Carbon")
bg_dna[bg_dna == 0] <- NA

c_dna_ag <- aggregate(c_dna[5:49], by=list(c_dna$L1), 
                  FUN=median, na.rm=TRUE)

c_dna_ag2 <- aggregate(c_dna[5:49], by=list(c_dna$L2), 
                                   FUN=median, na.rm=TRUE)

c_rna <- subset(bg_rna, Cycle =="Carbon")
bg_rna[bg_rna == 0] <- NA

c_rna_ag <- aggregate(c_rna[5:19], by=list(c_rna$L1), 
                   FUN=median, na.rm=TRUE)


m_c_dna <- melt(c_dna_ag)
m_c_rna <- melt(c_rna_ag)
m_c_dna_w_fact <- full_join(m_c_dna, dna_meta, by ="variable")
m_c_rna_w_fact <- full_join(m_c_rna, rna_meta, by ="variable")
m_c_dna_se <- summarySE(m_c_dna_w_fact, measurevar="value", groupvars=c("Group.1","Site","Depth"))



ggplot(m_c_dna_se, aes(x=Site, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  #coord_cartesian(ylim = c(1000, 2500)) +
  facet_wrap( ~Group.1) +
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

ggplot(m_c_rna_w_fact, aes(x=Site, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~Group.1, scales = "free_y") +
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



#totes

m_c_dna_tots <- summarySE(m_c_dna_w_fact, measurevar="value", groupvars=c("Group.1"))
m_c_rna_tots <- summarySE(m_c_rna_w_fact, measurevar="value", groupvars=c("Group.1"))

ggplot(m_c_dna_tots, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #coord_cartesian(ylim = c(0, 6000)) +
  coord_flip() +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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

ggplot(m_c_rna_tots, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #coord_cartesian(ylim = c(0, 6000)) +
  coord_flip() +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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

###################################################################
# subset of carbon


Methylotrophy <- subset(bg_dna, L1 =="15_Methylotrophy")
mty <- aggregate(Methylotrophy[5:49], by=list(Methylotrophy$L2), 
                      FUN=median, na.rm=TRUE)

Methanogenesis <- subset(bg_dna, L1 =="14_Methanogenesis")
mth <- aggregate(Methanogenesis[5:49], by=list(Methanogenesis$L2), 
                 FUN=median, na.rm=TRUE)
c_fix <- subset(bg_dna, L1 =="C_fix")
cfix <- aggregate(c_fix[5:49], by=list(c_fix$L2), 
                 FUN=median, na.rm=TRUE)

Methylotrophy_rna <- subset(bg_rna, L1 =="15_Methylotrophy")
mty_rna <- aggregate(Methylotrophy_rna[5:19], by=list(Methylotrophy_rna$L2), 
                 FUN=median, na.rm=TRUE)
Methanogenesis_rna <- subset(bg_rna, L1 =="14_Methanogenesis")
mth_rna <- aggregate(Methanogenesis_rna[5:19], by=list(Methanogenesis_rna$L2), 
                     FUN=median, na.rm=TRUE)
c_fix_rna <- subset(bg_rna, L1 =="C_fix")
cfix_rna <- aggregate(c_fix_rna[5:19], by=list(c_fix_rna$L2), 
                  FUN=median, na.rm=TRUE)

methyl_d <- melt(mty)
methyl_r <- melt(mty_rna)
methan_d <- melt(mth)
methan_r <- melt(mth_rna)
c_fix_d <- melt(cfix)
c_fix_r <- melt(cfix_rna)

methyl_dna_fact <- full_join(methyl_d, dna_meta, by ="variable")
methyl_rna_fact <- full_join(methyl_r, rna_meta, by ="variable")
methan_dna_fact <- full_join(methan_d, dna_meta, by ="variable")
methan_rna_fact <- full_join(methan_r, rna_meta, by ="variable")
cfix_dna_fact <- full_join(c_fix_d, dna_meta, by ="variable")
cfix_rna_fact <- full_join(c_fix_r, rna_meta, by ="variable")

d_methyl_stats <- compare_means(value ~ Depth, data = methyl_dna_fact, 
                              group.by = "Group.1", method = "wilcox.test", p.adjust.method= "fdr")
write.csv(d_methyl_stats, "methyl_d_pw_stats")
r_methyl_stats <- compare_means(value ~ Depth, data = methyl_rna_fact, 
                              group.by = "Group.1", method = "wilcox.test", p.adjust.method= "fdr")
write.csv(d_methyl_stats, "methyl_r_pw_stats")
d_methan_stats <- compare_means(value ~ Depth, data = methan_dna_fact, 
                                group.by = "Group.1", method = "wilcox.test", p.adjust.method= "fdr")
write.csv(d_methan_stats, "methan_d_pw_stats")
r_methan_stats <- compare_means(value ~ Depth, data = methan_rna_fact, 
                                group.by = "Group.1", method = "wilcox.test", p.adjust.method= "fdr")
write.csv(r_methan_stats, "methan_r_pw_stats")
cfix_dna_fact[is.na(cfix_dna_fact)] <- 0
cfix_rna_fact[is.na(cfix_rna_fact)] <- 0

d_cfix_stats <- compare_means(value ~ Depth, data = cfix_dna_fact, 
                                group.by = "Group.1", method = "wilcox.test", p.adjust.method= "fdr")
write.csv(d_cfix_stats, "cfix_d__pairwise_stats")

r_cfix_stats <- compare_means(value ~ Depth, data = cfix_rna_fact, 
                                group.by = "Group.1", method = "wilcox.test", p.adjust.method= "fdr")
write.csv(r_cfix_stats, "cfix_r_pairwise_stats")

make_box(methyl_dna_fact)
make_box(methyl_rna_fact)
make_box(methan_dna_fact)
make_box(methan_rna_fact)
make_box(cfix_dna_fact)
make_box(cfix_rna_fact)

make_box <- function(x){
  ggplot(x, aes(x=Depth, y = value, group = Depth , fill = Group.1)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  stat_compare_means(method = "kruskal.test",label= "p.format", label.y.npc = 0.9) +
  facet_wrap(~ Group.1, scales = "free_y") +
  theme_bw() +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none") }


methyl_dna_tots <- summarySE(methyl_dna_fact, measurevar="value", groupvars=c("Group.1"))
methyl_rna_tots <- summarySE(methyl_rna_fact, measurevar="value", groupvars=c("Group.1"))
methan_dna_tots <- summarySE(methan_dna_fact, measurevar="value", groupvars=c("Group.1"))
methan_rna_tots <- summarySE(methan_rna_fact, measurevar="value", groupvars=c("Group.1"))
cfix_dna_tots <- summarySE(cfix_dna_fact, measurevar="value", groupvars=c("Group.1"))
cfix_rna_tots <- summarySE(cfix_rna_fact, measurevar="value", groupvars=c("Group.1"))


make_bar(methyl_dna_tots)
make_bar(methyl_rna_tots)
make_bar(methan_dna_tots)
make_bar(methan_rna_tots)
make_bar(cfix_dna_tots)
make_bar(cfix_rna_tots)


make_bar = function(x){ggplot(x, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #coord_cartesian(ylim = c(0, 6000)) +
  coord_flip() +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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
}

ggplot(m_c_rna_tots, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #coord_cartesian(ylim = c(0, 6000)) +
  coord_flip() +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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


a <- compare_means(value ~ Depth, data = meth_fact, 
                   group.by = "L2", method = "anova", p.adjust.method= "fdr")

b <- compare_means(value ~ Depth, data = m_c_rna_w_fact, 
                   group.by = "Group.1", method = "kruskal.test", p.adjust.method= "fdr")

write.csv(b, "c_l1_adjustedpvals_rna")

ggplot(m_c_dna_w_fact, aes(x=Depth, y = value, group = Depth , fill = Group.1)) +
  geom_boxplot() +
  #ylim(5, 6.25) +
  scale_fill_brewer(palette = "Set1") +
  stat_compare_means(method = "kruskal.test",label= "p.format", label.y.npc = 0.9, p.adjust.method="fdr") +
  #coord_cartesian(ylim = c(1000, 2500)) +
  facet_wrap(~ Group.1, scales = "free_y") +
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

##### methanol taxa



head(dna_tmp)
meth_genes <- subset(dna_tmp, L1 == "15_Methylotrophy")
meth_genes_rna <- subset(rna_tmp, L1 == "15_Methylotrophy")

m <- as.data.frame(meth_genes$Gene)
write.csv(m, "meth_dna", quote =FALSE, row.names = FALSE)
m_rna <- as.data.frame(meth_genes_rna$Gene)
write.csv(m_rna, "meth_rna", quote =FALSE, row.names = FALSE)



## xoxf taxa 

x_genes <- read.csv("xoxf_taxa.csv", header = T, sep = "\t")

x_genes_dna <- full_join(dna_tmp, x_genes, by = "Gene")
x_genes_rna <- full_join(rna_tmp, x_genes, by = "Gene")

subset(dna_tmp, Gene == "10625556")

x_genes_class <- aggregate(x_genes_dna[5:49], by=list(x_genes_dna$Class), 
                             FUN=sum, na.rm = TRUE)

write.csv(x_genes_class, "x_genes_class", quote = FALSE, row.names = F)

x_genes_rna_class <- aggregate(x_genes_rna[5:22], by=list(x_genes_rna$Class), 
                           FUN=median, na.rm=TRUE)

write.csv(x_genes_rna, "x_genes_rna_class", quote = FALSE, row.names = F)

amm_genes_dna <- read.csv("amm_genes_class", header = T)
amm_genes_rna <- read.csv("amm_genes_rna_class", header = T)

melt_am_dna <- melt(amm_genes_dna)
j_dna_amm <- full_join(melt_am_dna, dna_meta, by ="variable")
final_dna_amm <- summarySE(j_dna_amm, measurevar="value", groupvars=c("Group.1","Site","Depth"))
melt_am_rna <- melt(amm_genes_rna)
final_rna_amm <- full_join(melt_am_rna, rna_meta, by ="variable")



ggplot(j_dna_amm, aes(x=Group.1, y = value, fill = Site, group = Depth)) +
  geom_bar(stat = "summary", fun.y = "median", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  #geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  coord_flip() +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")

ggplot(final_dna_amm, aes(x=Group.1, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Spectral") +
  coord_flip() +
  theme_bw() +
  facet_wrap( ~ Site) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")

ggplot(final_rna_amm, aes(x=Group.1, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  #geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Spectral") +
  #coord_cartesian(ylim = c(1000, 2500)) +
  theme_bw() +
  coord_flip() +
  facet_wrap( ~ Site) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")





#########################
# Nitrogen
#########################

n_dna <- subset(bg_dna, Cycle =="Nitrogen")

n_dna_ag <- aggregate(n_dna[5:49], by=list(n_dna$L2), 
                      FUN=median, na.rm=TRUE)

n_rna <- subset(bg_rna, Cycle =="Nitrogen")

n_rna_ag <- aggregate(n_rna[5:19], by=list(n_rna$L2), 
                      FUN=median, na.rm=TRUE)


m_n_dna <- melt(n_dna_ag)
m_n_rna <- melt(n_rna_ag)
m_n_dna_w_fact <- full_join(m_n_dna, dna_meta, by ="variable")
m_n_rna_w_fact <- full_join(m_n_rna, rna_meta, by ="variable")
m_n_dna_se <- summarySE(m_n_dna_w_fact, measurevar="value", groupvars=c("Group.1","Site","Depth"))

ggplot(m_n_dna_se, aes(x=Site, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~Group.1, scales = "free_y") +
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

ggplot(m_n_rna_w_fact, aes(x=Site, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~Group.1, scales = "free_y") +
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

m_n_dna_tots <- summarySE(m_n_dna_w_fact, measurevar="value", groupvars=c("Group.1"))

ggplot(m_n_dna_tots, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0, 6000)) +
  coord_flip() +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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

m_n_rna_tots <- summarySE(m_n_rna_w_fact, measurevar="value", groupvars=c("Group.1"))

ggplot(m_n_rna_tots, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #coord_cartesian(ylim = c(0, 6000)) +
  coord_flip() +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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


n_dna_stats <- compare_means(value ~ Depth, data = m_n_dna_w_fact, 
                              group.by = "Group.1", method = "kruskal.test", p.adjust.method= "fdr")
write.csv(n_dna_stats, "n_dna_stats")

n_rna_stats <- compare_means(value ~ Depth, data = m_n_rna_w_fact, 
                              group.by = "Group.1", method = "kruskal.test", p.adjust.method= "fdr")
write.csv(n_rna_stats, "n_rna_stats")

make_box(m_n_dna_w_fact)
make_box(m_n_rna_w_fact)


make_box <- function(x){
  ggplot(x, aes(x=Depth, y = value, group = Depth , fill = Group.1)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1") +
    stat_compare_means(method = "kruskal.test",label= "p.format", label.y.npc = 0.9) +
    facet_wrap(~ Group.1, scales = "free_y") +
    theme_bw() +
    theme(legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=14),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust =1),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 14),
          panel.border = element_rect(fill = NA, size=0.5),
          strip.background = element_rect(color = "black", size = 0.5),
          legend.position="none") }


#ammonia oxidation genes



amm_genes <- read.csv("amm_genes.fasta.kaiju.names.out", header = T, sep = "\t")
amm_genes_dna <- full_join(dna_tmp, amm_genes, by = "Gene")
amm_genes_rna <- full_join(rna_tmp, amm_genes, by = "Gene")


amm_genes_class <- aggregate(amm_genes_dna[5:49], by=list(amm_genes_dna$Class), 
                             FUN=median, na.rm=T)

write.csv(amm_genes_class, "amm_genes_class", quote = FALSE, row.names = F)

amm_genes_rna <- aggregate(amm_genes_rna[5:25], by=list(amm_genes_rna$Class), 
                           FUN=sum, na.rm=TRUE)

write.csv(amm_genes_rna, "amm_genes_rna_class", quote = FALSE, row.names = F)

amm_genes_dna <- read.csv("amm_genes_class", header = T)
amm_genes_rna <- read.csv("amm_genes_rna_class", header = T)

melt_am_dna <- melt(amm_genes_dna)
j_dna_amm <- full_join(melt_am_dna, dna_meta, by ="variable")
final_dna_amm <- summarySE(j_dna_amm, measurevar="value", groupvars=c("Group.1","Site","Depth"))
melt_am_rna <- melt(amm_genes_rna)
final_rna_amm <- full_join(melt_am_rna, rna_meta, by ="variable")



ggplot(j_dna_amm, aes(x=Group.1, y = value, fill = Site, group = Depth)) +
  geom_bar(stat = "summary", fun.y = "median", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  #geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  coord_flip() +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")

ggplot(final_dna_amm, aes(x=Group.1, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Spectral") +
  coord_flip() +
  theme_bw() +
  facet_wrap( ~ Site) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")

ggplot(final_rna_amm, aes(x=Group.1, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  #geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Spectral") +
  #coord_cartesian(ylim = c(1000, 2500)) +
  theme_bw() +
  coord_flip() +
  facet_wrap( ~ Site) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust =1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")


#######################
#sulfur
######################

s_dna <- subset(bg_dna, Cycle =="Sulfur")

s_dna_ag <- aggregate(s_dna[5:49], by=list(s_dna$L2), 
                      FUN=median, na.rm=TRUE)

s_rna <- subset(bg_rna, Cycle =="Sulfur")

s_rna_ag <- aggregate(s_rna[5:19], by=list(s_rna$L2), 
                      FUN=median, na.rm=TRUE)


m_s_dna <- melt(s_dna_ag)
m_s_rna <- melt(s_rna_ag)
m_s_dna_w_fact <- full_join(m_s_dna, dna_meta, by ="variable")
m_s_dna_se <- summarySE(m_s_dna_w_fact, measurevar="value", groupvars=c("Group.1","Site","Depth"))
m_s_dna_tots <- summarySE(m_s_dna_w_fact, measurevar="value", groupvars=c("Group.1"))
m_s_rna_tots <- summarySE(m_s_rna_w_fact, measurevar="value", groupvars=c("Group.1"))

ggplot(m_s_dna_se, aes(x=Site, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~Group.1, scales = "free_y") +
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

ggplot(m_s_rna_w_fact, aes(x=Site, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~Group.1, scales = "free_y") +
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

ggplot(m_s_dna_tots, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  ylim(0,6200) +
  coord_flip() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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

ggplot(m_s_rna_tots, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  ylim(0,7000) +
  coord_flip() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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


s_dna_stats <- compare_means(value ~ Depth, data = m_s_dna_w_fact, 
                             group.by = "Group.1", method = "wilcox.test", p.adjust.method= "fdr")
write.csv(s_dna_stats, "s_dna_stats")

s_rna_stats <- compare_means(value ~ Depth, data = m_s_rna_w_fact, 
                             group.by = "Group.1", method = "wilcox.test", p.adjust.method= "fdr")
write.csv(s_rna_stats, "s_rna_pw_stats")

make_box(m_s_dna_w_fact)
make_box(m_s_rna_w_fact)


#########################
# Phosphorus
#########################

p_dna <- subset(bg_dna, Cycle =="Phosphorus")

p_dna_ag <- aggregate(p_dna[5:49], by=list(p_dna$L2), 
                      FUN=median, na.rm=TRUE)

tail(bg_rna)
p_rna <- subset(bg_rna, Cycle =="Phosphate")

p_rna_ag <- aggregate(p_rna[5:19], by=list(p_rna$L2), 
                      FUN=median, na.rm=TRUE)


m_p_dna <- melt(p_dna_ag)
m_p_rna <- melt(p_rna_ag)
m_p_dna_w_fact <- full_join(m_p_dna, dna_meta, by ="variable")
m_p_rna_w_fact <- full_join(m_p_rna, rna_meta, by ="variable")
m_p_dna_se <- summarySE(m_p_dna_w_fact, measurevar="value", groupvars=c("Group.1","Site","Depth"))
m_p_dna_tots <- summarySE(m_p_dna_w_fact, measurevar="value", groupvars=c("Group.1"))
m_p_rna_tots <- summarySE(m_p_rna_w_fact, measurevar="value", groupvars=c("Group.1"))

ggplot(m_p_dna_se, aes(x=Site, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~Group.1, scales = "free_y") +
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

ggplot(m_p_rna_w_fact, aes(x=Site, y = value, group = Depth , fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~Group.1, scales = "free_y") +
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

ggplot(m_p_dna_tots, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  ylim(0,4000) +
  coord_flip() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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

ggplot(m_p_rna_tots, aes(x=reorder(Group.1, value, FUN=mean), y = value, fill = Group.1)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  ylim(0,4000) +
  coord_flip() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
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

p_dna_stats <- compare_means(value ~ Depth, data = m_p_dna_w_fact, 
                             group.by = "Group.1", method = "kruskal.test", p.adjust.method= "fdr")
write.csv(p_dna_stats, "p_dna_stats")

p_rna_stats <- compare_means(value ~ Depth, data = m_p_rna_w_fact, 
                             group.by = "Group.1", method = "kruskal.test", p.adjust.method= "fdr")
write.csv(p_rna_stats, "p_rna_pw_stats")

make_box(m_p_dna_w_fact)
make_box(m_p_rna_w_fact)

############################################
# genome
############################################

# import bin abundanaces
setwd("/home/jonny/Documents/data/anvi_2015/final_bins/bins_across_samples")
bin_covs <- read.csv("abundance.txt", sep = ",", header = T)
ktedono <- subset(bin_covs, bins =="Bin_2_7_1")

k_dna <- ktedono[1:46]
k_rna <- ktedono[47:61]

m_k_dna <- melt(k_dna)
m_k_rna <- melt(k_rna)
#write.csv(dna_factors, "dna_factors", row.names = F, quote = F)
dna_factors <- read.csv("dna_factors", header = T)
#write.csv(rna_factors_sub, "rna_factors", row.names = F, quote = F)
rna_factors <- read.csv("rna_factors", header = T)

final_k_dna <- full_join(m_k_dna, dna_factors, by = "variable")
final_dna_k_se <- summarySE(final_k_dna, measurevar="value", groupvars=c("Site", "Depth"))
final_k_rna <- full_join(m_k_rna, rna_factors, by = "variable")


ggplot(final_dna_k_se, aes(x=Site, y = value, group = Depth , fill = Site)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  #ylim(5, 6.25) +
  scale_fill_brewer(palette = "Set1") +
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

ggplot(final_k_rna, aes(x=Site, y = value, group = Depth , fill = Site)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge(0.9)) +
  #geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  #ylim(5, 6.25) +
  scale_fill_brewer(palette = "Set1") +
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


# import bins w/ foam genes
#bin gene covs

setwd("/home/jonny/Documents/data/anvi_2015/bin_genes")
# gene PA
gene_presence <- read.csv("bins_to_func_e50_no_methano.csv", header = T)
melted_gene_counts <- melt(gene_presence)
k_gene_counts <- subset(melted_gene_counts, variable =="Bin_2_7_1")
k_gene_counts_final <- subset(k_gene_counts, value > 0)

ggplot(k_gene_counts_final, aes(x =L2, y= value, fill = L1)) +
  geom_bar(stat = "identity", color = "black", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  coord_flip() +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")


# gene expression

bin_gene_covs <- read.csv("final_genes_w_anns_covs", sep="\t" , header = T)

ktedono_genes <- subset(bin_gene_covs, bin =="Bin_2_7_1")
ktedono_genes_rna <- ktedono_genes[51:71]
kteodono_gnees_rna <- cbind(ktedono_genes[1:5], ktedono_genes_rna)

total_activity <- colSums(ktedono_genes_rna)

melt(total_activity)
write.csv(total_activity, "tots_exp")
tots <- read.csv("tots_exp", header = T)

ggplot(tots, aes(x= Depth, y = total.gene.coverage, fill = Depth)) +
  scale_fill_brewer(palette = "Set1") +
  geom_bar(stat ="identity", color = "black") +
  facet_wrap( ~ Site) +
  theme_bw() +
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size=16),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")

write.csv(kteodono_gnees_rna, "k_genes_rna", row.names = FALSE, quote =FALSE)
k_gene_exp <- read.csv("k_genes_rna", header = TRUE)

#remove gene_names not in e50
genes_e50 <- read.csv("genes_e50", sep = " ", header = FALSE)
k_gene_exp_e50 <- subset(k_gene_exp, gene %in% genes_e50$V1)


k_gene_ko <- aggregate(k_gene_exp_e50[6:20], by=list(k_gene_exp_e50$KO), 
                         FUN=median, na.rm=TRUE)
k_gene_ko_cats <- k_gene_exp_e50[3:5]
colnames(k_gene_ko)[1] <- "KO"

k_gene_ko_w_cat <- full_join(k_gene_ko_cats,k_gene_ko, by = "KO")
m_ko_exp <- melt(k_gene_ko_w_cat)

final_ko_exp <- full_join(m_ko_exp, rna_factors, by = "variable")

ferm_ko_exp <- subset(final_ko_exp, L1 == "01_Fermentation")
hyd_ko_exp <- subset(final_ko_exp, L1 == "08_Hydrocarbon degradation")
carb_ko_exp <- subset(final_ko_exp, L1 == "09_Carbohydrate Active enzyme - CAZy")
nit_ko_exp <- subset(final_ko_exp, L1 == "11_Nitrogen cycle")
meth_ko_exp <- subset(final_ko_exp, L1 == "15_Methylotrophy")
sulf_ko_exp <- subset(final_ko_exp, L1 == "18_Sulfur compounds metabolism")

make_grids(ferm_ko_exp)
make_grids(nit_ko_exp)

make_grids2(final_ko_exp)

subset



make_grids = function(x){
  ggplot(x, aes(x=L2, y = value, fill = L2)) +
    geom_bar(stat = "summary", fun.y = "median", colour = "black", position = position_dodge(0.9)) +
    #ylim(5, 6.25) +
    #geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=14),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 14),
          panel.border = element_rect(fill = NA, size=0.5),
          strip.background = element_rect(color = "black", size = 0.5),
          legend.position="none")+
    facet_grid( Site ~ Depth)  
}

make_grids2 = function(x){
  ggplot(x, aes(x=L1, y = value, fill = L1)) +
    geom_bar(stat = "summary", fun.y = "median", colour = "black", position = position_dodge(0.9)) +
    #ylim(5, 6.25) +
    #geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    coord_flip() +
    theme(legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=14),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 14),
          panel.border = element_rect(fill = NA, size=0.5),
          strip.background = element_rect(color = "black", size = 0.5),
          legend.position="none")+
    facet_grid( Site ~ Depth)  
}

make_grids2(final_ko_exp)

ggplot(carb_ko_exp, aes(x=L2, y = value, fill = L2)) +
  geom_bar(stat = "summary", fun.y = "sum", colour = "black", position = position_dodge(0.9)) +
  #ylim(5, 6.25) +
  #geom_errorbar(aes(ymin=value-se, ymax=value+se), colour = "black", position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")+
  facet_grid( Site ~ Depth)

a <- subset(final_ko_exp, Site == "G1" )
b <- subset(a, Depth == "15-20" )
b <- read.csv("b", header = T)

ggplot(b, aes(x =L2, y= value, fill = L1)) +
  geom_bar(stat = "summary", fun.y = "median", color = "black", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  coord_flip() +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="none")





