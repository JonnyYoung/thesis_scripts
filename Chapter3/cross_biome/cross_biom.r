library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggvegan)
library(plyr)
library(factoextra)
library(data.table)
library(RColorBrewer)
############
# load data
############
##################################################################
# cross biom comparison (using diamond aligned data)
##################################################################
setwd("/home/jonny/Documents/data/2013_data/cross_biom_2/final_cross_biom") #### <--- set working directoy here
# data
phy <- read.csv("phylum.csv", sep=',', header=TRUE)
names(phy)[1] <- "Label"
class <- read.csv("class.csv", sep=',', header=TRUE)
names(class)[1] <- "Label"
ord <- read.csv("order.csv", sep=',', header=TRUE)
names(ord)[1] <- "Label"
fam <- read.csv("family.csv", sep=',', header=TRUE)
names(fam)[1] <- "Label"
gen <- read.csv("genus.csv", sep=',', header=TRUE)
names(gen)[1] <- "Label"
seed1 <- read.csv("seed-lvl1.csv", sep=',', header=TRUE)
names(seed1)[1] <- "Label"
seed2 <- read.csv("seed-lvl2.csv", sep=',', header=TRUE)
names(seed2)[1] <- "Label"
seed3 <- read.csv("seed-lvl3.csv", sep=',', header=TRUE)
names(seed3)[1] <- "Label"
# factors
fact <- read.csv("factors.csv", sep=',', header=TRUE)

#############################################################
# prepare data
############################################################


remove_low = function(x){
  dna3 <- subset(x, rowSums(x[,2:105])>100)
  tax_red2 <- dna3[,-1]
  rownames(tax_red2) <- dna3[,1]
  as.data.frame(t(tax_red2))
}

remove_low2 = function(x){
  dna3 <- x
  tax_red2 <- dna3[,-1]
  rownames(tax_red2) <- dna3[,1]
  as.data.frame(t(tax_red2))
}

phy_f <- remove_low2(phy)
class_f <- remove_low2(class)
ord_f <- remove_low2(ord)
fam_f <- remove_low2(fam)
gen_f <- remove_low2(gen)
seed1_f <- remove_low2(seed1)
seed2_f <- remove_low2(seed2)
seed3_f <- remove_low2(seed3)

fix_names <- read.csv("names.csv", sep=',', header=TRUE)
fix_names_nope <- read.csv("names_nope.csv", sep=',', header=TRUE)

remove_bad = function(x){
  row.names(x) <- fix_names$filename
   a <- subset(x, row.names(x) != "Agricultural Soil 9")
   b <- subset(a, row.names(a) != "Agricultural Soil 10")
   c <- subset(b, row.names(b) != "Agricultural Soil 35")
   d <- subset(c, row.names(c) != "Agricultural Soil 22")
   e <- subset(d, row.names(d) != "Agricultural Soil 22")
   f <- subset(e, row.names(e) != "Tropical forest 1")
   g <- subset(f, row.names(f) != "Tropical forest 2")
   h <- subset(g, row.names(g) != "Agricultural Soil 32")
   i <- subset(h, row.names(h) != "Agricultural Soil 33")
   j <- subset(i, row.names(i) != "Agricultural Soil 41")
   k <- subset(j, row.names(j) != "Agricultural Soil 42")
   k1 <- subset(k, row.names(k) != "Australian Grassland 10")
   k2 <- subset(k1, row.names(k1) != "Australian Grassland 11")
   k3 <- subset(k2, row.names(k2) != "Australian Grassland 12")
   k4 <- subset(k3, row.names(k3) != "Australian Grassland 13")
   k5 <- subset(k4, row.names(k4) != "Australian Grassland 14")
   k6 <- subset(k5, row.names(k5) != "Australian Grassland 15")
   k7 <- subset(k6, row.names(k6) != "Australian Grassland 16")
   k8 <- subset(k7, row.names(k7) != "Australian Grassland 17")
   k9 <- subset(k8, row.names(k8) != "Australian Grassland 18")
   subset(k9, row.names(k9) != "Australian Grassland 19")
 }

phy_clust <- remove_bad(phy_f)
class_clust <- remove_bad(class_f)
ord_clust <- remove_bad(ord_f)
fam_clust <- remove_bad(fam_f)
gen_clust <- remove_bad(gen_f)
gen_clust <- remove_bad(gen_f)
seed1_clust <- remove_bad(seed1_f)
seed2_clust <- remove_bad(seed2_f)
seed3_clust <- remove_bad(seed3_f)
 
##############################################################
# alpha div
#############################################################

#shannon
adiv_phy_sh <- diversity(phy_clust, index = "shannon")
adiv_class_sh <- diversity(class_clust, index = "shannon")
adiv_ord_sh <- diversity(ord_clust, index = "shannon")
adiv_fam_sh <- diversity(fam_clust, index = "shannon")
adiv_gen_sh <- diversity(gen_clust, index = "shannon")
adiv_seed1_sh <- diversity(seed1_clust, index = "shannon")
adiv_seed2_sh <- diversity(seed2_clust, index = "shannon")
adiv_seed3_sh <- diversity(seed3_clust, index = "shannon")

adiv_all_sh <- rbind(adiv_gen_sh, adiv_seed3_sh)
adiv_all_sh <- t(adiv_all_sh)
adiv2 <- setDT(as.data.frame(adiv_all_sh), keep.rownames = TRUE)[]
fact2 <- fact
colnames(fact2) <- c("name", "Group", "rn")
adiv_plot <- full_join(adiv2, fact2, by ="rn")
adiv_plot_f <- subset(adiv_plot, adiv_phy_sh != "NA")
adiv_plot_f_m <- melt(adiv_plot_f)


colourCount = 10
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
scale_color_manual(name = "Soil type", values = getPalette(colourCount))

b <- subset(adiv_plot_f_m, Group == "Alpine Paleosol")

ggplot(b, aes(x = name, y = value, colour = Group)) +
  geom_point() +
  facet_wrap( ~ variable, scales = "free_y")


ggplot(adiv_plot_f_m, aes(x = Group, y = value, colour = Group)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5) +
  scale_color_manual(values = getPalette(colourCount)) +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11),
  panel.border = element_rect(fill = NA, size=0.5),
  strip.background = element_rect(color = "black", size = 0.5)) +
  theme(legend.position="none") +
  facet_wrap( ~ variable, scales = "free_y")


############################################################
# clustering
############################################################

library("factoextra")
library("cluster")

# standardize (for clustering)
n_phy_f <- decostand(phy_clust, method = "hellinger")
phy_dist <- vegdist(n_phy_f, method = "euclidean")
n_class_f <- decostand(class_clust, method = "hellinger")
class_dist <- vegdist(n_class_f, method = "euclidean")
n_ord_f <- decostand(ord_clust, method = "hellinger")
ord_dist <- vegdist(n_ord_f, method = "euclidean")
n_fam_f <- decostand(fam_clust, method = "hellinger")
fam_dist <- vegdist(n_fam_f, method = "euclidean")
n_gen_f <- decostand(gen_clust, method = "hellinger")
gen_dist <- vegdist(n_gen_f, method = "euclidean")
n_seed1_f <- decostand(seed1_clust, method = "hellinger")
seed1_dist <- vegdist(n_seed1_f, method = "euclidean")
n_seed2_f <- decostand(seed2_clust, method = "hellinger")
seed1_dist <- vegdist(n_seed2_f, method = "euclidean")
n_seed3_f <- decostand(seed3_clust, method = "hellinger")
seed3_dist <- vegdist(n_seed3_f, method = "euclidean")


# use standardized values
# choose ideal k value using three methods:
# 'sihlouette' 'wss' or 'gap_stat'
# 'gap_stat' peforms better with standardised data



p <- fviz_nbclust(n_gen_f, kmeans, method = "silhouette", k.max = 15, nboot = 500) + ggtitle(NULL)
res.hc <- hclust(dist , method = "ward.D2")
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = 8, border = 2:5) # add rectangle)()

#### best dendogram_ fix sample names first
res <- hcut(gen_dist, k = 9, stand = TRUE)
vert_clust <- fviz_dend(res, cex = 1.0, lwd = 1.0, k_colors = "Dark2", horiz = TRUE, color_labels_by_k = TRUE, labels_track_height = 20, main = NULL) +
theme_void()

circ_clust <- fviz_dend(res, cex = 1.0, lwd = 1.0, k_colors = "npg", type = "circular", horiz = TRUE, color_labels_by_k = TRUE, labels_track_height = 20, main = NULL) +
  theme_void()

ggsave("vert2",
       plot = vert_clust,
       device = "svg",
       width = 5,
       height = 10,
       scale = 1.5,
       dpi = 300)

##########################################
# ordinations
#########################################

make_pcoa = function(x){
  helling_Taxa<-decostand(x,method="hellinger")
  Eco_Dist<-vegdist(helling_Taxa,dist="euclidean")
  pcoa <- capscale(Eco_Dist~1)
  fort <- fortify(pcoa)
  subset(fort, Score =="sites")
}

fps_p <- make_pcoa(phy_clust)
fps_c <- make_pcoa(class_clust)
fps_f <- make_pcoa(fam_clust)
fps_o <- make_pcoa(ord_clust)
fps_g <- make_pcoa(gen_clust)
fps_seed1 <- make_pcoa(seed1_clust)
fps_seed2 <- make_pcoa(seed2_clust)
fps_seed3 <- make_pcoa(seed3_clust)


fps_p$Rank <- "Phylum" 
fps_c$Rank <- "Class"
fps_o$Rank <- "Order"
fps_f$Rank <- "Family"
fps_g$Rank <- "Genus"
fps_seed1$Rank <- "Seed level 1"
fps_seed2$Rank <- "Seed level 2"
fps_seed3$Rank <- "Seed level 3"

allm <- rbind(fps_p, fps_c, fps_o, fps_f, fps_g, fps_seed1, fps_seed2, fps_seed3)
pca_final <- full_join(allm, fact, by = "Label")
pca_final$Rank_f = factor(pca_final$Rank, levels=c('Phylum','Class','Order','Family', 'Genus', 'Seed level 1', 'Seed level 2', 'Seed level 3'))

#remove weird samples

#filtered_pca <- subset(pca_final, Group !='Remove')
#filtered_pca2 <- subset(filtered_pca, Label != 'X2012_19_Corn_Yes')
#filtered_pca3 <- subset(filtered_pca2, Group !='Tropical forest')
# remove all read 1 samples from german soil
  
sub1 <- subset(pca_final, Rank_f == "Class")
sub2 <- subset(pca_final, Rank_f == "Genus")
sub2$Dim2 <- (sub2$Dim2*-1)
sub3 <- subset(pca_final, Rank_f == "Seed level 1")
sub4 <- subset(pca_final, Rank_f == "Seed level 3")
two <- rbind(sub1, sub2, sub3, sub4)



 #all together
ggplot(pca_final, aes(x = Dim1, y = Dim2, colour=Group)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set1") +
  #stat_ellipse(linetype = 2, level = 0.95) +
  #scale_shape_manual(values=1:nlevels(filtered_pca2$Group)) +
  facet_wrap( ~Rank_f, ncol = 3, scales = 'free') +
  xlab("PCoA 1") + ylab("PCoA 2")
  

#display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
#                   colorblindFriendly=FALSE)

library(RColorBrewer)

colourCount = 10
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
#scale_color_manual(name = "Soil type", values = getPalette(colourCount))


#just class, genus and seed1 --- use this
pcoa <- ggplot(two, aes(x = Dim1, y = Dim2, color = Group, shape = Group)) +
  geom_point(size = 2, stroke = 1.2) +
  scale_color_manual(name = "Soil type", values = getPalette(colourCount)) +
  scale_shape_manual(name = "Soil type", values=1:nlevels(two$Group)) +
    #theme(legend.position = c(0.6, 0.2)) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12)) +
  #stat_ellipse(linetype = 2, level = 0.95) +
  facet_wrap( ~Rank_f, ncol = 2, scales = 'free') +
  theme(strip.text.x = element_text(size = 12, face="bold")) +
  xlab("PCoA 1") + ylab("PCoA 2")

ggsave("ordinations.png",
       plot = pcoa,
       device = "png",
       width = 10,
       height = 7,
       dpi = 300)

#########################################
# get axis variation explained
#########################################

extract_importance = function(x){
  helling_Taxa<-decostand(x,method="hellinger")
  Eco_Dist<-vegdist(helling_Taxa,dist="euclidean")
  pcoa <- capscale(Eco_Dist~1)
  a <- (summary(eigenvals(pcoa)))
  b <- t(a$importance)
  b[1:2,]
}

ax_p <- as.data.frame(extract_importance(phy_clust))
ax_c <- as.data.frame(extract_importance(class_clust))
ax_o <- as.data.frame(extract_importance(ord_clust))
ax_f <- as.data.frame(extract_importance(fam_clust))
ax_g <- as.data.frame(extract_importance(gen_clust))
ax_sd1 <- as.data.frame(extract_importance(seed1_clust))
ax_sd2<- as.data.frame(extract_importance(seed2_clust))
ax_sd3 <- as.data.frame(extract_importance(seed3_clust))

ax_p$Rank <- "Phylum" 
ax_c$Rank <- "Class"
ax_o$Rank <- "Order"
ax_f$Rank <- "Family"
ax_g$Rank <- "Genus"
ax_s$Rank <- "Species"
ax_sd1$Rank <- "SEED level 1"
ax_sd2$Rank <- "SEED level 2"
ax_sd3$Rank <- "SEED level 3"

setDT(ax_p, keep.rownames = TRUE)[]
setDT(ax_c, keep.rownames = TRUE)[]
setDT(ax_o, keep.rownames = TRUE)[]
setDT(ax_f, keep.rownames = TRUE)[]
setDT(ax_g, keep.rownames = TRUE)[]
setDT(ax_s, keep.rownames = TRUE)[]
setDT(ax_sd1, keep.rownames = TRUE)[]
setDT(ax_sd2, keep.rownames = TRUE)[]
setDT(ax_sd3, keep.rownames = TRUE)[]

all_dna_ax <-rbind(ax_p, ax_c, ax_o, ax_f, ax_g, ax_sd1, ax_sd2, ax_sd3)

all_ax_melt <- melt(all_dna_ax)
all_ax_melt_prop <-subset(all_ax_melt, variable =="Proportion Explained")
final_ax <- dcast(all_ax_melt_prop, Rank ~ rn, value.var = "value")
write.csv(final_ax, file ="pcoa_ax_scores.csv", quote = FALSE, row.names = FALSE)



##########################################
#adonis / anosim / bdisper
########################################

remove_bad_fact = function(x){
  a <- subset(x, Label != "Agricultural Soil 9")
  b <- subset(a, Label != "Agricultural Soil 10")
  c <- subset(b, Label != "Agricultural Soil 35")
  d <- subset(c, Label != "Agricultural Soil 22")
  e <- subset(d, Label != "Agricultural Soil 22")
  f <- subset(e, Label != "Tropical forest 1")
  g <- subset(f, Label != "Tropical forest 2")
  h <- subset(g, Label != "Agricultural Soil 32")
  i <- subset(h, Label != "Agricultural Soil 33")
  j <- subset(i, Label != "Agricultural Soil 41")
  k <- subset(j, Label != "Agricultural Soil 42")
  k1 <- subset(k, Label != "Australian Grassland 10")
  k2 <- subset(k1, Label != "Australian Grassland 11")
  k3 <- subset(k2, Label != "Australian Grassland 12")
  k4 <- subset(k3, Label != "Australian Grassland 13")
  k5 <- subset(k4, Label != "Australian Grassland 14")
  k6 <- subset(k5, Label != "Australian Grassland 15")
  k7 <- subset(k6, Label != "Australian Grassland 16")
  k8 <- subset(k7, Label != "Australian Grassland 17")
  k9 <- subset(k8, Label != "Australian Grassland 18")
  k10 <- subset(k9, Label != "Australian Grassland 19")
  k11 <- subset(k10, Label != "Temperate grassland")
  k12 <- subset(k11, Label != "Boreal forest")
  k13 <- subset(k12, Label != "Temperate forest")
  subset(k13, Label != "Arctic tundra")
} 

remove_single_samp = function(x){
  a <- subset(x, row.names(x) != "Arctic tundra")
  b <- subset(a, row.names(a) != "Boreal forest")
  subset(b, row.names(b) != "Temperate grassland")
  }


red_fact <- remove_bad_fact(fact)
#row.names(red_fact) <- red_fact$Label
#levels(red_fact$Group)
#levels(red_fact$Group) <-(droplevels(red_fact$Group))
# just write red_fact to csv and reload
#write.csv(red_fact, file = "red_fact", quote =FALSE)
anosim_fact <- read.csv("red_fact", header = TRUE, row.names = 1)
row.names(anosim_fact) <- anosim_fact$Label

anosim_fact_nope <- read.csv("red_fact_no_polar", header = TRUE, row.names = 1)

bdis_prep = function(x){ 
  a <- remove_single_samp(x)
  b <- vegdist(a, method = "euclidean")
  with(anosim_fact, betadisper(b, Group, type = "centroid"))
  }


bdis_class <- bdis_prep(n_class_f)
bdis_gen <- bdis_prep(n_gen_f)
bdis_seed1 <- bdis_prep(n_seed1_f)
bdis_seed3 <- bdis_prep(n_seed3_f)

anova(bdis_class)
anova(bdis_gen)
anova(bdis_seed1)
anova(bdis_seed3)

bdis_boxplot <- function(x) {
  a <- cbind((as.data.frame(x$distances)), (as.data.frame((x$group))))
  colnames(a) <- c("Dist", "group")
  a
}

bdis_plot1 <- bdis_boxplot(bdis_class)
bdis_plot1$rank = "Class"
bdis_plot2 <- bdis_boxplot(bdis_gen)
bdis_plot2$rank = "Genus"
bdis_plot3 <- bdis_boxplot(bdis_seed1)
bdis_plot3$rank = "Seed level 1"
bdis_plot4 <- bdis_boxplot(bdis_seed3)
bdis_plot4$rank = "Seed level 3"

bdis_all <- rbind(bdis_plot1, bdis_plot2, bdis_plot3, bdis_plot4)

bids_plots <- ggplot(bdis_all, aes(x = group, y = Dist, colour = group)) +
  geom_boxplot() +
  xlab("") +
  ylab("Distance to centroid") +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none") +
  facet_wrap( ~ rank, scales = 'free')


ggsave("bdis_same_scale.png",
       plot = bids_plots,
       device = "png",
       width = 8,
       height = 8,
       dpi = 300)

##functions for simper, adonis and bdisper

##anosim###
do_anosim = function(x){ 
  a <- remove_single_samp(x)
  b <- vegdist(a, method = "euclidean")
  with(anosim_fact, anosim(b, Group, permutations = 9999,distance = "euclidean"))
}


class_ano <- do_anosim(n_class_f)
gen_ano <- do_anosim(n_gen_f)
seed1_ano <- do_anosim(n_seed1_f)
seed3_ano <- do_anosim(n_seed3_f)
class_ano
gen_ano
seed1_ano
seed3_ano
summary(class_ano)
summary(gen_ano)
summary(seed1_ano)

anosim_plots <- function(x){
  a <-cbind((as.data.frame(x$dis.rank)),(as.data.frame(x$class.vec)))
  colnames(a) <- c("rank", "cat")
  a
}

class_ano_p <- anosim_plots(class_ano)
gen_ano_p <- anosim_plots(gen_ano)
seed1_ano_p <- anosim_plots(seed1_ano)
seed3_ano_p <- anosim_plots(seed3_ano)

class_ano_p$class <- "Class"
gen_ano_p$class <- "Genus"
seed1_ano_p$class <- "Seed level 1"
seed3_ano_p$class <- "seed level 3"

ano_merged <- rbind(class_ano_p, gen_ano_p, seed1_ano_p, seed3_ano_p)


anosim_plots <- ggplot(ano_merged, aes(x = cat, y = rank, colour = cat)) +
  geom_boxplot() +
  xlab("") +
  ylab("Ranked dissimilarity") +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none") +
  facet_wrap( ~ class)

ggsave("anosim.png",
       plot = anosim_plots,
       device = "png",
       width = 8,
       height = 8,
       dpi = 300)

##print results to file

adon_class <- adonis(reduced_class ~ Group, data = anosim_fact, method = "euclidean", permutations = 9999)
#As to the strength of this effect? ~35% of the sums of squares can be explained by `location
adon_gen <- adonis(reduced_gen ~ Group, data = anosim_fact, method = "euclidean", permutations = 9999)
adon_seed1 <- adonis(reduced_seed1 ~ Group, data = anosim_fact, method = "euclidean", permutations = 9999)
adon_seed3 <- adonis(reduced_seed3 ~ Group, data = anosim_fact, method = "euclidean", permutations = 9999)


#######################################
#indicator species
#######################################

library(indicspecies)
library(stats)
reduced_class <- remove_single_samp(n_class_f)
reduced_gen <- remove_single_samp(n_gen_f)
reduced_seed1 <- remove_single_samp(n_seed1_f)
reduced_seed2 <- remove_single_samp(n_seed2_f)
reduced_seed3 <- remove_single_samp(n_seed3_f)

sign_groups <- as.vector(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,4,4,4,4,4,4,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7))

## Look for species whose abundance is significantly higher in one group
#other groups individually
class_signassoc <- signassoc(a, cluster=sign_groups, mode=1, control = how(nperm=9999))
## Look for species whose abundance is significantly higher in sites belonging 
## to one group as opposed to sites not belonging to it.
#all other groups
#signassoc(a, cluster=sign_groups, mode=0, control = how(nperm=999))

class_signassoc <- signassoc(reduced_class, cluster=sign_groups, mode=1, control = how(nperm=9999))
best_class_alps <- subset(class_signassoc, best == 1 & psidak <= 0.05)
#write.csv(best_class_alps, file = "best_class_alps.csv", quote = FALSE)
gen_signassoc <- signassoc(reduced_gen, cluster=sign_groups, mode=1, control = how(nperm=9999))
best_gen_alps <- subset(gen_signassoc, best == 1 $ psidak <= 0.01)
#write.csv(best_gen_alps, file = "best_gen_alps.csv", quote = FALSE)
seed1_signassoc <- signassoc(reduced_seed1, cluster=sign_groups, mode=1, control = how(nperm=9999))
best_seed1_alps <- subset(seed1_signassoc, best == 1 & psidak <=0.05)
#write.csv(best_seed1_alps, file = "best_seed1_alps.csv", quote = FALSE)
seed2_signassoc <- signassoc(reduced_seed2, cluster=sign_groups, mode=1, control = how(nperm=9999))
best_seed2_alps <- subset(seed2_signassoc, best == 1 & psidak <= 0.001)
#write.csv(best_seed2_alps, file = "best_seed2_alps.csv", quote = FALSE)
#seed3_signassoc <- signassoc(reduced_seed3, cluster=sign_groups, mode=1, control = how(nperm=9999))
#best_seed3_alps <- subset(seed3_signassoc, best == 1 & psidak <= 0.001)
#write.csv(best_seed3_alps, file = "best_seed2_alps.csv", quote = FALSE)



## make boxplots of these abudnances get from sig diff script

#first get names of sig hits for class and seed lvl1


alps_class_hits <- as.data.frame(rownames(best_class_alps))
best_gen_alps_sub <- subset(best_gen_alps, psidak <= 0.001)
alps_gen_hits <- as.data.frame(rownames(best_gen_alps_sub))
alps_seed_hits <- as.data.frame(rownames(best_seed_alps))
colnames(alps_class_hits) <- "Taxa"
colnames(alps_gen_hits) <- "Taxa"
colnames(alps_seed_hits) <- "Taxa"


a <- remove_single_samp(seed1_clust)
b <- decostand(a, method = "total")
data <- b
groups <- factor(anosim_fact$Group)

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
kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
pdf("KW_correction.pdf")
plot(kruskal.wallis.table$p.value,
     kruskal.wallis.table$E.value,
     main='Multitesting corrections',
     xlab='Nominal p-value',
     ylab='Multitesting-corrected statistics',
     log='xy',
     col='blue',
     panel.first=grid(col='#BBBBBB',lty='solid'))
lines(kruskal.wallis.table$p.value,
      kruskal.wallis.table$FWER,
      pch=20,col='darkgreen', type='p'
)
lines(kruskal.wallis.table$p.value,
      kruskal.wallis.table$q.value,
      pch='+',col='darkred', type='p'
)
abline(h=kruskal.wallis.alpha, col='red', lwd=2)
legend('topleft', legend=c('E-value', 'p-value', 'q-value'), col=c('blue', 'darkgreen','darkred'), lwd=2,bg='white',bty='o')
dev.off()

last.significant.element <- max(which(kruskal.wallis.table$q.value <= kruskal.wallis.alpha))
selected <- 1:last.significant.element
diff.cat.factor <- kruskal.wallis.table$id[selected]
diff.cat <- as.vector(diff.cat.factor)

print(kruskal.wallis.table[selected,])

#" q = ",round(kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"],5),sep=""
#Now we plot taxa significantly different between the categories
df<-NULL
for(i in diff.cat){
  tmp<-data.frame(data[,i],groups,rep(paste(i),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Type","Taxa")

df2 <- subset(df, Taxa %in% alps_seed_hits$Taxa)

write.csv(df2, file = "df2.csv", quote= FALSE)
df3 <- read.csv("df2.csv", header = TRUE)
n_pages <- ceiling(length(levels(df3$Taxa))/12)




alps_class_plot <- ggplot(df3,aes(Type,Value,colour=Type)) +
  ylab("Relative abundance (%)") +
  xlab("")+
  geom_boxplot()+geom_jitter(alpha = 0.5)+theme_bw() +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none") +
  facet_wrap( ~ Taxa, scales = 'free_y', ncol=3, nrow=3)


alps_seed_plot <- ggplot(df3,aes(Type,Value,colour=Type)) +
  ylab("Relative abundance (%)") +
  xlab("")+
  geom_boxplot()+geom_jitter(alpha = 0.5)+theme_bw() +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 10))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=10)) +
  theme(legend.position="none") +
  facet_wrap( ~ Taxa, scales = 'free_y', ncol=3, nrow=4)


ggsave("alps_seed_abund.svg",
       plot = alps_seed_plot,
       device = "svg",
       width = 10,
       height = 8,
       dpi = 300)


alps_gen_plot # paginated

gen_p1 <- ggplot(df3,aes(Type,Value,colour=Type)) +
  ylab("Relative abundance") +
  xlab("")+
  geom_boxplot()+geom_jitter(alpha = 0.5)+theme_bw() +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none") +
  facet_wrap_paginate( ~ Taxa , scales = 'free_y', ncol=4, nrow=3, page = 1)

gen_p2 <- ggplot(df3,aes(Type,Value,colour=Type)) +
  ylab("Relative abundance") +
  xlab("")+
  geom_boxplot()+geom_jitter(alpha = 0.5)+theme_bw() +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none") +
  facet_wrap_paginate( ~ Taxa , scales = 'free_y', ncol=4, nrow=3, page = 2)

##ggsave all##