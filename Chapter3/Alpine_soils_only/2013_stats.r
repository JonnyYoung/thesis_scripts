#################################################################################################################
# stats processing pipeline for 2013 data
#################################################################################################################
setwd("~/Documents/data/2013_data/alps") # <---- set working directy here
#################
# laod libraraies
#################
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggvegan)
library(plyr)
library(caret)
############
# load data
############
##################################################################
# cross biom comparison (using diamond aligned data)
##################################################################
#see other script
################################################################
# alps alone (use kaiju for taxa)
################################################################
# data
a_phy <- read.csv("2013.s45.phy.csv", sep=',', header=TRUE)
a_class <- read.csv("2013.s45.class.csv", sep=',', header=TRUE)
a_ord <- read.csv("2013.s45.ord.csv", sep=',', header=TRUE)
a_fam <- read.csv("2013.s45.fam.csv", sep=',', header=TRUE)
a_gen <- read.csv("2013.s45.gen.csv", sep=',', header=TRUE)
a_spec <- read.csv("2013.s45.spec.csv", sep=',', header=TRUE)
setwd("~/Documents/data/2013_data/alps")
a_seed1 <- read.csv("SEED-Level1.txt", sep='\t', header=TRUE)
names(a_seed1)[1] <- "Function"
a_seed2 <- read.csv("SEED-Level2.txt", sep='\t', header=TRUE)
names(a_seed2)[1] <- "Function"
a_seed3 <- read.csv("SEED-Level3.txt", sep='\t', header=TRUE)
names(a_seed3)[1] <- "Function"

# factors

alps_fact <- read.table("factors.txt", sep = "\t", header = T)
names(alps_fact)[1] <- "Label"

# functions to remove v.low abundance hits

remove_low_taxa = function(x){
  tax_red <- subset(x, Genus!="unclassified" & Genus!="cannot be assigned")
  dna3 <- subset(tax_red, rowSums(tax_red[,2:16])>10)
  tax_red2 <- dna3[,-1]
  rownames(tax_red2) <- dna3[,1]
  as.data.frame(t(tax_red2))
}

remove_low_fun = function(x){
  dna3 <- subset(x, rowSums(x[,2:16])>10)
  tax_red2 <- dna3[,-1]
  rownames(tax_red2) <- dna3[,1]
  as.data.frame(t(tax_red2))
}

a_phy_f <- remove_low_taxa(a_phy)
a_class_f <- remove_low_taxa(a_class)
a_ord_f <- remove_low_taxa(a_ord)
a_fam_f <- remove_low_taxa(a_fam)
a_gen_f <- remove_low_taxa(a_gen)
a_spec_f <- remove_low_taxa(a_spec)
a_seed1_f <- remove_low_fun(a_seed1)
a_seed2_f <- remove_low_fun(a_seed2)
a_seed3_f <- remove_low_fun(a_seed3)

###################
# rarefaction plots
###################

# rarefaction plots # 
rf_plot = function(x){
S <- specnumber(x) # observed number of species
raremax <- min(rowSums(x))
Srare <- rarefy(x, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(x, step = 50000, sample = raremax, col = "blue", cex = 0.6)
}

rf_phy <- rf_plot(a_phy_f)
rf_class <- rf_plot(a_class_f)
rf_ord <- rf_plot(a_ord_f)
rf_fam <- rf_plot(a_fam_f)
rf_gen <- rf_plot(a_gen_f)
rf_spec <- rf_plot(a_spec_f)
rf_s1 <- rf_plot(a_seed1_f)
rf_s2 <- rf_plot(a_seed2_f)
rf_s3 <- rf_plot(a_seed3_f)


#########################
# unsupervised clustering
#########################
library("factoextra")
library("cluster")

# standardize (for clustering)
n_a_phy_f <- decostand(a_phy_f, method = "hellinger")
n_a_class_f <- decostand(a_class_f, method = "hellinger")
n_a_ord_f <- decostand(a_ord_f, method = "hellinger")
n_a_fam_f <- decostand(a_fam_f, method = "hellinger")
n_a_gen_f <- decostand(a_gen_f, method = "hellinger")
n_a_spec_f <- decostand(a_spec_f, method = "hellinger")
n_a_seed1_f <- decostand(a_seed1_f, method = "hellinger")
n_a_seed2_f <- decostand(a_seed2_f, method = "hellinger")
n_a_seed3_f <- decostand(a_seed3_f, method = "hellinger")

# use standardized values
# choose ideal k value using three methods:
# 'sihlouette' 'wss' or 'gap_stat'
# 'gap_stat' peforms better with standardised data

fviz_nbclust(a_seed_f, kmeans, method = "silhouette", nboot = 50)

# make cluster dendrogram and highlight clustering
res.hc <- hclust(n_a_gen_f, method = "ward.D2")
#grp<- cutree(res.hc, k= 3)
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = 2, border = 2:5) # add rectangle)

# diffeent type of den
km.res <- kmeans(test, 2, nstart = 25)
fviz_cluster(km.res, data = test, frame.type = "convex")+
  theme_minimal()

res <- hcut(a_seed2_f, k = 4, stand = TRUE)
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF","#E7B800"))


p <- fviz_nbclust(gen_, kmeans, method = "silhouette") + ggtitle(NULL)
res.hc <- hclust(dist , method = "ward.D2")
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = 8, border = 2:5) # add rectangle)()

gen_dist <- vegdist(n_a_gen_f, method ="euclidean")

#### best dendogram_ fix sample names first
res <- hcut(gen_dist, stand = F, k = 4, hc_method = "average")

vert_clust <- fviz_dend(res, cex = 1.0, lwd = 1.0, 
                        k_colors =  c("#4DAF4A", "#E41A1C", "#377EB8", "#984EA3"),
                        horiz = TRUE,
                        rect = TRUE,
                        lower_rect = -0.2,
                        color_labels_by_k = TRUE,
                        labels_track_height = 0.1,
                        main = NULL) +  theme_void()

k_clus2 <- fviz_cluster(res, data=n_a_gen_f, repel = TRUE, palette = "Set1", show.clust.cent = FALSE, main = NULL) +
  theme_bw() +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size = 12))


ggsave("clus.svg",
       plot = vert_clust,
       device = "svg",
       width = 8,
       height = 6,
       scale = 1,
       dpi = 300)


ggsave("k_clus2.svg",
       plot = k_clus2,
       device = "svg",
       width = 8,
       height = 6,
       scale = 1,
       dpi = 300)


#############################################################
# testing different dissimilarity measures and transformations
#############################################################
#test on species matix + seed vlv2
#1. using standardized
std_euc <- vegdist(n_a_gen_f, method = "euclidean")
std_man <- vegdist(n_a_gen_f, method = "manhattan")
std_gow <- vegdist(n_a_gen_f, method = "gower")
std_altgor <- vegdist(n_a_gen_f, method = "altGower")
std_can <- vegdist(n_a_gen_f, method = "canberra")
std_bray <- vegdist(n_a_gen_f, method = "bray")
std_kul <- vegdist(n_a_gen_f, method = "kulczynski")
std_mor <- vegdist(n_a_gen_f, method = "morisita")
std_horn <- vegdist(n_a_gen_f, method = "binomial")
std_cao <- vegdist(n_a_gen_f, method = "cao")

#2. hellinger
hel_euc <- vegdist(n_a_gen_f, method = "euclidean")
hel_man <- vegdist(n_a_gen_f, method = "manhattan")
hel_gow <- vegdist(n_a_gen_f, method = "gower")
hel_altgor <- vegdist(n_a_gen_f, method = "altGower")
hel_can <- vegdist(n_a_gen_f, method = "canberra")
hel_bray <- vegdist(n_a_gen_f, method = "bray")
hel_kul <- vegdist(n_a_gen_f, method = "kulczynski")
hel_mor <- vegdist(n_a_gen_f, method = "morisits")
hel_horn <- vegdist(n_a_gen_f, method = "binomial")
hel_cao <- vegdist(n_a_gen_f, method = "cao")

#3. log
log_euc <- vegdist(n_a_gen_f, method = "euclidean")
log_man <- vegdist(n_a_gen_f, method = "manhattan")
log_gow <- vegdist(n_a_gen_f, method = "gower")
log_altgor <- vegdist(n_a_gen_f, method = "altGower")
log_can <- vegdist(n_a_gen_f, method = "canberra")
log_bray <- vegdist(n_a_gen_f, method = "bray")
log_kul <- vegdist(n_a_gen_f, method = "kulczynski")
log_mor <- vegdist(n_a_gen_f, method = "morisits")
log_horn <- vegdist(n_a_gen_f, method = "binomial")
log_cao <- vegdist(n_a_gen_f, method = "cao")


# export files


##################
# ordination plots
##################

# cross biom
# all taxonomic levels plus functional levels

# Alps alone
# all taxonomic levels plus functional levels

make_pcoa = function(x){
  helling_Taxa<-decostand(x,method="hellinger")
  Eco_Dist<-vegdist(helling_Taxa,dist="euclidean")
  capscale(Eco_Dist~1)
}

a_pcoa_p <- make_pcoa(a_phy_f)
a_pcoa_c <- make_pcoa(a_class_f)
a_pcoa_o <- make_pcoa(a_ord_f)
a_pcoa_f <- make_pcoa(a_fam_f)
a_pcoa_g <- make_pcoa(a_gen_f)
a_pcoa_s <- make_pcoa(a_spec_f)
a_pcoa_seed1 <- make_pcoa(a_seed1_f)
a_pcoa_seed2 <- make_pcoa(a_seed2_f)
a_pcoa_seed3 <- make_pcoa(a_seed3_f)


a_fp_p <- fortify(a_pcoa_p)
a_fp_c <- fortify(a_pcoa_c)
a_fp_o <- fortify(a_pcoa_o)
a_fp_f <- fortify(a_pcoa_f)
a_fp_g <- fortify(a_pcoa_g)
a_fp_s <- fortify(a_pcoa_s)
a_fp_seed1 <- fortify(a_pcoa_seed1)
a_fp_seed2 <- fortify(a_pcoa_seed2)
a_fp_seed3 <- fortify(a_pcoa_seed3)

##subset and add taxrank descriptor##

fps_p <- subset(a_fp_p, Score =="sites")
fps_c <- subset(a_fp_c, Score =="sites")
fps_o <- subset(a_fp_o, Score =="sites")
fps_f <- subset(a_fp_f, Score =="sites")
fps_g <- subset(a_fp_g, Score =="sites")
fps_s <- subset(a_fp_s, Score =="sites")
fps_seed1 <- subset(a_fp_seed1, Score =="sites")
fps_seed2 <- subset(a_fp_seed2, Score =="sites")
fps_seed3 <- subset(a_fp_seed3, Score =="sites")

fps_p$Rank <- "Phylum" 
fps_c$Rank <- "Class"
fps_o$Rank <- "Order"
fps_f$Rank <- "Family"
fps_f$Dim2 <- (fps_f$Dim2*-1)
fps_g$Rank <- "Genus"
fps_s$Rank <- "Species"
fps_seed1$Rank <- "SEED level 1"
fps_seed2$Rank <- "SEED level 2"
fps_seed3$Rank <- "SEED level 3"
fps_seed3$Dim2 <- (fps_seed3$Dim2*-1)

allm <- rbind(fps_p, fps_c, fps_o, fps_f, fps_g, fps_s, fps_seed1, fps_seed2, fps_seed3)
pca_final <- full_join(allm, alps_fact, by = "Label")
pca_final$Rank_f = factor(pca_final$Rank, levels=c('Phylum','Class','Order','Family', 'Genus', 'Species', 'SEED level 1', 'SEED level 2', 'SEED level 3'))

#find convex hulls for horizons

find_hull <- function(my_data) my_data[chull(my_data[,1], my_data[,2]), ]
h1 <- subset(pca_final, Rank_f == "Phylum")
h2 <- subset(pca_final, Rank_f == "Class")
h3 <- subset(pca_final, Rank_f == "Order")
h4 <- subset(pca_final, Rank_f == "Family")
h5 <- subset(pca_final, Rank_f == "Genus")
h6 <- subset(pca_final, Rank_f == "Species")
h7 <- subset(pca_final, Rank_f == "SEED level 1")
h8 <- subset(pca_final, Rank_f == "SEED level 2")
h9 <- subset(pca_final, Rank_f == "SEED level 3")
ch1 <- ddply(h1, "Horizon", find_hull)
ch2 <- ddply(h2, "Horizon", find_hull)
ch3 <- ddply(h3, "Horizon", find_hull)
ch4 <- ddply(h4, "Horizon", find_hull)
ch5 <- ddply(h5, "Horizon", find_hull)
ch6 <- ddply(h6, "Horizon", find_hull)
ch7 <- ddply(h7, "Horizon", find_hull)
ch8 <- ddply(h8, "Horizon", find_hull)
ch9 <- ddply(h9, "Horizon", find_hull)
hor_hulls <- rbind(ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8,ch9)

#plot just genus pcoa ############################
# fix inverted scales for some
#family and seed level 3



# plot by horizon "using hulls ################################################################
ord_all <- ggplot(pca_final, aes(x = Dim1, y = Dim2, colour= Horizon)) +
  geom_point(size = 2) +
  geom_polygon(data=hor_hulls, aes(x=Dim1, y=Dim2, group=Horizon, fill=Horizon), alpha=0.2) +
  scale_color_brewer(name= "Horizon", 
                     palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 11),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5)) +
  facet_wrap( ~Rank_f, ncol = 3, scales = 'free') +
  xlab("PCoA 1") + ylab("PCoA 2")
  
ggsave("all_ordination.svg",
       plot = ord_all,
       device = "svg",
       width = 8,
       height = 8,
       scale = 1,
       dpi = 300)




# plot by site
ggplot(pca_final, aes(x = Dim1, y = Dim2, colour= Site)) +
  geom_point(size = 2) +
  scale_color_brewer(name= "Site", 
                     palette = "Set1") +
  stat_ellipse(linetype = 2, level = 0.9) +
  facet_wrap( ~Rank_f, ncol = 3) +
  xlab("PCoA 1") + ylab("PCoA 2")


# get % variabtion explained by pcoA axes
extract_importance = function(x){
  a <- (summary(eigenvals(x)))
  b <- t(a$importance)
  b[1:2,]
}

ax_p <- as.data.frame(extract_importance(a_pcoa_p))
ax_c <- as.data.frame(extract_importance(a_pcoa_c))
ax_o <- as.data.frame(extract_importance(a_pcoa_o))
ax_f <- as.data.frame(extract_importance(a_pcoa_f))
ax_g <- as.data.frame(extract_importance(a_pcoa_g))
ax_s <- as.data.frame(extract_importance(a_pcoa_s))
ax_sd1 <- as.data.frame(extract_importance(a_pcoa_seed1))
ax_sd2<- as.data.frame(extract_importance(a_pcoa_seed2))
ax_sd3 <- as.data.frame(extract_importance(a_pcoa_seed3))

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

all_dna_ax <-rbind(ax_p, ax_c, ax_o, ax_f, ax_g, ax_s, ax_sd1, ax_sd2, ax_sd3)

all_ax_melt <- melt(all_dna_ax)
all_ax_melt_prop <-subset(all_ax_melt, variable =="Proportion Explained")
final_ax <- dcast(all_ax_melt_prop, Rank ~ rn, value.var = "value")
write.csv(final_ax, file ="pcoa_ax_scores.csv", quote = FALSE, row.names = FALSE)





# Env variables

##################
# Alpha diversirty
#################

#shannon
adiv_phy_sh <- diversity(a_phy_f, index = "shannon")
adiv_class_sh <- diversity(a_class_f, index = "shannon")
adiv_ord_sh <- diversity(a_ord_f, index = "shannon")
adiv_fam_sh <- diversity(a_fam_f, index = "shannon")
adiv_gen_sh <- diversity(a_gen_f, index = "shannon")
adiv_spec_sh <- diversity(a_spec_f, index = "shannon")
adiv_seed1_sh <- diversity(a_seed1_f, index = "shannon")
adiv_seed2_sh <- diversity(a_seed2_f, index = "shannon")
adiv_seed3_sh <- diversity(a_seed3_f, index = "shannon")


#invsimpson
adiv_phy_sim <- diversity(a_phy_f, index = "invsimpson")
adiv_class_sim <- diversity(a_class_f, index = "invsimpson")
adiv_ord_sim <- diversity(a_ord_f, index = "invsimpson")
adiv_fam_sim <- diversity(a_fam_f, index = "invsimpson")
adiv_gen_sim <- diversity(a_gen_f, index = "invsimpson")
adiv_spec_sim <- diversity(a_spec_f, index = "invsimpson")
adiv_seed1_sim <- diversity(a_seed1_f, index = "invsimpson")
adiv_seed2_sim <- diversity(a_seed2_f, index = "invsimpson")
adiv_seed3_sim <- diversity(a_seed3_f, index = "invsimpson")

# merge dfs and plot

#adiv_all_sh <- rbind(adiv_phy_sim, adiv_class_sim, adiv_ord_sim, adiv_fam_sim, adiv_gen_sim, adiv_spec_sim, adiv_seed1_sim, adiv_seed2_sim, adiv_seed3_sim)
adiv_all_sh <- rbind(adiv_phy_sh, adiv_class_sh, adiv_ord_sh, adiv_fam_sh, adiv_gen_sh, adiv_spec_sh, adiv_seed1_sh, adiv_seed2_sh, adiv_seed3_sh)
adiv_all_sh <- t(adiv_all_sh)
adiv_all_sh <- setDT(as.data.frame(adiv_all_sh), keep.rownames = TRUE)[]
names(adiv_all_sh)[1] <- "Label"
adiv_all_w_fact <- full_join(adiv_all_sh, alps_fact, by = "Label")
colnames(adiv_all_w_fact) <- c("Label",
                               "Phylum: R2=0.51 P=2.7E-03   **",
                               "Class: R2=0.47 P=4.8E-03    **",
                               "Order: R2=0.55 P=1.7E-03    **",
                               "Family: R2=0.34 P=2.2E-02   **",
                               "Genus: R2=0.58 P=9.1E-04    ***",
                               "Species: R2=0.074 P=3.3E-01",
                               "SEED l1: R2=0.081 P=3.1E-01",
                               "SEED l2: R2=0.57 P=1.1E-03    **",
                               "SEED l3: R2=0.022 P=6.0E-01",
                               "Site",
                               "Horizon",
                               "Depth")
adiv_sh_final <- melt(adiv_all_w_fact, id.vars = c("Label", "Site", "Horizon", "Depth"))
adiv_gen_sub <- subset(adiv_sh_final, variable =="adiv_gen_sh")
adiv_seed2_sub <- subset(adiv_sh_final, variable =="adiv_seed2_sh")
final_adiv_plot <- rbind(adiv_gen_sub, adiv_seed2_sub)

adiv_p_sub <- subset(adiv_sh_final, variable =="Phylum")
adiv_c_sub <- subset(adiv_sh_final, variable =="Class")
adiv_o_sub <- subset(adiv_sh_final, variable =="Order")
adiv_f_sub <- subset(adiv_sh_final, variable =="Family")
adiv_g_sub <- subset(adiv_sh_final, variable =="Genus")
adiv_s_sub <- subset(adiv_sh_final, variable =="Species")
adiv_sd1_sub <- subset(adiv_sh_final, variable =="SEED level 1")
adiv_sd2_sub <- subset(adiv_sh_final, variable =="SEED level 2")
adiv_sd3_sub <- subset(adiv_sh_final, variable =="SEED level 3")

lm_eqn = function(df){
  m = lm(Depth ~ value, df);
  eq <- substitute(r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 2)))
  r <- as.character(as.expression(eq))
  f <- summary(m)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  p2 <- format(p, digits = 2, scientific = TRUE) 
  d <- cbind(r, p2)
  d
}

p1 <- lm_eqn(adiv_p_sub)
p2 <-lm_eqn(adiv_c_sub)
p3 <- lm_eqn(adiv_f_sub)
p4 <- lm_eqn(adiv_o_sub)
p5 <- lm_eqn(adiv_g_sub)
p6 <- lm_eqn(adiv_s_sub)
p7 <- lm_eqn(adiv_sd1_sub)
p8 <- lm_eqn(adiv_sd2_sub)
p9 <- lm_eqn(adiv_sd3_sub)

finalP <- rbind(p1, p2, p3, p4, p5, p6, p7, p8, p9)
write.csv(finalP, file ="simp_div_pvals")

write.csv(adiv_sh_final,
          file = "adiv_no_funct",
          quote = F)
plot2 <- read.csv(file = "adiv_no_funct", row.names = 1, header = T)

a <- ggplot(plot2, aes(x= Depth, y= value))+
  geom_point(shape = 1) +
  geom_smooth(method = 'lm') +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 11),
        panel.border = element_rect(fill = NA, size=0.5),
        strip.background = element_rect(color = "black", size = 0.5)) +
  facet_wrap( ~ variable, scales = 'free') +
    xlab("Depth (cm)") + ylab("Shannon diversity")



#########################
# observed beta diversity
#########################

a_phy_fb <- a_phy_f
a_class_fb <- a_class_f
a_ord_fb <- a_ord_f 
a_fam_fb <- a_fam_f 
a_gen_fb <- a_gen_f 
a_spec_fb <- a_spec_f 
a_seed1_fb <- a_seed1_f
a_seed2_fb <- a_seed2_f
a_seed3_fb <- a_seed3_f

bdiv_hor = function(x){ 
  bdiv_phy <-setDT(x, keep.rownames = TRUE)[]
  names(bdiv_phy)[1] <- "Label"
  phy_readya <- left_join(alps_fact, bdiv_phy, by = "Label")
  #add line here to subset by site first
  ah <- subset(phy_readya, Horizon=="Ah")
  bw <-subset(phy_readya,Horizon=="Bw")
  cox <-subset(phy_readya,Horizon=="Cox")
  ah[1:3] <- list(NULL)
  bw[1:3] <- list(NULL)
  cox[1:3] <- list(NULL)
  ###you need to load the beta.div function from anohter script!
  Ob_Beta1<-beta.div(ah,method="hellinger",nperm=100)
  Ob_Beta2<-beta.div(bw,method="hellinger",nperm=100)
  Ob_Beta3<-beta.div(cox,method="hellinger",nperm=100)
  ah_bdiv <- Ob_Beta1$SStotal_BDtotal
  bw_bdiv <- Ob_Beta2$SStotal_BDtotal
  cox_bdiv <- Ob_Beta3$SStotal_BDtotal
  bdiv <- as.data.frame(cbind(ah_bdiv,bw_bdiv,cox_bdiv))
  row.names(bdiv) <- c("ss", "bd")
  bdiv
}
phy_bd <- bdiv_hor(a_phy_fb)
class_bd <- bdiv_hor(a_class_fb)
ord_bd <- bdiv_hor(a_ord_fb)
fam_bd <- bdiv_hor(a_fam_fb)
gen_bd <- bdiv_hor(a_gen_fb)
spec_bd <- bdiv_hor(a_spec_fb)
seed1_bd <- bdiv_hor(a_seed1_fb)
seed2_bd <- bdiv_hor(a_seed2_fb)
seed3_bd <- bdiv_hor(a_seed3_fb)

phy_bd$Rank <- "Phylum"
class_bd$Rank <- "Class"
ord_bd$Rank <- "Order"
fam_bd$Rank <- "Family"
gen_bd$Rank <- "Genus"
spec_bd$Rank <- "Species"
seed1_bd$Rank <- "Seed Level 1"
seed2_bd$Rank <- "Seed Level 2"
seed3_bd$Rank <- "Seed Level 3"

test <- rbind(phy_bd[2,], class_bd[2,], ord_bd[2,], fam_bd[2,], gen_bd[2,], spec_bd[2,], seed1_bd[2,], seed2_bd[2,], seed3_bd[2,])
final_bd_table <- melt(test)

lm

ggplot(final_bd_table, aes(x= variable, y= value)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Rank, scales = 'free') 
  
################
# test for multivariate dispersions
# if significantly different then adonis/ anosim resutls are suspect
# ie significance means one group vaires more than others and differences may be attributable to this
################

setwd("~/Documents/data/2013_data/alps")
alps_fact_anosim <- read.table("factors.txt", sep = "\t", header = T, row.names = 1)

do_bdisp = function(x){
  dist <- make_dist(x)
  with(alps_fact_anosim, betadisper(dist, Horizon, type = "median", bias.adjust = TRUE))
}

do_bdisp_site = function(x){
  dist <- make_dist(x)
  with(alps_fact_anosim, betadisper(dist, Site, type = "median", bias.adjust = TRUE))
}

######by horizon########

bd_p <- do_bdisp(n_a_phy_f)
bd_c <- do_bdisp(n_a_class_f)
bd_o <- do_bdisp(n_a_ord_f)
bd_f <- do_bdisp(n_a_fam_f)
bd_g <- do_bdisp(n_a_gen_f)
bd_s <- do_bdisp(n_a_gen_f)
bd_s1 <- do_bdisp(n_a_seed1_f)
bd_s2 <- do_bdisp(n_a_seed2_f)
bd_s3 <- do_bdisp(n_a_seed3_f)

anova(bd_p)
anova(bd_c)
anova(bd_o)
anova(bd_f)
anova(bd_g)
anova(bd_s)
anova(bd_s1)
anova(bd_s2)
anova(bd_s3)

######by site ########

bd_p_s <- do_bdisp_site(n_a_phy_f)
bd_c_s <- do_bdisp_site(n_a_class_f)
bd_o_s <- do_bdisp_site(n_a_ord_f)
bd_f_s <- do_bdisp_site(n_a_fam_f)
bd_g_s <- do_bdisp_site(n_a_gen_f)
bd_s_s <- do_bdisp_site(n_a_gen_f)
bd_s1_s <- do_bdisp_site(n_a_seed1_f)
bd_s2_s <- do_bdisp_site(n_a_seed2_f)
bd_s3_s <- do_bdisp_site(n_a_seed3_f)

anova(bd_p_s)
anova(bd_c_s)
anova(bd_o_s)
anova(bd_f_s)
anova(bd_g_s)
anova(bd_s_s)
anova(bd_s1_s)
anova(bd_s2_s)
anova(bd_s3_s)


plot(bdis)
mod.HSD <- TukeyHSD(bdis)
plot(mod.HSD)
permutest(bdis, pairwise = TRUE, permutations = 99)
boxplot(bd_g)


bdis_plot <- bdis_boxplot(bd_p_s)

bdis_boxplot <- function(x) {
  a <- cbind((as.data.frame(x$distances)), (as.data.frame((x$group))))
  colnames(a) <- c("Dist", "group")
  a
}

ggplot(bdis_plot, aes(x = group, y = Dist, colour = group)) +
  geom_boxplot() +
  xlab("") +
  ylab("Distance to centroid") +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none")






#########################
# Anosim / adonis
#########################


#########
# ANOSIM
#########

setwd("~/Documents/data/2013_data/alps")
alps_fact_anosim <- read.table("factors.txt", sep = "\t", header = T, row.names = 1)

make_dist = function(x){
  helling_Taxa<-decostand(x,method="hellinger")
  vegdist(helling_Taxa,dist="euclidean")
}

do_anosim = function(x){
  dist <- make_dist(x)
  anosim(dist, alps_fact_anosim$Site, permutations = 9999)
}

anosim_plots <- function(x){
  a <-cbind((as.data.frame(x$dis.rank)),(as.data.frame(x$class.vec)))
  colnames(a) <- c("rank", "cat")
  a
}


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

ggsave("ansoim_site_300dpi.png",
       plot = anosim_ploty,
       device = "png",
       width = 8,
       height = 6,
       scale = 1,
       dpi = 300)


#########
# Adonis
#########

adon_hor = function(x){
  dist <- make_dist(x)
  adonis2(dist ~ Horizon , data=alps_fact_anosim, permutations = 99999)
}

adon_site = function(x){
  dist <- make_dist(x)
  adonis2(dist ~ Site , data=alps_fact_anosim, permutations = 99999)
}

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

ad_phy_hor
ad_class_hor
ad_ord_hor
ad_fam_hor
ad_gen_hor
ad_spec_hor
ad_seed1_hor
ad_seed2_hor
ad_seed3_hor
ad_phy_site
ad_class_site
ad_ord_site
ad_fam_site
ad_gen_site
ad_spec_site
ad_seed1_site
ad_seed2_site
ad_seed3_site


##########################
# SIMPER
##########################

do_simp = function(x){
simp <- with(alps_fact_anosim, simper(x, Horizon, permutations = 9999))
summary(simp, ordered = TRUE)
}

summary

simp_class <- do_simp(n_a_phy_f)
ab <- simp_class$Ah_Bw
ac <- simp_class$Ah_Cox
bc <- simp_class$Bw_Cox
write.csv(ab, file = "phy_a_b_simper", quote = FALSE)
write.csv(ac, file = "phy_a_c_simper", quote = FALSE)
write.csv(bc, file = "phy_b_c_simper", quote = FALSE)


#then get individual tables and sort by p values, take top 10 from each


#######################
# indispecies
#######################
library(indicspecies)
library(stats)

sign_groups <- as.vector(c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3))

phy_signassoc <- signassoc(n_a_phy_f, cluster=sign_groups, mode=1, alternative = "greater", control = how(nperm=9999))
best_phy_alps <- subset(phy_signassoc, psidak <= 0.05)
write.csv(best_class_alps, file = "signassoc_phy", quote = FALSE)

class_signassoc <- signassoc(n_a_class_f, cluster=sign_groups, mode=1, alternative = "greater", control = how(nperm=9999))
best_class_alps <- subset(class_signassoc, psidak <= 0.05)
write.csv(best_class_alps, file = "signassoc_class", quote = FALSE)

gen_signassoc <- signassoc(n_a_gen_f, cluster=sign_groups, mode=1, alternative = "greater", control = how(nperm=9999))
best_gen_alps <- subset(gen_signassoc, psidak <= 0.05)
write.csv(best_gen_alps, file = "signassoc_gen", quote = FALSE)

seed1_signassoc <- signassoc(n_a_seed1_f, cluster=sign_groups, mode=1, alternative = "greater", control = how(nperm=9999))
best_seed1_alps <- subset(seed1_signassoc, psidak <= 0.05)
write.csv(best_seed1_alps, file = "signassoc_seed1", quote = FALSE)


## make boxplots of these abudnances get from sig diff script

#first get names of sig hits for class and seed lvl1

alps_phy_hits <- as.data.frame(rownames(best_phy_alps))
alps_class_hits <- as.data.frame(rownames(best_class_alps))
alps_gen_hits <- as.data.frame(rownames(best_gen_alps))
alps_seed1_hits <- as.data.frame(rownames(best_seed1_alps))
alps_gen_hits <- as.data.frame(rownames(best_gen_alps_sub))
alps_seed1_hits <- as.data.frame(rownames(best_seed1_alps))
colnames(alps_phy_hits) <- "Taxa"
colnames(alps_class_hits) <- "Taxa"
colnames(alps_gen_hits) <- "Taxa"
colnames(alps_seed1_hits) <- "Taxa"


##for plotting selec variables which are alo above a certain % abundance

abund_filter = function(x){
  tmp1 <- decostand(x, method = "total")
  tmp2 <- as.data.frame((colSums(tmp1))/15)
  colnames(tmp2) <- "abund"
  subset(tmp2, abund >= 1.0e-3)
}

abund_filter(n_a_phy_f)

##x = raw raw abundance, y = alps_x_hits
filt_tophits = function(x, y) {
  a  <- abund_filter(x)
  subset(y, Taxa %in% (row.names(a)))  
  }

phyhits_0.01 <- filt_tophits(a_phy_f, alps_phy_hits)
classhits_0.01 <- filt_tophits(a_class_f, alps_class_hits)
genhits_0.01 <- filt_tophits(a_gen_f, alps_gen_hits)
seed1hits_0.01 <- filt_tophits(a_seed1_f, alps_seed1_hits)
seed2hits_0.01 <- filt_tophits(a_seed2_f, alps_seed2_hits)

###########################################
# Sig differences charts (run other script)
###########################################

deco_seed1_f <-decostand(a_seed1_f, method = "total")
deco_seed2_f <-decostand(a_seed2_f, method = "total")
deco_class_f <-decostand(a_class_f, method = "total")

data <- deco_seed2_f
groups <- factor(alps_fact_anosim$Horizon)

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

df2 <- subset(df, Taxa %in% seed2hits_0.01$Taxa)
write.csv(df2, file = "df2.csv", quote= FALSE)
df3 <- read.csv("df2.csv", header = TRUE)
n_pages <- ceiling(length(levels(df3$Taxa))/12)




seed_0.01 <- ggplot(df3,aes(Type,Value,colour=Type)) +
  ylab("Relative abundance") +
  xlab("")+
  geom_boxplot()+geom_jitter(alpha = 0.5)+theme_bw() +
  scale_color_brewer(name = "", palette = "Set1") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none") +
  facet_wrap( ~ Taxa, scales = 'free_y', ncol=3, nrow=3)




alps_seed_plot <- ggplot(df3,aes(Type,Value,colour=Type)) +
  ylab("Relative abundance") +
  xlab("")+
  geom_boxplot()+geom_jitter(alpha = 0.5)+theme_bw() +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 10))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=10)) +
  theme(legend.position="none") +
  facet_wrap( ~ Taxa, scales = 'free_y', ncol=4, nrow=4)


ggsave("seed_a0001_q0.1.svg",
       plot = seed_0.01,
       device = "svg",
       width = 10,
       height = 8,
       dpi = 300)


alps_gen_plot # paginated

library(ggforce)

seed_p1 <- ggplot(df3,aes(Type,Value,colour=Type)) +
  ylab("Relative abundance") +
  xlab("")+
  geom_boxplot()+geom_jitter(alpha = 0.5)+theme_bw() +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none") +
  facet_wrap_paginate( ~ Taxa , ncol=3, 'free_y', nrow=3, page = 2)

seed_p2 <- ggplot(df3,aes(Type,Value,colour=Type)) +
  ylab("Relative abundance") +
  xlab("")+
  geom_boxplot()+geom_jitter(alpha = 0.5)+theme_bw() +
  scale_color_brewer(name = "", palette = "Dark2") +
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 11))+
  theme(axis.title=element_text(size=12, face = "bold")) +
  theme(strip.text = element_text(size=11)) +
  theme(legend.position="none") +
  facet_wrap_paginate( ~ Taxa , scales = 'free_y', ncol=3, nrow=3, page = 2)

##ggsave all##
ggsave("seed1_0.01.svg",
       plot = seed3_0.01,
       device = "svg",
       width = 10,
       height = 8,
       dpi = 300)