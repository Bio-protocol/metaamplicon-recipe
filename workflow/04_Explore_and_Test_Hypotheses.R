# -----------------------------------------------------------------------------#
# Meta-amplicon analysis recipe
# Exploring cleaned data using phyloseq and corncob packages
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     Biostrings v 2.54.0
#                     corncob v 0.1.0
#                     vegan v 2.6.0
#                     patchwork v 1.0.1
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(corncob); packageVersion("corncob")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")

source("./lib/bbdml_helper.R")

#################################################################################
#                               Main workflow                                   #
#  Explore alpha and beta diversity, visualize data set, test hypotheses,       #
#  search for differentially abundant taxa                                      #
#                                                                               #
#################################################################################


# Load cleaned phyloseq object ####

ps <- readRDS("./cache/clean_phyloseq_object.RDS")


# Getting to know your phyloseq data ####

# number of taxa
ntaxa(ps)

# number of samples
nsamples(ps)

# sample names
sample_names(ps)
rank_names(ps)
# taxa names
taxa_names(ps)

# overview of taxonomy
tax_table(ps)
tax_table(ps)[,1] %>% table()
tax_table(ps)[,2] %>% table()
tax_table(ps)[,3] %>% table()
tax_table(ps)[,4] %>% table()
tax_table(ps)[,5] %>% table()
tax_table(ps)[,6] %>% table()

# sample metadata
sample_data(ps) %>% as_tibble()

# ASV table
otu_table(ps)

# how many sequences observed in each sample?
otu_table(ps) %>% rowSums()

# how many times was each taxon observed?
otu_table(ps) %>% colSums()

# how many different samples was each taxon found in?
asv <- otu_table(ps) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame
asv[asv>0] <- 1 # convert to presence/absence
colSums(asv) # sum of presence (present = 1; absent = 0)

# what was the most widespread taxon (not abundance)
widespread_seq <- names(asv)[which(colSums(asv) == max(colSums(asv)))] # this gives long sequence
tax_table(ps)[widespread_seq,] # this pull that row from ASV table

# what was most abundant (raw counts) taxon?
abund_seq <- which(otu_table(ps) %>% colSums() == max(otu_table(ps) %>% colSums()))
tax_table(ps)[abund_seq,]
otu_table(ps)[,abund_seq]

# access the phylogenetic tree
phy_tree(ps)
plot_tree(ps,color="Phylum")
# this tree needs some refining, but it's will work for this tutorial



# Alpha diversity metrics ####

# since we have no singletons (consequence of DADA2), we cannot use Chao1 estimate!
# plot alpha diversity for every sample
plot_richness(ps, 
              measures = c("Observed","Shannon","Simpson"), 
              color = "Status", 
              sortby = "Observed") +
  theme(axis.text.x = element_blank())

# plot, grouped by colony color with added boxplot
plot_richness(ps, 
              x = "Status",
              measures = c("Observed","Shannon","Simpson"), 
              color = "Status", 
              sortby = "Observed") +
  geom_boxplot(alpha = .5) +
  theme_minimal()
ggsave("./graphs/alpha_diversity_boxplot.png",dpi=300,height = 4,width = 6)

# transform raw counts to relative abundance ####
ps_ra <- transform_sample_counts(ps, fun = function(x){x/sum(x)})

# Beta-diversity ####

# Ordination
dca <- ordinate(ps_ra)
plot_ordination(ps_ra,dca,color = "Status") 

(
  ord1 <- plot_ordination(ps_ra,dca,color = "Status",shape="Island") +
  geom_point(size=4)  + theme_minimal() +
    theme(legend.position = "top") +
    labs(title = "DCA - Bray")
)

# try another ordination method
nmds <- ordinate(ps_ra,method = "NMDS")

(
ord2 <- plot_ordination(ps_ra,nmds,color = "Status",shape="Island") +
  geom_point(size=4)  + theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "NMDS - Bray")
)


# also try with unifrac distance, which takes phylogeny into account
unifrac.dist <- UniFrac(ps_ra)
unifrac <- ordinate(ps_ra,method = "NMDS",distance = unifrac.dist)

(
ord3 <- plot_ordination(ps_ra,unifrac,color = "Status",shape="Island") +
  geom_point(size=4) + theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "NMDS - Unifrac")
)

# combine all plots into one figure for comparison
ord1 / ord2 / ord3
ggsave("./graphs/ordinations.png",dpi=500,width = 6,height = 8)

# permanova ####
# pull out components
asv <- otu_table(ps_ra) %>% as("matrix") %>% as.data.frame()
meta <- sample_data(ps_ra) %>% as.data.frame()

# run permanova model with colony_color and Island as predictors (with interaction term included)
permanova.bray <- vegan::adonis(asv ~ meta$Status * meta$Island,method = "bray")
permanova.bray

# save output to a file
sink("./cache/permANOVA_Table.txt")
permanova.bray
sink(NULL)

# try with jaccard distance as well
permanova.jaccard <- vegan::adonis(asv ~ meta$Status * meta$Island,method = "jaccard")
permanova.jaccard

# Island appears to be a significant indicator of community structure in either case...

# Differential abundance/dispersion tests ####

# Here, we use the corncob package to test for taxa that have differential abundance/dispersion between groups


# use non-transformed data!
set.seed(123)
da_analysis_island <- differentialTest(formula = ~ Island, #abundance
                                             phi.formula = ~ 1, #dispersion
                                             formula_null = ~ 1, #mean
                                             phi.formula_null = ~ 1,
                                             test = "Wald", boot = FALSE,
                                             data = ps,
                                             fdr_cutoff = 0.05,
                                             full_output = TRUE)
plot(da_analysis_island)


# find the significant taxa
da_analysis_island$significant_taxa
da_analysis_island$significant_taxa %>% otu_to_taxonomy(data=ps)

# This is a helper function for plotting corncob results. It's found in "scripts/bbdml_helper.R" 
bbdml_obj <- multi_bbdml(da_analysis_island,
                         ps_object = ps,
                         mu_predictor = "Island",
                         phi_predictor = "Island",
                         taxlevels = 6)

# another helper function found in the same file
plot_multi_bbdml(bbdml_obj,
                 color="Island", 
                 pointsize = 3)


# This saves a plot for each significant taxon in your environment called bbdml_plot_N (1 through number of sig. taxa)

# view all the plots together using patchwork package
(bbdml_plot_1 / bbdml_plot_2 / bbdml_plot_3)
ggsave("./graphs/bbdml_plot_all_sig_taxa.png",height = 6, width = 8)

