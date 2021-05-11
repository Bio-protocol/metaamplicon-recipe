# -----------------------------------------------------------------------------#
# Meta-amplicon analysis recipe
# Building and adding a phylogeny to the cleaned phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     vegan v 2.5.6
#                     phangorn v 2.2.5
#                     phyloseq v 1.30.0
#                     msa v 1.18.0
#                     ape v 5.4
#                     seqinr v 3.6.1
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Perform multiple sequence alignment of all ASVs, build distance matrix,       # 
# construct and refine a phylogenetic tree, add the tree to the phyloseq object #
#           With larger data sets, this can be a long process...                #
# Further, proper phylogenetics is beyond the scope of this tutorial.           #
#################################################################################

# Packages and functions ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
library(msa); packageVersion("msa")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")

# Read in phyloseq object from first script output ####
ps <- readRDS("./cache/ps_not-cleaned.RDS")
sample_data(ps)

# simplify ASV names
seqs <- rownames(tax_table(ps))
names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree

# Multiple sequence alignment  ####
alignment <- msa(seqs,method = "Muscle", type = "dna",verbose = TRUE,order = "input",maxiters = 10)

# save progress 
saveRDS(alignment,"./cache/16S_dna_alignment_muscle.RDS")

# Convert to phangorn format
phang.align = as.phyDat(alignment, type = "DNA")

# Model testing
mt <- modelTest(phang.align)
mt$Model

# distance - maximum likelihood ####
dm <- dist.ml(phang.align,model = "JC69")


# save progress
saveRDS(dm,"./cache/16S_ML_Distance.RDS")

# Initial neighbor-joining and upgma trees ####
treeNJ <- NJ(dm) # Note, tip order != sequence order
treeUPGMA <- upgma(dm)

# find parsimony scores for each tree
parsimony(treeNJ, phang.align)
parsimony(treeUPGMA, phang.align)

# search through treespace with NNI and SPR algorithm
optim <- optim.parsimony(treeNJ,phang.align)


# save progress
saveRDS(treeNJ, "./cache/16S_treeNJ.RDS")
saveRDS(optim, "./cache/16S_tree_optim.RDS")



# Maximum likelihood ####
fit = pml(treeNJ, data=phang.align)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")


# edit tip labels
fit$tree$tip.label <- seqs


# save trees
saveRDS(fit,"./cache/16S_fit_treeNJ.RDS")



# add tree to phyloseq object ####
ps2 <- phyloseq(tax_table(tax_table(ps)),
                otu_table(otu_table(ps)),
                sample_data(sample_data(ps)),
                phy_tree(fit$tree))


# Save updated phyloseq object with tree
saveRDS(ps2, "./cache/ps_not-cleaned_w_tree.RDS")


