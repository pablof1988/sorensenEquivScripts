library(goSorensen)
library(equivStandardTest)
source("adjSignifPvals.R")
data(humanEntrezIDs)


# package goSorensen authomatically charges object "allOncoGeneLists" 
allOncoGeneLists
sapply(allOncoGeneLists, length)

# Formerly, the original set of gene lists was reduced to those with almost 100 genes:
# allOncoGeneLists <- allOncoGeneLists[sapply(allOncoGeneLists, length) >= 100]
# sapply(allOncoGeneLists, length)

# The generic function "equivTestSorensen" in package "goSorensen" implements the equivalence test based
# on the Sorensen-Dice distance.
# There are many methods of this function, for different classes of objects passed as arguments.
# Providing two gene lists (essentially, to "character" vectors of gene identifiers) as arguments, it returns
# an object of class "equivSDhtest" inheriting from "htest".
# Providing an object of class "list" of length k with each of its elements representing a gene list (i.e., a "list"
# of k "character" vectors), all possible k*(k - 1)/2 pairwise tests are performed and an object of class "equivSDhtestList"
# (in fact, a "list") of "htest" objects is generated.
# There are also methods for data already summarized in form of 2x2 contingency tables of joint enrichment.

# Examples:
# Equivalence test between gene lists 'waldman' and 'atlas', in dataset 'cancerGeneLists', 
# at level 4 of the BP ontology:
waldman_atlas.BP.4 <- equivTestSorensen(allOncoGeneLists[["waldman"]], allOncoGeneLists[["atlas"]], 
                                        geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                                        onto = "BP", GOLevel = 4, listNames = c("waldman", "atlas"))
getTable(waldman_atlas.BP.4)
getPvalue(waldman_atlas.BP.4)
getSE(waldman_atlas.BP.4)

# Bootstrap approach:
boot.waldman_atlas.BP.4 <- equivTestSorensen(allOncoGeneLists[["waldman"]], allOncoGeneLists[["atlas"]], 
                                             boot = TRUE,
                                             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                                             onto = "BP", GOLevel = 4, listNames = c("waldman", "atlas"))
getTable(boot.waldman_atlas.BP.4)
getPvalue(boot.waldman_atlas.BP.4)
getSE(boot.waldman_atlas.BP.4)

# A faster way to obtain boot.waldman_atlas.BP.4:
boot.waldman_atlas.BP.4 <- upgrade(waldman_atlas.BP.4, boot = TRUE)
getPvalue(boot.waldman_atlas.BP.4)

waldman_atlas.BP.4_strict <- upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
waldman_atlas.BP.4_strict

class(waldman_atlas.BP.4_strict)
upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25), conf.level = 0.99)

# All pairwise equivalence tests at level 4 of the BP ontology (quite time consuming, you can jump
# this sentence and use the dataset 'BP.4' which is directly charged with package 'goSorensen')
BP.4 <- equivTestSorensen(allOncoGeneLists,
                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                          onto = "BP", GOLevel = 4)

getPvalue(BP.4)
getPvalue(BP.4, simplify = FALSE)

getDissimilarity(BP.4)
getDissimilarity(BP.4, simplify = FALSE)

getUpper(BP.4)
getUpper(BP.4, simplify = FALSE)

getSE(BP.4)
getSE(BP.4, simplify = FALSE)

getTable(BP.4)

BP.4_strict <- upgrade(BP.4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
BP.4_strict
getPvalue(BP.4_strict)

# (Very time consuming. Alternatively, you may use the dataset 'cancerEquivSorensen' directly,
# it is automatically charged with the package 'goSorensen'),
# By default, the tests are iterated over all GO ontologies and for levels 3 to 10:
cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists, 
                                            geneUniverse = humanEntrezIDs, 
                                            orgPackg = "org.Hs.eg.db")
set.seed(123)
boot.cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists,
                                                 boot = TRUE,
                                                 geneUniverse = humanEntrezIDs, 
                                                 orgPackg = "org.Hs.eg.db")
# Or, much more faster:
set.seed(123)
boot.cancerEquivSorensen <- upgrade(cancerEquivSorensen, boot = TRUE)

# # (Also very time consuming.) It is not required to iterate over all ontologies,
# # "allEquivTestSorensen" iterates these procedures over the specified GO ontologies and levels, 
# # e.g., to iterate only over the GO ontologies MF and BP and levels 5 and 6:
# allEquivTestSorensen(allOncoGeneLists, ontos = c("MF", "BP"), GOLevels = 5:6,
#                      geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# # With a more strict equivalence limit d0:
# allEquivTestSorensen(allOncoGeneLists, ontos = c("MF", "BP"), GOLevels = 5:6, d0 = 1 / (1 + 2 * (10/9)), #d0 =0.3103
#                      geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

cancerEquivSorensen

# From 21 possible p-values (21 = 7 * (7 - 1) / 2) and after the Holm's adjustment for testing multiplicity, 
# identify those who are <= 0.05? (Excluding NA values)
# BUT JUMP TO LINE 210 FOR A MORE INFORMATIVE OUTPUT, INCLUDING THE 2x2 CONTINGENCY TABLES OF ENRICHMENT
# ---------------------------------------------------------------------------------------------------------
# For ontology BP:
# Under d0 = 0.4444
sapply(getPvalue(cancerEquivSorensen, onto = "BP")[["BP"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# Under a more restrictive d0 = 0.2857
sapply(getPvalue(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)), onto = "BP")[["BP"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})
# ************** Highly stable results along all GO levels for BP ontology **************************
# ************** a set of equivalent lists seems to be identified          **************************

# ---------------------------------------------------------------------------------------------------------
# For ontology CC:
# Under d0 = 0.4444
sapply(getPvalue(cancerEquivSorensen, onto = "CC")[["CC"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# Under a more restrictive d0 = 0.2857
sapply(getPvalue(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)), onto = "CC")[["CC"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})
# **** Stable but possibly less interesting (less similarity between lists) results for CC ontology ****

# ---------------------------------------------------------------------------------------------------------
# For ontology MF:
# Under d0 = 0.4444
sapply(getPvalue(cancerEquivSorensen, onto = "MF")[["MF"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# Under a more restrictive d0 = 0.2857
sapply(getPvalue(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)), onto = "MF")[["MF"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# **** Stable but possibly less interesting (less similarity between lists) results for MF ontology ****
# ---------------------------------------------------------------------------------------------------------

# 2x2 contingecy tables of joint enrichment:
getTable(cancerEquivSorensen)
getTable(cancerEquivSorensen, GOLevel = "level 6")
getTable(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getTable(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getTable(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getTable(cancerEquivSorensen$BP$`level 4`)

# p-values:
getPvalue(cancerEquivSorensen)
getPvalue(cancerEquivSorensen, simplify = FALSE)
getPvalue(cancerEquivSorensen, GOLevel = "level 6")
getPvalue(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
getPvalue(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getPvalue(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getPvalue(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
getPvalue(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getPvalue(cancerEquivSorensen$BP$`level 4`)

# Sorensen-Dice dissimilarity:
getDissimilarity(cancerEquivSorensen)
getDissimilarity(cancerEquivSorensen, simplify = FALSE)
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6")
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getDissimilarity(cancerEquivSorensen$BP$`level 4`)

# Upper confidence limits:
getUpper(cancerEquivSorensen)
getUpper(cancerEquivSorensen, simplify = FALSE)
getUpper(cancerEquivSorensen, GOLevel = "level 6")
getUpper(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
getUpper(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getUpper(cancerEquivSorensen$BP$`level 4`)

# Standard error of the Sorensen-Dice dissimilarity estimate:
getSE(cancerEquivSorensen)
getSE(cancerEquivSorensen, simplify = FALSE)
getSE(cancerEquivSorensen, GOLevel = "level 6")
getSE(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
getSE(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getSE(cancerEquivSorensen$BP$`level 4`)

cancerEquivSorensen2 <- upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
cancerEquivSorensen2

# CAUTION! Some of these "significant" results may have a very low reliability if the joint enrichment
# frequencies are extremely low.
# This function returns the adjusted significant p-values jointly with the enrichment contingency tables,
# in order to put these values in an adequate context (e.g., those who are not very credible):
signifPvals_d0_0.4444 <- adjSignifPvals(cancerEquivSorensen)

signifPvals_d0_0.4444$BP
signifPvals_d0_0.4444$CC
signifPvals_d0_0.4444$MF

# For a more restrictive d0 = 0.2857:
signifPvals_d0_0.2857 <- adjSignifPvals(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)))

signifPvals_d0_0.2857$BP
signifPvals_d0_0.2857$CC
signifPvals_d0_0.2857$MF

# *********************************************************************************
# Bootstrap approach, tends to be conservative under low enrichment frequencies so the positive
# results seem more reliable (the number of valid bootstrap replicates, over 10000, is also displayed):
boot.signifPvals_d0_0.4444 <- adjSignifPvals(boot.cancerEquivSorensen)

boot.signifPvals_d0_0.4444$BP
boot.signifPvals_d0_0.4444$CC
boot.signifPvals_d0_0.4444$MF

# For a more restrictive d0 = 0.2857:
set.seed(123)
boot.signifPvals_d0_0.2857 <- adjSignifPvals(upgrade(boot.cancerEquivSorensen, d0 = 1/(1 + 2*1.25), boot = TRUE))

boot.signifPvals_d0_0.2857$BP
boot.signifPvals_d0_0.2857$CC
boot.signifPvals_d0_0.2857$MF

