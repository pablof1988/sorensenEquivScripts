library(goSorensen)
source("..\\adjSignifPvals.R")
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

# 2x2 contingency table of jointly enriched BP GO terms at level 4 between gene lists waldman and atlas:
getTable(waldman_atlas.BP.4)
# Equivalence test p-value (normal approximation):
getPvalue(waldman_atlas.BP.4)
# Sorensen-Dice dissimilarity standard error:
getSE(waldman_atlas.BP.4)

# Bootstrap approach:
boot.waldman_atlas.BP.4 <- equivTestSorensen(allOncoGeneLists[["waldman"]], allOncoGeneLists[["atlas"]], 
                                             boot = TRUE,
                                             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                                             onto = "BP", GOLevel = 4, listNames = c("waldman", "atlas"))
getTable(boot.waldman_atlas.BP.4)
getPvalue(boot.waldman_atlas.BP.4)
getSE(boot.waldman_atlas.BP.4)

# Building the GO terms mutual enrichment contingency table is the slowest step.
# A faster way to obtain boot.waldman_atlas.BP.4 is to upgrade waldman_atlas.BP.4 to the
# bootstrap approach (but without unnecessarily building again the contingency table):
boot.waldman_atlas.BP.4 <- upgrade(waldman_atlas.BP.4, boot = TRUE)
getPvalue(boot.waldman_atlas.BP.4)

# Similarly, the test can be performed directly from the contingency table:
boot.waldman_atlas.BP.4 <- equivTestSorensen(getTable(waldman_atlas.BP.4), boot = TRUE)
getTable(boot.waldman_atlas.BP.4)
getPvalue(boot.waldman_atlas.BP.4)
getSE(boot.waldman_atlas.BP.4)

# Package goSorensen was designed under the object-oriented programming paradigm. All operations,
# like 'equivTestSorensen', have method functions to perform the same operation from adequate classes
# of objects, like performing the equivalence test from two gene lists, from scratch, or to perform
# the equivalence test from the enrichment contingency table associated to them.

# Using an stricter equivalence limit:
waldman_atlas.BP.4_strict <- upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
waldman_atlas.BP.4_strict

class(waldman_atlas.BP.4_strict)
upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25), conf.level = 0.99)

# All pairwise equivalence tests at level 4 of the BP ontology (quite time consuming, you can jump
# this sentence and use the dataset 'BP.4' which is directly charged with package 'goSorensen')
BP.4 <- equivTestSorensen(allOncoGeneLists,
                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                          onto = "BP", GOLevel = 4)

# All p-values in vector form:
getPvalue(BP.4)
# The symmetric matrix of all p-values:
getPvalue(BP.4, simplify = FALSE)

# All Sorensen-Dice dissimiliraties:
getDissimilarity(BP.4)
getDissimilarity(BP.4, simplify = FALSE)

# 95% upper limits of all one-sided confidence intervals for the 
# Sorensen-Dice dissimilarity:
getUpper(BP.4)
getUpper(BP.4, simplify = FALSE)

getSE(BP.4)
getSE(BP.4, simplify = FALSE)

getTable(BP.4)

BP.4_strict <- upgrade(BP.4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
BP.4_strict
getPvalue(BP.4_strict)

# (Extremely time consuming, for the same reason as before: building all
# contingency tables. 
# Alternatively, you may use the dataset 'cancerEquivSorensen' directly,
# it is automatically charged with the package 'goSorensen'),
# By default, the tests are iterated over all GO ontologies and for levels 3 to 10:
cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists, 
                                            geneUniverse = humanEntrezIDs, 
                                            orgPackg = "org.Hs.eg.db")

# The same but with the bootstrap approach. Even more time consuming and clearly
# unnecessary, see below for a faster approach:
# set.seed(123)
# boot.cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists,
#                                                  boot = TRUE,
#                                                  geneUniverse = humanEntrezIDs, 
#                                                  orgPackg = "org.Hs.eg.db")

# It takes its time, but it is much more faster:
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

# The above multiple equivalence tests are performed without any adjustment for
# testing multiplicity, which is the responsibility of the user. Some suggestions
# follow on how to perform this task.

# From 21 possible p-values (21 = 7 * (7 - 1) / 2) and after the Holm's adjustment for 
# testing multiplicity, identify those who are <= 0.05:
# Function adjSignifPvals in script "adjSignifPvals.R" returns the adjusted significant p-values
# jointly with the enrichment contingency tables, in order to put these values in an adequate
# context (e.g., those who are not very credible due to low table frequencies):
signifPvals_d0_0.4444 <- adjSignifPvals(cancerEquivSorensen)

signifPvals_d0_0.4444$BP
signifPvals_d0_0.4444$CC
signifPvals_d0_0.4444$MF

# For a more restrictive d0 = 0.2857:
signifPvals_d0_0.2857 <- adjSignifPvals(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)))

signifPvals_d0_0.2857$BP
signifPvals_d0_0.2857$CC
signifPvals_d0_0.2857$MF

# The bootstrap version of these tests tends to be more exact (less danger of commiting a type I
# error) and, in any case, tends to be conservative under low enrichment frequencies. So the positive
# results seem more reliable. The number of valid bootstrap replicates, over 10000, is also 
# displayed. Under low table frequencies, some generated bootstrap tables are not adequate for
# Sorensen-Dice computations, but this induces a conservative tendency in the test:
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

# High level of consistency between normal and bootstrap results is a general trend, 
# with greater p-values in the bootstrap case, as expected. For example:
signifPvals_d0_0.2857$BP$`level 5`
boot.signifPvals_d0_0.2857$BP$`level 5`
