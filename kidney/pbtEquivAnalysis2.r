# ******************************************************************************
#             PBT GENE LISTS JOINT ENRICHMENT EQUIVALENCE ANALYSIS
# ******************************************************************************

library(goSorensen)
source("..\\adjSignifPvals.R")
# data(humanEntrezIDs)

load("pbtGeneLists2.rda")

# pbtGeneLists2
sapply(pbtGeneLists2, length)

# All pairwise equivalence tests for all levels from 3 to 10 and for all GO ontologies
# (Extremely time consuming. To save time jump, uncomment and run 'load("pbtEquivSorensen2.rda")'
# above (line 65):
pbtAllOntosAndLevels2 <- allEquivTestSorensen(pbtGeneLists2, 
                                              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# save(pbtAllOntosAndLevels2, file = "pbtAllOntosAndLevels2.rda")

boot.pbtAllOntosAndLevels2 <- upgrade(pbtAllOntosAndLevels2, boot = TRUE)
# save(pbtAllOntosAndLevels2, file = "pbtAllOntosAndLevels2.rda")

# load("pbtEquivSorensen2.rda")
# load("boot.pbtEquivSorensen2.rda")


# From 91 possible p-values (91 = 14 * (14 - 1) / 2) and after the Holm's adjustment for testing multiplicity, 
# identify those who are <= 0.05?, excluding NA values.
# Some of these "significant" results may have low reliability due low enrichment
# frequencies (poor accuracy of the underlying asymptotic theory):
# The function "adjSignifPvals" returns the adjusted significant p-values jointly with the enrichment
# contingency tables, in order to put these values in an adequate context (e.g., those who are not very
# credible):
signifPvals_d0_0.4444 <- adjSignifPvals(pbtAllOntosAndLevels2)


signifPvals_d0_0.4444$BP
signifPvals_d0_0.4444$CC
signifPvals_d0_0.4444$MF

# For a more restrictive d0 = 0.2857:
signifPvals_d0_0.2857 <- adjSignifPvals(upgrade(pbtAllOntosAndLevels2, d0 = 1/(1 + 2*1.25)))
signifPvals_d0_0.2857$BP
signifPvals_d0_0.2857$CC
signifPvals_d0_0.2857$MF


# Bootstrap approach, tends to be conservative under low enrichment frequencies so the positive
# results seem more reliable (the number of valid bootstrap replicates, over 10000, is also displayed):
boot.signifPvals_d0_0.4444 <- adjSignifPvals(boot.pbtAllOntosAndLevels2)

boot.signifPvals_d0_0.4444$BP
boot.signifPvals_d0_0.4444$CC
boot.signifPvals_d0_0.4444$MF

# For a more restrictive d0 = 0.2857:
set.seed(123)
boot.signifPvals_d0_0.2857 <- adjSignifPvals(upgrade(boot.pbtAllOntosAndLevels2, d0 = 1/(1 + 2*1.25), boot = TRUE))
boot.signifPvals_d0_0.2857$BP
boot.signifPvals_d0_0.2857$CC
boot.signifPvals_d0_0.2857$MF

