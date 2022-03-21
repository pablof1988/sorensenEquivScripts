date()
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                         onto = "MF", GOLevel = 3, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                         onto = "MF", GOLevel = 4, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                         onto = "MF", GOLevel = 5, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                         onto = "MF", GOLevel = 6, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                  onto = "MF", GOLevel = 7, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                  onto = "MF", GOLevel = 8, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", 
                  onto = "MF", GOLevel = 9, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
set.seed(123)
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", boot = TRUE, 
                         onto = "MF", GOLevel = 3, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
set.seed(123)
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", boot = TRUE, 
                         onto = "MF", GOLevel = 4, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
set.seed(123)
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", boot = TRUE, 
                         onto = "MF", GOLevel = 5, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
set.seed(123)
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", boot = TRUE, 
                         onto = "MF", GOLevel = 6, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
set.seed(123)
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", boot = TRUE, 
                         onto = "MF", GOLevel = 7, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
set.seed(123)
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", boot = TRUE, 
                         onto = "MF", GOLevel = 8, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()

date()
set.seed(123)
tst <- equivTestSorensen(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["sanger"]], 
                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", boot = TRUE, 
                         onto = "MF", GOLevel = 9, listNames = c("Vogelstein", "sanger"))
tst
getTable(tst)
date()