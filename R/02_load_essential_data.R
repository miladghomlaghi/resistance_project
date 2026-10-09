

load("./Data/matching/rob8_core_to_match")
load("./Data/topologies/frequency_matlab_210220.RData")
load("./Data/matching/merged_Signor_PathwaysPlus.RData")
load("./Data/databases_proteins/targetsSignor.RData")
load(paste0("./Data/matching/", "Targets.RData", collapse = ""))
load(paste0("./Data/matching/", "SIGNOR_NETSCAN.RData", collapse = ""))

topos_NL<- topos[, !1:9 %in% c(1,5,9)]
