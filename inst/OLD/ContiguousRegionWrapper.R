
CpGsOrdered_char <- OrderCpGsByLocation(CpGs_char = row.names (betaValues_mtx),
                                        CpGlocations_df,
                                        Output="CpGs")

betaValuesOrdered_mtx <- betaValues_mtx [CpGsOrdered_char,]

betaCluster_mtx <- t(betaValuesOrdered_mtx)

comethylatedSubregionsCpGs_df <- FindComethylatedRegions(betaCluster_mtx,
                                                         minCpGs_int=3,
                                                         CpGlocations_df,
                                                         threshold_r_num=0.5)



