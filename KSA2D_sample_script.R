
###############

### Data preparation


devtools::install_github("ginnyintifa/KSA2D")
library(KSA2D)





protData_p = protData_prep(protData_filename = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_cleanData/noNorm_unique_all_annotate_prot.tsv",
                           fudge_factor = 0.01)


psiteData_p =   psiteData_prep(psiteData_filename = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_cleanData/noNorm_unique_all_annotate_phospho.tsv",
                               uniprot_gn_filename = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_webSiteFiles/sp_20200513_uniprot_gn.tsv",
                               fudge_factor = 0.01,
                               protData_p = protData_p)



### Get the column names 

library(data.table)
grp = fread("/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_webSiteFiles/KHU_grouping20200521.tsv",
            stringsAsFactors = F)

co_id = grp$id_no[which(grp$group == "CO")]
cl_id = grp$id_no[which(grp$group == "CL")]
sl_id = grp$id_no[which(grp$group == "SL")]


allg_id = c(co_id, cl_id, sl_id)


coln_5min =  paste0("p_5min_",allg_id)
coln_0min =  paste0("p_0min_",allg_id)
coln_30min =  paste0("p_30min_",allg_id)


set.seed(123)
seeds30 = matrix(sample(c(1:30000),30000,replace = F), nrow = 1000, ncol = 30)

### ID analysis 

all5_prot = comparison_time_points_1d(s1_col_name = coln_5min ,
                                      s2_col_name = coln_0min,
                                      d1_data = protData_p,
                                      nna_cutoff = 6,
                                      all_seeds = seeds30,
                                      permute_times = 500,
                                      working_dir = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_sampleRun/",
                                      compare_name = "test_all_5min_prot")




all5_psite = comparison_time_points_1d( s1_col_name =  coln_5min,
                                        s2_col_name =  coln_0min,
                                        d1_data = psiteData_p,
                                        nna_cutoff = 6,
                                        all_seeds = seeds30,
                                        permute_times = 500,
                                        working_dir = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_sampleRun/",
                                        compare_name = "test_all_5min_psite")

### Find kS relationships



ks_list = ks_map(protData_filename = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_cleanData/noNorm_unique_all_annotate_prot.tsv",
                 psiteData_filename = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_cleanData/noNorm_unique_all_annotate_phospho.tsv",
                 uniprot_gn_filename = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_webSiteFiles/sp_20200513_uniprot_gn.tsv",
                 fudge_factor = 0.01,
                 ksNetwork_filename = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_webSiteFiles/mergeNetwork.tsv",
                 working_dir = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_sampleRun/",
                 ks_outputName = "test_KS_network.tsv")

### 2D analysis 

comparison_time_points_2dZ_5_30(s1_d1_col_name = coln_5min,
                                s2_d1_col_name = coln_0min,
                                s1_d2_col_name = coln_30min,
                                s2_d2_col_name = coln_0min,
                                d1_data = ks_list[[1]],
                                d2_data = ks_list[[2]],
                                nna_cutoff = 6,
                                all_seeds = seeds30,
                                permute_times = 500,
                                working_dir = "/Users/ginny/Google Drive/KSA2D_github/KSA2Dpackage_sampleRun/",
                                compare_name = "test_all_5prot_30psite")




