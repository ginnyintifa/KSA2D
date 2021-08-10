

#' Find kinase-substrate relationships 

#' @param protData_filename file name of protein data frame  
#' @param psiteData_filename file name of phosphorylation data frame  
#' @param uniprot_gn_filename file name of the uniprot accession number and gene name matching info, downloadable from our website
#' @param fudge_factor a factor between 0 and 1, default to 0.01 to be added to the entire data for robustness
#' @param ksNetwork_filename file name of known kinase-substrate relationships, downloadable from our website
#' @param working_dir the directory output files will be deposited in 
#' @param ks_outputName a label for the mapped relationships
#' @import data.table dplyr magrittr 
#' @keywords map
#' @export
#' @examples
#' ks_map



ks_map = function(protData_filename,
                  psiteData_filename,
                  uniprot_gn_filename,
                  fudge_factor,
                  ksNetwork_filename,
                  sub_norm = T,
                  working_dir,
                  ks_outputName)


{
  
sel_prot= fread(protData_filename,
                stringsAsFactors = F,
                data.table = F)

p_psite = fread(psiteData_filename,
                stringsAsFactors = F,
                data.table = F)

ks_network = fread(ksNetwork_filename,
                   stringsAsFactors = F,
                   data.table = F)

uniprot_gn= fread(uniprot_gn_filename,
                  stringsAsFactors = F,
                  data.table = F)


parse_psite = p_psite%>%
  dplyr::left_join(uniprot_gn, by  = "accession")%>%
  dplyr::select(accession, geneName,res, phosphoPos, resPos, everything())

#### here I think it is better the convert every thing to log2 scale first 



### add a fudge factor 
fprot = quantile(sel_prot[,-c(1:3)],probs = fudge_factor, na.rm = T)

sel_prot_log_f = data.frame(sel_prot[,c(1:3)],log2(sel_prot[,-c(1:3)]+fprot), stringsAsFactors = F)


fpsite = quantile(parse_psite[,-c(1:5)],probs = fudge_factor, na.rm = T)

parse_psite_log_f = data.frame(parse_psite[,c(1:5)],log2(parse_psite[,-c(1:5)]+fpsite),  stringsAsFactors = F)




data_ks = rbindlist(lapply(1:nrow(ks_network), function(x) {
  
  
  if(x%%100 ==0)
    cat(x, "\n")
  kin = ks_network$kinase[x]
  sub = ks_network$substrate[x]
  subsite = ks_network$substrate_site[x]
  source  = paste0(ks_network$source1[x], "_",ks_network$source2[x],"_", ks_network$source3[x])
  kin_prot = sel_prot_log_f%>%
    dplyr::filter(geneName == kin)
  
  kin_prot_data = kin_prot[,-c(1:3)]
  
  kin_prot_patients = colnames(kin_prot_data)[which(!is.na(kin_prot_data))]
  
  
  
  if(nrow(kin_prot)>0)
  {
    
    
    prot_df = data.frame(barcode = colnames(kin_prot_data),  
                         prot = as.numeric(kin_prot_data),
                         stringsAsFactors = F)
    
    
    sub_psite = parse_psite_log_f%>%
      dplyr::filter(geneName == sub)
    
    
    sub_psite_data = sub_psite[,-c(1:5)]
    
    
    if(nrow(sub_psite)>0)
    {
      ### find the protein data for the substrate 
      
      
      substrate_prot = sel_prot_log_f%>%
        dplyr::filter(geneName == sub)
      
      substrate_prot_data = substrate_prot[,-c(1:3)]
      
      substrate_prot_df = data.frame(barcode = colnames(substrate_prot_data),  
                                     substrate_prot = as.numeric(substrate_prot_data),
                                     stringsAsFactors = F)
      
      
      ### exact site?
      subsite_psite = sub_psite%>%
        dplyr::filter(resPos == subsite)
      
      
      if(nrow(subsite_psite)>0)
      {
        
        
        subsite_psite_data = subsite_psite[,-c(1:5)]
        
        
        psite_df = data.frame(barcode = colnames(subsite_psite_data),  
                              psite = as.numeric(subsite_psite_data),
                              stringsAsFactors = F)
        
        prot_psite_join = prot_df %>%
          dplyr::left_join(psite_df, by = "barcode")%>%
          dplyr::left_join(substrate_prot_df, by = "barcode")%>%
          dplyr::mutate(subProt_psite = psite - substrate_prot)
        
        #### can I transpose this?
        
        t_join_data = t(prot_psite_join[,-1])
        
        df_name = paste(kin,sub, subsite, source, rownames(t_join_data), sep = "_")
        
        nn_df_name = gsub("_NA","",df_name )
        
        
        rownames(t_join_data) = NULL
        
        t_join_df = data.frame(name = nn_df_name,
                               t_join_data,
                               stringsAsFactors = F)
        
        
        colnames(t_join_df) = c("name", prot_psite_join$barcode)
        
        
        
        return(t_join_df)
        
        
      }
      
    }
  }
  
  
  
}))


write.table(data_ks, 
            paste0(working_dir, ks_outputName),
            quote = F, row.names = F, sep = "\t")




ks_pair = unique(data_ks)

n_pair = nrow(ks_pair)/4
pair_seq = c(1:n_pair)
kprot_row = (pair_seq-1)*4+1
supsite_row = (pair_seq-1)*4+4
pair_names = gsub("_prot","", ks_pair$name[kprot_row])
kprot_data = ks_pair[kprot_row,]
supsite_data = ks_pair[supsite_row,]



if(sub_norm = F)
{
  site_row = (pair_seq-1)*4 + 2
  supsite_data = ks_pair[site_row,]
  
}


pair_data_list = list(kprot_data, supsite_data)

return(pair_data_list)

}


