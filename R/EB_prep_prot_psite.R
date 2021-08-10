
#' convert to log2 scale 
#'
#' @param protData_filename file name of protein data frame  
#' @param fudge_factor a factor between 0 and 1, default to 0.01 to be added to the entire data for robustness
#
#' @keywords prepare data
#' @import dplyr data.table magrittr 
#' @export
#' @examples
#' protData_prep


protData_prep = function(protData_filename,
                         fudge_factor = 0.01)
{

  sel_prot= fread(protData_filename,
                  stringsAsFactors = F,
                  data.table = F)
  
  ### prepare a 1 col name for each here 
  
  prot_names = paste(sel_prot$accession, sel_prot$geneName, sep = "_")
  
  ### add a fudge factor 
  fprot = quantile(sel_prot[,-c(1:3)],probs = fudge_factor, na.rm = T)
  
  prot_log_f = data.frame(prot_names,log2(sel_prot[,-c(1:3)]+fprot), stringsAsFactors = F)
  
  
  return(prot_log_f)
  
}


#' convert to log2 scale and normalise psites with mother protein abundance 
#'
#' @param psiteData_filename file name of phosphorylation data frame  
#' @param uniprot_gn_filename file name of the uniprot accession number and gene name matching info, downloadable from our website
#' @param fudge_factor a factor between 0 and 1, default to 0.01 to be added to the entire data for robustness
#' @param protData_p processed protein data, the output object of protData_prep function 
#' @keywords prepare data
#' @export
#' @examples
#' protData_prep

psiteData_prep = function(psiteData_filename,
                          uniprot_gn_filename,
                          sub_norm = T,
                          fudge_factor,
                          protData_p)
{
  
  
  p_psite = fread(psiteData_filename,
                  stringsAsFactors = F,
                  data.table = F)
  
  uniprot_gn = fread(uniprot_gn_filename,
                     stringsAsFactors = F, 
                     data.table = F)
  prot_log_f = protData_p
  
  parse_psite = p_psite%>%
    dplyr::left_join(uniprot_gn, by  = "accession")%>%
    dplyr::select(accession, geneName,res, phosphoPos, resPos, everything())
  
  psite_names = paste(parse_psite$accession, parse_psite$geneName, 
                      parse_psite$resPos,
                      sep = "_")
  fpsite = quantile(parse_psite[,-c(1:5)],probs = fudge_factor, na.rm = T)
  
  psite_log_f = data.frame(psite_names,log2(parse_psite[,-c(1:5)]+fpsite),  stringsAsFactors = F)
  
  norm_psite_log_f = psite_log_f
  
  if(sub_norm == T)
  {
    
  
  norm_psite_log_f = rbindlist(lapply(1:nrow(psite_log_f), function(i) {
    
    
    if(i%%1000 ==0)
      cat(i, "\n")
    this_name = psite_log_f$psite_names[i]
    
    this_prot = paste(unlist(strsplit(this_name,split = "_"))[1:2], collapse = "_")
    
    fil_prot = prot_log_f%>%
      dplyr::filter(prot_names == this_prot)
    if(nrow(fil_prot)>0)
    {
      norm_psite_data = psite_log_f[i, -1] - fil_prot[,-1]
      
      norm_psite_df = data.frame(psite_names = this_name, 
                                 norm_psite_data,
                                 stringsAsFactors = F)
      
      return(norm_psite_df)
    }
    
    
  }))
  
  norm_psite_log_f = data.frame(norm_psite_log_f, stringsAsFactors = F)
  }
  
  
  
  return(norm_psite_log_f)
  
  
  
}


