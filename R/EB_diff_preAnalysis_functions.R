### this is a file to record all the functions used in the first step of the program 
### cleaning 



get_phospho_parsed = function(pho)
{
  
  

  parse_site = rbindlist(lapply(1:nrow(pho), function(x){
    #parse_site = rbindlist(lapply(12000:13000, function(x){
    
    
    
    this_mod = pho$Modifications[x]
    
    #### extract phosphorylation from here 
    
    ps = unlist(strsplit(this_mod, split = "Phospho "))[2]
    resp = unlist(strsplit(ps, split = ";"))
    
    this_pho = rbindlist(lapply(1:length(resp),function(i){
      
      this_site = resp[i]
      noP = gsub("\\(\\d*\\)","",this_site)
      noB = gsub("\\[", "", noP)
      noBB =gsub("\\]","",noB)
      noBBB = gsub(" ","",noBB)
      
      get_split = unlist(strsplit(noBBB, split = ""))
      get_res = get_split[1]
      get_pos = as.integer(paste(get_split[-1],collapse = ""))
      this_site_df = data.frame(res = get_res, pos = get_pos,
                                stringsAsFactors = F)
      
      return(this_site_df)
      
    }))
    
    
    
    if(nchar(pho$Accession[x])>1)
    {
      all_acc = unlist(strsplit(pho$Accessions[x], split = "; "))
      
      
      total_df = rbindlist(lapply(1:length(all_acc),function(k) {
        get_info = data.frame(accession = all_acc[k],
                              sequence = pho$Sequence[x],
                              res = this_pho$res,
                              pos = this_pho$pos,
                              stringsAsFactors = F)
        
        this_data = pho[x,8:23]
        get_df = cbind(get_info, this_data[rep(1,nrow(this_pho)),])
        
        rownames(get_df) = NULL
        
        return(get_df)
      }))
      
      if(x%%1000 ==0 )
        cat(x, "\n")
      
      return(total_df)
      
    }
    
    
  }))
  
  
  
  
  
  return(parse_site)
}




f = function(parse_pho)
{
  #parse_pho = parse_pho1
  res_na = which(is.na(parse_pho$res))
  
  pos_na = which(is.na(parse_pho$pos))
  
  get_na = union(res_na, pos_na)
  ppho = parse_pho[-get_na,]
  
  return(ppho)
}









