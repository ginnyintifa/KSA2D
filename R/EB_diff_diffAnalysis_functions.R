#EB_diff_diffAnalysis_functions



noNA = function(x)
{
  return(length(which(!is.na(x))))
}


shuffle_time_points = function(sample_size, numPermu, seeds)
{

  these_shuffle = lapply(1:numPermu, function(x) {
    
    this_row = rep(0,sample_size)
    for(i in 1:sample_size)
    {
      set.seed(seeds[x, i])
      this_row[i] = sample(c(1:2),  1, replace = F)
      
    }
    
    return(this_row)
  })
  
  
  return(these_shuffle)
  
  
}



s1_s2_filter = function(s1_d1,
                        s1_d2,
                        s2_d1,
                        s2_d2,
                        pair_names = pair_names,
                        nna_cutoff = 6)
{
  
   # 
  nas1_d1 = s1_d1
  nas1_d2 = s1_d2
  nas2_d1 = s2_d1
  nas2_d2 = s2_d2
  
  
  ### keep the ones with No NAs 
  for(i in 1:nrow(s1_d1))
  {
    s1_d1_na = which(is.na(s1_d1[i,]))
    s1_d2_na = which(is.na(s1_d2[i,]))
    s1_na = union(s1_d1_na, s1_d2_na)
    
    s2_d1_na = which(is.na(s2_d1[i,]))
    s2_d2_na = which(is.na(s2_d2[i,]))
    s2_na = union(s2_d1_na, s2_d2_na)
    
    s12_na = union(s1_na, s2_na)
    
      # 
    if(length(s12_na)>0)
    {
      
      nas1_d1[i,s12_na] = NA
      nas1_d2[i,s12_na] = NA
      
      nas2_d1[i,s12_na] = NA
      nas2_d2[i,s12_na] = NA
      
    }
  }
  
  
  nas1_d1 = as.matrix(nas1_d1)
  nas1_d2 = as.matrix(nas1_d2)
  nas2_d1 = as.matrix(nas2_d1)
  nas2_d2 = as.matrix(nas2_d2)
  
  
  nas1_nna = apply(nas1_d1, 1, noNA)  ### same as   nas2_nna = apply(nas2_d1, 1, noNA)
  
  
  keep_gene = which(nas1_nna>=nna_cutoff)
  keep_gene_names = pair_names[keep_gene]
  
  nas1_d1_use = nas1_d1[keep_gene,]
  nas1_d2_use = nas1_d2[keep_gene,]
  nas2_d1_use = nas2_d1[keep_gene,]  
  nas2_d2_use = nas2_d2[keep_gene,]
  
  
  use_list = list(nas1_d1_use, nas1_d2_use, nas2_d1_use, nas2_d2_use, keep_gene_names)
  
  return(use_list)
  
  
}




s1_s2_filter_1d = function(s1_d1,
                           s2_d1,
                           names,
                           nna_cutoff)
{
  # 
  
  nas1_d1 = s1_d1
  nas2_d1 = s2_d1
  
  
  ### keep the ones with No NAs 
  for(i in 1:nrow(s1_d1))
  {
    s1_na = which(is.na(s1_d1[i,]))
    
    s2_na = which(is.na(s2_d1[i,]))
    
    s12_na = union(s1_na, s2_na)
    
    if(length(s12_na)>0)
    {
      
      nas1_d1[i,s12_na] = NA
      
      nas2_d1[i,s12_na] = NA
      
    }
  }
  
  
  nas1_d1 = as.matrix(nas1_d1)
  nas2_d1 = as.matrix(nas2_d1)
  
  
  nas1_nna = apply(nas1_d1, 1, noNA)  ### same as   nas2_nna = apply(nas2_d1, 1, noNA)
  
  
  keep_gene = which(nas1_nna>=nna_cutoff)
  keep_gene_names = names[keep_gene]
  
  nas1_d1_use = nas1_d1[keep_gene,]
  nas2_d1_use = nas2_d1[keep_gene,]  
  
  
  use_list = list(nas1_d1_use, nas2_d1_use, keep_gene_names)
  
  return(use_list)
  
  
}



#### add an exchangeability factor s0 



generate_oneSample_tstat_1d = function(s1_df1, s2_df1)
{
  
  
  # s1_df1 = s1_d1_use
  # s2_df1 = s2_d1_use
  # 
  s12_df1 = s1_df1 - s2_df1
  
  s12_noNA = apply(s12_df1, 1, noNA)
  
  s12_mean = apply(s12_df1, 1, mean, na.rm = T)
  s12_sd = apply(s12_df1, 1, sd, na.rm = T)
  
  s0 = quantile(s12_sd, probs = 0.1, na.rm = T)
  s12_sd_adjust = s12_sd + s0 
  
  
  
  oneSample_tstat = s12_mean/(s12_sd_adjust/sqrt(s12_noNA))#***
  
  return(oneSample_tstat)
  
  
}



generate_oneSample_tstat_1d_90 = function(s1_df1, s2_df1)
{
  
  
  # s1_df1 = s1_d1_use
  # s2_df1 = s2_d1_use
  # 
  s12_df1 = s1_df1 - s2_df1
  
  s12_noNA = apply(s12_df1, 1, noNA)
  
  s12_mean = apply(s12_df1, 1, mean, na.rm = T)
  s12_sd = apply(s12_df1, 1, sd, na.rm = T)
  
  s0 = quantile(s12_sd, probs = 0.9, na.rm = T)
  s12_sd_adjust = s12_sd + s0 
  
  
  
  oneSample_tstat = s12_mean/(s12_sd_adjust/sqrt(s12_noNA))#***
  
  return(oneSample_tstat)
  
  
}



compute_log2foldchage = function(s1_df1, s2_df1)
{
  
  # 
  # s1_df1 = s1_d1_use
  # s2_df1 = s2_d1_use
  # # 
  # 
  s12_df1 = s1_df1 - s2_df1
  
  s12_mean = apply(s12_df1, 1, mean, na.rm = T)
  
  return(s12_mean)
  
  
}






##### functions for producing the density plots in 2d scenario 
### 0704  change the orientation 


density_persp = function(pdf_name,
                         x,
                         y,
                         z,
                         x_margin_dens_x,
                         x_margin_dens_y,
                         y_margin_dens_x,
                         y_margin_dens_y,
                         my_xlab,
                         my_ylab,
                         my_zlab,
                         my_main,
                         my_theta,
                         my_phi,
                         my_shade)
  
{
  

  pdf(pdf_name, useDingbats = F)
  
  
  zl= c(0, 1.5*max(z))
  
  
  col.pal = colorRampPalette(c("blue","red"))
  colors = col.pal(100)
  z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
  z.facet.range<-cut(z.facet.center, 100)
  
  
  
  trmat = persp(x,y,z, 
                theta=my_theta,
                phi = my_phi,
                box = T,
                border = NA,
                zlim = zl, 
                xlab = my_xlab,
                ylab = my_ylab,
                zlab = my_zlab,
                main = my_main,
                col = colors[z.facet.range],
                shade=my_shade, 
                lwd = 0.01,
                ticktype = "detailed")
  

  ### calculate a scale factor for the two margnial distributions 

  max_z = zl[2]
  max_xy = max(c(max(x_margin_dens_y),max(y_margin_dens_y)))
  ms = (max_z/max_xy)*0.9
  
  
  
  lines( trans3d( x_margin_dens_x, max(y), ms*x_margin_dens_y,
                  trmat), lwd=1, col=rgb(0.2,0.6,1,0.3) )
  
  
  lines( trans3d( c(0,0), c(max(y),max(y)),c(0,zl[2]),
                  trmat), lwd=0.5,lty = 2)
  

  lines( trans3d( min(x), y_margin_dens_x, ms*y_margin_dens_y,
                  trmat), lwd=1, col=rgb(0.2,0.6,1,0.3) )
  
  lines( trans3d( c(min(x),min(x)), c(0,0),c(0,zl[2]),
                  trmat), lwd=0.5,lty = 2)
  
  
  
  dev.off()
  
  
}





generate_twoSample_tstat_1d = function(s1_df1, s2_df1)
{
  
  # s1_df1 = g1_common_fc
  # s2_df1 = g2_common_fc  
  # 
  s1_noNA = apply(s1_df1, 1, noNA)
  s2_noNA = apply(s2_df1, 1, noNA)
  
  
  s1_mean = apply(s1_df1, 1, mean, na.rm = T)
  s2_mean = apply(s2_df1, 1, mean, na.rm = T)
  s2_s1_diff = s1_mean-s2_mean   
  
  s1_sd = apply(s1_df1, 1, sd, na.rm = T)
  s2_sd = apply(s2_df1, 1, sd, na.rm = T)
  
  s1_sd_s0 = quantile(s1_sd, probs = 0.1, na.rm = T)
  s2_sd_s0 = quantile(s2_sd, probs = 0.1, na.rm = T)
  
  
  
  s1_sd_adjust = s1_sd + s1_sd_s0
  s2_sd_adjust = s2_sd + s2_sd_s0
  
  twoSample_tstat = s2_s1_diff/(sqrt(s1_sd_adjust^2/s1_noNA + s2_sd_adjust^2/s2_noNA)) #***
  
  return(twoSample_tstat)
  
  
}





generate_permute_cols = function(s1_size, s2_size, numPermu, seeds)
{
  #   s1_size = 10
  #   s2_size = 10
  #   numPermu = 20
  #   will be comparison between c1 and c2 in next stage 
  
  get_cols_s = matrix(0, nrow = numPermu, ncol = s1_size/2)
  
  
  ### always assume s1 and s2 have same sizes
  
  
  for(i in 1:numPermu)
  {
    set.seed(seeds[i])
    
    this_s = sample(c(1:s1_size), floor(s1_size/2), replace = F)
    
    get_cols_s[i,] = this_s
    
    
  }
  
  
  
  return(get_cols_s)
  
  
}















