



#' Compare between 2 conditions (paired) in 1D
#'
#' @param s1_col_name names of the columns in the first condition 
#' @param s2_col_name names of the columns in the second condition
#' @param d1_data dataframe for the data, rows are genes, columns are patients in 2 conditions
#' @param nna_cutoff genes with values in at least the number of samples are included in comparison
#' @param all_seeds seeds for the permutation
#' @param permute_times number of permutations 
#' @param working_dir the directory output files will be deposited in 
#' @param compare_name a label for the anlaysis 
#' @import dplyr data.table magrittr splines stringr MASS 
#' @keywords 1D
#' @export
#' @examples
#' camparison_time_points_1d


comparison_time_points_1d = function(s1_col_name,
                                     s2_col_name,
                                     d1_data,
                                     nna_cutoff,
                                     all_seeds,
                                     permute_times,
                                     working_dir,
                                     compare_name)
  
{
  
  s1_d1 = d1_data[,s1_col_name]
  s2_d1 = d1_data[,s2_col_name]
  
  #### filter data
  s1_s2_use = s1_s2_filter_1d(
    s1_d1 = s1_d1,
    s2_d1 = s2_d1,
    names = d1_data[,1],
    nna_cutoff = nna_cutoff
  )
  

  cat("data filtered.", "\n")
  #### generate Z
  s1_d1_use = s1_s2_use[[1]]
  s2_d1_use = s1_s2_use[[2]]
  
  s1_s2_Z_use = generate_oneSample_tstat_1d(s1_df1 = s1_d1_use,
                                            s2_df1 = s2_d1_use) 
  

  s1_s2_lfc = compute_log2foldchage(s1_df1 = s1_d1_use,
                                    s2_df1 = s2_d1_use)
  
  t_Z =unlist(lapply(1:nrow(s1_d1_use), function(x)
    
  {
    this_t   = t.test(s1_d1_use[x,], s2_d1_use[x,], paired = T)
    
    this_p = this_t$p.value
    
    return(this_p)
  })) 
  
  
  
  
  Z_pdf_name = paste0(working_dir,compare_name, "_Z_hist.pdf")
  
  pdf(Z_pdf_name,useDingbats = F)
  hist(s1_s2_Z_use, breaks  = 100)
  dev.off()
  cat("Z generated.", "\n")
  
  #### generate null z
  
  null_permute = 
    shuffle_time_points( sample_size = ncol(s1_d1_use),
                         numPermu = permute_times,
                        # numPermu = 100,
                         seeds = all_seeds)
  
  ### reduce the null_permute to the unique number 
  
  if(length(null_permute)>2^length(s1_col_name))
  {
    null_permute = unique(null_permute)
    permute_times = 2^length(s1_col_name)
  }
  
  
  s1_s2_z_null_use = lapply(1:permute_times, function(x)
  {
    # x = 15
    this_assign = null_permute[[x]]
    
    #first = which(this_assign ==1)
    which_flip = which(this_assign == 2)
    
    this_c1_d1 = s1_d1_use
    this_c2_d1 = s2_d1_use
    
    this_c1_d1[,which_flip] = s2_d1_use[,which_flip]
    this_c2_d1[,which_flip] = s1_d1_use[,which_flip]
    
    
    
    permute_use = s1_s2_filter_1d(s1_d1 = this_c1_d1,
                                  s2_d1 = this_c2_d1,
                                  names = s1_s2_use[[3]],
                                  nna_cutoff = nna_cutoff)
    
    
    this_z= generate_oneSample_tstat_1d(s1_df1 = permute_use[[1]],
                                        s2_df1 = permute_use[[2]])
    if(x%%100 == 0)
      cat(x, "\n")
    
    return(this_z)
    
  })
  
  
  
  s1_s2_z_null_use_cat =  Reduce(c, s1_s2_z_null_use)
  
  z_null_pdf_name = paste0(working_dir,compare_name, "_z_null_hist.pdf")
  pdf(z_null_pdf_name,useDingbats = F)
  
  hist(s1_s2_z_null_use_cat, breaks  = 100)
  dev.off()
  

  
    
  
  Z_dens = density(s1_s2_Z_use, n = 2000)
  z_null_dens = density(s1_s2_z_null_use_cat, n = 2000)
  
  
  
  
  Zz_dens_pdf_name = paste0(working_dir,compare_name, "_Zz_dens.pdf")
  ylimit = max(z_null_dens$y)
  xmin = min(c(z_null_dens$x, Z_dens$x))
  xmax = max(c(z_null_dens$x, Z_dens$x))
  pdf(Zz_dens_pdf_name,useDingbats = F)
  plot(z_null_dens, col = "green", ylim = c(0,ylimit), xlim = c(xmin,xmax))
  lines(Z_dens)
  dev.off()
  
  cat("z null generated.", "\n")
  
  
  #############################################
  #############################################
  #############################################
  
  
  
  
  Z_zero_which = which(abs(Z_dens$x-0) == min(abs(Z_dens$x-0)))[1]
  Z_zero_dens = Z_dens$y[Z_zero_which]
  
  z_null_zero_which = which(abs(z_null_dens$x-0) == min(abs(z_null_dens$x-0)))[1]
  z_null_zero_dens = z_null_dens$y[z_null_zero_which]
  
  
  s1_s2_p0_use = Z_zero_dens/z_null_zero_dens
  ### calculate f/f0 for each Z 
  
 
  f02f = rep(0, length(s1_s2_Z_use))
  
  for(i in 1:length(s1_s2_Z_use))
  {
    this_value = s1_s2_Z_use[i]
    
    ### find the most similar in z_null_dens 
    dis_this = abs(this_value-Z_dens$x)
    nearest = which(dis_this == min(dis_this))[1]
    this_Z_f = Z_dens$y[nearest]
    
    dis_this_null = abs(this_value-z_null_dens$x)
    nearest_null = which(dis_this_null == min(dis_this_null))[1]
    this_Z_f0 = z_null_dens$y[nearest_null]
    
    f02f[i] = this_Z_f0/this_Z_f
    
  }
  
  Z_p0 = s1_s2_p0_use*f02f
  Z_p1 = 1-Z_p0
  
  p1_corrected = Z_p1
  p1_corrected[which(Z_p1<0)]=0
  lfdr = 1-p1_corrected
  
  
  pc = rep(rgb(0.1,0,0,0.5), length(s1_s2_Z_use))
  pc[which( p1_corrected>0.9)] = rgb(0.5,0,0,0.5)
  
  pc[which(s1_s2_Z_use>0 & p1_corrected>0.95)] = rgb(1,0,0,0.5)
  pc[which(s1_s2_Z_use<0 & p1_corrected>0.95)] = rgb(1,0.5,0,0.5)
  
  ### new overlay 
  posterior_pdf_name = paste0(working_dir,compare_name, "_Z_posterior.pdf")
  
  ### when plot remove the outlier ones 
  
  pdf(posterior_pdf_name,useDingbats = F)
  
  plot(x = s1_s2_Z_use,
       y = p1_corrected,
       ylim = c(0,1),
       type = "p",
       col = pc,
       pch = 16)
  dev.off()
  

  ## a scatter plot
  sanity_pdf_name = paste0(working_dir,compare_name, "_ttest_check.pdf")

  ### when plot remove the outlier ones

  #leftones = which(s1_s2_Z_use< 0)

  pdf(sanity_pdf_name,useDingbats = F)

  plot(x = t_Z,
       y = lfdr,
       xlim = c(0,1),
       ylim = c(0,1),
       type = "p",
       cex = 0.5,
       xlab = "ttest p value",
       ylab = "local fdr",
       pch = 16)
  abline(a= 0, b = 1,lty = 2)

  plot(x = s1_s2_Z_use,
       y = t_Z,
       xlim = c(-15,15),
       ylim = c(0,1),
       type = "p",
       cex = 0.5,
       xlab = "Z",
       ylab = "ttest pvalue/lfdr",
       pch = 16)
  lines(x = s1_s2_Z_use,
        y = lfdr,
        type = "p",
        cex = 0.5,
        col = "red",
        pch = 16)
  abline(v =0,lty = 2)


  dev.off()

  # 
  result_df = data.frame(name =s1_s2_use[[3]],
                         fc = s1_s2_lfc,
                         Z = s1_s2_Z_use,
                         p1 = Z_p1,
                         p1_corrected = p1_corrected,
                         fdr = 1-p1_corrected,
                         stringsAsFactors = F)%>%
    dplyr::arrange(desc(p1))
  
  result_df_name = paste0(working_dir, compare_name,"_Z_result.tsv")
  
  write.table(result_df, result_df_name,
              quote = F, row.names = F, sep = "\t")
  
  
  f1_dens  = rep(0, length(Z_dens$x))
  
  for(i in 1:length(Z_dens$x))
  {
    
    this_value = Z_dens$x[i]
    this_dens = Z_dens$y[i]
    
    ### find the most similar in z_null_dens 
    dis_this = abs(this_value-z_null_dens$x)
    nearest = which(dis_this == min(dis_this))[1]
    this_null_dens = z_null_dens$y[nearest]
    
    f1_dens[i] = (this_dens-s1_s2_p0_use*this_null_dens)/(1-s1_s2_p0_use)  ### here is problematic in the 2d scenario, check this 
    
  }
 
    f1_dens[which(f1_dens<0)] = 0 
  
  
  overlay_dens_pdf_name =  paste0(working_dir, compare_name,"_dens_overlay.pdf")
  ylimit = max(c(max(z_null_dens$y),Z_dens$y, f1_dens))
  
  pdf(overlay_dens_pdf_name,useDingbats = F)
  plot(x = Z_dens$x, y = f1_dens, type = "l",ylim = c(0,ylimit), xlim = c(-10,10), col = "red")
  lines(z_null_dens, col = "green")
  lines(Z_dens)
  dev.off()
  
  
  #### new overlay
  pdf( paste0(working_dir, compare_name,"_dens_overlay_prop.pdf"))
  
  this_y_lim = max(c(Z_dens$y, s1_s2_p0_use * z_null_dens$y,(1-s1_s2_p0_use)*f1_dens))
  
  plot(x = Z_dens$x, y = Z_dens$y,ylim = c(0,this_y_lim), xlim = c(-10,10), type = "l")
  lines(x = z_null_dens$x, y= s1_s2_p0_use * z_null_dens$y, col = "green")
  lines(x = Z_dens$x, y = (1-s1_s2_p0_use)*f1_dens, type = "l", col = "red")
  
  dev.off()
  
  
  
  cat("posterior calculated.", "\n")
  
  
  return(s1_s2_p0_use)
  
  
  
  
}



####

comparison_time_points_1d_90 = function(s1_col_name,
                                     s2_col_name,
                                     d1_data,
                                     nna_cutoff,
                                     all_seeds,
                                     permute_times,
                                     compare_name,
                                     working_dir)
  
{
 
  s1_d1 = d1_data[,s1_col_name]
  s2_d1 = d1_data[,s2_col_name]
  
  #### filter data
  s1_s2_use = s1_s2_filter_1d(
    s1_d1 = s1_d1,
    s2_d1 = s2_d1,
    names = d1_data[,1],
    nna_cutoff = 6
  )
  
  ### this function is not very correct 
  
  cat("data filtered.", "\n")
  #### generate Z
  s1_d1_use = s1_s2_use[[1]]
  s2_d1_use = s1_s2_use[[2]]
  
  s1_s2_Z_use = generate_oneSample_tstat_1d_90(s1_df1 = s1_d1_use,
                                            s2_df1 = s2_d1_use) 
  
  ### produce the original log2 fold change 
  
  s1_s2_lfc = compute_log2foldchage(s1_df1 = s1_d1_use,
                                    s2_df1 = s2_d1_use)
  
  t_Z =unlist(lapply(1:nrow(s1_d1_use), function(x)
    
  {
    this_t   = t.test(s1_d1_use[x,], s2_d1_use[x,], paired = T)
    
    this_p = this_t$p.value
    
    return(this_p)
  })) 
  
  
  
  
  Z_pdf_name = paste0(working_dir,compare_name, "_Z_hist.pdf")
  
  pdf(Z_pdf_name,useDingbats = F)
  hist(s1_s2_Z_use, breaks  = 100)
  dev.off()
  cat("Z generated.", "\n")
  
  #### generate null z
  
  null_permute = 
    shuffle_time_points( sample_size = ncol(s1_d1_use),
                         numPermu = permute_times,
                         seeds = all_seeds)
  
  
  
  s1_s2_z_null_use = lapply(1:permute_times, function(x)
  {
    # x = 15
    this_assign = null_permute[[x]]
    
    #first = which(this_assign ==1)
    which_flip = which(this_assign == 2)
    
    this_c1_d1 = s1_d1_use
    this_c2_d1 = s2_d1_use
    
    this_c1_d1[,which_flip] = s2_d1_use[,which_flip]
    this_c2_d1[,which_flip] = s1_d1_use[,which_flip]
    
    
    
    permute_use = s1_s2_filter_1d(s1_d1 = this_c1_d1,
                                  s2_d1 = this_c2_d1,
                                  names = s1_s2_use[[3]],
                                  nna_cutoff = 6)
    
    
    this_z= generate_oneSample_tstat_1d_90(s1_df1 = permute_use[[1]],
                                        s2_df1 = permute_use[[2]])
    if(x%%100 == 0)
      cat(x, "\n")
    
    return(this_z)
    
  })
  
  
  
  s1_s2_z_null_use_cat =  Reduce(c, s1_s2_z_null_use)
  
  z_null_pdf_name = paste0(working_dir,compare_name, "_z_null_hist.pdf")
  pdf(z_null_pdf_name,useDingbats = F)
  
  hist(s1_s2_z_null_use_cat, breaks  = 100)
  dev.off()
  
  
  
  Z_dens = density(s1_s2_Z_use, n = 2000)
  z_null_dens = density(s1_s2_z_null_use_cat, n = 2000)
  
  
  Zz_dens_pdf_name = paste0(working_dir,compare_name, "_Zz_dens_more.pdf")
  ylimit = max(z_null_dens$y)
  pdf(Zz_dens_pdf_name,useDingbats = F)
  plot(z_null_dens, col = "green", ylim = c(0,ylimit), xlim = c(-10,10))
  lines(Z_dens)
  dev.off()
  
  cat("z null generated.", "\n")
  
  
  #############################################
  #############################################
  #############################################
  
  
  
  
  Z_zero_which = which(abs(Z_dens$x-0) == min(abs(Z_dens$x-0)))[1]
  Z_zero_dens = Z_dens$y[Z_zero_which]
  
  z_null_zero_which = which(abs(z_null_dens$x-0) == min(abs(z_null_dens$x-0)))[1]
  z_null_zero_dens = z_null_dens$y[z_null_zero_which]
  
  
  s1_s2_p0_use = Z_zero_dens/z_null_zero_dens
  ### calculate f/f0 for each Z 
  
  
  f02f = rep(0, length(s1_s2_Z_use))
  
  for(i in 1:length(s1_s2_Z_use))
  {
    this_value = s1_s2_Z_use[i]
    
    ### find the most similar in z_null_dens 
    dis_this = abs(this_value-Z_dens$x)
    nearest = which(dis_this == min(dis_this))[1]
    this_Z_f = Z_dens$y[nearest]
    
    dis_this_null = abs(this_value-z_null_dens$x)
    nearest_null = which(dis_this_null == min(dis_this_null))[1]
    this_Z_f0 = z_null_dens$y[nearest_null]
    
    f02f[i] = this_Z_f0/this_Z_f
    
  }
  
  Z_p0 = s1_s2_p0_use*f02f
  Z_p1 = 1-Z_p0
  
  
  ### new overlay 
  posterior_pdf_name = paste0(working_dir,compare_name, "_Z_posterior.pdf")
  
  ### when plot remove the outlier ones 
  
  pdf(posterior_pdf_name,useDingbats = F)
  
  plot(x = s1_s2_Z_use,
       y = Z_p1,
       ylim = c(-1,1),
       type = "p",
       pch = 16)
  dev.off()
  
  p1_corrected = Z_p1
  p1_corrected[which(Z_p1<0)]=0
  lfdr = 1-p1_corrected
  ### a scatter plot 
  sanity_pdf_name = paste0(working_dir,compare_name, "_ttest_check.pdf")
  
 
  pdf(sanity_pdf_name,useDingbats = F)
  
  plot(x = t_Z,
       y = lfdr,
       xlim = c(0,1),
       ylim = c(0,1),
       type = "p",
       cex = 0.5,
       xlab = "ttest p value",
       ylab = "local fdr",
       pch = 16)
  abline(a= 0, b = 1,lty = 2)
  
  plot(x = s1_s2_Z_use,
       y = t_Z,
       xlim = c(-15,15),
       ylim = c(0,1),
       type = "p",
       cex = 0.5,
       xlab = "Z",
       ylab = "ttest pvalue/lfdr",
       pch = 16)
  lines(x = s1_s2_Z_use,
        y = lfdr,
        type = "p",
        cex = 0.5,
        col = "red",
        pch = 16)
  abline(v =0,lty = 2)
  
  
  dev.off()
  
  
  result_df = data.frame(name =s1_s2_use[[3]],
                         fc = s1_s2_lfc,
                         Z = s1_s2_Z_use,
                         p1 = Z_p1,
                         p1_corrected = p1_corrected,
                         fdr = 1-p1_corrected,
                         stringsAsFactors = F)%>%
    dplyr::arrange(desc(p1))
  
  result_df_name = paste0(working_dir, compare_name,"_Z_result.tsv")
  
  write.table(result_df, result_df_name,
              quote = F, row.names = F, sep = "\t")
  
  
  f1_dens  = rep(0, length(Z_dens$x))
  
  for(i in 1:length(Z_dens$x))
  {
    
    this_value = Z_dens$x[i]
    this_dens = Z_dens$y[i]
    
    ### find the most similar in z_null_dens 
    dis_this = abs(this_value-z_null_dens$x)
    nearest = which(dis_this == min(dis_this))[1]
    this_null_dens = z_null_dens$y[nearest]
    
    f1_dens[i] = (this_dens-s1_s2_p0_use*this_null_dens)/(1-s1_s2_p0_use)  ### here is problematic in the 2d scenario, check this 
    
  }
   f1_dens[which(f1_dens<0)] = 0 
  
  
  overlay_dens_pdf_name =  paste0(working_dir, compare_name,"_dens_overlay.pdf")
  ylimit = max(c(max(z_null_dens$y),Z_dens$y, f1_dens))
  
  pdf(overlay_dens_pdf_name,useDingbats = F)
  plot(x = Z_dens$x, y = f1_dens, type = "l",ylim = c(0,ylimit), xlim = c(-10,10), col = "red")
  lines(z_null_dens, col = "green")
  lines(Z_dens)
  dev.off()
  
  
  #### new overlay
  pdf( paste0(working_dir, compare_name,"_dens_overlay_prop.pdf"))
  
  this_y_lim = max(c(Z_dens$y, s1_s2_p0_use * z_null_dens$y,(1-s1_s2_p0_use)*f1_dens))
  
  plot(x = Z_dens$x, y = Z_dens$y,ylim = c(0,this_y_lim), xlim = c(-10,10), type = "l")
  lines(x = z_null_dens$x, y= s1_s2_p0_use * z_null_dens$y, col = "green")
  lines(x = Z_dens$x, y = (1-s1_s2_p0_use)*f1_dens, type = "l", col = "red")
  
  dev.off()
  
  
  
  cat("posterior calculated.", "\n")
  
  
  return(s1_s2_p0_use)
  
  
  
  
}






comparison_groups_afterTimeLog = function(s1_g1_col_name,
                                          s2_g1_col_name,
                                          s1_g2_col_name,
                                          s2_g2_col_name,
                                          d1_data,
                                          nna_cutoff,
                                          all_seeds,
                                          permute_times,
                                          working_dir,
                                          compare_name)
  
   

{
  
  

  
  s1_g1_d1 = d1_data[,s1_g1_col_name]
  s2_g1_d1 = d1_data[,s2_g1_col_name]
  
  s1_g2_d1 = d1_data[,s1_g2_col_name]
  s2_g2_d1 = d1_data[,s2_g2_col_name]
  
  ##### figure out the complete NA for g1 and g2 
  s1_g1_na = apply(s1_g1_d1, 1, noNA)
  s2_g1_na = apply(s2_g1_d1, 1, noNA)
  
  g1_na_name = d1_data[which(s1_g1_na ==0 & s2_g1_na == 0),1]
  
  s1_g2_na = apply(s1_g2_d1, 1, noNA)
  s2_g2_na = apply(s2_g2_d1, 1, noNA)
  
  g2_na_name = d1_data[which(s1_g2_na ==0 & s2_g2_na == 0),1]
  
  
  
  
  #### filter data
  g1_s1_s2_use = s1_s2_filter_1d(
    s1_d1 = s1_g1_d1,
    s2_d1 = s2_g1_d1,
    names = d1_data[,1],
    nna_cutoff = nna_cutoff
  )
  
  
  
  g1_s1_d1_use = g1_s1_s2_use[[1]]
  g1_s2_d1_use = g1_s1_s2_use[[2]]
  
  
  ### this is to be compared between g1 and g2 
  
  g1_fc_d1_use = g1_s1_d1_use - g1_s2_d1_use
  ################
  
  
  
  g2_s1_s2_use = s1_s2_filter_1d(
    s1_d1 = s1_g2_d1,
    s2_d1 = s2_g2_d1,
    names = d1_data[,1],
    nna_cutoff = nna_cutoff
  )
  
  g2_s1_d1_use = g2_s1_s2_use[[1]]
  g2_s2_d1_use = g2_s1_s2_use[[2]]
  
  g2_fc_d1_use = g2_s1_d1_use - g2_s2_d1_use
  
  ### get the names in common
  g1_name = data.frame(name =g1_s1_s2_use[[3]], g1_rows = c(1:nrow(g1_fc_d1_use)), stringsAsFactors = F ) 
  g2_name = data.frame(name =g2_s1_s2_use[[3]], g2_rows = c(1:nrow(g2_fc_d1_use)), stringsAsFactors = F ) 
  
  g1_not_g2 = intersect(g1_s1_s2_use[[3]], g2_na_name)
  
  write.table(g1_not_g2, paste0(working_dir, compare_name,"_not_appear_g2.tsv"),
              sep = "\t", quote = F, col.names = F, row.names = F)
  
  g2_not_g1 = intersect(g2_s1_s2_use[[3]], g1_na_name)
  

  write.table(g2_not_g1, paste0(working_dir, compare_name,"_not_appear_g1.tsv"),
              sep = "\t", quote = F, col.names = F, row.names = F)
  
  
  common_names = g1_name%>%
    dplyr::left_join(g2_name, by = "name")%>%
    na.omit()
  
  
  g1_common_fc = g1_fc_d1_use[common_names$g1_rows,]
  g2_common_fc = g2_fc_d1_use[common_names$g2_rows,]
  
  ##################################################################################################
  ##################################################################################################
  ##################################################################################################
  
  g1_g2_t = rep(1, nrow(common_names))
  g1_g2_fc = rep(0, nrow(common_names))
  
  
  for(i in 1:nrow(common_names))
  {
    
    this_test = t.test(g1_common_fc[i,], g2_common_fc[i,])
    
    g1_g2_fc[i] = mean(g1_common_fc[i,], na.rm = T)- mean(g2_common_fc[i,], na.rm = T)
    g1_g2_t[i] = this_test$p.value
    
    
    
  }
  
  ### do BH correction
  p_adjust = p.adjust(g1_g2_t, method = "BH")
  
  g1_g2_ttest_df = data.frame(name = common_names$name, 
                              fc = g1_g2_fc,
                              p = g1_g2_t,
                              p_adjust = p_adjust,
                              stringsAsFactors = F)
  
  write.table(g1_g2_ttest_df, paste0(working_dir, compare_name,"_ttest_output.tsv"),
              quote = F, row.names = F, sep = "\t")
  
  pdf(paste0(working_dir, compare_name,"_ttest_vol.pdf"), useDingbats = F)
  plot(g1_g2_ttest_df$fc, -log10(g1_g2_ttest_df$p),
       type = "p", pch = 16, 
       cex = 0.5)
  dev.off()
  ##################################################################################################
  ##################################################################################################
  ##################################################################################################
  
  #### get other functions ready 
  
  
  s1_s2_Z_use = generate_twoSample_tstat_1d(s1_df1 = g1_common_fc,
                                            s2_df1 = g2_common_fc)
  
  s1_s2_lfc = compute_log2foldchage(s1_df1 = g1_common_fc,
                                    s2_df1 = g2_common_fc)
  
  
  
  Z_pdf_name = paste0(working_dir,compare_name, "_Z_hist.pdf")
  
  pdf(Z_pdf_name,useDingbats = F)
  hist(s1_s2_Z_use, breaks  = 100)
  dev.off()
  cat("Z generated.", "\n")
  
  
  
  null_permute = generate_permute_cols(s1_size = ncol(g1_common_fc),
                                       s2_size = ncol(g2_common_fc),
                                       numPermu = permute_times,
                                       seeds = all_seeds)
  
  
  
  
  s1_s2_z_null_use = lapply(1:permute_times, function(x)
  {
    
    
    
    this_shift = null_permute[x,]
    
    
    this_g1_d1 = g1_common_fc
    this_g2_d1 = g2_common_fc
    
    this_g1_d1[,this_shift] = g2_common_fc[,this_shift]
    this_g2_d1[,this_shift] = g1_common_fc[,this_shift]
    
    
    
    permute_use = s1_s2_filter_1d(s1_d1 = this_g1_d1,
                                  s2_d1 = this_g2_d1,
                                  names = common_names$name,
                                  nna_cutoff = 6)
    
    
    this_z = generate_twoSample_tstat_1d(s1_df1 = permute_use[[1]],
                                         s2_df1 = permute_use[[2]])
    if(x%%100 == 0)
      cat(x, "\n")
    
    return(this_z)
    
  })
  
  
  
  
  s1_s2_z_null_use_cat =  Reduce(c, s1_s2_z_null_use)
  
  z_null_pdf_name = paste0(working_dir,compare_name, "_z_null_hist.pdf")
  pdf(z_null_pdf_name,useDingbats = F)
  
  hist(s1_s2_z_null_use_cat, breaks  = 100)
  dev.off()
  
  
  
  Z_dens = density(s1_s2_Z_use, n = 2000)
  z_null_dens = density(s1_s2_z_null_use_cat, n = 2000)
  
  
  Zz_dens_pdf_name = paste0(working_dir,compare_name, "_Zz_dens_more.pdf")
  ylimit = max(z_null_dens$y)
  pdf(Zz_dens_pdf_name,useDingbats = F)
  plot(z_null_dens, col = "green", ylim = c(0,ylimit), xlim = c(-10,10))
  lines(Z_dens)
  dev.off()
  
  cat("z null generated.", "\n")
  
  
  
  
  
  Z_zero_which = which(abs(Z_dens$x-0) == min(abs(Z_dens$x-0)))[1]
  Z_zero_dens = Z_dens$y[Z_zero_which]
  
  z_null_zero_which = which(abs(z_null_dens$x-0) == min(abs(z_null_dens$x-0)))[1]
  z_null_zero_dens = z_null_dens$y[z_null_zero_which]
  
  
  s1_s2_p0_use = Z_zero_dens/z_null_zero_dens
  ### calculate f/f0 for each Z 
  
  
  f02f = rep(0, length(s1_s2_Z_use))
  
  for(i in 1:length(s1_s2_Z_use))
  {
    this_value = s1_s2_Z_use[i]
    
    ### find the most similar in z_null_dens 
    dis_this = abs(this_value-Z_dens$x)
    nearest = which(dis_this == min(dis_this))[1]
    this_Z_f = Z_dens$y[nearest]
    
    dis_this_null = abs(this_value-z_null_dens$x)
    nearest_null = which(dis_this_null == min(dis_this_null))[1]
    this_Z_f0 = z_null_dens$y[nearest_null]
    
    f02f[i] = this_Z_f0/this_Z_f
    
  }
  
  Z_p0 = s1_s2_p0_use*f02f
  Z_p1 = 1-Z_p0
  
  p1_corrected = Z_p1
  p1_corrected[which(Z_p1<0)]=0
  lfdr = 1-p1_corrected
  
  ### new overlay 
  
  pc = rep(rgb(0.1,0,0,0.5), length(s1_s2_Z_use))
  pc[which( p1_corrected>0.9)] = rgb(0.5,0,0,0.5)
  
  pc[which(s1_s2_Z_use>0 & p1_corrected>0.95)] = rgb(1,0,0,0.5)
  pc[which(s1_s2_Z_use<0 & p1_corrected>0.95)] = rgb(1,0.5,0,0.5)
  
  posterior_pdf_name = paste0(working_dir,compare_name, "_Z_posterior.pdf")
  
  pdf(posterior_pdf_name,useDingbats = F)
  
  plot(x = s1_s2_Z_use,
       y = p1_corrected,
       ylim = c(0,1),
       type = "p",
       col = pc,
       pch = 16)
  dev.off()
  

  ### a scatter plot 
  sanity_pdf_name = paste0(working_dir,compare_name, "_ttest_check.pdf")
  
  ### when plot remove the outlier ones 
  
  #leftones = which(s1_s2_Z_use< 0)
  
  pdf(sanity_pdf_name,useDingbats = F)
  
  plot(x = g1_g2_ttest_df$p,
       y = lfdr,
       xlim = c(0,1),
       ylim = c(0,1),
       type = "p",
       cex = 0.5,
       xlab = "ttest p value",
       ylab = "local fdr",
       pch = 16)
  abline(a= 0, b = 1,lty = 2)
  
  plot(x = s1_s2_Z_use,
       y = g1_g2_ttest_df$p,
       xlim = c(-15,15),
       ylim = c(0,1),
       type = "p",
       cex = 0.5,
       xlab = "Z",
       ylab = "ttest pvalue/lfdr",
       pch = 16)
  lines(x = s1_s2_Z_use,
        y = lfdr,
        type = "p",
        cex = 0.5,
        col = "red",
        pch = 16)
  abline(v =0,lty = 2)
  
  
  dev.off()
  
  
  result_df = data.frame(name =common_names$name,
                         fc = s1_s2_lfc,
                         pvalue = g1_g2_ttest_df$p,
                         Z = s1_s2_Z_use,
                         p1 = Z_p1,
                         p1_corrected = p1_corrected,
                         fdr = 1-p1_corrected,
                         stringsAsFactors = F)%>%
    dplyr::arrange(desc(p1))
  
  result_df_name = paste0(working_dir, compare_name,"_Z_result.tsv")
  
  write.table(result_df, result_df_name,
              quote = F, row.names = F, sep = "\t")
  
  
  f1_dens  = rep(0, length(Z_dens$x))
  
  for(i in 1:length(Z_dens$x))
  {
    
    this_value = Z_dens$x[i]
    this_dens = Z_dens$y[i]
    
    ### find the most similar in z_null_dens 
    dis_this = abs(this_value-z_null_dens$x)
    nearest = which(dis_this == min(dis_this))[1]
    this_null_dens = z_null_dens$y[nearest]
    
    f1_dens[i] = (this_dens-s1_s2_p0_use*this_null_dens)/(1-s1_s2_p0_use)  ### here is problematic in the 2d scenario, check this 
    
  }
   f1_dens[which(f1_dens<0)] = 0 
  
  
  overlay_dens_pdf_name =  paste0(working_dir, compare_name,"_dens_overlay.pdf")
  ylimit = max(c(max(z_null_dens$y),Z_dens$y, f1_dens))
  
  pdf(overlay_dens_pdf_name,useDingbats = F)
  plot(x = Z_dens$x, y = f1_dens, type = "l",ylim = c(0,ylimit), xlim = c(-10,10), col = "red")
  lines(z_null_dens, col = "green")
  lines(Z_dens)
  dev.off()
  
  
  #### new overlay
  pdf( paste0(working_dir, compare_name,"_dens_overlay_prop.pdf"))
  
  this_y_lim = max(c(Z_dens$y, s1_s2_p0_use * z_null_dens$y,(1-s1_s2_p0_use)*f1_dens))
  
  plot(x = Z_dens$x, y = Z_dens$y,ylim = c(0,this_y_lim), xlim = c(-10,10), type = "l")
  lines(x = z_null_dens$x, y= s1_s2_p0_use * z_null_dens$y, col = "green")
  lines(x = Z_dens$x, y = (1-s1_s2_p0_use)*f1_dens, type = "l", col = "red")
  
  dev.off()
  
  
  
  cat("posterior calculated.", "\n")
  
  return(s1_s2_p0_use)
  
}








