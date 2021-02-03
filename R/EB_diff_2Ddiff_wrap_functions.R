



#' Compare between 2 conditions (paired) in 2D
#'
#' @param s1_col_name names of the columns in the first condition 
#' @param s2_col_name names of the columns in the second condition
#' @param d1_data dataframe for the data of the first dimension, rows are genes, columns are patients in 2 conditions
#' @param d2_data dataframe for the data of the second dimension, rows are genes, columns are patients in 2 conditions
#' @param nna_cutoff genes with values in at least the number of samples are included in comparison
#' @param all_seeds seeds for the permutation
#' @param permute_times number of permutations 
#' @param working_dir the directory output files will be deposited in 
#' @param compare_name a label for the anlaysis 
#' @import dplyr data.table magrittr splines stringr MASS scatterplot3d 
#' @keywords 2D
#' @export
#' @examples
#' camparison_time_points_2dZ



comparison_time_points_2dZ = function(s1_col_name,
                                      s2_col_name,
                                      d1_data,
                                      d2_data,
                                      nna_cutoff,
                                      all_seeds,
                                      permute_times,
                                      working_dir,
                                      compare_name)
{
  
  d1_data = data.frame(d1_data)
  d2_data = data.frame(d2_data)
  
  
  ### get the data for comparison
  s1_d1 = d1_data[,s1_col_name]
  s1_d2 = d2_data[,s1_col_name]
  s2_d1 = d1_data[,s2_col_name]
  s2_d2 = d2_data[,s2_col_name]
  
  #### filter data
  
  s1_s2_use = s1_s2_filter(s1_d1 = s1_d1,
                           s1_d2 = s1_d2,
                           s2_d1 = s2_d1,
                           s2_d2 = s2_d2,
                           pair_names = d1_data[,1],
                           nna_cutoff = nna_cutoff)
  
  
  cat("pairs filtered.", "\n")
  
  
  s1_d1_use = s1_s2_use[[1]]
  s1_d2_use = s1_s2_use[[2]]
  s2_d1_use = s1_s2_use[[3]]
  s2_d2_use = s1_s2_use[[4]]
  
  
  ### now for each dimension build the Z score 
  
  
  s1_s2_Z_use_d1 = generate_oneSample_tstat_1d(s1_df1 = s1_d1_use,
                                               s2_df1 = s2_d1_use) 
  
  s1_s2_Z_use_d2 = generate_oneSample_tstat_1d(s1_df1 = s1_d2_use,
                                               s2_df1 = s2_d2_use) 
  
  s1_s2_lfc_d1 = compute_log2foldchage(s1_df1 = s1_d1_use,
                                       s2_df1 = s2_d1_use)
  
  s1_s2_lfc_d2 = compute_log2foldchage(s1_df1 = s1_d2_use,
                                       s2_df1 = s2_d2_use)

  
  null_permute =  shuffle_time_points( sample_size = ncol(s1_d1_use),
                         numPermu = permute_times,
                         seeds = all_seeds)
  

  s1_s2_z_null_use = lapply(1:permute_times, function(x)
  {
    
    this_assign = null_permute[[x]]
    
    which_flip = which(this_assign == 2)
    
    this_c1_d1 = s1_d1_use
    this_c1_d2 = s1_d2_use
    this_c2_d1 = s2_d1_use
    this_c2_d2 = s2_d2_use
    
    this_c1_d1[,which_flip] = s2_d1_use[,which_flip]
    this_c2_d1[,which_flip] = s1_d1_use[,which_flip]
    
    this_c1_d2[,which_flip] = s2_d2_use[,which_flip]
    this_c2_d2[,which_flip] = s1_d2_use[,which_flip]
    
    
    permute_use = s1_s2_filter(s1_d1 = this_c1_d1,
                               s1_d2 = this_c1_d2,
                               s2_d1 = this_c2_d1,
                               s2_d2 = this_c2_d2,
                               pair_names = s1_s2_use[[5]],
                               nna_cutoff = nna_cutoff)
    
    
    this_z_d1 = generate_oneSample_tstat_1d(s1_df1 = permute_use[[1]],
                                            s2_df1 = permute_use[[3]]) 
    
    this_z_d2 = generate_oneSample_tstat_1d(s1_df1 = permute_use[[2]],
                                            s2_df1 = permute_use[[4]]) 
    
    
    this_z= cbind(this_z_d1, this_z_d2)
    
    if(x%%100 == 0)
      cat(x, "\n")
    
    return(this_z)
    
  })
  
  

  
  s1_s2_z_null_use_cat =  Reduce(rbind, s1_s2_z_null_use)
  
 
  
  Z_dens = kde2d(s1_s2_Z_use_d1,
                 s1_s2_Z_use_d2,
                 n = 200)
  
  Z_x_dens =  density(s1_s2_Z_use_d1)
  Z_y_dens =  density(s1_s2_Z_use_d2)
  
  density_persp( pdf_name =  paste0(working_dir, compare_name,"_2d_Z_dens_function.pdf"),
                 x = Z_dens$x,
                 y = Z_dens$y,
                 z = Z_dens$z,
                 x_margin_dens_x = Z_x_dens$x,
                 x_margin_dens_y = Z_x_dens$y,
                 y_margin_dens_x = Z_y_dens$x,
                 y_margin_dens_y = Z_y_dens$y,
                 my_xlab = "kinase",
                 my_ylab = "substrate",
                 my_zlab = "density",
                 my_main = "Z",
                 my_theta = 120,
                 my_phi = 15,
                 my_shade = 0.5)
  

  #########
  
  
  z_null_dens = kde2d(s1_s2_z_null_use_cat[,1],
                      s1_s2_z_null_use_cat[,2],
                      n = 200) 
  
  z_null_x_dens =  density(s1_s2_z_null_use_cat[,1])
  z_null_y_dens =  density(s1_s2_z_null_use_cat[,2])
  
  density_persp( pdf_name =  paste0(working_dir, compare_name,"_2d_z_null_dens_function.pdf"),
                 x = z_null_dens$x,
                 y = z_null_dens$y,
                 z = z_null_dens$z,
                 x_margin_dens_x = z_null_x_dens$x,
                 x_margin_dens_y = z_null_x_dens$y,
                 y_margin_dens_x= z_null_y_dens$x,
                 y_margin_dens_y = z_null_y_dens$y,
                 my_xlab = "kinase",
                 my_ylab = "substrate",
                 my_zlab = "density",
                 my_main = "z",
                 my_theta = 120,
                 my_phi = 15,
                 my_shade = 0.5)
  
  

  ###########
  
  ### find the square which is closest to the 0 box 
  
  
  
  Z_zero_which_x = which(abs(Z_dens$x-0) == min(abs(Z_dens$x-0)))[1]
  Z_zero_which_y = which(abs(Z_dens$y-0) == min(abs(Z_dens$y-0)))[1]
  
  Z_zero_dens = Z_dens$z[Z_zero_which_x, Z_zero_which_y]
  
  z_null_zero_which_x = which(abs(z_null_dens$x-0) == min(abs(z_null_dens$x-0)))[1]
  z_null_zero_which_y = which(abs(z_null_dens$y-0) == min(abs(z_null_dens$y-0)))[1]
  
  z_null_zero_dens = z_null_dens$z[z_null_zero_which_x, z_null_zero_which_y]
  
  
  s1_s2_p0_use = Z_zero_dens/z_null_zero_dens
  ### calculate f/f0 for each Z 
  
  
  f02f = rep(0, length(s1_s2_Z_use_d1))
  
  for(i in 1:length(s1_s2_Z_use_d1))
  {
    this_value_x = s1_s2_Z_use_d1[i]
    this_value_y = s1_s2_Z_use_d2[i]
    
    ### find the most similar in z_null_dens 
    dis_this_x = abs(this_value_x-Z_dens$x)
    nearest_x = which(dis_this_x == min(dis_this_x))[1]
    dis_this_y = abs(this_value_y-Z_dens$y)
    nearest_y = which(dis_this_y == min(dis_this_y))[1]
    
    this_Z_f = Z_dens$z[nearest_x, nearest_y]
    
    dis_this_null_x = abs(this_value_x-z_null_dens$x)
    nearest_null_x = which(dis_this_null_x == min(dis_this_null_x))[1]
    dis_this_null_y = abs(this_value_y-z_null_dens$y)
    nearest_null_y = which(dis_this_null_y == min(dis_this_null_y))[1]
    
    
    this_Z_f0 = z_null_dens$z[nearest_null_x, nearest_null_y]
    
    f02f[i] = this_Z_f0/this_Z_f
    
  }
  
  
  
  Z_p0 = s1_s2_p0_use*f02f
  Z_p1 = 1-Z_p0
  
  
  
  
  p1_corrected = Z_p1
  p1_corrected[which(Z_p1<0)]=0
  lfdr = 1-p1_corrected
  
  
  posterior_pdf_name = paste0(working_dir,compare_name, "_2d_Z_posterior.pdf")
  
  ### when plot remove the outlier ones 
  pdf(posterior_pdf_name,useDingbats = F)
  
  get_col = rep(rgb(0.1,0,0,0.5), length(p1_corrected))
  get_col[which(p1_corrected>0.9)] = rgb(0.5,0,0,0.5)
  get_col[which(p1_corrected>0.99)] = rgb(1,0.5,0,0.5)
  
  get_col[which(p1_corrected>0.99 & s1_s2_Z_use_d1>0 & s1_s2_Z_use_d2>0)] = rgb(1,0,0,0.5)
  
  #pw = which(Z_p1>0.99)
  
  
  scatterplot3d(x = s1_s2_Z_use_d1,
                y = s1_s2_Z_use_d2,
                z = p1_corrected,
                # type = "h",
                color = get_col,
                # box = F,
                pch = 16,
                xlab = "kinase",
                ylab = "subtrate",
                zlab = "posterior")
  
  dev.off()
  
  
  
  

  
  result_df = data.frame(name =s1_s2_use[[5]],
                         Z_kinase = s1_s2_Z_use_d1,
                         Z_substrate = s1_s2_Z_use_d2,
                         fc_kinase = s1_s2_lfc_d1,
                         fc_substrate = s1_s2_lfc_d2,
                         p1 = Z_p1,
                         p1_corrected = p1_corrected,
                         fdr = 1-p1_corrected,
                         stringsAsFactors = F)%>%
    dplyr::arrange(desc(p1))
  
  result_df_name = paste0(working_dir, compare_name,"_Z_result.tsv")
  
  write.table(result_df, result_df_name,
              quote = F, row.names = F, sep = "\t")
  #
  
  
  
  f1_dens  = Z_dens$z
  
  for(i in 1:length(Z_dens$x))
  {
    for(j in 1:length(Z_dens$y))
    {
      this_value_x = Z_dens$x[i]
      this_value_y  = Z_dens$y[j]
      this_dens = Z_dens$z[i,j]
      
      ### find the most similar in z_null_dens 
      dis_this_x = abs(this_value_x-z_null_dens$x)
      nearest_x = which(dis_this_x == min(dis_this_x))[1]
      dis_this_y = abs(this_value_y-z_null_dens$y)
      nearest_y = which(dis_this_y == min(dis_this_y))[1]
      
      this_null_dens = z_null_dens$z[nearest_x, nearest_y]
      
      
      f1_dens[i,j] = (this_dens-s1_s2_p0_use*this_null_dens)/(1-s1_s2_p0_use) 
    }
    
    
  }
 
    f1_dens[which(f1_dens<0)] = 0 
  
  
  f1_x_marginal = rowSums(f1_dens)/sum( rowSums(f1_dens))
  f1_y_marginal = colSums(f1_dens)/sum(colSums(f1_dens))
  
  
  density_persp( pdf_name =  paste0(working_dir, compare_name,"_2d_f1_function.pdf"),
                 x = Z_dens$x,
                 y = Z_dens$y,
                 z = f1_dens,
                 x_margin_dens_x = Z_dens$x,
                 x_margin_dens_y = f1_x_marginal,
                 y_margin_dens_x= Z_dens$y,
                 y_margin_dens_y = f1_y_marginal,
                 my_xlab = "kinase",
                 my_ylab = "substrate",
                 my_zlab = "density",
                 my_main = "f1",
                 my_theta = 120,
                 my_phi = 15,
                 my_shade = 0.5)
  
  

  cat("posterior calculated.", "\n")
  
  
  return(s1_s2_p0_use)
  
  
  
}



####



#' Compare between 2 conditions (paired) in 2D different conditions in different dimemsions 
#'
#' @param s1_d1_col_name names of the columns in the first condition in the first dimension data  
#' @param s2_d1_col_name names of the columns in the second condition in the first dimension data
#' @param s1_d2_col_name names of the columns in the first condition in the second dimension data  
#' @param s2_d2_col_name names of the columns in the second condition in the second dimension data
#' @param d1_data dataframe for the data of the first dimension, rows are genes, columns are patients in 2 conditions
#' @param d2_data dataframe for the data of the second dimension, rows are genes, columns are patients in 2 conditions
#' @param nna_cutoff genes with values in at least the number of samples are included in comparison
#' @param all_seeds seeds for the permutation
#' @param permute_times number of permutations 
#' @param working_dir the directory output files will be deposited in 
#' @param compare_name a label for the anlaysis 
#' @import dplyr data.table magrittr splines stringr MASS scatterplot3d 
#' @keywords 2D no condition/dimension requirement 
#' @export
#' @examples
#' camparison_time_points_2dZ_5_30

comparison_time_points_2dZ_5_30 = function(s1_d1_col_name,
                                           s2_d1_col_name,
                                           s1_d2_col_name,
                                           s2_d2_col_name,
                                           d1_data,
                                           d2_data,
                                           nna_cutoff,
                                           all_seeds,
                                           permute_times,
                                           compare_name,
                                           working_dir)
{
  
  
  d1_data = data.frame(d1_data)
  d2_data = data.frame(d2_data)
  
  
  ### get the data for comparison
  s1_d1 = d1_data[,s1_d1_col_name]
  s2_d1 = d1_data[,s2_d1_col_name]
  s1_d2 = d2_data[,s1_d2_col_name]
  s2_d2 = d2_data[,s2_d2_col_name]
  
  #### filter data
  
  s1_s2_use = s1_s2_filter(s1_d1 = s1_d1,
                           s1_d2 = s1_d2,
                           s2_d1 = s2_d1,
                           s2_d2 = s2_d2,
                           pair_names = d1_data[,1],
                           nna_cutoff = nna_cutoff)
  
  
  cat("pairs filtered.", "\n")
  
  
  s1_d1_use = s1_s2_use[[1]]
  s1_d2_use = s1_s2_use[[2]]
  s2_d1_use = s1_s2_use[[3]]
  s2_d2_use = s1_s2_use[[4]]
  
  
  ### now for each dimension build the Z score 
  
  
  s1_s2_Z_use_d1 = generate_oneSample_tstat_1d(s1_df1 = s1_d1_use,
                                               s2_df1 = s2_d1_use) 
  
  s1_s2_Z_use_d2 = generate_oneSample_tstat_1d(s1_df1 = s1_d2_use,
                                               s2_df1 = s2_d2_use) 
  
  s1_s2_lfc_d1 = compute_log2foldchage(s1_df1 = s1_d1_use,
                                       s2_df1 = s2_d1_use)
  
  s1_s2_lfc_d2 = compute_log2foldchage(s1_df1 = s1_d2_use,
                                       s2_df1 = s2_d2_use)
  
 
  null_permute = 
    shuffle_time_points( sample_size = ncol(s1_d1_use),
                         numPermu = permute_times,
                         seeds = all_seeds)
  
  
  
  s1_s2_z_null_use = lapply(1:permute_times, function(x)
  {
    
    # x = 1  
    this_assign = null_permute[[x]]
    
    #first = which(this_assign ==1)
    which_flip = which(this_assign == 2)
    
    this_c1_d1 = s1_d1_use
    this_c1_d2 = s1_d2_use
    this_c2_d1 = s2_d1_use
    this_c2_d2 = s2_d2_use
    
    this_c1_d1[,which_flip] = s2_d1_use[,which_flip]
    this_c2_d1[,which_flip] = s1_d1_use[,which_flip]
    
    this_c1_d2[,which_flip] = s2_d2_use[,which_flip]
    this_c2_d2[,which_flip] = s1_d2_use[,which_flip]
    
    
    permute_use = s1_s2_filter(s1_d1 = this_c1_d1,
                               s1_d2 = this_c1_d2,
                               s2_d1 = this_c2_d1,
                               s2_d2 = this_c2_d2,
                               pair_names = s1_s2_use[[5]],
                               nna_cutoff = nna_cutoff)
    
    
    this_z_d1 = generate_oneSample_tstat_1d(s1_df1 = permute_use[[1]],
                                            s2_df1 = permute_use[[3]]) 
    
    this_z_d2 = generate_oneSample_tstat_1d(s1_df1 = permute_use[[2]],
                                            s2_df1 = permute_use[[4]]) 
    
    
    this_z= cbind(this_z_d1, this_z_d2)
    
    if(x%%100 == 0)
      cat(x, "\n")
    
    return(this_z)
    
  })
  
  
  ### combine all the zs
  
  
  
  
  s1_s2_z_null_use_cat =  Reduce(rbind, s1_s2_z_null_use)
  
  
  Z_dens = kde2d(s1_s2_Z_use_d1,
                 s1_s2_Z_use_d2,
                 n = 200)
  
  Z_x_dens =  density(s1_s2_Z_use_d1)
  Z_y_dens =  density(s1_s2_Z_use_d2)
  
  density_persp( pdf_name =  paste0(working_dir, compare_name,"_2d_Z_dens_function.pdf"),
                 x = Z_dens$x,
                 y = Z_dens$y,
                 z = Z_dens$z,
                 x_margin_dens_x = Z_x_dens$x,
                 x_margin_dens_y = Z_x_dens$y,
                 y_margin_dens_x = Z_y_dens$x,
                 y_margin_dens_y = Z_y_dens$y,
                 my_xlab = "kinase",
                 my_ylab = "substrate",
                 my_zlab = "density",
                 my_main = "Z",
                 my_theta = 20,
                 my_phi = 15,
                 my_shade = 0.5)
  
  
  #########
  
  
  z_null_dens = kde2d(s1_s2_z_null_use_cat[,1],
                      s1_s2_z_null_use_cat[,2],
                      n = 200) 
  
  z_null_x_dens =  density(s1_s2_z_null_use_cat[,1])
  z_null_y_dens =  density(s1_s2_z_null_use_cat[,2])
  
  density_persp( pdf_name =  paste0(working_dir, compare_name,"_2d_z_null_dens_function.pdf"),
                 x = z_null_dens$x,
                 y = z_null_dens$y,
                 z = z_null_dens$z,
                 x_margin_dens_x = z_null_x_dens$x,
                 x_margin_dens_y = z_null_x_dens$y,
                 y_margin_dens_x= z_null_y_dens$x,
                 y_margin_dens_y = z_null_y_dens$y,
                 my_xlab = "kinase",
                 my_ylab = "substrate",
                 my_zlab = "density",
                 my_main = "z",
                 my_theta = 20,
                 my_phi = 15,
                 my_shade = 0.5)
  
  
  
  ###########
  
  ### find the square which is closest to the 0 box 
  
  
  
  Z_zero_which_x = which(abs(Z_dens$x-0) == min(abs(Z_dens$x-0)))[1]
  Z_zero_which_y = which(abs(Z_dens$y-0) == min(abs(Z_dens$y-0)))[1]
  
  Z_zero_dens = Z_dens$z[Z_zero_which_x, Z_zero_which_y]
  
  z_null_zero_which_x = which(abs(z_null_dens$x-0) == min(abs(z_null_dens$x-0)))[1]
  z_null_zero_which_y = which(abs(z_null_dens$y-0) == min(abs(z_null_dens$y-0)))[1]
  
  z_null_zero_dens = z_null_dens$z[z_null_zero_which_x, z_null_zero_which_y]
  
  
  s1_s2_p0_use = Z_zero_dens/z_null_zero_dens
  ### calculate f/f0 for each Z 
  
  
  f02f = rep(0, length(s1_s2_Z_use_d1))
  
  for(i in 1:length(s1_s2_Z_use_d1))
  {
    this_value_x = s1_s2_Z_use_d1[i]
    this_value_y = s1_s2_Z_use_d2[i]
    
    ### find the most similar in z_null_dens 
    dis_this_x = abs(this_value_x-Z_dens$x)
    nearest_x = which(dis_this_x == min(dis_this_x))[1]
    dis_this_y = abs(this_value_y-Z_dens$y)
    nearest_y = which(dis_this_y == min(dis_this_y))[1]
    
    this_Z_f = Z_dens$z[nearest_x, nearest_y]
    
    dis_this_null_x = abs(this_value_x-z_null_dens$x)
    nearest_null_x = which(dis_this_null_x == min(dis_this_null_x))[1]
    dis_this_null_y = abs(this_value_y-z_null_dens$y)
    nearest_null_y = which(dis_this_null_y == min(dis_this_null_y))[1]
    
    
    this_Z_f0 = z_null_dens$z[nearest_null_x, nearest_null_y]
    
    f02f[i] = this_Z_f0/this_Z_f
    
  }
  
  
  
  Z_p0 = s1_s2_p0_use*f02f
  Z_p1 = 1-Z_p0
  
  
  
  
  p1_corrected = Z_p1
  p1_corrected[which(Z_p1<0)]=0
  lfdr = 1-p1_corrected
  
  
  posterior_pdf_name = paste0(working_dir,compare_name, "_2d_Z_posterior.pdf")
  
  ### when plot remove the outlier ones 
  pdf(posterior_pdf_name,useDingbats = F)
  
  get_col = rep(rgb(0.1,0,0,0.5), length(p1_corrected))
  get_col[which(p1_corrected>0.9)] = rgb(0.5,0,0,0.5)
  get_col[which(p1_corrected>0.99)] = rgb(1,0.5,0,0.5)
  
  get_col[which(p1_corrected>0.99 & s1_s2_Z_use_d1>0 & s1_s2_Z_use_d2>0)] = rgb(1,0,0,0.5)
  
  #pw = which(Z_p1>0.99)
  
  
  scatterplot3d(x = s1_s2_Z_use_d1,
                y = s1_s2_Z_use_d2,
                z = p1_corrected,
                # type = "h",
                color = get_col,
                #box = F,
                pch = 16,
                angle = 45,
                xlab = "kinase",
                ylab = "subtrate",
                zlab = "posterior")
  
   #   
  dev.off()
  
  
  
  
  
  
  result_df = data.frame(name =s1_s2_use[[5]],
                         Z_kinase = s1_s2_Z_use_d1,
                         Z_substrate = s1_s2_Z_use_d2,
                         fc_kinase = s1_s2_lfc_d1,
                         fc_substrate = s1_s2_lfc_d2,
                         p1 = Z_p1,
                         p1_corrected = p1_corrected,
                         fdr = 1-p1_corrected,
                         stringsAsFactors = F)%>%
    dplyr::arrange(desc(p1))
  
  result_df_name = paste0(working_dir, compare_name,"_Z_result.tsv")
  
  write.table(result_df, result_df_name,
              quote = F, row.names = F, sep = "\t")
  #
  
  
  
  f1_dens  = Z_dens$z
  
  for(i in 1:length(Z_dens$x))
  {
    for(j in 1:length(Z_dens$y))
    {
      this_value_x = Z_dens$x[i]
      this_value_y  = Z_dens$y[j]
      this_dens = Z_dens$z[i,j]
      
      ### find the most similar in z_null_dens 
      dis_this_x = abs(this_value_x-z_null_dens$x)
      nearest_x = which(dis_this_x == min(dis_this_x))[1]
      dis_this_y = abs(this_value_y-z_null_dens$y)
      nearest_y = which(dis_this_y == min(dis_this_y))[1]
      
      this_null_dens = z_null_dens$z[nearest_x, nearest_y]
      
      
      f1_dens[i,j] = (this_dens-s1_s2_p0_use*this_null_dens)/(1-s1_s2_p0_use) 
    }
    
    
  }

    f1_dens[which(f1_dens==0)] = 0
  
  
  
  f1_x_marginal = rowSums(f1_dens)/sum( rowSums(f1_dens))
  f1_y_marginal = colSums(f1_dens)/sum(colSums(f1_dens))
  
  
  density_persp( pdf_name =  paste0(working_dir, compare_name,"_2d_f1_function.pdf"),
                 x = Z_dens$x,
                 y = Z_dens$y,
                 z = f1_dens,
                 x_margin_dens_x = Z_dens$x,
                 x_margin_dens_y = f1_x_marginal,
                 y_margin_dens_x= Z_dens$y,
                 y_margin_dens_y = f1_y_marginal,
                 my_xlab = "kinase",
                 my_ylab = "substrate",
                 my_zlab = "density",
                 my_main = "f1",
                 my_theta = 20,
                 my_phi = 15,
                 my_shade = 0.5)
  
  
  
  
  
  
  cat("posterior calculated.", "\n")
  
  
  return(s1_s2_p0_use)
  
  
  
}

