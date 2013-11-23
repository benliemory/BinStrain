Single_est <-
function (Counts_Sample, SNP_Pattern, N = nrow(SNP_Pattern),filename =NA)
  {
    Counts_sample_dealed = Preprocessing(Counts_Sample, SNP_Pattern, N)
    single_row = which(rowSums(Counts_sample_dealed$Mat_Zs) == 1)
    single_all_info = Counts_sample_dealed$Mat_Zs[single_row,]
    single_num_read = as.numeric(as.character(Counts_sample_dealed$Sample_num_reads[single_row] ))
    single_num_allele = as.numeric(Counts_sample_dealed$Sample_num_allele[single_row] )
    single_vari_name = sapply(1: nrow(single_all_info), function(x){which(single_all_info[x,]==1)})
    pi_est = single_num_allele/single_num_read
    single_est = data.frame(vari_name = single_vari_name, p_est = pi_est)
    single_names = single_vari_name[!duplicated(single_vari_name)]
    single_names = single_names[order(single_names)]
    
    para_est_mean = sapply(single_names,function(x)
    {mean(single_est[single_est$vari_name == x,]$p_est)})
    para_est_sd = sapply(single_names,function(x)
    {sd(single_est[single_est$vari_name == x,]$p_est)})
    
    
    para_details = lapply(single_names,function(x)
    {Counts_sample_dealed$Sample[as.numeric(row.names(single_est[single_est$vari_name == x,])),]})
    
    para_means = lapply(single_names,function(x)
    {single_est[as.numeric(row.names(single_est[single_est$vari_name == x,])),]})
    
    if (is.na(filename) == 0)
    { filename = as.character(paste(filename,".pdf",sep=""))
      pdf(filename)
      plot(single_est$vari_name,single_est$p_est,main = filename,
           xlab = "beta_i",ylab="beta_i_estimate")
      dev.off()
    }

    results=list(para_means = para_means ,
                 para_details = para_details,
                 para_est_mean = data.frame(vari_name = single_names,para_est_mean = para_est_mean) ,
                 para_est_sd =  data.frame(vari_name = single_names,para_est_sd = para_est_sd) 
                 
    )
  }
