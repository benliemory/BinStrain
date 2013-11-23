Whole_est <-
function(results,data)
{
  require(quadprog)
  #results = Mix1_results
  #data = Mix1_dealed
  beta = rep(-1,dim(data$Mat_Zs)[2])
  for (i in 1:length(results$para_means))
  { 
    if (max(results$para_means[[i]]$p_est)<0.05)
    {beta[unique(results$para_means[[i]]$vari_name)]=0}
    
  }
  
  #plot(results$para_means[[1]]$vari_name,
  #     results$para_means[[1]]$p_est,xlim=c(1,32),ylim=c(0,1))

  
  Mat_Zs = as.data.frame(data$Mat_Zs)
  Z_new0 = Mat_Zs[,which(beta==-1),drop = F]
  Z_new = Z_new0[which(rowSums(Z_new0) ==1),,drop=F]
  
  index = as.numeric(row.names(Z_new))
  n = as.numeric(as.character((data$Sample_num_reads)))
  x = as.numeric(data$Sample_num_allele)
  pi = x/n
  pi = pi[index]
  
  all.names = names(Z_new)
  v.names = sapply(1:dim(Z_new)[1], function(i) all.names[which(Z_new[i,] ==1)])
  aa = strsplit(v.names,split="V")
  v.names = sapply(1:length(aa),function(i) aa[[i]][2]) 
  
  est_means = data.frame(var.name = v.names, pi.est = pi)
  
  
  para_means = lapply(unique(est_means$var.name),function(x)
    est_means[which(est_means$var.name == x),])
  
  for (i in 1:length(para_means))
  { 
    if (max(para_means[[i]]$pi.est)<0.05)
    { index0 = as.numeric(as.character(unique(para_means[[i]]$var.name)))
      beta[index0]=0}
    
  }
  
  
  X = Mat_Zs[,which(beta==-1),drop=F]
  index11 = !rowSums(X)==0
  X = X[index11,,drop=F]
  pi = x/n
  pi = pi[index11]
  aa = lm.restricted(pi,as.matrix(X))
  
  beta[which(beta== -1)] =aa
  
  beta
}
