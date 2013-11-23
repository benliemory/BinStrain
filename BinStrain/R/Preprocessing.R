Preprocessing <-
function (Counts_Sample, SNP_Pattern, N = nrow(SNP_Pattern))
  {
    IDs = SNP_Pattern$V1
    
    Sample = Counts_Sample[which(Counts_Sample$Value %in% IDs), ]
    
    Sample_num_reads = Sample$Coverage
    
    
    ### if have no "coverage"
    #Sample_num_reads = rowSums(Sample[,4:9])
    
    names(Sample)[4:7]=c("A","G","C","T")
    
    
    
    
    
    Mat_Zs = matrix(0, nrow = dim(Sample)[1], ncol = dim(SNP.Pattern)[2]-3)
    
    if (dim(Sample)[1] == N) 
    {   Missing = NA
        Zs = SNP_Pattern[,4:dim(SNP_Pattern)[2]]
        for (i in 1:N)
        { index = as.numeric(Zs[i,])
          index = index[!is.na(index)]
          Mat_Zs[i,index-1] = 1
        }
    }  else
      
    {
      Missing = IDs[which(!(IDs %in% Sample$Value ) == T)]
      Zs = SNP_Pattern[-which(SNP_Pattern$V1 %in% Missing),4:dim(SNP_Pattern)[2]]
      for (i in 1:dim(Sample)[1])
      { index = as.numeric(Zs[i,])
        index = index[!is.na(index)]
        Mat_Zs[i,index-1] = 1
      }                          
    }
    
    Sample_num_allele = numeric(dim(Sample)[1])
    if (all(is.na(Missing)) == FALSE)
    {SNP_Pattern_del_missing = SNP_Pattern[-which(SNP_Pattern$V1 %in% Missing),]
     Sample_num_allele = sapply(1:dim(Sample)[1], function(x) 
     {Sample[[as.character(SNP_Pattern_del_missing$V3)[x]]][x]})
    } else
      Sample_num_allele = sapply(1:dim(Sample)[1], function(x) 
      {Sample[[as.character(SNP_Pattern$V3)[x]]][x]})
    
    
    
    
    list(Sample = Sample, Sample_num_reads = Sample_num_reads,
         Sample_num_allele = Sample_num_allele, Mat_Zs = Mat_Zs,
         Missing = Missing   
    )
    
  }
