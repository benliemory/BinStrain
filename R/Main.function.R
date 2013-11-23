Main.function <-
function(Sample,SNP.Pattern)
{
  Sample2 = Preprocessing(Sample, SNP.Pattern)
  Result2 = Single_est(Sample, SNP.Pattern)
  betas = Whole_est(Result2,Sample2)
  betas[betas<1e-10] = 0
  betas = betas/sum(betas,na.rm= T)
  write.table(betas, file="./betas.txt",col.names= F,quote= F)
}
