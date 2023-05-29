library(tibble)
library(metan)

#mmn data
KEM_MMN<- read.csv("KEM_MMN.csv")
head(KEM_MMN)
KEM_MMN<- column_to_rownames(KEM_MMN, var = "child_no")

#spearman
#KEM_MMN_S<- KEM_MMN[,-c(18,19)]
#MS<- corr_coef(KEM_MMN_S, type = "linear", method = "spearman", use = "everything")
#plot(MS)

#KEM_MMN_S1<- KEM_MMN[,c(6,14,15)]
#MS1<- corr_coef(KEM_MMN_S1, type = "linear", method = "spearman", use = "everything")
#plot(MS1)

#pearson
KEM_MMN_P<- KEM_MMN[,-c(3,8,10,11,12,17,18,19)]
MP<- corr_coef(KEM_MMN_P, type = "linear", method = "pearson", use = "everything")
plot(MP)

KEM_MMN_P1<- KEM_MMN[,c(6,10,15,17)]
MP1<- corr_coef(KEM_MMN_P1, type = "linear", method = "pearson", use = "everything")
plot(MP1)

#pheno data
KEM_pheno<- read.csv("KEM_pheno.csv")
head(KEM_pheno)
KEM_pheno<- column_to_rownames(KEM_pheno, var = "child_no")

#spearman
#KEM_pheno_S<- KEM_pheno[,-c(1,2)]
#PS<- corr_coef(KEM_pheno_S, type = "linear", method = "spearman", use = "everything")
#plot(PS)

#KEM_pheno_S1<- KEM_pheno[,c(?)]
#PS1<- corr_coef(KEM_pheno_S1, type = "linear", method = "spearman", use = "everything")
#plot(PS1)

#pearson
KEM_pheno_P<- KEM_pheno[,-c(1,2)]
PP<- corr_coef(KEM_pheno_P, type = "linear", method = "pearson", use = "everything")
plot(PP)

KEM_pheno_P1<- KEM_pheno[,c(5,14,18)]
PP1<- corr_coef(KEM_pheno_P1, type = "linear", method = "pearson", use = "everything")
plot(PP1)
