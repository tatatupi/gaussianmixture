rm(list=ls(all=TRUE))
library(MASS)
source("functions.R")

# Tamanho amostral:
n = 500

# ParÃ¢metros do modelo para a distribuiÃ§Ã£o amostral:
mu = 11; s2 = 0.64; nu = 0.2

# HiperparÃ¢metros das distribuiÃ§Ãµes a priori (somente pa-
# ra mu e sigma2):
m = 11; V = 1; a = 7; d = 4

# gerado amostras
sam = sampling(n, nu, mu, s2)


#definindo os vetores de média e matriz de covariância das variáveis
mu_vector = c(11, 0.6, 0.2)
cov_matrix = diag(1, 3, 3)


