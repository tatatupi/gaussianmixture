rm(list=ls(all=TRUE))
source("ScaMixNorFun.R")

# Tamanho amostral:
n = 500

# Parâmetros do modelo para a distribuição amostral:
mu = 11; s2 = 0.64; nu = 0.2

# Hiperparâmetros das distribuições a priori (somente pa-
# ra mu e sigma2):
m = 11; V = 1; a = 7; d = 4

# Gerando uma amostra de tamanho n hierarquicamente:

set.seed(122019)
u = U_pdf(n, nu)
sam = rnorm(n, mu, sqrt(s2*u))
hist(sam, breaks=100, prob=T)

# Em geral, se trabalha com o logaritmo da verossimilhan-
# ça, uma vez que o produtório é sempre muito baixo.

logA(sam, mu, s2, nu)
hist(logA(sam, mu, s2, nu), prob=T, breaks=50)
sum(logA(sam, mu, s2, nu))

logh(n, mu, s2, nu, m, V, a, d)

logkpost(sam, mu, s2, nu, m, V, a, d) #Soma as log-veros-
                                      #similhanças.

# Método SIR:

# Primeiramente, como não sabemos gerar da posteriori, i-
# remos considerar uma densidade conhecida, da qual sabe-
# mos gerar valores aleatórios (amostras).

# Escolhemos como densidade conhecida a normal multivari-
# ada, com k=3. Desta forma, como nem todos os argumentos
# da posteriori variam na reta real, devemos reparametri-
# zá-los de tal forma que todos os suportes dos novos pa-
# râmetros coincidam com os das médias na normal multiva-
# riada. Logo, a reparametrização será feita tomando fun-
# ções que levem sigma e nu de seus domínios para imagens
# na reta real (funções logarítmica e logística).

logkpost_re(sam, mu, s2, nu, m, V, a, d)

# Tamanhos para a amostragem da posteriori via SIR:

am1 = 500; am2 = 5000; am3 = 50000

# Geremos agora da distrib. normal multivariada com média
# de cada componente dada pelo respectivo valor no modelo
# e variância dada pela matriz identidade de ordem 3:

library(mvtnorm) # Pacote com funções densidade e gerado-
                 # ra da distrib. normal multivariada.

set.seed(122019)
r1_q = rmvnorm(am1, mean=c(11, log(0.64), log(0.2/0.8)),
               sigma=diag(3))
set.seed(122019)
r2_q = rmvnorm(am2, mean=c(11, log(0.64), log(0.2/0.8)),
               sigma=diag(3))
set.seed(122019)
r3_q = rmvnorm(am3, mean=c(11, log(0.64), log(0.2/0.8)),
               sigma=diag(3))

summary(r1_q[,1])
summary(r2_q[,1])
summary(r3_q[,1])
summary(exp(r1_q[,2]))
summary(exp(r2_q[,2]))
summary(exp(r3_q[,2]))
summary(1/(1 + exp(-r1_q[,3])))
summary(1/(1 + exp(-r2_q[,3])))
summary(1/(1 + exp(-r3_q[,3])))

d1_q = dmvnorm(r1_q, mean=c(11, log(0.64), log(0.2/0.8)),
               sigma=diag(3))
d2_q = dmvnorm(r2_q, mean=c(11, log(0.64), log(0.2/0.8)),
               sigma=diag(3))
d3_q = dmvnorm(r3_q, mean=c(11, log(0.64), log(0.2/0.8)),
               sigma=diag(3))

# Com amostras aleatórias e sua densidades podemos calcu-
# lar os pesos (probabilidade de inclusão) para cada pon-
# to amostral, os avaliando no núcleo da posteriori dese-
# jada para geração:

w1 = numeric(am1); w2 = numeric(am2); w3 = numeric(am3)
aux11 = numeric(am1); aux21 = numeric(am1)
aux12 = numeric(am2); aux22 = numeric(am2)
aux13 = numeric(am3); aux23 = numeric(am3)

for(i in 1:am1) {
  aux11[i] = logkpost_re(
    X = sam, mu = r1_q[i,1], sigma2 = exp(r1_q[i,2]),
    nu = 1/(1 + exp(-r1_q[i,3])), m, V, a, d)
}
aux21 = aux11/d1_q
w1 = aux21/sum(aux21); sum(w1)

for(i in 1:am2) {
  aux12[i] = logkpost_re(
    X = sam, mu = r2_q[i,1], sigma2 = exp(r2_q[i,2]),
    nu = 1/(1 + exp(-r2_q[i,3])), m, V, a, d)
}
aux22 = aux12/d2_q
w2 = aux22/sum(aux22); sum(w2)

for(i in 1:am3) {
  aux13[i] = logkpost_re(
    X = sam, mu = r3_q[i,1], sigma2 = exp(r3_q[i,2]),
    nu = 1/(1 + exp(-r3_q[i,3])), m, V, a, d)
}
aux23 = aux13/d3_q
w3 = aux23/sum(aux23); sum(w3)

set.seed(122019)
am1_mu_p = sample(x=r1_q[,1], size=am1,
                  replace=T, prob=w1)
am1_s2_p = sample(x=exp(r1_q[,2]), size=am1,
                  replace=T, prob=w1)
am1_nu_p = sample(x=1/(1+exp(-r1_q[,3])), size=am1,
                  replace=T, prob=w1)

set.seed(122019)
am2_mu_p = sample(x=r2_q[,1], size=am2,
                  replace=T, prob=w2)
am2_s2_p = sample(x=exp(r2_q[,2]), size=am2,
                  replace=T, prob=w2)
am2_nu_p = sample(x=1/(1+exp(-r2_q[,3])), size=am2,
                  replace=T, prob=w2)

set.seed(122019)
am3_mu_p = sample(x=r3_q[,1], size=am3,
                  replace=T, prob=w3)
am3_s2_p = sample(x=exp(r3_q[,2]), size=am3,
                  replace=T, prob=w3)
am3_nu_p = sample(x=1/(1+exp(-r3_q[,3])), size=am3,
                  replace=T, prob=w3)

library(moments) # Contém funções para extrair assimetri-
                 # a e curtose amostrais.

st_sam_post = function(Xpar) {
  medi = mean(Xpar)
  vari = var(Xpar)
  assi = skewness(Xpar)
  curt = kurtosis(Xpar)
  return(list(media=medi, variancia=vari,
              assimetria=assi, curtose=curt))
}

st_sam_post(am1_mu_p); hist(am1_mu_p, breaks=50, prob=T)
st_sam_post(am2_mu_p); hist(am2_mu_p, breaks=50, prob=T)
st_sam_post(am3_mu_p); hist(am3_mu_p, breaks=50, prob=T)

st_sam_post(am1_s2_p); hist(am1_s2_p, breaks=50, prob=T)
st_sam_post(am2_s2_p); hist(am2_s2_p, breaks=50, prob=T)
st_sam_post(am3_s2_p); hist(am3_s2_p, breaks=50, prob=T)

st_sam_post(am1_nu_p); hist(am1_nu_p, breaks=50, prob=T)
st_sam_post(am2_nu_p); hist(am2_nu_p, breaks=50, prob=T)
st_sam_post(am3_nu_p); hist(am3_nu_p, breaks=50, prob=T)
