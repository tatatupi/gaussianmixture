rm(list=ls(all=TRUE))
source("ScaMixNorFun.R")

# Tamanho amostral:
n = 500

# Parâmetros do modelo para a distribuição amostral:
mu = 11; s2 = 0.64; nu = 0.2

# Hiperparâmetros das distribuições a priori (somente pa-
# ra mu e sigma2):
m = 11; V = 1; a = 7; d = 4

# Distribuição de mu a priori (dados s2, m e V):
set.seed(122019)
s_mu = rnorm(n=5000, mean=m, sd=sqrt(V*s2))
hist(s_mu, prob=T, breaks=100)
median(s_mu); mean(s_mu); min(s_mu); max(s_mu)
# Massa probabilística concentrada quase totalmente entre
# 9 e 13 e aproximadamente centrada em 11.

# Distribuição de sigma2 a priori (dados a e d):
set.seed(122019)
s_s2 = 1/rgamma(n=5000, shape=a, rate=d)
hist(s_s2, prob=T, breaks=100)
median(s_s2); mean(s_s2); min(s_s2); max(s_s2)
# Massa probabilística concentrada quase totalmente entre
# 0.2 e 1.5 e aproximadamente centrada em 0.64.

# A priori, nu está distribuído uniformemente em (0,1).

# Gerando uma amostra de tamanho n hierarquicamente:

set.seed(122019)
u = U_pdf(n, nu)
sam = rnorm(n, mu, sqrt(s2*u))
par(mar = c(5,5,3,2))
hist(sam, breaks=100, prob=T, main="",
     xlim = c(-10, 40), ylim = c(0, 0.5),
     xlab = "Amostra do modelo",
     ylab = "Densidade empírica")

library(moments)
mean(sam); var(sam); skewness(sam); kurtosis(sam)

# Em geral, se trabalha com o logaritmo da verossimilhan-
# ça, uma vez que o produtório é sempre muito baixo.

logA(sam, mu, s2, nu)
hist(logA(sam, mu, s2, nu), prob=T, breaks=50)
sum(logA(sam, mu, s2, nu))

logh(n, mu, s2, nu, m, V, a, d)

logkpost(sam, mu, s2, nu, m, V, a, d)

# Método da Quadratura de Riemann:

L1 = 15  # Número de intervalos.
L2 = 50
L3 = 100

# Note que, a posteriori, não temos nenhuma ideia de como
# cada parâmetro está distribuído no seu respectivo espa-
# paço paramétrico. Conjuntamente, não conseguimos extra-
# ir amostras a posteriori dos 3 parâmetros.

# Para escolher o intervalo de integração associado a ca-
# da parâmetro, fixemos para os outros dois os seus valo-
# res correspondentes no modelo e calculemos o exponenci-
# al do logaritmo da posteriori (de modo a evitar proble-
# mas numéricos) como uma função do parâmetro para o qual
# queremos definir o intervalo. Na região em que a função
# tiver maior maior massa, é onde construiremos nosso in-
# tervalo de intergração.

# Começando por mu, temos que:

mu_sup = seq(9, 13, 0.01) # A massa posteriori pode estar
t_mu = length(mu_sup)     # concentrada ao menos no mesmo
kp_mu = numeric(t_mu)     # intervalo que a priori.

for(i in 1:t_mu) {
  kp_mu[i] = exp(logkpost(X=sam,
    mu=mu_sup[i], sigma2=s2, nu=nu, m=m, V=V, a=a, d=d))
}
plot(mu_sup, kp_mu, type="l")

mu_sup = seq(10.8, 11.2, 0.001) # Redefinindo.
t_mu = length(mu_sup)
kp_mu = numeric(t_mu)

for(i in 1:t_mu) {
  kp_mu[i] = exp(logkpost(X=sam,
    mu=mu_sup[i], sigma2=s2, nu=nu, m=m, V=V, a=a, d=d))
}
plot(mu_sup, kp_mu, type="l", main = "",
     xlab = expression(paste(mu)), ylab = "")

# Para mu o intervalo de integ. será de 10.85 a 11.13.

# Para sigma2, temos que:

s2_sup = seq(0.2, 1.5, 0.001)
t_s2 = length(s2_sup)
kp_s2 = numeric(t_s2)

for(i in 1:t_s2) {
  kp_s2[i] = exp(logkpost(X=sam,
    mu=mu, sigma2=s2_sup[i], nu=nu, m=m, V=V, a=a, d=d))
}
plot(s2_sup, kp_s2, type="l")

s2_sup = seq(0.4, 0.8, 0.001) # Redefinindo.
t_s2 = length(s2_sup)
kp_s2 = numeric(t_s2)

for(i in 1:t_s2) {
  kp_s2[i] = exp(logkpost(X=sam,
    mu=mu, sigma2=s2_sup[i], nu=nu, m=m, V=V, a=a, d=d))
}
plot(s2_sup, kp_s2, type="l", main = "",
     xlab = expression(paste(sigma^2)), ylab = "")

# Para sigma2 o intervalo de integ. será de 0.48 a 0.78.

# Por fim, para nu temos que:

nu_sup = seq(0, 1, 0.001)
t_nu = length(nu_sup)
kp_nu = numeric(t_nu)

for(i in 1:t_nu) {
  kp_nu[i] = exp(logkpost(X=sam,
    mu=mu, sigma2=s2, nu=nu_sup[i], m=m, V=V, a=a, d=d))
}
plot(nu_sup, kp_nu, type="l")

nu_sup = seq(0.1, 0.3, 0.001) #Redefinindo.
t_nu = length(nu_sup)
kp_nu = numeric(t_nu)

for(i in 1:t_nu) {
  kp_nu[i] = exp(logkpost(X=sam,
    mu=mu, sigma2=s2, nu=nu_sup[i], m=m, V=V, a=a, d=d))
}
plot(nu_sup, kp_nu, type="l", main = "",
     xlab = expression(paste(nu)), ylab = "")

# Para nu o intervalo de integ. será de 0.13 a 0.26.

# Observe que nenhum dos três intervalos escolhidos ficou
# centrado no mesmo valor assumido para o modelo. Mais um
# sinal da importância de considerar tanto a informação a
# priori quanto a dada pela amostra.

# Temos então os intervalos:

mu_step1 = (11.13 - 10.85)/L1
mu_step2 = (11.13 - 10.85)/L2
mu_step3 = (11.13 - 10.85)/L3
mu_gr1 = seq(10.85, 11.13, mu_step1)
mu_gr2 = seq(10.85, 11.13, mu_step2)
mu_gr3 = seq(10.85, 11.13, mu_step3)
mu_gr1; mu_gr2; mu_gr3

s2_step1 = (0.78 - 0.48)/L1
s2_step2 = (0.78 - 0.48)/L2
s2_step3 = (0.78 - 0.48)/L3
s2_gr1 = seq(0.48, 0.78, s2_step1)
s2_gr2 = seq(0.48, 0.78, s2_step2)
s2_gr3 = seq(0.48, 0.78, s2_step3)
s2_gr1; s2_gr2; s2_gr3

nu_step1 = (0.26 - 0.13)/L1
nu_step2 = (0.26 - 0.13)/L2
nu_step3 = (0.26 - 0.13)/L3
nu_gr1 = seq(0.13, 0.26, nu_step1)
nu_gr2 = seq(0.13, 0.26, nu_step2)
nu_gr3 = seq(0.13, 0.26, nu_step3)
nu_gr1; nu_gr2; nu_gr3

# E a grade:

grid_tri1 = cbind(mu_gr1, s2_gr1, nu_gr1); grid_tri1
grid_tri2 = cbind(mu_gr2, s2_gr2, nu_gr2); grid_tri2
grid_tri3 = cbind(mu_gr3, s2_gr3, nu_gr3); grid_tri3
l1 = nrow(grid_tri1); l1
l2 = nrow(grid_tri2); l2
l3 = nrow(grid_tri3); l3

# Como os tamanhos dos passos são iguais em todas as três
# dimensões, podemos calcular os produtos antes de inici-
# ar a rotina iterativa:

prod_step1 = mu_step1*s2_step1*nu_step1; prod_step1
prod_step2 = mu_step2*s2_step2*nu_step2; prod_step2
prod_step3 = mu_step3*s2_step3*nu_step3; prod_step3

# Função para calcular a constante de proporcionalidade:

cprop = function(l, X, mgr, s2gr, ngr, prst, m, V, a, d){
  c = 0
  ar.aux = array(dim = rep(l, 3))
  for (i in 1:l) {
    for (j in 1:l) {
      for (k in 1:l) {
        aux = logkpost(X, mgr[i], s2gr[j], ngr[k],
                       m, V, a, d)
        ar.aux[i,j,k] = aux
        c = c + exp(aux)*prst
      }
    }
  }
  hist(ar.aux, prob=T, breaks=1000)
  return(list(c, min(ar.aux), max(ar.aux),
              median(ar.aux), mean(ar.aux)))
}

c1 = cprop(l=l1, X=sam,
           mgr=mu_gr1, s2gr=s2_gr1, ngr=nu_gr1,
           prst=prod_step1, m=m, V=V, a=a, d=d)
c2 = cprop(l=l2, X=sam,
           mgr=mu_gr2, s2gr=s2_gr2, ngr=nu_gr2,
           prst=prod_step2, m=m, V=V, a=a, d=d)
c3 = cprop(l=l3, X=sam,
           mgr=mu_gr3, s2gr=s2_gr3, ngr=nu_gr3,
           prst=prod_step3, m=m, V=V, a=a, d=d)
c1; c2; c3

# Vamos precisar dos produtos dos grides 2 a 2:

pr1_dup12 = mu_step1*s2_step1
pr1_dup13 = mu_step1*nu_step1
pr1_dup23 = s2_step1*nu_step1
pr1_dup12; pr1_dup13; pr1_dup23

pr2_dup12 = mu_step2*s2_step2
pr2_dup13 = mu_step2*nu_step2
pr2_dup23 = s2_step2*nu_step2
pr2_dup12; pr2_dup13; pr2_dup23

pr3_dup12 = mu_step3*s2_step3
pr3_dup13 = mu_step3*nu_step3
pr3_dup23 = s2_step3*nu_step3
pr3_dup12; pr3_dup13; pr3_dup23

# Para calcular densidades marginais a posteriori, consi-
# dere as funções:

postmu_quarie = function(l, X, mgr, s2gr, ngr, prst,
                         m, V, a, d, c) {
    postmu = numeric(l)
    for (i in 1:l) {
      postconj = 0
      for (j in 1:l) {
        for(k in 1:l) {
          aux = logkpost(X, mgr[i], s2gr[j], ngr[k],
                         m, V, a, d)
          postconj = postconj + exp(aux)*prst
        }
      }
      postmu[i] = postconj/c
    }
    return(postmu)
}

pmq1 = postmu_quarie(
  l=l1, X=sam, mgr=mu_gr1, s2gr=s2_gr1, ngr=nu_gr1,
  prst=pr1_dup23, m=m, V=V, a=a, d=d, c=c1[[1]])
pmq2 = postmu_quarie(
  l=l2, X=sam, mgr=mu_gr2, s2gr=s2_gr2, ngr=nu_gr2,
  prst=pr2_dup23, m=m, V=V, a=a, d=d, c=c2[[1]])
pmq3 = postmu_quarie(
  l=l3, X=sam, mgr=mu_gr3, s2gr=s2_gr3, ngr=nu_gr3,
  prst=pr3_dup23, m=m, V=V, a=a, d=d, c=c3[[1]])

pmq1; pmq2; pmq3
hist(pmq1); hist(pmq2); hist(pmq3)

plot(mu_gr1, pmq1, type = "l", main = "",
     xlab = expression(paste(mu)), ylab = "")
plot(mu_gr2, pmq2, type = "l", main = "",
     xlab = expression(paste(mu)), ylab = "")
plot(mu_gr3, pmq3, type = "l", main = "",
     xlab = expression(paste(mu)), ylab = "")

posts2_quarie = function(l, X, mgr, s2gr, ngr, prst,
                         m, V, a, d, c) {
  posts2 = numeric(l)
  for (i in 1:l) {
    postconj = 0
    for (j in 1:l) {
      for(k in 1:l) {
        aux = logkpost(X, mgr[j], s2gr[i], ngr[k],
                       m, V, a, d)
        postconj = postconj + exp(aux)*prst
      }
    }
    posts2[i] = postconj/c
  }
  return(posts2)
}

psq1 = posts2_quarie(
  l=l1, X=sam, mgr=mu_gr1, s2gr=s2_gr1, ngr=nu_gr1,
  prst=pr1_dup13, m=m, V=V, a=a, d=d, c=c1[[1]])
psq2 = posts2_quarie(
  l=l2, X=sam, mgr=mu_gr2, s2gr=s2_gr2, ngr=nu_gr2,
  prst=pr2_dup13, m=m, V=V, a=a, d=d, c=c2[[1]])
psq3 = posts2_quarie(
  l=l3, X=sam, mgr=mu_gr3, s2gr=s2_gr3, ngr=nu_gr3,
  prst=pr3_dup13, m=m, V=V, a=a, d=d, c=c3[[1]])

psq1; psq2; psq3
hist(psq1); hist(psq2); hist(psq3)

plot(s2_gr1, psq1, type = "l", main = "",
     xlab = expression(paste(sigma^2)), ylab = "")
plot(s2_gr2, psq2, type = "l", main = "",
     xlab = expression(paste(sigma^2)), ylab = "")
plot(s2_gr3, psq3, type = "l", main = "",
     xlab = expression(paste(sigma^2)), ylab = "")

postnu_quarie = function(l, X, mgr, s2gr, ngr, prst,
                         m, V, a, d, c) {
  postnu = numeric(l)
  for (i in 1:l) {
    postconj = 0
    for (j in 1:l) {
      for(k in 1:l) {
        aux = logkpost(X, mgr[j], s2gr[k], ngr[i],
                       m, V, a, d)
        postconj = postconj + exp(aux)*prst
      }
    }
    postnu[i] = postconj/c
  }
  return(postnu)
}

pnq1 = postnu_quarie(
  l=l1, X=sam, mgr=mu_gr1, s2gr=s2_gr1, ngr=nu_gr1,
  prst=pr1_dup12, m=m, V=V, a=a, d=d, c=c1[[1]])
pnq2 = postnu_quarie(
  l=l2, X=sam, mgr=mu_gr2, s2gr=s2_gr2, ngr=nu_gr2,
  prst=pr2_dup12, m=m, V=V, a=a, d=d, c=c2[[1]])
pnq3 = postnu_quarie(
  l=l3, X=sam, mgr=mu_gr3, s2gr=s2_gr3, ngr=nu_gr3,
  prst=pr3_dup12, m=m, V=V, a=a, d=d, c=c3[[1]])

pnq1; pnq2; pnq3
hist(pnq1); hist(pnq2); hist(pnq3)

plot(nu_gr1, pnq1, type = "l", main = "",
     xlab = expression(paste(nu)), ylab = "")
plot(nu_gr2, pnq2, type = "l", main = "",
     xlab = expression(paste(nu)), ylab = "")
plot(nu_gr3, pnq3, type = "l", main = "",
     xlab = expression(paste(nu)), ylab = "")

# Para calcular média, variância, assimetria e curtose de
# cada parâmetro a posteriori, considere as funções a se-
# guir:

stat_post = function(gr, marg, prst) {
  media = 0
  var = 0
  assim = 0
  cur = 0
  l = length(gr)

  aux1 = sum(gr*marg*prst)     #Aproxima 1º momento.
  aux2 = sum((gr^2)*marg*prst) #Aproxima 2º momento.
  aux3 = sum((gr^3)*marg*prst) #Aproxima 3º momento.
  aux4 = sum((gr^4)*marg*prst) #Aproxima 4º momento.

  media = aux1
  var = aux2 - (media)^2
  assim = (aux3 - 3*media*var - media^3)/(var^(3/2))
  cur = (aux4 - 4*media*aux3 +
         6*media^2*aux2 - 3*(media^4))/(var^2)
  return(list(media, var, assim, cur))
}

stat_post(mu_gr1, pmq1, mu_step1)
stat_post(mu_gr2, pmq2, mu_step2)
stat_post(mu_gr3, pmq3, mu_step3)

stat_post(s2_gr1, psq1, s2_step1)
stat_post(s2_gr2, psq2, s2_step2)
stat_post(s2_gr3, psq3, s2_step3)

stat_post(nu_gr1, pnq1, nu_step1)
stat_post(nu_gr2, pnq2, nu_step2)
stat_post(nu_gr3, pnq3, nu_step3)
