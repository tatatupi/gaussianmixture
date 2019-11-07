rm(list=ls(all=TRUE))

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

# Gerando da variável auxiliar U:

U_pdf = function(n, nu) {
  u = sample(c(100, 1), prob=c(nu, 1-nu),
             size=n, replace=T)
}

#table(U_pdf(100, nu))
#table(U_pdf(1000, nu))
#table(U_pdf(10000, nu))

# Gerando uma amostra de tamanho n hierarquicamente:

set.seed(122019)
u = U_pdf(n, nu)
sam = rnorm(n, mu, sqrt(s2*u))
hist(sam, breaks=100, prob=T)

# Definindo a distribuição *a posteriori*:

# Em geral, se trabalha com o logaritmo da verossimilhan-
# ça, uma vez que o produtório é muito baixo mesmo quando
# sigma2 é pequeno. Na função a seguir calculamos o loga-
# ritmo do produtório apenas da soma envolvendo os 2 com-
# ponentes da mistura:

logA = function(X, mu, sigma2, nu) {
  n = length(X)
  k1 = (nu/10)*exp(-(X-mu)^2/(2*100*sigma2))
  k2 = (1-nu)*exp(-(X-mu)^2/(2*sigma2))
  k = log(k1 + k2)
  return(k)
}

logA(sam, mu, s2, nu)
hist(logA(sam, mu, s2, nu), prob=T, breaks=50)
sum(logA(sam, mu, s2, nu))

# Função para calcular logaritmo do produto de termos en-
# volvendo apenas parâmetros e hiperparâmetros livres dos
# valores na amostra (presentes nos componentes da mistu-
# ra):

logh = function(n, mu, sigma2, nu, m, V, a, d) {
  k1 = -((n + 1)/2 + a + 1)*log(sigma2)
  k2 = -((mu - m)^2/(2*V) + d)/sigma2
  k = k1 + k2
  return(k)
}

logh(n, mu, s2, nu, m, V, a, d)

# Por fim, para o logaritmo (do núcleo) da distribuição a
# posteriori, temos que:

logkpost = function(X, mu, sigma2, nu, m, V, a, d) {
  n = length(X)
  lA = logA(X, mu, sigma2, nu)
  lh = logh(n, mu, sigma2, nu, m, V, a, d)
  lkp = sum(lA) + lh
  return(lkp)
}

# Método da Quadratura de Riemann:

L = 50 # Número de intervalos.

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
plot(mu_sup, kp_mu, type="l")

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
plot(s2_sup, kp_s2, type="l")

# Para sigma2 o intervalo de integ. será de 0.48 a 0.76.

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
plot(nu_sup, kp_nu, type="l")

# Para nu o intervalo de integ. será de 0.13 a 0.26.

# Observe que nenhum dos três intervalos escolhidos ficou
# centrado no mesmo valor assumido para o modelo. Mais um
# sinal da importância de considerar tanto a informação a
# priori quanto a dada pela amostra.

# Temos então os intervalos:

mu_step = (11.13 - 10.85)/L
mu_gr = seq(10.85, 11.13, mu_step)
mu_gr

s2_step = (0.76 - 0.48)/L
s2_gr = seq(0.48, 0.76, s2_step)
s2_gr

nu_step = (0.26 - 0.13)/L
nu_gr = seq(0.13, 0.26, nu_step)
nu_gr

# E a grade:

grid_tri = cbind(mu_gr, s2_gr, nu_gr); grid_tri
l = nrow(grid_tri); l

# Como os tamanhos dos passos são iguais em todas as três
# dimensões, podemos calcular os produtos antes de inici-
# ar a rotina iterativa:

prod_step = mu_step*s2_step*nu_step; prod_step

# Função para calcular a constante de proporcionalidade:

cprop = function(l, X, mgr, s2gr, ngr, prst, m, V, a, d){
  c = 0
  n = length(X)
  ar.aux = array(dim = rep(l, 3))
  for (i in 1:l) {
    for (j in 1:l) {
      for (k in 1:l) {
        aux1 = sum(logA(X, mgr[i], s2gr[j], ngr[k]))
        aux2 = logh(n, mgr[i], s2gr[j], ngr[k],
                    m, V, a, d)
        aux = aux1 + aux2
        ar.aux[i,j,k] = aux
        c = c + exp(aux)*prst
      }
    }
  }
  hist(ar.aux, prob=T, breaks=1000)
  return(list(c, min(ar.aux), max(ar.aux),
              median(ar.aux), mean(ar.aux)))
}

c = cprop(l=l, X=sam, mgr=mu_gr, s2gr=s2_gr, ngr=nu_gr,
           prst=prod_step, m=m, V=V, a=a, d=d)
c

# Pacote para cálculo de assimetria e curtose amostrais:

#install.packages("moments")
library(moments)

# Vamos precisar dos produtos dos grides 2 a 2:

pr_dup12 = mu_step*s2_step
pr_dup13 = mu_step*nu_step
pr_dup23 = s2_step*nu_step
pr_dup12; pr_dup13; pr_dup23

# Para calcular densidades marginais a posteriori, consi-
# dere as funções:

postmu_quarie = function(l, X, mgr, s2gr, ngr, prst,
                         m, V, a, d, c) {
    postmu = numeric(l)
    n = length(X)
    for (i in 1:l) {
      postconj = 0
      for (j in 1:l) {
        for(k in 1:l) {
          aux1 = sum(logA(X, mgr[i], s2gr[j], ngr[k]))
          aux2 = logh(n, mgr[i], s2gr[j], ngr[k],
                      m, V, a, d)
          aux = aux1 + aux2
          postconj = postconj + exp(aux)*prst
        }
      }
      postmu[i] = postconj/c
    }

    return(list(postmu, mean(postmu), var(postmu),
                median(postmu),
                min(postmu), max(postmu),
                skewness(postmu), kurtosis(postmu)))
}

pmq = postmu_quarie(l=l, X=sam, mgr=mu_gr, s2gr=s2_gr,
                    ngr=nu_gr, prst=pr_dup23, m=m, V=V,
                    a=a_1, d=d_1, c=c[[1]])
pmq
hist(pmq[[1]])

posts2_quarie = function(l, X, mgr, s2gr, ngr, prst,
                         m, V, a, d, c) {
  posts2 = numeric(l)
  n = length(X)
  for (i in 1:l) {
    postconj = 0
    for (j in 1:l) {
      for(k in 1:l) {
        aux1 = sum(logA(X, mgr[j], s2gr[i], ngr[k]))
        aux2 = logh(n, mgr[j], s2gr[i], ngr[k],
                    m, V, a, d)
        aux = aux1 + aux2
        postconj = postconj + exp(aux)*prst
      }
    }
    posts2[i] = postconj/c
  }

  return(list(posts2, mean(posts2), var(posts2),
              median(posts2),
              min(posts2), max(posts2),
              skewness(posts2), kurtosis(posts2)))
}

psq = posts2_quarie(l=l, X=sample2_1,
  mgr=mu_gr_1, s2gr=sigma2_gr_1, ngr=nu_gr,
  prst=pr1_dup13, m=m, V=V, a=a_1, d=d_1, c=c1[[1]])
psq1
hist(psq1[[1]])

psq2 = posts2_quarie(l=l, X=sample2_2,
  mgr=mu_gr_2, s2gr=sigma2_gr_2, ngr=nu_gr,
  prst=pr2_dup13, m=m, V=V, a=a_2, d=d_2, c=c2[[1]])
psq2
hist(psq2[[1]])

psq3 = posts2_quarie(l=l, X=sample2_3,
  mgr=mu_gr_3, s2gr=sigma2_gr_3, ngr=nu_gr,
  prst=pr3_dup13, m=m, V=V, a=a_3, d=d_3, c=c3[[1]])
psq3
hist(psq3[[1]])

# As estimativas para a média de sigma2, no primeiro e no
# terceiro casos, foram bem menores ao comparar com o va-
# lor de sigma2 no modelo.

postnu_quarie = function(l, X, mgr, s2gr, ngr, prst,
                         m, V, a, d, c) {
  postnu = numeric(l)
  n = length(X)
  for (i in 1:l) {
    postconj = 0
    for (j in 1:l) {
      for(k in 1:l) {
        aux1 = sum(logA(X, mgr[j], s2gr[k], ngr[i]))
        aux2 = logh(n, mgr[j], s2gr[k], ngr[i],
                    m, V, a, d)
        aux = aux1 + aux2
        postconj = postconj + exp(aux)*prst
      }
    }
    postnu[i] = postconj/c
  }

  return(list(postnu, mean(postnu), var(postnu),
              median(postnu),
              min(postnu), max(postnu),
              skewness(postnu), kurtosis(postnu)))
}

pnq1 = postnu_quarie(l=l, X=sample2_1,
  mgr=mu_gr_1, s2gr=sigma2_gr_1, ngr=nu_gr,
  prst=pr1_dup12, m=m, V=V, a=a_1, d=d_1, c=c1[[1]])
pnq1
hist(pnq1[[1]])

pnq2 = postnu_quarie(l=l, X=sample2_2,
  mgr=mu_gr_2, s2gr=sigma2_gr_2, ngr=nu_gr,
  prst=pr2_dup12, m=m, V=V, a=a_2, d=d_2, c=c2[[1]])
pnq2
hist(pnq2[[1]])

pnq3 = postnu_quarie(l=l, X=sample2_3,
  mgr=mu_gr_3, s2gr=sigma2_gr_3, ngr=nu_gr,
  prst=pr3_dup12, m=m, V=V, a=a_3, d=d_3, c=c3[[1]])
pnq3
hist(pnq3[[1]])

# Distribuição totalmente assimétrica, com valor estimado
# para a média muito longe do especificado no modelo. Pe-
# lo menos pode-se concluir que a integração numérica pa-
# ra nu não é alterada por variações em sigma_2.
