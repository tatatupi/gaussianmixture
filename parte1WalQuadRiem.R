rm(list=ls(all=TRUE))

# Tamanho amostral:
n = 500

# Parâmetros do modelo para a distribuição amostral:
mu = 11
sigma2_1 = 4
sigma2_2 = 0.16
sigma2_3 = 1
nu = 0.2

# Hiperparâmetros das distribuições a priori:
m = 11
V = 1
a_1 = 7; d_1 = 26
a_2 = 7; d_2 = 1
a_3 = 7; d_3 = 6.5

# Para ter uma ideia da dispersão de mu a priori:
set.seed(122019)
hist(rnorm(n=10000, mean=m, sd=sqrt(V*sigma2_1)),
     prob=T, breaks=100)
hist(rnorm(n=10000, mean=m, sd=sqrt(V*sigma2_2)),
     prob=T, breaks=100)
hist(rnorm(n=10000, mean=m, sd=sqrt(V*sigma2_3)),
     prob=T, breaks=100)
# Massa probabilística concentrada quase que totalmente
# entre 5 e 17 (caso 1); 10 e 12 (caso 2) e 8 e 14 (ca-
# so 3) e todas aproximadamente centradas em 11.

# Para ter uma ideia da dispersão de sigma2 a priori:
set.seed(122019)
s1 = 1/rgamma(n=10000, shape=a_1, rate=d_1)
hist(s1, prob=T, breaks=100)
median(s1); mean(s1); min(s1); max(s1)
s2 = 1/rgamma(n=10000, shape=a_2, rate=d_2)
hist(s2, prob=T, breaks=100)
median(s2); mean(s2); min(s2); max(s2)
s3 = 1/rgamma(n=10000, shape=a_3, rate=d_3)
hist(s3, prob=T, breaks=100)
median(s3); mean(s3); min(s3); max(s3)
# Massa probabilística concentrada quase que totalmente
# entre 1.5 e 10.5 e aproximadamente centrada em 4 (ca-
# so 1), concent. quase que totalmente entre .05 e .0.4
# e aproximadamente centrada em 0.16 (caso 2) e concen-
# trada quase que totalmente entre 0.4 e 3.4, aproxima-
# damente centrada em 1 (caso 3).

# A priori, o parâmetro nu é uniforme em (0,1).

# Como gerar indicadoras da mistura dado nu:
#Inu_pdf = function(n, nu, a = 0, b = 1) {
#  unif = runif(n)
#  #(unif < nu)*a + (unif >= nu)*b # Taiguara (complem.).
#  (unif >= nu)*a + (unif < nu)*b  # Walmir (val. orig.).
#}

#mean(Inu_pdf(100, nu))
#mean(Inu_pdf(1000, nu))
#mean(Inu_pdf(10000, nu))

# Como gerar da variável auxiliar U (Dani):

U_pdf = function(n, nu) {
  u = sample(c(100, 1), prob=c(nu, 1-nu),
             size=n, replace=T)
}

#table(U_pdf(100, nu))
#table(U_pdf(1000, nu))
#table(U_pdf(10000, nu))

# Gerando uma amostra de tamanho $n$:

# 1ª forma: gerando da própria densidade da mistura:

#nu_s = nu_pdf(n, nu)
#sample = nu_s*(mu + sigma2*100*rnorm(n)) +
#         (1 - nu_s)*(mu + sigma2*rnorm(n))

# Note que devemos usar todos os parâmetros do modelo,
# como se os conhecêssemos:

#set.seed(122019)
#Inu = Inu_pdf(n, nu)
#set.seed(122019)
#sample1 = Inu*(mu + sqrt(sigma2/100)*rnorm(n)) +
#          (1 - Inu)*(mu + sqrt(sigma2)*rnorm(n))
#hist(sample, breaks=50, prob=T)
#sample1_1 = Inu*(rnorm(n, mu, sqrt(sigma2_1/100))) +
#            (1 - Inu)*(rnorm(n, mu, sqrt(sigma2_1)))
#hist(sample1_1, breaks=50, prob=T)
#sample1_2 = Inu*(rnorm(n, mu, sqrt(sigma2_2/100))) +
#  (1 - Inu)*(rnorm(n, mu, sqrt(sigma2_2)))
#hist(sample1_2, breaks=50, prob=T)
#sample1_3 = Inu*(rnorm(n, mu, sqrt(sigma2_3/100))) +
#  (1 - Inu)*(rnorm(n, mu, sqrt(sigma2_3)))
#hist(sample1_3, breaks=50, prob=T)

# 2ª forma: gerando de forma hierárquica (Dani):

set.seed(122019)
u = U_pdf(n, nu)
sample2_1 = rnorm(n, mu, sqrt(sigma2_1/u))
hist(sample2_1, breaks=100, prob=T)
sample2_2 = rnorm(n, mu, sqrt(sigma2_2/u))
hist(sample2_2, breaks=100, prob=T)
sample2_3 = rnorm(n, mu, sqrt(sigma2_3/u))
hist(sample2_3, breaks=100, prob=T)

# Dúvida: devemos mesmo dividir a variância da normal por
# 100 como diz o enunciado na geração hierárquica? O cer-
# to não seria multiplicar? Pelo o que os artigos indicam
# devemos sim dividir.

# Definindo a distribuição *a posteriori*:

A = function(X, mu, sigma2, nu) {
  n = length(X)
  k1 = (nu/10)*exp(-(X-mu)^2/(2*100*sigma2))
  k2 = (1-nu)*exp(-(X-mu)^2/(2*sigma2))
  k = k1 + k2
  return(k)
}

#A(sample2_1, mu, sigma2_1, nu)
#A(sample2_2, mu, sigma2_2, nu)
#A(sample2_3, mu, sigma2_3, nu)
#hist(A(sample2_1, mu, sigma2_1, nu), prob=T, breaks=50)
#hist(A(sample2_2, mu, sigma2_2, nu), prob=T, breaks=50)
#hist(A(sample2_3, mu, sigma2_3, nu), prob=T, breaks=50)
#prod(A(sample2_1, mu, sigma2_1, nu))
#prod(A(sample2_2, mu, sigma2_2, nu))
#prod(A(sample2_3, mu, sigma2_3, nu))

# Claramente, devemos trabalhar com o logaritmo de A, uma
# vez que o produtório é muito baixo como esperado, mesmo
# quando sigma2 é pequeno.

logA = function(X, mu, sigma2, nu) {
  n = length(X)
  k1 = (nu/10)*exp(-(X-mu)^2/(2*100*sigma2))
  k2 = (1-nu)*exp(-(X-mu)^2/(2*sigma2))
  k = log(k1 + k2)
  return(k)
  #return(sum(k))
}

logA(sample2_1, mu, sigma2_1, nu)
logA(sample2_2, mu, sigma2_2, nu)
logA(sample2_3, mu, sigma2_3, nu)
hist(logA(sample2_1, mu, sigma2_1, nu),
     prob=T, breaks=50)
hist(logA(sample2_2, mu, sigma2_2, nu),
     prob=T, breaks=50)
hist(logA(sample2_3, mu, sigma2_3, nu),
     prob=T, breaks=50)
sum(logA(sample2_1, mu, sigma2_1, nu))
sum(logA(sample2_2, mu, sigma2_2, nu))
sum(logA(sample2_3, mu, sigma2_3, nu))

# logA é negativa e tratável. A diferença é pequena mesmo
# se o parâmetro sigma2 variar bastante.

h = function(n, mu, sigma2, nu, m, V, a, d) {
  k1 = (1/sigma2)^((n + 1)/2 + a + 1)
  k2 = exp(-((mu - m)^2/(2*V) + d)/sigma2)
  k = k1*k2
  return(k)
}

#h(n, mu, sigma2_1, nu, m, V, a_1, d_1)
#h(n, mu, sigma2_2, nu, m, V, a_2, d_2)
#h(n, mu, sigma2_3, nu, m, V, a_3, d_3)
# Muito pequeno se sigma2 é grande;
# Muito grande se sigma2 é pequeno.

logh = function(n, mu, sigma2, nu, m, V, a, d) {
  k1 = -((n + 1)/2 + a + 1)*log(sigma2)
  k2 = -((mu - m)^2/(2*V) + d)/sigma2
  k = k1 + k2
  return(k)
  #return(sum(k))
}

logh(n, mu, sigma2_1, nu, m, V, a_1, d_1)
logh(n, mu, sigma2_2, nu, m, V, a_2, d_2)
logh(n, mu, sigma2_3, nu, m, V, a_3, d_3)
# Negativo se sigma2 é grande, mas tratável;
# Positivo se sigma2 é pequeno, idem.

# Para encontrar o máximo da distribuição a posteriori,
# será melhor trabalhar com o logaritmo do seu núcleo:

logkerpost_1 = logA(sample2_1, mu, sigma2_1, nu) +
               logh(sample2_1, mu, sigma2_1, nu,
                    m, V, a_1, d_1)
logkerpost_2 = logA(sample2_2, mu, sigma2_2, nu) +
               logh(sample2_2, mu, sigma2_2, nu,
                    m, V, a_2, d_2)
logkerpost_3 = logA(sample2_3, mu, sigma2_3, nu) +
               logh(sample2_3, mu, sigma2_3, nu,
                    m, V, a_3, d_3)
hist(logkerpost_1, prob=T, breaks = 50)
hist(logkerpost_2, prob=T, breaks = 50)
hist(logkerpost_3, prob=T, breaks = 50)

kerpost_1 = exp(logkerpost_1)
kerpost_2 = exp(logkerpost_2)
kerpost_3 = exp(logkerpost_3)
hist(kerpost_1, prob=T, breaks = 50)
hist(kerpost_2, prob=T, breaks = 50)
hist(kerpost_3, prob=T, breaks = 50)

# Método da Quadratura de Riemann:

L = 50 # Número de intervalos.

# Note que, a posteriori, não temos ideia de como cada um
# dos parâmetros estão distribuídos dado o respectivo es-
# paço paramétrico. Assim, precisaríamos amostrar valores
# de cada um deles, de suas densidades marginais a poste-
# riori correspondentes, para ter uma ideia do respectivo
# tamanho do intervalo de integração.

# Como a amostra tem tamanho 500, pelo Teorema Central do
# Limite é possível aproximar a variância a posteriori de
# cada parâmetro por 1/500 da variância a priori.

# Para mu, aproximamos sua distribuição a posteriori por:

mu = 11
sigma2_1 = 4
sigma2_2 = 0.16
sigma2_3 = 1
nu = 0.2

m = 11
V = 1
# Variância teórica de mu a priori: 4; 0.16; 1
# Variância aproximada de mu a posteriori: 4/n; 0.16/n e
# 1/n
# Para se ter uma ideia da dispersão de mu a posteriori:
set.seed(122019)
hist(rnorm(n=10000, mean=m, sd=sqrt(V*sigma2_1/n)),
     prob=T, breaks=100)
hist(rnorm(n=10000, mean=m, sd=sqrt(V*sigma2_2/n)),
     prob=T, breaks=100)
hist(rnorm(n=10000, mean=m, sd=sqrt(V*sigma2_3/n)),
     prob=T, breaks=100)
# Massa probabilística concentrada quase que totalmente
# entre 10.7 e 11.3 (caso 1), 10.95 e 11.05 (caso 2) ou
# 10.9 e 11.1 (caso 3), todas aproximadamente centradas
# em 11.

# Temos os intervalos:

mu_step_1 = (11.3-10.7)/L
mu_gr_1 = seq(10.7, 11.3, mu_step_1)
mu_gr_1

mu_step_2 = (11.05-10.95)/L
mu_gr_2 = seq(10.95, 11.05, mu_step_2)
mu_gr_2

mu_step_3 = (11.1-10.9)/L
mu_gr_3 = seq(10.9, 11.1, mu_step_3)
mu_gr_3

# Para sigma2, aproximamos sua distrib. a posteriori por:

a_1 = 7; d_1 = 26
a_2 = 7; d_2 = 1
a_3 = 7; d_3 = 6.5
# Variância teórica de sigma2 a priori:
26^2/((7-1)^2*(7-2))  #3.7556
1^2/((7-1)^2*(7-2))   #0.0056
6.5^2/((7-1)^2*(7-2)) #0.2347
# Variância aproximada de mu a posteriori:
3.7556/n; 0.0056/n; 0.2347/n
# Para ter uma ideia da dispersão de sigma2 a posteriori:
set.seed(122019)
hist(rnorm(n=10000, mean=sigma2_1, sd=sqrt(3.7556/n)),
     prob=T, breaks=100)
hist(rnorm(n=10000, mean=sigma2_2, sd=sqrt(0.0056/n)),
     prob=T, breaks=100)
hist(rnorm(n=10000, mean=sigma2_3, sd=sqrt(0.2347/n)),
     prob=T, breaks=100)
# Massa probabilística concentrada quase que totalmente
# entre 3.8 e 4.2 e aproximadamente centrada em 4 (caso
# 1), concentrada quase que totalmente entre .15 e .17,
# e aproximadamente centrada em 0.16 (caso 2) e concen-
# trada quase que totalmente entre 0.94 e 1.06, aproxi-
# madamente centrada em 1 (caso 3).

# Temos os intervalos:

sigma2_step_1 = (4.2-3.8)/L
sigma2_gr_1 = seq(3.8, 4.2, sigma2_step_1)
sigma2_gr_1

sigma2_step_2 = (0.17-0.15)/L
sigma2_gr_2 = seq(0.15, 0.17, sigma2_step_2)
sigma2_gr_2

sigma2_step_3 = (1.06-0.94)/L
sigma2_gr_3 = seq(0.94, 1.06, sigma2_step_3)
sigma2_gr_3

# Para nu, sua variância a priori é igual a 1/12. Se con-
# siderarmos que mu=0.2 no modelo, a distribuição de nu a
# posteriori pode ser aproximada pelo TCL por:

set.seed(122019)
hist(rnorm(n=10000, mean=nu, sd=sqrt(1/(12*500))),
     prob=T, breaks=100)

# Temos o intervalo:

nu_step = (0.24-0.16)/L
nu_gr = seq(0.16, 0.24, nu_step)
nu_gr

# E por fim, as grades:

grid_tri_1 = cbind(mu_gr_1, sigma2_gr_1, nu_gr)
grid_tri_2 = cbind(mu_gr_2, sigma2_gr_2, nu_gr)
grid_tri_3 = cbind(mu_gr_3, sigma2_gr_3, nu_gr)
grid_tri_1
grid_tri_2
grid_tri_3

# Observe que no segundo grid temos variações bem menores
# para mu e sigma2, dado que sigma2 é pequeno.

l = nrow(grid_tri_1); l

# Como os tamanhos dos passos são iguais em todas as três
# dimensões, podemos calcular os produtos antes de inici-
# ar a rotina iterativa:

prod_step_1 = mu_step_1*sigma2_step_1*nu_step
prod_step_2 = mu_step_2*sigma2_step_2*nu_step
prod_step_3 = mu_step_3*sigma2_step_3*nu_step
prod_step_1
prod_step_2
prod_step_3

# Função para calcular a constante de proporcionalidade:

cprop = function(l, X, mgr, s2gr, ngr, prst,
                  m, V, a, d, const=0) {
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
        c = c + exp(aux + const)*prst
      }
    }
  }
  c_ver = c/exp(const)
  hist(ar.aux, prob=T, breaks=1000)
  return(list(c, c_ver, min(ar.aux), max(ar.aux),
              median(ar.aux), mean(ar.aux)))
}

# No primeiro caso (sigma2 = 4), temos que:

c1 = cprop(l=l, X=sample2_1, mgr=mu_gr_1,
           s2gr=sigma2_gr_1, ngr=nu_gr,
           prst=prod_step_1, m=m, V=V, a=a_1, d=d_1)
c1

# A soma de c será na prática nula se a variância do mo-
# delo (sigma2) for alta (a precisão do R não basta para
# detectar que cada parcela é diferente de zero). Assim,
# não podemos estimar c corretamente sem usar algum ter-
# mo multiplicativo. Refazendo com este termo (exp(640))
# temos:

cprop(l=l, X=sample2_1, mgr=mu_gr_1, s2gr=sigma2_gr_1,
      ngr=nu_gr, prst=prod_step_1, m=m, V=V,
      a=a_1, d=d_1, const=640)

# Mesmo multiplicando e dividindo, temos o mesmo resulta-
# do de antes.

# No segundo caso (sigma2 = 0.16), temos que:

c2 = cprop(l=l, X=sample2_2, mgr=mu_gr_2,
           s2gr=sigma2_gr_2, ngr=nu_gr,
           prst=prod_step_2, m=m, V=V, a=a_2, d=d_2)
c2

# Em ambos os casos (variância grande ou pequena) a qua-
# dratura é demorada com muitas iterações, e mesmo assim
# não se obtém um número facilmente tangível. E se a va-
# riância for unitária:

c3 = cprop(l=l, X=sample2_3, mgr=mu_gr_3,
           s2gr=sigma2_gr_3, ngr=nu_gr,
           prst=prod_step_3, m=m, V=V, a=a_3, d=d_3)
c3

# Mesmo com variância unitária, c tem um valor muito pe-
# queno.

# Pacote para cálculo de assimetria e curtose amostrais:

#install.packages("moments")
library(moments)

# Vamos precisar dos produtos dos grides 2 a 2:

pr1_dup12 = mu_step_1*sigma2_step_1
pr1_dup13 = mu_step_1*nu_step
pr1_dup23 = sigma2_step_1*nu_step
pr1_dup12; pr1_dup13; pr1_dup23

pr2_dup12 = mu_step_2*sigma2_step_2
pr2_dup13 = mu_step_2*nu_step
pr2_dup23 = sigma2_step_2*nu_step
pr2_dup12; pr2_dup13; pr2_dup23

pr3_dup12 = mu_step_3*sigma2_step_3
pr3_dup13 = mu_step_3*nu_step
pr3_dup23 = sigma2_step_3*nu_step
pr3_dup12; pr3_dup13; pr3_dup23

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

pmq1 = postmu_quarie(l=l, X=sample2_1,
  mgr=mu_gr_1, s2gr=sigma2_gr_1, ngr=nu_gr,
  prst=pr1_dup23, m=m, V=V, a=a_1, d=d_1, c=c1[[1]])
pmq1
hist(pmq1[[1]])

pmq2 = postmu_quarie(l=l, X=sample2_2,
  mgr=mu_gr_2, s2gr=sigma2_gr_2, ngr=nu_gr,
  prst=pr2_dup23, m=m, V=V, a=a_2, d=d_2, c=c2[[1]])
pmq2
hist(pmq2[[1]])

pmq3 = postmu_quarie(l=l, X=sample2_3,
  mgr=mu_gr_3, s2gr=sigma2_gr_3, ngr=nu_gr,
  prst=pr3_dup23, m=m, V=V, a=a_3, d=d_3, c=c3[[1]])
pmq3
hist(pmq3[[1]])

# Todas as estimativas foram ou muito pequenas se a vari-
# ância era maior do que 1, ou muito grandes se a variân-
# cia era menor do que ou igual a 1 (mas um menos discre-
# pantes neste último caso).

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

psq1 = posts2_quarie(l=l, X=sample2_1,
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
