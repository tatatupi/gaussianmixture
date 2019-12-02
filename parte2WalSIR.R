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

# Escolhemos como densidade conhecida a distribuição nor-
# mal multivariada com dimensão 3.

# Como nem todo argumento da posteriori varia na reta re-
# al, vamos reparametrizá-los de modo que todos os supor-
# tes dos novos parâmetros coincidam com os das médias na
# normal multivariada.

# Logo, a reparametrização será feita tomando funções que
# levem sigma e nu dos seus domínios para imagens na reta
# real (funções logarítmica e logística).

logkpost_re(sam, mu, s2, nu, m, V, a, d)

# Tamanhos para o total de pesos utilizados e de amostras
# da posteriori via SIR:

aw1 = 5000; aw2 = 50000; aw3 = 500000
am1 = aw1/10; am2 = aw2/10; am3 = aw3/10

# Na escolha dos parâmetros da distribuição normal multi-
# variada, o vetor de médias será dado pelos valores ver-
# dadeiros de cada parâmetro do modelo gerador dos dados.

# Já na matriz de covariância, sua diagonal principal se-
# rá dada pelo quadrado de um sexto do valor da distância
# entre os pontos de concentração da distribuição de cada
# parâmetro, dados os demais. Vimos na quadratura de Rie-
# mann que:

# mu está concentrado entre 10.85 e 11.13;
# sigma2 está concentrado entre 0.48 e 0.78;
# nu está concentrado entre 0.13 e 0.26.

# Reparametrizando, temos que:

# par1 está concentrado entre 10.85 e 11.13;
# par2 está concentrado entre log(0.48) e log(0.78);
# par3 está concentrado entre log(0.13/0.87) e
# log(0.26/0.74);

((11.13-10.85)/6)^2
((log(0.78)-log(0.48))/6)^2
((log(0.26/0.74)-log(0.13/0.87))/6)^2

# Logo, diag(Matriz) = c(0.0022, 0.0065, 0.0203).

# Embora cada uma das distribuições a posteriori não seja
# exatamente dada pela normal univariada, a ideia é obter
# uma variância tal que 99.7% da massa probabilística es-
# tivesse contida (3 desvios padrões), na prática 100%, e
# ao mesmo tempo cobrir o mínimo possível de regiões cuja
# massa probabilística seja praticamente nula.

# Para os demais elementos da matriz de covariância, ire-
# mos supor que a correlação entre os parâmetros é nula.

# Geremos agora da distribuição normal multivariada:

library(mvtnorm) # Pacote com funções densidade e gerado-
                 # ra da distrib. normal multivariada.

set.seed(122019)
r1_q = rmvnorm(aw1, mean=c(11, log(0.64), log(0.2/0.8)),
               sigma=diag(c(0.0022, 0.0065, 0.0203)))
set.seed(122019)
r2_q = rmvnorm(aw2, mean=c(11, log(0.64), log(0.2/0.8)),
               sigma=diag(c(0.0022, 0.0065, 0.0203)))
set.seed(122019)
r3_q = rmvnorm(aw3, mean=c(11, log(0.64), log(0.2/0.8)),
               sigma=diag(c(0.0022, 0.0065, 0.0203)))

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

weiSIR = function(X, r_q, m, V, a, d, d_q) {
  k = nrow(r_q)
  aux1 = aux2 = w = numeric(k)
  for(i in 1:k) {
    aux1[i] = logkpost_re(X=X, mu=r_q[i,1],
      sigma2=exp(r_q[i,2]), nu=1/(1+exp(-r_q[i,3])),
      m=m, V=V, a=a, d=d)
  }
  aux2 = aux1/d_q
  w = aux2/sum(aux2)
  return(w)
}

w1 = weiSIR(sam, r1_q, m, V, a, d, d1_q)
w2 = weiSIR(sam, r2_q, m, V, a, d, d2_q)
w3 = weiSIR(sam, r3_q, m, V, a, d, d3_q)

set.seed(122019)
am1_mu_p = sample(x=r1_q[,1], size=am1,
                  replace=T, prob=w1)
set.seed(122019)
am1_s2_p = sample(x=exp(r1_q[,2]), size=am1,
                  replace=T, prob=w1)
set.seed(122019)
am1_nu_p = sample(x=1/(1+exp(-r1_q[,3])), size=am1,
                  replace=T, prob=w1)

set.seed(122019)
am2_mu_p = sample(x=r2_q[,1], size=am2,
                  replace=T, prob=w2)
set.seed(122019)
am2_s2_p = sample(x=exp(r2_q[,2]), size=am2,
                  replace=T, prob=w2)
set.seed(122019)
am2_nu_p = sample(x=1/(1+exp(-r2_q[,3])), size=am2,
                  replace=T, prob=w2)

set.seed(122019)
am3_mu_p = sample(x=r3_q[,1], size=am3,
                  replace=T, prob=w3)
set.seed(122019)
am3_s2_p = sample(x=exp(r3_q[,2]), size=am3,
                  replace=T, prob=w3)
set.seed(122019)
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

st_sam_post(am1_mu_p)
hist(am1_mu_p, breaks=50, prob=T, main="",
     xlab=expression(paste(mu)), ylab="")
st_sam_post(am2_mu_p)
hist(am2_mu_p, breaks=50, prob=T, main="",
     xlab=expression(paste(mu)), ylab="")
st_sam_post(am3_mu_p)
hist(am3_mu_p, breaks=50, prob=T, main="",
     xlab=expression(paste(mu)), ylab="")

st_sam_post(am1_s2_p);
hist(am1_s2_p, breaks=50, prob=T, main="",
     xlab=expression(paste(sigma^2)), ylab="")
st_sam_post(am2_s2_p)
hist(am2_s2_p, breaks=50, prob=T, main="",
     xlab=expression(paste(sigma^2)), ylab="")
st_sam_post(am3_s2_p)
hist(am3_s2_p, breaks=50, prob=T, main="",
     xlab=expression(paste(sigma^2)), ylab="")

st_sam_post(am1_nu_p)
hist(am1_nu_p, breaks=50, prob=T, main="",
     xlab=expression(paste(nu)), ylab="")
st_sam_post(am2_nu_p)
hist(am2_nu_p, breaks=50, prob=T, main="",
     xlab=expression(paste(nu)), ylab="")
st_sam_post(am3_nu_p)
hist(am3_nu_p, breaks=50, prob=T, main="",
     xlab=expression(paste(nu)), ylab="")
