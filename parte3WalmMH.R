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

# Método MCMC (Metropolis-Hastings ou MH):

# Primeiramente, como não sabemos gerar da posteriori, i-
# remos considerar uma densidade conhecida, da qual sabe-
# mos gerar valores aleatórios (amostras), assim como foi
# feito para o método SIR.

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

# Tamanhos para a amostragem da posteriori via MH:

am1 = 500; am2 = 5000; am3 = 50000

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

# Embora não seja obrigatório fazer isso para o algoritmo
# Metropolis-Hastings, a eficiência do mesmo aumenta con-
# sideravelmente quando a distribuição desejada é repara-
# metrizada de modo que os suportes das quantidades alea-
# tórias a serem amostradas (parâmetros na distribuição a
# posteriori) sejam similares aos das quantidades aleató-
# rias da distribuição proposta (as três componentes para
# a distribuição normal multivariada).

# Devido à simetria da distribuição proposta (normal mul-
# tivariada), a razão de teste do algoritmo MH é reduzida
# às densidades da distribuição desejada nos pontos gera-
# dos consecutivamente da distribuição proposta.

# Para o burn-in na geração da cadeia em cada cenário, i-
# remos desprezar os 100 valores iniciais no primeiro ce-
# nário, os 1000 no segundo e os 10000 no terceiro.

library(mvtnorm) # Pacote com funções densidade e gerado-
                 # ra da distrib. normal multivariada.

# Geremos de uma uniforme padrão para controlar a aceita-
# ção de cada amostra gerada da distribuição proposta:

set.seed(122019); up1 = runif(am1)
set.seed(122019); up2 = runif(am2)
set.seed(122019); up3 = runif(am3)

# Temos então o seguinte vetor de valores iniciais para o
# algoritmo MH, os quais estão distantes dos valores ver-
# dadeiros no modelo da distribuição desejada, embora com
# probabilidade positiva:

x = c(10.86, log(0.50), log(0.14/0.86))
#x = c(10.9, log(0.55), log(0.17/0.83))
#x = c(10.95, log(0.60), log(0.19/0.81))
# Não importa a distância do ponto de maior massa da dis-
# tribuição desejada. Se a proposta tiver alta variância,
# será difícil para o algoritmo MH transitar entre os es-
# tados.

exp(logkpost_re(sam, x[1], exp(x[2]), # Massa probabilís-
                1/(1+exp(-x[3])),     # é bem pequena mas
                m, V, a, d))          # não nula.

# Vetores para receber amostras da distribuição desejada:

am1_mu_p = am1_s2_p = am1_nu_p = numeric(am1)
am2_mu_p = am2_s2_p = am2_nu_p = numeric(am2)
am3_mu_p = am3_s2_p = am3_nu_p = numeric(am3)

# Função para o algoritmo MCMC de Metropolis-Hastings:

XparMH = function(X, m, V, a, d, x, var_d,
                  up, seed) {
  k = length(up)
  am_mu_p = am_s2_p = am_nu_p = numeric(k)
  set.seed(122019)
  for(i in 1:k) {
    y = rmvnorm(1, mean=c(x[1], x[2], x[3]),
                sigma=diag(var_d))
    auy = exp(logkpost_re(sam, y[1], exp(y[2]),
      1/(1 + exp(-y[3])), m, V, a, d))
    aux = exp(logkpost_re(sam, x[1], exp(x[2]),
      1/(1 + exp(-x[3])), m, V, a, d))
    acep = min(c(1, auy/aux))
    if(up[i] <= acep) {
      am_mu_p[i] = y[1]
      am_s2_p[i] = exp(y[2])
      am_nu_p[i] = 1/(1+exp(-y[3]))
      x = y
    }
    else {
      am_mu_p[i] = x[1]
      am_s2_p[i] = exp(x[2])
      am_nu_p[i] = 1/(1+exp(-x[3]))
    }
  }
  
  return(list(XparMH_mu=am_mu_p,
              XparMH_s2=am_s2_p,
              XparMH_nu=am_nu_p))
}

am1_p = XparMH(X=sam, m=m, V=V, a=a, d=d, x=x,
               var_d=c(0.0022, 0.0065, 0.0203),
               up=up1, seed=122019)
am2_p = XparMH(X=sam, m=m, V=V, a=a, d=d, x=x,
               var_d=c(0.0022, 0.0065, 0.0203),
               up=up2, seed=122019)
am3_p = XparMH(X=sam, m=m, V=V, a=a, d=d, x=x,
               var_d=c(0.0022, 0.0065, 0.0203),
               up=up3, seed=122019)

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

st_sam_post(am1_p$XparMH_mu)
hist(am1_p$XparMH_mu, breaks=50, prob=T)
st_sam_post(am2_p$XparMH_mu)
hist(am2_p$XparMH_mu, breaks=50, prob=T)
st_sam_post(am3_p$XparMH_mu)
hist(am3_p$XparMH_mu, breaks=50, prob=T)

st_sam_post(am1_p$XparMH_s2)
hist(am1_p$XparMH_s2, breaks=50, prob=T)
st_sam_post(am2_p$XparMH_s2)
hist(am2_p$XparMH_s2, breaks=50, prob=T)
st_sam_post(am3_p$XparMH_s2)
hist(am3_p$XparMH_s2, breaks=50, prob=T)

st_sam_post(am1_p$XparMH_nu)
hist(am1_p$XparMH_nu, breaks=50, prob=T)
st_sam_post(am2_p$XparMH_nu)
hist(am2_p$XparMH_nu, breaks=50, prob=T)
st_sam_post(am3_p$XparMH_nu)
hist(am3_p$XparMH_nu, breaks=50, prob=T)

# Verificando a convergência e autocorrelação das cadeias
# geradas para cada parâmetro...