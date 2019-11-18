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

# Escolhemos como densidade conhecida a normal multivari-
# ada, com k=3. Desta forma, como nem todos os argumentos
# da posteriori variam na reta real, podemos reparametri-
# zá-los de tal forma que todos os suportes dos novos pa-
# râmetros coincidam com os das médias na normal multiva-
# riada. Logo, a reparametrização será feita tomando fun-
# ções que levem sigma e nu de seus domínios para imagens
# na reta real (funções logarítmica e logística).

logkpost_re(sam, mu, s2, nu, m, V, a, d)

# Embora não seja obrigatório fazer isso para o algoritmo
# Metropolis-Hastings, a eficiência do mesmo aumenta con-
# sideravelmente quando a distribuição desejada é repara-
# metrizada de modo que os suportes das quantidades alea-
# tórias a serem amostradas (parâmetros na distribuição a
# posteriori) sejam similares aos das quantidades aleató-
# rias da distribuição proposta (as três componentes para
# a distribuição normal multivariada).

# Tamanhos para a amostragem da posteriori via MCMC/MH:

am1 = 500; am2 = 5000; am3 = 50000

# Devido à simetria da distribuição proposta (normal mul-
# tivariada), a razão de teste do algoritmo MH é reduzida
# às densidades da distribuição desejada nos pontos gera-
# dos consecutivamente da distribuição proposta.

# Para o burn-in na geração da cadeia em cada cenário, i-
# remos desprezar os 100 valores iniciais no primeiro ce-
# nário, os 1000 no segundo e os 10000 no terceiro.

# Agora, geremos a partir da distribuição normal multiva-
# riada com média de cada componente dada pelo respectivo
# valor no modelo e variância dada pela matriz identidade
# de ordem 3:

library(mvtnorm) # Pacote com funções densidade e gerado-
                 # ra da distrib. normal multivariada.

# Geremos de uma uniforme padrão para controlar a aceita-
# ção de cada amostra gerada da distribuição proposta:

set.seed(122019); up1 = runif(am1)
set.seed(122019); up2 = runif(am2)
set.seed(122019); up3 = runif(am3)

# Vimos no método da quadratura de Riemann que, fixados 2
# parâmetros:

# A distrib. de mu concentra-se entre 10.85 e 11.13;
# A distrib. de sigma2 concentra-se entre 0.48 e 0.78;
# A distrib. de nu concentra-se entre 0.13 e 0.26.

# Temos então o seguinte vetor de valores iniciais para o
# algoritmo MH, os quais estão distantes dos valores ver-
# dadeiros no modelo da distribuição desejada, embora com
# probabilidade positiva:

x = c(10.86, log(0.50), log(0.14/0.86))

exp(logkpost_re(sam, x[1], exp(x[2]), # Massa probabilís-
                1/(1+exp(-x[3])),     # é bem pequena mas
                m, V, a, d))          # não nula.

# Vetores para receber amostras da distribuição desejada:

am1_mu_p = am1_s2_p = am1_nu_p = numeric(am1)
am2_mu_p = am2_s2_p = am2_nu_p = numeric(am2)
am3_mu_p = am3_s2_p = am3_nu_p = numeric(am3)

# O algoritmo MCMC de Metropolis-Hastings:

set.seed(122019)
for(i in 1:am1) {
  y = rmvnorm(1, mean=c(x[1],x[2],x[3]), sigma=diag(3))
  auy = exp(logkpost_re(sam, y[1], exp(y[2]),
                        1/(1 + exp(-y[3])), m, V, a, d))
  aux = exp(logkpost_re(sam, x[1], exp(x[2]),
                        1/(1 + exp(-x[3])), m, V, a, d))
  acep = min(c(1, auy/aux))
  if(acep >= u[i]) {
    am1_mu_p[i] = y[1]
    am1_s2_p[i] = exp(y[2])
    am1_nu_p[i] = 1/(1+exp(-y[3]))
    x = y
  }
  else {
    am1_mu_p[i] = x[1]
    am1_s2_p[i] = exp(x[2])
    am1_nu_p[i] = 1/(1+exp(-x[3]))
  }
}

# Teste e implemente depois como função...