# Tamanho amostral:
n = 500

# Parâmetros do modelo para a distribuição amostral:
mu = 11
sigma2_1 = 9
sigma2_2 = 0.09
nu = 0.2

# Hiperparâmetros das distribuições a priori:
m = 11
V = 1
a_1 = 5; d_1 = 50
a_2 = 5; d_2 = 0.5

# Como gerar indicadoras da mistura dado nu:
Inu_pdf = function(n, nu, a = 0, b = 1) {
  unif = runif(n)
  #(unif < nu)*a + (unif >= nu)*b # Taiguara (complem.).
  (unif >= nu)*a + (unif < nu)*b  # Walmir (val. orig.).
}

mean(Inu_pdf(100, nu))
mean(Inu_pdf(1000, nu))
mean(Inu_pdf(10000, nu))

# Como gerar da variável auxiliar U (Dani):

U_pdf = function(n, nu) {
  u = sample(c(100, 1), prob=c(nu, 1-nu),
             size=n, replace=T)
}

table(U_pdf(100, nu))
table(U_pdf(1000, nu))
table(U_pdf(10000, nu))

# Gerando uma amostra de tamanho $n$:

# 1ª forma: gerando da própria densidade da mistura:

#nu_s = nu_pdf(n, nu)
#sample = nu_s*(mu + sigma2*100*rnorm(n)) +
#         (1 - nu_s)*(mu + sigma2*rnorm(n))

# Note que devemos usar todos os parâmetros do modelo,
# como se os conhecêssemos:

Inu = Inu_pdf(n, nu)
#sample1 = Inu*(mu + sqrt(sigma2/100)*rnorm(n)) +
#          (1 - Inu)*(mu + sqrt(sigma2)*rnorm(n))
#hist(sample, breaks=50, prob=T)
sample1_1 = Inu*(rnorm(n, mu, sqrt(sigma2_1/100))) +
            (1 - Inu)*(rnorm(n, mu, sqrt(sigma2_1)))
hist(sample1_1, breaks=50, prob=T)
sample1_2 = Inu*(rnorm(n, mu, sqrt(sigma2_2/100))) +
  (1 - Inu)*(rnorm(n, mu, sqrt(sigma2_2)))
hist(sample1_2, breaks=50, prob=T)

# 2ª forma: gerando de forma hierárquica (Dani):

u = U_pdf(n, nu)
sample2_1 = rnorm(n, mu, sqrt(sigma2_1/u))
hist(sample2_1, breaks=50, prob=T)
sample2_2 = rnorm(n, mu, sqrt(sigma2_2/u))
hist(sample2_2, breaks=50, prob=T)

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

A(sample2_1, mu, sigma2_1, nu)
A(sample2_2, mu, sigma2_2, nu)
hist(A(sample2_1, mu, sigma2_1, nu), prob=T, breaks=50)
hist(A(sample2_2, mu, sigma2_2, nu), prob=T, breaks=50)
prod(A(sample2_1, mu, sigma2_1, nu))
prod(A(sample2_2, mu, sigma2_2, nu))
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
hist(logA(sample2_1, mu, sigma2_1, nu),
     prob=T, breaks=50)
hist(logA(sample2_2, mu, sigma2_2, nu),
     prob=T, breaks=50)
sum(logA(sample2_1, mu, sigma2_1, nu))
sum(logA(sample2_2, mu, sigma2_2, nu))
# logA é negativa e tratável. A diferença é pequena mesmo
# se o parâmetro sigma2 variar bastante.

h = function(n, mu, sigma2, nu, m, V, a, d) {
  k1 = (1/sigma2)^((n + 1)/2 + a + 1)
  k2 = exp(-((mu - m)^2/(2*V) + d)/sigma2)
  k = k1*k2
  return(k)
}

h(n, mu, sigma2_1, nu, m, V, a_1, d_1)
h(n, mu, sigma2_2, nu, m, V, a_2, d_2)
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
hist(logkerpost_1, prob=T, breaks = 50)
hist(logkerpost_2, prob=T, breaks = 50)

kerpost_1 = exp(logkerpost_1)
kerpost_2 = exp(logkerpost_2)
hist(kerpost_1, prob=T, breaks = 50)
hist(kerpost_2, prob=T, breaks = 50)

# Método 1.1: Quadratura de Riemann:

L = 100 # Mesmo reduzindo, vê-se ao final que a aproxi-
        # mação é péssima.

# Para ter uma ideia da dispersão de mu:
hist(rnorm(n=10000, mean=m, sd=sqrt(V*sigma2_1)),
     prob=T, breaks=100)
hist(rnorm(n=10000, mean=m, sd=sqrt(V*sigma2_2)),
     prob=T, breaks=100)
# Massa probabilística concentrada quase que totalmente
# entre 1 e 21 (caso 1) e entre 10 e 12 (caso 2), ambas
# aproximadamente centradas em 11.

mu_step_1 = (21-1)/L
mu_gr_1 = seq(1, 21, mu_step_1)
mu_gr_1

mu_step_2 = (12-10)/L
mu_gr_2 = seq(10, 12, mu_step_2)
mu_gr_2

# Para ter uma ideia da dispersão de sigma2, se grande ou
# pequeno:
hist(1/rgamma(n=10000, shape=a_1, rate=d_1),
     prob=T, breaks=100)
hist(1/rgamma(n=10000, shape=a_2, rate=d_2),
     prob=T, breaks=100)
# Massa probabilística concentrada quase que totalmente
# entre 0 e 60 e aproximadamente centrada em 9 (caso 1)
# concentrada quase que totalmente entre 0 e 0.6 e apro-
# ximadamente centrada em 0.09 (caso 2).

sigma2_step_1 = 60/L
sigma2_gr_1 = seq(0.1, 60.1, sigma2_step_1)
sigma2_gr_1

sigma2_step_2 = 0.6/L
sigma2_gr_2 = seq(0.001, 0.601, sigma2_step_2)
sigma2_gr_2

# Como a priori de nu é U(0,1), é suficiente que partici-
# onemos este intervalo:

nu_step = 1/L
nu_gr = seq(0, 1, nu_step)
nu_gr

grid_tri_1 = cbind(mu_gr_1, sigma2_gr_1, nu_gr)
grid_tri_2 = cbind(mu_gr_2, sigma2_gr_2, nu_gr)
grid_tri_1
grid_tri_2

# Observe que no segundo grid temos variações bem menores
# para mu e sigma2, dado que sigma2 é pequeno.

l = nrow(grid_tri_1)
l

# Como os tamanhos dos passos são iguais em todas as três
# dimensões, podemos calcular o produto antes de começar
# a rotina iterativa:

prod_step_1 = mu_step_1*sigma2_step_1*nu_step
prod_step_2 = mu_step_2*sigma2_step_2*nu_step
prod_step_1
prod_step_2

# No caso trivariado, este produto será muito baixo mesmo
# se o número de intervalos for pequeno. O tamanho de ca-
# da intervalo necessariamente será menor que 1, indepen-
# dentemente da dimensão.

inv_k = 0 # Inicializando a soma que aproximará o inverso
          # da constante de proporcionalidade.

for (i in 1:l) {
  for (j in 1:l) {
    for (k in 1:l) {
      aux1=sum(logA(sample2_1, mu_gr_1[i],
                    sigma2_gr_1[j], nu_gr[k]))
      aux2=logh(n, mu_gr_1[i], sigma2_gr_1[j], nu_gr[k],
                m, V, a_1, d_1)
      aux = aux1 + aux2
      inv_k = inv_k + exp(aux)*prod_step_1
    }
  }
}

inv_k

# Evidentemente, a soma de inv_k será na prática nula se
# a variância do modelo (sigma2) for alta (a precisão do
# R não basta para detectar que cada parcela é diferente
# de zero). Neste caso, não podemos estimar inv_k corre-
# tamente.

inv_k = 0 # Inicializando a soma que aproximará o inverso
# da constante de proporcionalidade.

for (i in 1:l) {
  for (j in 1:l) {
    for (k in 1:l) {
      aux1=sum(logA(sample2_2, mu_gr_2[i],
                    sigma2_gr_2[j], nu_gr[k]))
      aux2=logh(n, mu_gr_2[i], sigma2_gr_2[j], nu_gr[k],
                m, V, a_2, d_2)
      aux = aux1 + aux2
      inv_k = inv_k + exp(aux)*prod_step_2
    }
  }
}

inv_k

# Agora, a soma de inv_k é na prática infinita, quando o
# valor de sigma2 for baixo. Logo, aqui também não pode-
# mos estimar inv_k corretamente.

# Em ambos os casos (variância grande ou pequena) a qua-
# dratura é demorada com muitas iterações, e mesmo assim
# não se obtém uma boa aproximação.

# Método 1.2: Quadratura do Trapézio:

