# Gerando da variável auxiliar U:

U_pdf = function(n, nu) {
  u = sample(c(100, 1), prob=c(nu, 1-nu),
             size=n, replace=T)
}

#table(U_pdf(100, nu))
#table(U_pdf(1000, nu))
#table(U_pdf(10000, nu))

# Definindo a distribuição *a posteriori*:

# Na função a seguir, calculamos o logaritmo do produtório
# apenas da soma envolvendo os 2 componentes da mistura:

logA = function(X, mu, sigma2, nu) {
  n = length(X)
  k1 = (nu/10)*exp(-(X-mu)^2/(2*100*sigma2))
  k2 = (1-nu)*exp(-(X-mu)^2/(2*sigma2))
  k = log(k1 + k2)
  return(k)
}

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

# Função para calcular logaritmo do produto de termos ad-
# vindos do determinante do jacobiano para reparametriza-
# ção da distribuição a posteriori:

logr = function(sigma2, nu) {
  k = log(sigma2) + 3*log(nu) - log(1-nu)
  return(k)
}

# Logaritmo (do núcleo) da distribuição a posteriori:

logkpost = function(X, mu, sigma2, nu, m, V, a, d) {
  n = length(X)
  lA = logA(X, mu, sigma2, nu)
  lh = logh(n, mu, sigma2, nu, m, V, a, d)
  lkp = sum(lA) + lh
  return(lkp)
}

# Logaritmo (do núcleo) da distribuição a posteriori, com
# reparametrização:

logkpost_re = function(X, mu, sigma2, nu, m, V, a, d) {
  n = length(X)
  lA = logA(X, mu, sigma2, nu)
  lh = logh(n, mu, sigma2, nu, m, V, a, d)
  lr = logr(sigma2, nu)
  lkrp = sum(lA) + lh + lr
  return(lkrp)
}