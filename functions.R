U_pdf = function(n, nu) {
  u = sample(c(100, 1), prob=c(nu, 1-nu),
             size=n, replace=T)
}

sampling = function(n, nu, mu, s2) {
  u = U_pdf(n, nu)
  sam = rnorm(n, mu, sqrt(s2*u))
}


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

# Por fim, para o logaritmo (do núcleo) da distribuição a
# posteriori, temos que:

logkpost = function(X, mu, sigma2, nu, m, V, a, d) {
  n = length(X)
  lA = logA(X, mu, sigma2, nu)
  lh = logh(n, mu, sigma2, nu, m, V, a, d)
  lkp = sum(lA) + lh
  return(lkp)
}


#calculo dos pesos da SIR, considerando multivariate

weight_calc = function(X, mu, sigma2, nu, m, V, a, d, mu_vec, sigma_vec) {
  
  
  
}
