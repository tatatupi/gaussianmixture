# Trabalho - Inferência Bayesiana


Considere uma amostra **X** = (X<sub>1</sub>, ..., X<sub>N</sub>) independente e identicamente distribuída da seguinte mistura finita de distribuições normais:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;f(x|\mu,\sigma^2,\nu)=\nu\phi(x|\mu,100\sigma²)+(1-\nu)\phi(x|\mu,\sigma²)x\in\mathcal{R}"/>

Onde <img src="https://latex.codecogs.com/svg.latex?\Large&space;\phi(x|M,V)"/> denota a função densidade de probabilidade da distribuição normal com média M e variância V avaliada no 
ponto <img src="https://latex.codecogs.com/svg.latex?\Large&space;x,\mu\in\mathcal{R},\sigma^2\in\mathcal{R^+}"/> e <img src="https://latex.codecogs.com/svg.latex?\Large&space;\nu\in(0,1)"/>

1. Assuma que, *a priori*, <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mu|\sigma²\sim{N}(m,V\sigma²),\sigma²\sim{GI}(a,d)"/> e <img src="https://latex.codecogs.com/svg.latex?\Large&space;\nu\sim{U}(0,1)"/>. Encontre a expressão da densidade *a posteriori* de (μ,σ²,ν).
1. Escolha valores apropriados para (μ,σ²,ν). Fixado estes parâmetros, gere uma amostra de tamanho n = 500 da distribuição amostral. Armazene estes valores para análises posteriores. Faça uma apresentação gráfica apropriada da amostra gerada. Determine a média, variância, curtose e assimetria da amostra gerada.


DICA: pode-se usar, sem provar, que a mistura finita de distribuições normais pode ser hierarquicamente representada por 

<img src="https://latex.codecogs.com/svg.latex?\Large&space;X_i|\mu,\sigma^2,\U_i=u_i)\sim{N}(\mu,\sigma^2u_i^{-1})"/> e <img src="https://latex.codecogs.com/svg.latex?\Large&space;U_i|\mu\sim{discreta}(1,100)"/>, com <img src="https://latex.codecogs.com/svg.latex?\Large&space;P(U_i=100)=\nu"/>
