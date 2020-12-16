library(regressionMultiple)
# Un premier exemple
#Génération des données 
#Soit un modèle linéaire de la forme $Y = a_1X_1 + a_2X_2 + a_3X_3 + ...+ a_nX_n + \epsilon$
#On tire epsilon:
mu <- 0
sigma <- 22 #bruit
n= 1000 #nombre de tirages
eps <- rnorm(n, mu,sigma) #vecteur des n epsilons

#Initialisation des paramètres pour l'exemple
p<- 3 #nombre de paramètres
a1 <- 5 #valeurs du premier paramètre a1, les ai sont des multiples de a1
a_i <- a1*(1:p) 
coeff <- a_i #p coefficients

# On initialise les variables X1, X2, X3,...Xn
X <- matrix(0, ncol=p, nrow=n)
for(i in 1:p){
  X[,i] <- runif(n)
}

# Générer un vecteur Y tel que Y = a1X1 + a2X2 + a3X3 +...+ anXn + eps

Y <- (X %*% coeff) + eps
#Y

# Estimation

hat <- solve(t(X)%*%X)%*%t(X)%*%Y
hat
coeff

# Comparaison estimation - coefficients réels

seuil <-2 #choix du seuil
identite(coeff,hat,seuil)
pourcentage <- identite(coeff,hat,seuil)/length(coeff) #diff,coeff et hat ont le même nombre de lignes
pourcentage

mean(order(hat) == order(coeff))

# Génération d'un modèle linéaire lm pour notre exemple

model<-lm(Y~X)
summary(model)

# Simulations

# Initialisation des paramètres fixes

# Clear the console & the R environment
cat("\014")
rm(list = ls())

N<-20 #nombre de simulations
n= 100 #nombre de tirages par simulations
sigma <- 22 #valeur de sigma choisie pour minimiser le bruit
seuil <- 2 #valeur de seuil pour tester l'exactitude de l'estimation des paramètres

p<- 5 #nombre de coefficients
a1 <- 5 #valeurs du premiers coefficient a1, les ai sont des multiples de a1
a_i <- a1*(1:p) 
coeff <- a_i #p coefficients



# On effectue des simulations
res<- replicate(N, simulation(sigma,coeff,n,seuil))
res

# Affichage de ces résultats

resultat_simu(res)

# Pourquoi ce résultat?

M <- 200
res<- replicate(M, histog_simu(sigma,coeff,n))
res_histo <- do.call(rbind, res)
library(reshape2)
library(ggplot2)
gg <- melt(res_histo)
ggplot(gg, aes(x=value, fill=Var2)) +
  geom_histogram(binwidth=1)+
  facet_grid(Var2~.)

# Faire varier sigma ( choix d'un sigma inférieur à la valeur précédente)

M <- 200
sigma_inf <- 3
res<- replicate(M, histog_simu(sigma_inf,coeff,n))
res_histo <- do.call(rbind, res)
library(reshape2)
library(ggplot2)
gg <- melt(res_histo)
ggplot(gg, aes(x=value, fill=Var2)) +
  geom_histogram(binwidth=0.3)+
  facet_grid(Var2~.)

# Application avec les données abalone

library("AppliedPredictiveModeling")
data(abalone)
matrice <- scale(as.matrix(abalone[,2:8]))
modele_lineaire <- lm(abalone$Rings ~ matrice)
modele_lineaire <- summary(modele_lineaire)
print(modele_lineaire)

# Calcul/estimation de epsilon = Y - X(XtX)-1XtY

#Y_abalone <- abalone$Rings #non centré
#Xi_abalone <- as.matrix(abalone[,2:8]) #les valeurs ne sont pas centrées
Y_abalone <- scale(abalone$Rings) #centré
Xi_abalone <- as.matrix(scale(abalone[,2:8])) #les valeurs sont centrées
eps_abalone <- Y_abalone - Xi_abalone%*%solve(t(Xi_abalone)%*%Xi_abalone) %*%t(Xi_abalone)%*%Y_abalone 
sigma_abalone <- sd(eps_abalone)
sigma_abalone


coeff1_abalone <- unname(modele_lineaire$coefficients[,1]) #extraction des valeurs "estimate"
coeff_abalone <- coeff1_abalone[2:length(coeff1_abalone)] #extraction des valeurs "estimate" en enlevant intercept
#print(coeff_abalone)
faux_coeff <- seq(min(coeff_abalone), max(coeff_abalone), length.out = length(coeff_abalone))
M <- 200
res<- replicate(M, histog_simu(sigma_abalone,faux_coeff,dim(abalone)[1]))
res_histo <- do.call(rbind, res)
library(reshape2)
library(ggplot2)
gg <- melt(res_histo)
ggplot(gg, aes(x=value, fill=Var2)) +
  geom_histogram(binwidth=0.3)+
  facet_grid(Var2~.)


# Comparaison avec les vraies valeurs coeff_abalone

M <- 200
res<- replicate(M, histog_simu(sigma_abalone,coeff_abalone,dim(abalone)[1]))
res_histo <- do.call(rbind, res)
library(reshape2)
library(ggplot2)
gg <- melt(res_histo)
ggplot(gg, aes(x=value, fill=Var2)) +
  geom_histogram(binwidth=0.3)+
  facet_grid(Var2~.)

