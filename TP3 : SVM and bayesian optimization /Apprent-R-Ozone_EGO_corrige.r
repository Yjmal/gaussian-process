rm(list = ls())
source("EGOutils.R")
ozone <- read.table("depSeuil.dat", sep = ",", header = TRUE)
# Vérification du contenu
summary(ozone)

# Changement du type de la variable jour en facteur
ozone[, "JOUR"] <- as.factor(ozone[, "JOUR"])

par(mfrow = c(1, 2))
options(repr.plot.width = 6, repr.plot.height = 3)
hist(ozone[, "O3obs"])
hist(ozone[, "NO2"])

# Même chose pour les autres variables
# hist(ozone[,"MOCAGE"]);hist(ozone[,"TEMPE"]);hist(ozone[,"RMH2O"])
# hist(ozone[,"NO"]);hist(ozone[,"VentMOD"]);hist(ozone[,"VentANG"])

ozone[, "SRMH2O"] <- sqrt(ozone[, "RMH2O"])
ozone[, "LNO2"] <- log(ozone[, "NO2"])
ozone[, "LNO"] <- log(ozone[, "NO"])

ozone <- ozone[, c(1:4, 8:13)]
summary(ozone)
pairs(ozone[, c(3, 4, 6:10)])


set.seed(111) # initialisation du générateur
# Extraction des échantillons
test.ratio <- .2   # part de l'échantillon test
npop <- nrow(ozone) # nombre de lignes dans les données
nvar <- ncol(ozone) # nombre de colonnes
# taille de l'échantillon test
ntest <- ceiling(npop * test.ratio) 
# indices de l'échantillon test
testi <- sample(1:npop, ntest)
# indices de l'échantillon d'apprentissage
appri <- setdiff(1:npop, testi) 

# construction de l'échantillon d'apprentissage
datappr <- ozone[appri, -11] 
# construction de l'échantillon test
datestr <- ozone[testi, -11] 
summary(datappr) # vérification

# construction de l'échantillon d'apprentissage
datappq <- ozone[appri,-2]
# construction de l'échantillon test 
datestq <- ozone[testi,-2] 
summary(datappq) # vérification

# --------------------------------------------------------------

library(e1071)
# SVM avec reglages par defaut
svmDEF <- svm(O3obs ~ ., data = datappr)
svmDEF

# ----------
# QUESTION 1
# ----------
kFoldError <- function(param){
  set.seed(0)   # pour figer les "folds" dans la validation croisee
  obj <- tune.svm(O3obs ~ ., data = datappr,
                  gamma = param[1], 
                  cost = param[2], 
                  epsilon = param[3])
  return(obj$best.performance)
}

parDEF <- c(svmDEF$gamma, svmDEF$cost, svmDEF$epsilon)
cat("\n(k-fold) error (default model):", 
    round(kFoldError(parDEF), 3),
    "\nCorresponding parameter values:", round(parDEF, 2))

# ----------
# QUESTION 2
# ----------
# Representer la fonction cout contre chaque variable
# En deduire un domaine de recherche pour cost, gamma et epsilon
# Valeur de la fonction cout pour (costOpt, gammaOpt, epsilonOpt) ?

svmTuneGamma <- tune.svm(O3obs ~ ., data = datappr,
                         gamma = seq(from = 0.01, by = 0.03, length.out = 6))
svmTuneCost <- tune.svm(O3obs ~ ., data = datappr,
                        cost = seq(from = 1, by = 0.5, length.out = 6))
svmTuneEpsilon <- tune.svm(O3obs ~ ., data = datappr,
                           epsilon = seq(from = 0, by = 0.1, length.out = 6))
gammaOpt <- as.numeric(svmTuneGamma$best.parameters)
costOpt <- as.numeric(svmTuneCost$best.parameters)
epsilonOpt <- as.numeric(svmTuneEpsilon$best.parameters)

par(mfrow = c(1, 3))
plot(svmTuneGamma)
plot(svmTuneCost)
plot(svmTuneEpsilon)
parMAR <- c(gammaOpt, costOpt, epsilonOpt)

cat("\n(k-fold) error (marginal optmization):", 
    round(kFoldError(parMAR), 3),
    "\nCorresponding parameter values:", round(parMAR, 2))


# ----------
# QUESTION 3
# ----------

n0 <- 10   # budget initial

library(DiceDesign)
X <- lhsDesign(n = n0, dimension = 3)$design
Xopt <- maximinESE_LHS(X, T0 = 0.005 * phiP(X), 
                       inner_it = 100, J = 50, it = 10)
X0 <- Xopt$design

lower <- c(0.01, 0.05, 0)
upper <- c(0.5, 3, 0.5)

for (i in 1:ncol(X0)){
  X0[, i] <- lower[i] + X0[, i] * (upper[i] - lower[i]) 
}
pairs(X0)

y0 <- apply(X0, 1, kFoldError)  
bestInitPar <- X0[which.min(y0), ]
cat("\nSmallest (k-fold) error found on X0:", round(min(y0), 3),
    "\nCorresponding parameter values:", round(bestInitPar, 2))



library(DiceOptim)
km0 <- km(~ 1, design = X0, response = y0)
km0         
plot(km0)  

options(warn = - 1)
ego <- EGO.nsteps(km0, kFoldError, nsteps = 20, 
                  lower = lower, upper = upper)

EGOpar <- ego$par
EGOvalue <- ego$value
bestPoint <- which.min(EGOvalue)
bestPar <- EGOpar[bestPoint, ]

cat("\nSmallest (k-fold) error found on additional EGO steps:", 
    round(EGOvalue[bestPoint], 3),
    "\nCorresponding parameter values:", round(EGOpar[bestPoint, ], 2))

# visualisation des points obtenus par EGO
visualizeEGO(initDesign = X0, initValues = y0,
             EGOpoints = EGOpar, EGOvalues = EGOvalue)
convergenceEGO(initValues = y0, EGOvalues = EGOvalue)


# ----------
# QUESTION 5 
# ----------
# a) 
cat("\n(k-fold) error (default model):", 
    round(kFoldError(parDEF), 3),
    "\nCorresponding parameter values:", round(parDEF, 2))
cat("\n(k-fold) error (marginal optmization):", 
    round(kFoldError(parMAR), 3),
    "\nCorresponding parameter values:", round(parMAR, 2))
cat("\nSmallest (k-fold) error found on X0:", round(min(y0), 3),
    "\nCorresponding parameter values:", round(X0[which.min(y0), ], 2))
cat("\nSmallest (k-fold) error found on additional EGO steps:", round(EGOvalue[bestPoint], 3),
    "\nCorresponding parameter values:", round(EGOpar[bestPoint, ], 2))


# b)
# Si l'optimisation globale revient a une optimisation/ chaque parametre,
# c'est le signe qu'il n'y a pas d'interaction : fonction additive de la forme
# J(x1, x2, x3) = J1(x1) + J2(x2) + J3(x3)


# ----------
# QUESTION 6 
# ----------
# Que donnent les 4 modeles regles en prediction ? (sur l'ensemble test)
# Que donnent les 4 modeles pour predire le depassement de seuil (> 150) ?

# Creation des modeles
svmMAR <- svm(O3obs ~ ., data = datappr, 
              gamma = parMAR[1], 
              cost = parMAR[2], 
              epsilon = parMAR[3])
svmLHS <- svm(O3obs ~ ., data = datappr, 
              gamma = bestInitPar[1], 
              cost = bestInitPar[2], 
              epsilon = bestInitPar[3])
svmEGO <- svm(O3obs ~ ., data = datappr, 
              gamma = bestPar[1], 
              cost = bestPar[2], 
              epsilon = bestPar[3])

# Erreur quadratique moyenne de prévision
pred.svmr <- predict(svmDEF, newdata = datestr)
RMSEsvmDEF <- mean((predict(svmDEF, newdata = datestr) - datestr[, "O3obs"])^2)
pred.svmMAR <- predict(svmMAR, newdata = datestr)
RMSEsvmPAR <- mean((predict(svmMAR, newdata = datestr) - datestr[, "O3obs"])^2)
pred.svmLHS <- predict(svmLHS, newdata = datestr)
RMSEsvmLHS <- mean((predict(svmLHS, newdata = datestr) - datestr[, "O3obs"])^2)
pred.svmEGO <- predict(svmEGO, newdata = datestr)
RMSEsvmEGO <- mean((predict(svmEGO, newdata = datestr) - datestr[, "O3obs"])^2)
cat("\nRMSE with SVM (default parameters):", RMSEsvmDEF)
cat("\nRMSE with SVM (marginal optimization):", RMSEsvmPAR)
cat("\nRMSE with SVM (LHS design):", RMSEsvmLHS)
cat("\nRMSE with SVM (additional EGO):", RMSEsvmEGO)

# Matrice de confusion pour la prévision du dépassement de seuil (régression)
print(table(pred.svmr > 150, datestr[, "O3obs"] > 150))
print(table(pred.svmMAR > 150, datestr[, "O3obs"] > 150))
print(table(pred.svmLHS > 150, datestr[, "O3obs"] > 150))
print(table(pred.svmEGO > 150, datestr[, "O3obs"] > 150))

