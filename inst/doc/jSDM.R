## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
fig.align = "center",
fig.retina = 1,
fig.width = 6, fig.height = 6,
cache = FALSE,
collapse = TRUE,
comment = "#>",
highlight = TRUE
)

## ----plot-DAG, echo=FALSE, out.width=500, fig.retina=1, fig.cap="(ref:cap-DAG)", eval=TRUE----
knitr::include_graphics("figures/DAG-jSDM-rand.png")

## ----libraries----------------------------------------------------------------
# Load libraries
library(jSDM)

## ----frog-picture, echo=FALSE, out.width=400, out.height=300, fig.retina=1, fig.cap="(ref:cap-frog)"----
knitr::include_graphics("figures/Litoria_ewingii.jpg")

## ----frogs-data---------------------------------------------------------------
# frogs data
data(frogs, package="jSDM")
head(frogs)

## ----arranging-frogs-data-----------------------------------------------------
# data.obs
PA_frogs <- frogs[,4:12]

# Normalized continuous variables
Env_frogs <- cbind(scale(frogs[,1]),frogs[,2],scale(frogs[,3]))
colnames(Env_frogs) <- colnames(frogs[,1:3])

## ----jSDM-probit--------------------------------------------------------------
mod_frogs_jSDM_probit <- jSDM_binomial_probit(
  # Chains
  burnin=1000, mcmc=1000, thin=1,
  # Response variable 
  presence_data = PA_frogs, 
  # Explanatory variables 
  site_formula = ~.,   
  site_data = Env_frogs,
  # Model specification 
  n_latent=2, site_effect="random",
  # Starting values
  alpha_start=0, beta_start=0,
  lambda_start=0, W_start=0,
  V_alpha=1, 
  # Priors
  shape_Valpha=0.1,
  rate_Valpha=0.1,
  mu_beta=0, V_beta=1,
  mu_lambda=0, V_lambda=1,
  # Various 
  seed=1234, verbose=1)

## ----plot-results-probit------------------------------------------------------
np <- nrow(mod_frogs_jSDM_probit$model_spec$beta_start)

## beta_j of the first two species
par(mfrow=c(2,2))
for (j in 1:2) {
  for (p in 1:np) {
    coda::traceplot(coda::as.mcmc(mod_frogs_jSDM_probit$mcmc.sp[[j]][,p]))
    coda::densplot(coda::as.mcmc(mod_frogs_jSDM_probit$mcmc.sp[[j]][,p]), 
                   main = paste(colnames(mod_frogs_jSDM_probit$mcmc.sp[[j]])[p],
                                ", species : ",j), cex.main=0.9)
  }
}

## lambda_j of the first two species
n_latent <- mod_frogs_jSDM_probit$model_spec$n_latent
par(mfrow=c(2,2))
for (j in 1:2) {
  for (l in 1:n_latent) {
    coda::traceplot(coda::as.mcmc(mod_frogs_jSDM_probit$mcmc.sp[[j]][,np+l]))
    coda::densplot(coda::as.mcmc(mod_frogs_jSDM_probit$mcmc.sp[[j]][,np+l]), 
                   main = paste(colnames(mod_frogs_jSDM_probit$mcmc.sp[[j]])
                                [np+l],", species : ",j), cex.main=0.9)
  }
}

## Latent variables W_i for the first two sites
par(mfrow=c(2,2))
for (l in 1:n_latent) {
  for (i in 1:2) {
  coda::traceplot(mod_frogs_jSDM_probit$mcmc.latent[[paste0("lv_",l)]][,i],
                  main = paste0("Latent variable W_", l, ", site ", i),
                  cex.main=0.9)
  coda::densplot(mod_frogs_jSDM_probit$mcmc.latent[[paste0("lv_",l)]][,i],
                 main = paste0("Latent variable W_", l, ", site ", i),
                 cex.main=0.9)
  }
}

## alpha_i of the first two sites
plot(coda::as.mcmc(mod_frogs_jSDM_probit$mcmc.alpha[,1:2]))

## V_alpha
par(mfrow=c(2,2))
coda::traceplot(mod_frogs_jSDM_probit$mcmc.V_alpha)
coda::densplot(mod_frogs_jSDM_probit$mcmc.V_alpha)
## Deviance
coda::traceplot(mod_frogs_jSDM_probit$mcmc.Deviance)
coda::densplot(mod_frogs_jSDM_probit$mcmc.Deviance)

## probit_theta
par (mfrow=c(1,2))
hist(mod_frogs_jSDM_probit$probit_theta_latent,
     main = "Predicted probit theta",
     xlab ="predicted probit theta")
hist(mod_frogs_jSDM_probit$theta_latent,
     main = "Predicted theta", 
     xlab ="predicted theta")

## ----correlation-matrix-probit------------------------------------------------
plot_residual_cor(mod_frogs_jSDM_probit)

## ----predictions-probit-------------------------------------------------------
# Sites and species concerned by predictions :
## 50 sites among the 104
Id_sites <- sample.int(nrow(PA_frogs), 50)
## All species 
Id_species <- colnames(PA_frogs)
# Simulate new observations of covariates on those sites 
simdata <- matrix(nrow=50, ncol = ncol(mod_frogs_jSDM_probit$model_spec$site_data))
colnames(simdata) <- colnames(mod_frogs_jSDM_probit$model_spec$site_data)
rownames(simdata) <- Id_sites
simdata <- as.data.frame(simdata)
simdata$Covariate_1 <- rnorm(50)
simdata$Covariate_3 <- rnorm(50)
simdata$Covariate_2 <- rbinom(50,1,0.5)

# Predictions 
theta_pred <- predict(mod_frogs_jSDM_probit, newdata=simdata, Id_species=Id_species,
                      Id_sites=Id_sites, type="mean")
hist(theta_pred, main="Predicted theta with simulated data", xlab="predicted theta")

