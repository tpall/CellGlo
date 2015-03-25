library(rjags)

set.seed(1337)
y <- rnorm(n = 20, mean = 10, sd = 5)
mean(y)
sd(y)
# The model specification
model_string <- "model{
for(i in 1:length(y)) {
y[i] ~ dnorm(mu, tau)
}
mu ~ dnorm(0, 0.0001)
sigma ~ dlnorm(0, 0.0625)
tau <- 1/pow(sigma, 2)
}"

# Running the model
model <- jags.model(textConnection(model_string), data = list(y = y), n.chains = 3, n.adapt= 10000)
update(model, 10000); # Burnin for 10000 samples
mcmc_samples <- coda.samples(model, variable.names=c("mu", "sigma"), n.iter=20000)
plot(mcmc_samples)
summary(mcmc_samples)

# 4 parameter curve
x <- 0:20
y <- 20 + (2 - 20)/(1 + (x/10)^5) + rnorm(21, sd=2)
dataList = list(y = y, x = x)
plot(y~x)
model_string <- "
model {
  for( i in 1:length(y) ) {
    y[i] ~ dnorm( mu[i] , tau )
    mu[i] <- upAsym + (y0 - upAsym)/ (1 + pow(x[i]/midPoint, slope))
  }
  tau ~ dgamma(sG ,rG )
  sG <- pow(m,2)/pow(d,2)
  rG <- m/pow(d,2)
  m ~ dgamma(1, 0.01)
  d ~ dgamma(1, 0.01)

  midPoint ~ dnorm(10, 0.0001) T(0,21)
  slope    ~ dnorm(5, 0.0001) T(0,)
  upAsym   ~ dnorm(30, 0.0001) T(0,40)
  y0       ~ dnorm(0, 0.0001) T(-20,20)
}"

model <- jags.model(textConnection(model_string), data = dataList, n.chains = 3, n.adapt= 10000)
update(model, 10000); # Burnin for 10000 samples
mcmc_samples <- coda.samples(model, variable.names=c("mu", "tau"), n.iter=20000)
plot(mcmc_samples)
summary(mcmc_samples)

# let's try my data
tec %>%
  filter(!treatment=="media") %>%
  filter(doses<1.50001e-05) %>%
  group_by(date,plate,GF) %>% 
  mutate(value=(value-min(value))/(range(value)%>%diff)) %>%
  filter(treatment%in%c("rhIgG-Fc","3MUT-Fc")) %>%
  mutate(value=value/median(value[treatment=="rhIgG-Fc"])) %>%
  group_by(date,GF,doses,treatment,treat2) %>%
  summarise(value=mean(value)) %>%
  filter(GF=="bFGF")

model_string <- "
model {
  for( i in 1:length(y) ) {
    y[i] ~ dnorm( mu[i] , tau )
    mu[i] <- upAsym + (y0 - upAsym)/ (1 + pow(x[i]/midPoint, slope))
  }
  tau ~ dgamma(sG ,rG )
  sG <- pow(m,2)/pow(d,2)
  rG <- m/pow(d,2)
  m ~ dgamma(1, 0.01)
  d ~ dgamma(1, 0.01)

  midPoint ~ dnorm(10, 0.0001) T(0,21)
  slope    ~ dnorm(5, 0.0001) T(0,)
  upAsym   ~ dnorm(30, 0.0001) T(0,40)
  y0       ~ dnorm(0, 0.0001) T(-20,20)
}"

model <- jags.model(textConnection(model_string), data = dataList, n.chains = 3, n.adapt= 10000)
update(model, 10000); # Burnin for 10000 samples
mcmc_samples <- coda.samples(model, variable.names=c("mu", "tau"), n.iter=20000)
plot(mcmc_samples)
summary(mcmc_samples)