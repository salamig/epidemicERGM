# EpidemicERGM

library(statnet)
library(EpiModel)
library(ndtv)

# 1) Network with homophily only by race
# initialize empty network with desired attributes
net <- network.initialize(n = 30, directed = FALSE)
net <- set.vertex.attribute(net, "race", rep(0:2, each = 10))

# Setting model parameters
# NOTE: this is where we can add attributes that will define structure in the way we like based on an ERGM formula
formation <- ~edges + nodefactor("race") + nodematch("race")
target.stats <- c(30, 15, 30, 26) # correspond to formation formula above

# because networks generated are dynamic, we also set parameters to determine the rate of dissolution of ties
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)

# Fitting model with desired parameters
mod1 <- netest(net, formation, target.stats, coef.diss, edapprox = TRUE)

# Model diagnostics
dx <- netdx(mod1, nsims = 5, nsteps = 500,
            nwstats.formula = ~edges + nodefactor("race", base = 0) +
              nodematch("race"))
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(dx)


# Simulate epidemic
param <- param.net(inf.prob = 0.03, act.rate = 1)
status.vector <- rbinom(30, 1, 0.02)
status.vector <- ifelse(status.vector == 1, "i", "s")
init <- init.net(status.vector = status.vector)

# Structural model features
control <- control.net(type = "SI", nsteps = 500, nsims = 10, epi.by = "race", ncores=2)

# simulate!
sim1 <- netsim(mod1, param, init, control)

# Plots of how the epidemic spread
plot(sim1, mean.line = FALSE, qnts = FALSE, sim.lines = TRUE)

# Plots of how the simulated dynamic networks looked at different timepoints
par(mfrow = c(1,2), mar = c(0,0,1,0))
plot(sim1, type = "network", at = 1, col.status = TRUE,
     main = "Prevalence at t1")
plot(sim1, type = "network", at = 50, col.status = TRUE,
     main = "Prevalence at t500")


# plot networks at specific points showing race
nw <- get_network(sim1, sim = 1)
out <- network.extract(nw, at=1)
plot(out, vertex.col="race")

# Make an animated plot of the networks over a specific duration
slice.par<-list(start=1, end=20, interval=1, aggregate.dur=1,rule="latest")
compute.animation(nw,slice.par=slice.par)

render.d3movie(nw,
               edge.col="darkgray",displaylabels=TRUE,
               label.cex=.6,label.col="blue",
               output.mode = 'htmlWidget',
               vertex.col = "race")



#### Non-temporal networks (to figure out good ways to parameterize different social processes)
par(mfrow=c(3,3))
for(i in seq(0,5, length.out = 9)){
  net <- network.initialize(n = 45, directed = FALSE)
  net <- set.vertex.attribute(net, "race", rep(0:1, each = 25))
  g.sim <- simulate(net ~ edges + nodematch("race"), coef = c(-5, i))
  plot(g.sim, main=paste("Homophily = ", i, sep=""))
}

net <- network.initialize(n = 50, directed = FALSE)
net <- set.vertex.attribute(net, "race", rep(0:1, each = 25))

g.sim <- simulate(net ~ edges , coef = c(-3))

plot(g.sim, vertex.col="race")












### Repeat basic model but on many nodes to see how diagnostics and epidemic spreads are functioning

# 1) Network with homophily only by race
# initialize empty network with desired attributes
net <- network.initialize(n = 1500, directed = FALSE)
net <- set.vertex.attribute(net, "race", rep(0:2, each = 500))

# Setting model parameters
formation <- ~edges + nodefactor("race") + nodematch("race")
target.stats <- c(2250, 375, 500, 1500) # correspond to formation formula above

coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)

# Fitting model with desired parameters
mod1 <- netest(net, formation, target.stats, coef.diss, edapprox = TRUE)

# Model diagnostics
dx <- netdx(mod1, nsims = 5, nsteps = 200,
            nwstats.formula = ~edges + nodefactor("race", base = 0) +
              nodematch("race"))
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(dx)


# Simulate epidemic
param <- param.net(inf.prob = 0.02, act.rate = 1)
status.vector <- rbinom(1500, 1, 0.02)
status.vector <- ifelse(status.vector == 1, "i", "s")
init <- init.net(status.vector = status.vector)

# Structural model features
control <- control.net(type = "SI", nsteps = 100, nsims = 10, epi.by = "race", ncores=2)

# simulate!
sim1 <- netsim(mod1, param, init, control)

plot(sim1, mean.line = FALSE, qnts = FALSE, sim.lines = TRUE)
par(mfrow = c(1,2), mar = c(0,0,1,0))
plot(sim1, type = "network", at = 1, col.status = TRUE,
     main = "Prevalence at t1")
plot(sim1, type = "network", at = 50, col.status = TRUE,
     main = "Prevalence at t500")


# plot single networks
nw <- get_network(sim1, sim = 1)
out <- network.extract(nw, at=1)
plot(out, vertex.col="race")

nw <- get_network(sim1, sim = 1)

slice.par<-list(start=1,end=20,interval=1,
                aggregate.dur=1,rule="latest")
compute.animation(nw,slice.par=slice.par)

render.d3movie(nw,
               edge.col="darkgray",displaylabels=TRUE,
               label.cex=.6,label.col="blue",
               output.mode = 'htmlWidget',
               vertex.col = "race")

plot(out, vertex.col="race")


par(mfrow=c(3,3))
for(i in seq(0,5, length.out = 9)){
  net <- network.initialize(n = 45, directed = FALSE)
  net <- set.vertex.attribute(net, "race", rep(0:1, each = 25))
  g.sim <- simulate(net ~ edges + nodematch("race"), coef = c(-5, i))
  plot(g.sim, main=paste("Homophily = ", i, sep=""))
}

net <- network.initialize(n = 50, directed = FALSE)
net <- set.vertex.attribute(net, "race", rep(0:1, each = 25))

g.sim <- simulate(net ~ edges , coef = c(-3))

plot(g.sim, vertex.col="race")
