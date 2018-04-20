## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(fig.width = 7)

## ---- message = FALSE----------------------------------------------------
library(disclap)
library(disclapmix)

## ------------------------------------------------------------------------
set.seed(1)

## ---- fig.keep = 'high', fig.cap="The probability mass function, $f(X = x; p, y)$, for the discrete Laplace distribution with dispersion parameter $p=0.3$ and location parameter $y=13$ from $x=8$ to $x=18$."----
p <- 0.3
y <- 13L # L for integer type
x <- seq(8L, 18L, by = 1L)
barplot(ddisclap(x - y, p), names = x, xlab = "x, e.g. Y-STR allele", 
  ylab = paste("Probability mass, f(X = x; ", p, ", ", y, ")", sep = ""))

## ------------------------------------------------------------------------
sum(ddisclap(x - y, p))

## ------------------------------------------------------------------------
set.seed(1) # Makes it possible to reproduce the simulation results
p <- 0.3 # Dispersion parameter 
y <- 13 # Location parameter 
x <- rdisclap(100, p) + y # Generate a sample using the rdisclap function

y_hat <- median(x)
y_hat
mu_hat <- mean(abs(x - y_hat))
mu_hat
p_hat <- mu_hat^(-1) * (sqrt(mu_hat^2 + 1) - 1)
p_hat # We expect 0.3

# The observed distribution of d's
tab <- prop.table(table(x))
tab

## ---- fig.keep = 'high', fig.cap = "Observed frequencies of the $x$'s compared to a discrete Laplace distribution with parameters estimated from the sample."----
plot(1L:length(tab), ddisclap(as.integer(names(tab)) - y_hat, p_hat), 
  type = "h", col = "#999999", lend = "butt", lwd = 50, 
  xlab = "x, e.g. Y-STR allele", ylab = "Probability mass", axes = FALSE)
axis(1, at = 1L:length(tab), labels = names(tab))
axis(2)
points(1L:length(tab), tab, type = "h", col = "#000000", 
  lend = "butt", lwd = 25)
legend("topright", c("Estimated distribution", "Observations"), 
  pch = 15, col = c("#999999", "#000000"))

## ----message = FALSE-----------------------------------------------------
library(disclapmix)

## ----message = FALSE, warning = FALSE------------------------------------
set.seed(1)
n <- 100 # number of individuals

# Locus 1
p1 <- 0.3 # Dispersion parameter 
m1 <- 13L # Location parameter (L means integer)
d1 <- rdisclap(n, p1) + m1 # Generate a sampling using the rdisclap function

# Locus 2
p2 <- 0.4
m2 <- 14L
d2 <- rdisclap(n, p2) + m2

# Locus 3
p3 <- 0.5
m3 <- 15L
d3 <- rdisclap(n, p3) + m3

db <- cbind(d1, d2, d3)
head(db)

fit <- disclapmix(db, clusters = 1L) # L means integer type

## ------------------------------------------------------------------------
fit$y

## ------------------------------------------------------------------------
fit$disclap_parameters

## ----message = FALSE-----------------------------------------------------
library(fwsim)

## ----cache = TRUE, tidy = FALSE------------------------------------------
set.seed(1)
generations <- 100L
pop_size <- 1e5L
loci <- 7L
mutation_rates <- seq(0.001, 0.01, length.out = loci)
mutation_rates
sim <- fwsim(G = generations, H0 = rep(0L, loci), N0 = pop_size, 
  mutmodel = mutation_rates, progress = FALSE)
summary(sim)
pop <- sim$population

## ------------------------------------------------------------------------
y <- c(14L, 12L, 28L, 22L, 10L, 11L, 13L)
for (i in 1L:loci) {
  pop[, i] <- pop[, i] + y[i]
}
head(pop)

## ------------------------------------------------------------------------
pop$PopFreq <- pop$N / sum(pop$N)

## ----warning = FALSE-----------------------------------------------------
set.seed(1)
n <- 500 # Data set size
types <- sample(x = 1L:nrow(pop), size = n, replace = TRUE, prob = pop$N)
types_table <- table(types)

alpha <- sum(types_table == 1) 
alpha / n # Singleton proportion

dataset <- pop[as.integer(names(types_table)), ]
dataset$Ndb <- types_table
head(dataset)

db <- as.matrix(pop[types, 1L:loci])
head(db)

## ----cache = TRUE--------------------------------------------------------
fit <- disclapmix(db, clusters = 1L)

# Estimated location parameters
fit$y 

# Estimated dispersion parameters
fit$disclap_parameters

## ----fig.keep = 'high', fig.height = 5, fig.width = 5--------------------
plot(mutation_rates, fit$disclap_parameters, xlab = "Mutation rate", ylab = "Estimated dispersion parameter")

## ----fig.keep = 'high', fig.height = 5, fig.width = 5, tidy = FALSE------
pred_popfreqs <- predict(fit, newdata = as.matrix(dataset[, 1L:loci]))
plot(dataset$PopFreq, pred_popfreqs, log = "xy", 
  xlab = "True population frequency", 
  ylab = "Estimated population frequency")
abline(a = 0, b = 1, lty = 1)
legend("bottomright", "y = x (predicted = true)", lty = 1)

## ----cache = TRUE--------------------------------------------------------
set.seed(1)

# Common parameters
generations <- 100L
pop_size <- 1e5L
loci <- 7L

mu1 <- seq(0.001, 0.005, length.out = loci)
sim1 <- fwsim(G = generations, H0 = y, N0 = pop_size, mutmodel = mu1, progress = FALSE)
pop1 <- sim1$population

mu2 <- seq(0.005, 0.01, length.out = loci)
sim2 <- fwsim(G = generations, H0 = c(14L, 13L, 29L, 23L, 11L, 13L, 13L), N0 = pop_size, mutmodel = mu2, progress = FALSE)
pop2 <- sim2$population

## ----warning = FALSE-----------------------------------------------------
set.seed(1)
n <- 500L # Data set size

n1 <- rbinom(1, n, 0.2)
c(n1, n1 / n)

n2 <- n - n1
c(n2, n2 / n)

types1 <- sample(x = 1L:nrow(pop1), size = n1, replace = TRUE, prob = pop1$N)
db1 <- pop1[types1, 1L:loci]

types2 <- sample(x = 1L:nrow(pop2), size = n2, replace = TRUE, prob = pop2$N)
db2 <- pop2[types2, 1L:loci]

db <- as.matrix(rbind(db1, db2))

# Singleton proportion
sum(table(apply(db, 1, paste, collapse = ";")) == 1L) / n

## ----message = FALSE, warning = FALSE, cache = TRUE----------------------
fits <- lapply(1L:5L, function(clusters) disclapmix(db, clusters = clusters, iterations = 100L))

## ------------------------------------------------------------------------
BIC <- sapply(fits, function(fit) fit$BIC_marginal)

## ------------------------------------------------------------------------
best_fit <- fits[[which.min(BIC)]]
best_fit

# Estimated a priori probability of originating from each subpopulation
best_fit$tau 

# Estimated location parameters
best_fit$y 

# Estimated dispersion parameters for each subpopulation
best_fit$disclap_parameters

## ------------------------------------------------------------------------
pop1$PopFreq <- pop1$N / sum(pop1$N)
pop2$PopFreq <- pop2$N / sum(pop2$N)

types1.table <- table(types1)
types2.table <- table(types2)

dataset1 <- pop1[as.integer(names(types1.table)), ]
dataset1$Ndb <- types1.table
sum(dataset1$Ndb)

dataset2 <- pop2[as.integer(names(types2.table)), ]
dataset2$Ndb <- types2.table
sum(dataset2$Ndb)

dataset <- merge(x = dataset1, y = dataset2, by = colnames(db), all = TRUE)
dataset[is.na(dataset)] <- 0

dataset$MixPopFreq <- (n1/n) * dataset$PopFreq.x + (n2/n) * dataset$PopFreq.y

dataset$Type <- "Only from pop1"
dataset$Type[dataset$Ndb.y > 0] <- "Only from pop2"
dataset$Type[dataset$Ndb.x > 0 & dataset$Ndb.y > 0] <- "Occurred in both"
dataset$Type <- factor(dataset$Type)

## ----fig.keep = 'high', fig.height = 5, fig.width = 5, tidy = FALSE------
pred_popfreqs <- predict(best_fit, newdata = as.matrix(dataset[, 1L:loci]))
plot(dataset$MixPopFreq, pred_popfreqs, log = "xy", 
  col = c("black", "#666666", "#cccccc")[dataset$Type], 
  pch = c(3, 0, 1)[dataset$Type],
  xlab = "True population frequency", 
  ylab = "Estimated population frequency")
abline(a = 0, b = 1, lty = 1)
legend("bottomright", c("y = x (predicted = true)", levels(dataset$Type)), 
  lty = c(1, rep(-1, 3)), 
  col = c("black", "black", "#666666", "#cccccc"), 
  pch = c(-1, 3, 0, 1))

## ------------------------------------------------------------------------
best_fit

## ----tidy = FALSE--------------------------------------------------------
set.seed(1)
sim_dbs <- lapply(1L:100L, function(i) {
  simulate(best_fit, nsim = nrow(db))
})
head(sim_dbs[[1]])

num_singletons <- function(d) {
  sum(table(apply(d, 1, paste, collapse = ";")) == 1L)
}

sim_dbs_singletons <- unlist(lapply(sim_dbs, function(sim_db) { 
  num_singletons(sim_db)
}))

summary(sim_dbs_singletons)
num_singletons(db)

## ----out.width = "\\textwidth", fig.width = 12, fig.height = 5-----------
mutpars.locus1 <- c(0.149,   2.08,    18.3,   0.149,   0.374,   27.4) # DYS19
mutpars.locus2 <- c(0.500,   1.18,    18.0,   0.500,   0.0183,  349)  # DYS389I
mutpars.locus3 <- c(0.0163,  17.7,    11.1,   0.0163,  0.592,   14.1)  # DYS391
mutpars <- matrix(c(mutpars.locus1, mutpars.locus2, mutpars.locus3), ncol = 3)
colnames(mutpars) <- c("DYS19", "DYS389I", "DYS391")
mutmodel <- init_mutmodel(modeltype = 2L, mutpars = mutpars)

mutmodel_not_mut(mutmodel, locus = 1L, alleles = 5L:20L)
mutmodel_dw_mut(mutmodel, locus = 1L, alleles = 5L:20L)
mutmodel_up_mut(mutmodel, locus = 1L, alleles = 5L:20L)

statdists <- approx_stationary_dist(mutmodel, alleles = 5L:20L)
bp <- barplot(statdists, beside = TRUE, ylim = c(0, 1.1*max(statdists)))
axis(3, bp, rep(rownames(statdists), ncol(mutmodel$mutpars)), cex.axis = 1, las = 3)

## ----tidy = FALSE, cache = TRUE------------------------------------------
rare_haplotype <- c(12L, 9L, 9L) # 0/126931 at yhrd.org r46
common_haplotype <- c(14L, 12L, 10L) # 8541/126931 at yhrd.org r46
H0 <- matrix(c(rare_haplotype, common_haplotype), 2, 3, byrow = TRUE)
set.seed(1)
sim <- fwsim(G = 5000L, H0 = H0, N0 = c(1e6L, 1e6L), 
  mutmodel = mutmodel, progress = FALSE)

## ----tidy = FALSE--------------------------------------------------------
sim$population[order(sim$population$N, decreasing = TRUE)[1L:10L], ]

## ----tidy = FALSE--------------------------------------------------------
1-sapply(1L:ncol(mutmodel$mutpars), function(locus) {
  mutmodel_not_mut(mutmodel, locus = locus, alleles = rare_haplotype[locus])
})

1-sapply(1L:ncol(mutmodel$mutpars), function(locus) {
  mutmodel_not_mut(mutmodel, locus = locus, alleles = common_haplotype[locus])
})

## ------------------------------------------------------------------------
sapply(1L:ncol(mutmodel$mutpars), function(locus) mutmodel_dw_mut(mutmodel, locus = locus, alleles = rare_haplotype[locus]))
sapply(1L:ncol(mutmodel$mutpars), function(locus) mutmodel_up_mut(mutmodel, locus = locus, alleles = rare_haplotype[locus]))

## ----tidy = FALSE--------------------------------------------------------
dists <- apply(sim$population[, -ncol(sim$population)], 1, function(h) {
  sum(abs(h - rare_haplotype))
})
dists_order <- order(dists, decreasing = FALSE)[1L:10L]
data.frame(sim$population[dists_order, ], 
  DistanceToRare = dists[dists_order])

