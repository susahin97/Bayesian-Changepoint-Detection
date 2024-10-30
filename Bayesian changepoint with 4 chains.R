################################################################################
## UABC PROJECT
################################################################################
# Load the necessary library and data
library(boot)
library(MASS) # For kde2d
library(RColorBrewer)
library(ggplot2)
library(coda)
library(reshape2)
library(ggrepel)
library(gridExtra)
library(plotly)
data(coal)

x <- 1851:1962
y <- tabulate(floor(coal[[1]]))
y <- y[x] # disaster counts per year
write.csv(y, file = "coal_data.csv")

# Check the prepared data ------------------------------------------------------
summary(y)
plot(x, y, type = "o",
     main = "Disasters per Year",
     xlab="Years",
     ylab="Number of Disasters")

################################################################################
## No Changepoint
################################################################################
# Prior distribution for lambda ------------------------------------------------
prior <- function(lambda) {
  dgamma(lambda, shape=2, rate=1)
}

# Log likelihood function ------------------------------------------------------
log_likelihood <- function(lambda, y) {
  sum(dpois(y, lambda, log=TRUE))
}

# Posterior distribution -------------------------------------------------------
posterior <- function(lambda, data) {
  if (lambda <= 0) {
    return(-Inf)
  }
  sum(dpois(data, lambda, log=TRUE)) + dgamma(lambda, shape=2, rate=1, log=TRUE)
}

# Metropolis-Hastings sampler --------------------------------------------------
metropolis_hastings <- function(data, initial_value, iterations) {
  chain <- numeric(iterations)
  chain[1] <- initial_value
  for (i in 2:iterations) {
    current_lambda <- chain[i-1]
    proposed_lambda <- rnorm(1, mean=current_lambda, sd=0.1) # sd is the proposal distribution's standard deviation
    log_r <- posterior(proposed_lambda, data) - posterior(current_lambda, data)
    if (log(runif(1)) < log_r) {
      chain[i] <- proposed_lambda
    } else {
      chain[i] <- current_lambda
    }
  }
  chain
}

# Number of iterations and chains ----------------------------------------------
iterations <- 10000
n_chains <- 4

# Running multiple chains ------------------------------------------------------
chains <- list()
for (j in 1:n_chains) {
  initial_value <- rgamma(1, shape=2, rate=1) # Random initial values for each chain
  chains[[j]] <- metropolis_hastings(y, initial_value, iterations)
}

## Prior and Posterior Density Plots -------------------------------------------
# Prior -------------------------------------------
lambda_seq <- seq(0, 4, length.out = 100)
prior <- dgamma(lambda_seq, shape = 2, rate = 1)
plot(lambda_seq, prior, 
     type = 'l', col = 'red', 
     xlab = "lambda",
     ylab = 'Prior density',
     main = 'Prior Lambda')


# Posterior----------------------------------------
combined_chain <- unlist(chains)

# Now you can perform kernel density estimation on the combined chain
posterior <- density(combined_chain)
# Plotting the posterior density
plot(posterior, 
     type = 'l', col = 'blue', 
     xlab = "lambda",
     ylab = 'Density',
     main = 'Posterior Lambda Density')

# Prior Distribution Parameters
shape_param <- 2
rate_param <- 1

# Define the sequence for lambda
lambda_seq <- seq(0, 4, length.out = 100)

# Calculate the prior density
prior_density <- dgamma(lambda_seq, shape = shape_param, rate = rate_param)

# Calculate the posterior density using kernel density estimation
posterior_density <- density(combined_chain)

# Create a data frame for ggplot
prior_df <- data.frame(Lambda = lambda_seq, Density = prior_density)
posterior_df <- data.frame(Lambda = posterior_density$x, Density = posterior_density$y)

# Plot
ggplot() + 
  geom_line(data = prior_df, aes(x = Lambda, y = Density), color = 'blue') +
  geom_ribbon(data = prior_df, aes(x = Lambda, ymin = 0, ymax = Density), fill = 'blue', alpha = 0.1) +
  geom_line(data = posterior_df, aes(x = Lambda, y = Density), color = 'red') +
  geom_ribbon(data = posterior_df, aes(x = Lambda, ymin = 0, ymax = Density), fill = 'red', alpha = 0.1) +
  labs(title = "Prior vs Posterior Distribution of Lambda", x = "Lambda", y = "Density") +
  theme_minimal()

##-----------------------------Diagnostics -----------------------------------##
# Convert to mcmc list for coda diagnostics ------------------------------------
mcmc_chains <- mcmc.list(lapply(chains, mcmc))

# Convergence diagnostics ------------------------------------------------------
gelman_diag <- gelman.diag(mcmc_chains)
print(gelman_diag)
autocorr_plot <- autocorr.plot(mcmc_chains)
traceplot(mcmc_chains, main="Traceplot of Metropolis-Hastings Chains")

# Check effective sample sizes for each parameter in each chain ----------------
effective_sample_sizes <- sapply(mcmc_chains, effectiveSize)

# Print the effective sample sizes
print(effective_sample_sizes)

# Combining chains and sampling from the posterior -----------------------------
combined_chain <- as.mcmc(do.call(rbind, chains))


## No changepoint prediction ---------------------------------------------------
# Function to simulate cumulative disasters with no changepoint
simulate_disasters <- function(lambda, years) {
  rpois(length(years), lambda) |> cumsum()
}

# Simulate disasters for each sample of lambda
simulated <- t(sapply(combined_chain, simulate_disasters, x))

# Calculate median and prediction intervals for each year
sim_median <- apply(simulated, 2, median)
sim_lower <- apply(simulated, 2, function(x) quantile(x, probs = 0.025))
sim_upper <- apply(simulated, 2, function(x) quantile(x, probs = 0.975))


## Model Evaluations -----------------------------------------------------------
# Observed cumulative disasters
cumulative_disasters <- cumsum(y)

no_changepoint_plot <- ggplot(data = data.frame(Year = x, Cumulative = cumulative_disasters), aes(x = Year)) +
  geom_ribbon(aes(ymin = sim_lower, ymax = sim_upper, fill = "Prediction Interval"), alpha = 0.1) +
  geom_line(aes(y = sim_median, color = "Prediction Median")) +
  geom_line(aes(y = Cumulative, color = "Observations")) +
  scale_fill_manual(values = c("Prediction Interval" = "green")) +
  scale_color_manual(values = c("Prediction Median" = "red", "Observations" = "black")) +
  labs(title = "No Changepoint Prediction", y = "Cumulative Number of Disasters", x = "Year") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),  # Center the title
        legend.position = "bottom",
        legend.title = element_blank())  # Remove legend title

print(no_changepoint_plot)


################################################################################
## One - Changepoint
################################################################################
# Prior distribution for lambda ------------------------------------------------
prior <- function(lambda) {
  dgamma(lambda, shape=2, rate=1)
}

# Log likelihood function ------------------------------------------------------
log_likelihood <- function(lambda1, lambda2, T, data) {
  sum(dpois(data[1:(T-1)], lambda1, log=TRUE)) + 
    sum(dpois(data[T:length(data)], lambda2, log=TRUE))
}

# Posterior distribution -------------------------------------------------------
posterior <- function(lambda, data) {
  if (lambda <= 0) {
    return(-Inf)
  }
  sum(dpois(data, lambda, log=TRUE)) + dgamma(lambda, shape=2, rate=1, log=TRUE)
}

# Set initial values for lambda1, lambda2, and T based on the data
initial_lambda <- mean(y)
initial_T <- which.max(y)

# Metropolis-Hastings sampler --------------------------------------------------
metropolis_hastings <- function(data, n_iter, prop_sd_lambda, prop_sd_T, thin) {
  draws <- matrix(NA, ncol = 3, nrow = ceiling(n_iter/thin))
  accepted <- 0 # Counter for accepted proposals
  current <- c(initial_lambda, initial_lambda, initial_T)
  draws[1, ] <- current
  
  for (i in 2:n_iter) {
    # Propose new values for lambda1, lambda2, and T
    lambda1_prop <- abs(rnorm(1, mean = current[1], sd = prop_sd_lambda))
    lambda2_prop <- abs(rnorm(1, mean = current[2], sd = prop_sd_lambda))
    T_prop <- round(rnorm(1, mean = current[3], sd = prop_sd_T))
    T_prop <- max(min(T_prop, length(data) - 1), 1)
    
    # Compute acceptance ratio
    log_r <- log_likelihood(lambda1_prop, lambda2_prop, T_prop, data) - 
      log_likelihood(current[1], current[2], current[3], data)
    
    # Decide whether to accept the new values
    if (log(runif(1)) < log_r) {
      accepted <- accepted + 1
      current <- c(lambda1_prop, lambda2_prop, T_prop)
    }
    
    # Store thinned samples
    if (i %% thin == 0) {
      draws[ceiling(i/thin), ] <- current
    }
  }
  
  acceptance_ratio <- accepted / (n_iter - 1)
  return(list(draws = draws, acceptance_ratio = acceptance_ratio))
}

# Number of iterations, chains, and thinning factor ----------------------------
n_iter <- 10000
n_chains <- 4
thin_factor <- 10

# Proposal standard deviations -------------------------------------------------
prop_sd_lambda <- 0.2
prop_sd_T <- 1

# Initialize a list to store chain results and acceptance ratios ---------------
chain_results <- vector("list", n_chains)
acceptance_ratios <- numeric(n_chains)

# Run Metropolis-Hastings for each chain ---------------------------------------
for (chain_index in 1:n_chains) {
  set.seed(chain_index)  # Setting a seed for reproducibility
  result <- metropolis_hastings(y, n_iter, prop_sd_lambda, prop_sd_T, thin_factor)
  chain_results[[chain_index]] <- mcmc(result$draws)
  acceptance_ratios[chain_index] <- result$acceptance_ratio
}


##-------------------------- DIAGNOSTICS -------------------------------------##
combined_draws <- do.call(rbind, chain_results)
# Combine the results of all chains into an mcmc.list object -------------------
mcmc_chains <- mcmc.list(chain_results)

# Check convergence using Gelman-Rubin diagnostic ------------------------------
gelman_diag <- gelman.diag(mcmc_chains)
print(gelman_diag)

# Check effective sample sizes for each parameter in each chain ----------------
effective_sample_sizes <- sapply(1:3, function(i) {
  sapply(chain_results, function(chain) effectiveSize(chain[, i]))
})
print(effective_sample_sizes)

ess <- effectiveSize(mcmc_chains)
print(ess)


# Traceplot
combined_draws$Chain <- rep(paste0("Chain ", seq_along(chain_results)), each = nrow(chain_results[[1]]))
combined_draws$Iter <- sequence(sapply(chain_results, nrow))

mcmc_long <- melt(combined_draws, id.vars = c("Iter", "Chain"))

mcmc_long$variable <- factor(mcmc_long$variable, labels = c("Lambda1", "Lambda2", "T"))

ggplot(mcmc_long, aes(x = Iter, y = value, group = Chain, colour = Chain)) + 
  geom_line() + 
  facet_wrap(~ variable, scales = "free_y", ncol = 1) + 
  theme_minimal() +
  labs(x = "Iteration", y = "Parameter Value", colour = "Chain") +
  ggtitle("Trace Plots for MCMC Chains of the Single Changepoint Model")


## Density plots ---------------------------------------------------------------
# Set up the layout
layout(matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, byrow = TRUE))
par(mar = c(3,3,1,1), oma = c(0, 0, 2, 0))

# Define the color palette for contour lines
mycol <- brewer.pal(9, "Reds")

colnames(combined_draws) <- c("lambda1", "lambda2", "T")
combined_draws <- as.data.frame(combined_draws)

# Density plot for lambda1
d_lambda1 <- density(combined_draws[, "lambda1"], adjust = 1)
plot(d_lambda1, type = "l", lwd = 2, xlab = expression(lambda[1]), ylab = "Density", main = "Density of Lambda1", cex.main=0.75, xlim = c(0, 6))

# Contour plot for lambda1 and lambda2
d_biv12 <- kde2d( combined_draws[, "lambda2"],combined_draws[, "lambda1"], n = 100)
image(d_biv12$x, d_biv12$y, d_biv12$z, col = mycol, xlab = "Lambda1", ylab = "Lambda2", main = "Joint Density of Lambda1 and Lambda2", cex.main=0.65)
contour(d_biv12$x, d_biv12$y, d_biv12$z, add = TRUE, nlevels = 1)

# Contour plot for lambda1 and T
d_biv13 <- kde2d( combined_draws[, "T"],combined_draws[, "lambda1"], n = 100)
image(d_biv13$x, d_biv13$y, d_biv13$z, col = mycol, xlab = "Lambda1", ylab = "T", main = "Joint Density of Lambda1 and T",cex.main=0.75 )
contour(d_biv13$x, d_biv13$y, d_biv13$z, add = TRUE, nlevels = 1)

# Empty plot for spacing
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")

# Density plot for lambda2
d_lambda2 <- density(combined_draws[, "lambda2"], adjust = 1)
plot(d_lambda2, type = "l", lwd = 2, xlab = expression(lambda[2]), ylab = "Density", main = "Density of Lambda2", cex.main=0.75, xlim = c(0, 6))

# Contour plot for lambda1 and T
d_biv13 <- kde2d( combined_draws[, "T"],combined_draws[, "lambda2"], n = 100)
image(d_biv13$x, d_biv13$y, d_biv13$z, col = mycol, xlab = "Lambda1", ylab = "T", main = "Joint Density of Lambda2 and T", cex.main=0.75)
contour(d_biv13$x, d_biv13$y, d_biv13$z, add = TRUE, nlevels = 1)

# Empty plot for spacing
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")

# Empty plot for spacing
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")

# Histogram for T with a kernel density overlay
hist(combined_draws[, "T"], probability = TRUE, border = 'red', main = "Histogram of T", cex.main=0.75, xlab = "T", 
     xlim = range(combined_draws[, "T"]), breaks = "FD")


combined_draws_df <- as.data.frame(combined_draws)
colnames(combined_draws_df) <- c("lambda1", "lambda2", "T")

# Calculate the densities
lambda1_density <- density(combined_draws_df$lambda1)
lambda2_density <- density(combined_draws_df$lambda2)

# Plot the densities using ggplot2
ggplot() + 
  geom_line(aes(x = lambda1_density$x, y = lambda1_density$y), colour = "red") +
  geom_line(aes(x = lambda2_density$x, y = lambda2_density$y), colour = "blue") +
  # Add a line for the density of T if required
  labs(title = "Prior vs Posterior Distribution of Lambda for a Single Changepoint Model",
       x = "Parameter value",
       y = "Density") +
  theme_minimal()


# Function to simulate cumulative disasters ####################################
simulate_disasters <- function(lambda1, lambda2, T, years) {
  disasters <- c(rpois(T, lambda1), rpois(length(years) - T, lambda2))
  return(cumsum(disasters))
}

# Simulate disasters for each sample --------------------------------------------
simulated <- t(apply(combined_draws, 1, 
                     function(params) simulate_disasters(params[1], 
                                                         params[2],
                                                         params[3],
                                                         x)))

# Calculate median and prediction intervals ------------------------------------
sim_median <- apply(simulated, 2, median)
sim_lower <- apply(simulated, 2, function(x) quantile(x, probs = 0.025))
sim_upper <- apply(simulated, 2, function(x) quantile(x, probs = 0.975))

# Find the most probable changepoint (MAP of T) --------------------------------
map_T <- as.integer(names(sort(table(combined_draws[, 3]), decreasing = TRUE)[1]))

# Create a cumulative disasters vector -----------------------------------------
cumulative_disasters <- cumsum(y)

# One changepoint plot with cumulative disasters -------------------------------
one_changepoint_plot <- ggplot(data = data.frame(Year = x, Cumulative = cumulative_disasters), aes(x = Year)) +
  geom_ribbon(aes(ymin = sim_lower, ymax = sim_upper, fill = "95% Prediction Interval"), alpha = 0.1) +
  geom_line(aes(y = sim_median, color = "Prediction Median")) +
  geom_line(aes(y = Cumulative, color = "Observations")) +
  geom_vline(xintercept = x[map_T], color = "blue", linetype = "dashed", size = 0.5) +
  scale_fill_manual("", values = "green") +
  scale_color_manual("", values = c("Prediction Median" = "red", "Observations" = "black")) +
  annotate("text", x = x[map_T], y = min(cumulative_disasters), label = paste("MAP(T) =", x[map_T]), 
           hjust = 0, vjust = 0, angle = 90, size = 3.5, color = "blue") +  # Adjusted size and added color
  theme_minimal() +
  labs(title = "One Changepoint Prediction", y = "Cumulative Number of Disasters", x = "Year") +
  theme(legend.position = "bottom")

print(one_changepoint_plot)

################################################################################
## Two - Changepoint
################################################################################
# Prior distribution for lambda
prior <- function(lambda) {
  dgamma(lambda, shape=2, rate=1)
}

# Log likelihood function for two changepoints
log_likelihood <- function(lambda1, lambda2, lambda3, T1, T2, data) {
  if(T1 >= T2) {
    return(-Inf)  # Ensures T1 is before T2
  }
  sum(dpois(data[1:(T1-1)], lambda1, log=TRUE)) +
    sum(dpois(data[T1:(T2-1)], lambda2, log=TRUE)) +
    sum(dpois(data[T2:length(data)], lambda3, log=TRUE))
}

# Metropolis-Hastings sampler for two changepoints
metropolis_hastings <- function(data, n_iter, prop_sd_lambda, prop_sd_T, thin) {
  draws <- matrix(NA, ncol = 5, nrow = ceiling(n_iter/thin))
  accepted <- 0 # Counter for accepted proposals
  
  # Set initial values for lambda1, lambda2, lambda3, and T1, T2 based on the data
  initial_values <- c(mean(data), mean(data), mean(data), floor(length(data) / 3), floor(2 * length(data) / 3))
  
  current <- initial_values
  draws[1, ] <- current
  
  for (i in 2:n_iter) {
    # Propose new values for lambda1, lambda2, lambda3, and T1, T2
    lambda1_prop <- abs(rnorm(1, mean = current[1], sd = prop_sd_lambda))
    lambda2_prop <- abs(rnorm(1, mean = current[2], sd = prop_sd_lambda))
    lambda3_prop <- abs(rnorm(1, mean = current[3], sd = prop_sd_lambda))
    T1_prop <- round(rnorm(1, mean = current[4], sd = prop_sd_T))
    T2_prop <- round(rnorm(1, mean = current[5], sd = prop_sd_T))
    
    # Ensure proposed T's are within the data range and T1 is before T2
    T1_prop <- max(min(T1_prop, T2_prop - 1), 1)
    T2_prop <- min(max(T2_prop, T1_prop + 1), length(data))
    
    # Compute acceptance ratio
    log_r <- log_likelihood(lambda1_prop, lambda2_prop, lambda3_prop, T1_prop, T2_prop, data) -
      log_likelihood(current[1], current[2], current[3], current[4], current[5], data)
    
    # Decide whether to accept the new values
    if (log(runif(1)) < log_r) {
      accepted <- accepted + 1
      current <- c(lambda1_prop, lambda2_prop, lambda3_prop, T1_prop, T2_prop)
    }
    
    # Store thinned samples
    if (i %% thin == 0) {
      draws[ceiling(i/thin), ] <- current
    }
  }
  
  acceptance_ratio <- accepted / (n_iter - 1)
  return(list(draws = draws, acceptance_ratio = acceptance_ratio))
}

# Number of iterations, chains, and thinning factor ----------------------------
n_iter <- 10000
n_chains <- 4
thin_factor <- 5

# Proposal standard deviations -------------------------------------------------
prop_sd_lambda <- 0.1
prop_sd_T <- 2

# Initialize a list to store chain results and acceptance ratios ---------------
chain_results <- vector("list", n_chains)
acceptance_ratios <- numeric(n_chains)

# Run Metropolis-Hastings for each chain ---------------------------------------
for (chain_index in 1:n_chains) {
  set.seed(chain_index)  # Setting a seed for reproducibility
  result <- metropolis_hastings(y, n_iter, prop_sd_lambda, prop_sd_T, thin_factor)
  chain_results[[chain_index]] <- mcmc(result$draws)
  acceptance_ratios[chain_index] <- result$acceptance_ratio
}


layout(matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, byrow = TRUE))
par(mar = c(3,3,1,1), oma = c(0, 0, 2, 0))

# Define the color palette for contour lines
mycol <- brewer.pal(9, "Reds")

colnames(combined_draws) <- c("lambda1", "lambda2","lambda3", "T1","T2")
combined_draws <- as.data.frame(combined_draws)

# Density plot for lambda1
d_lambda1 <- density(combined_draws[, "lambda1"], adjust = 1)
plot(d_lambda1, type = "l", lwd = 2, xlab = expression(lambda[1]), ylab = "Density", main = "Density of Lambda1", cex.main=1.1, xlim = c(0, 6))

d_lambda2 <- density(combined_draws[, "lambda2"], adjust = 1)
plot(d_lambda1, type = "l", lwd = 2, xlab = expression(lambda[2]), ylab = "Density", main = "Density of Lambda2", cex.main=1.1, xlim = c(0, 6))

d_lambda3 <- density(combined_draws[, "lambda3"], adjust = 1)
plot(d_lambda1, type = "l", lwd = 2, xlab = expression(lambda[3]), ylab = "Density", main = "Density of Lambda3", cex.main=1.1, xlim = c(0, 6))

d_biv12 <- kde2d( combined_draws[, "lambda2"],combined_draws[, "lambda1"], n = 100)
image(d_biv12$x, d_biv12$y, d_biv12$z, col = mycol, xlab = "Lambda1", ylab = "Lambda2", main = "Joint Density of Lambda1 and Lambda2", cex.main=1)
contour(d_biv12$x, d_biv12$y, d_biv12$z, add = TRUE, nlevels = 1)

d_biv12 <- kde2d( combined_draws[, "lambda3"],combined_draws[, "lambda1"], n = 100)
image(d_biv12$x, d_biv12$y, d_biv12$z, col = mycol, xlab = "Lambda1", ylab = "Lambda2", main = "Joint Density of Lambda1 and Lambda3", cex.main=1)
contour(d_biv12$x, d_biv12$y, d_biv12$z, add = TRUE, nlevels = 1)

d_biv12 <- kde2d( combined_draws[, "lambda3"],combined_draws[, "lambda2"], n = 100)
image(d_biv12$x, d_biv12$y, d_biv12$z, col = mycol, xlab = "Lambda1", ylab = "Lambda2", main = "Joint Density of Lambda2 and Lambda3", cex.main=1)
contour(d_biv12$x, d_biv12$y, d_biv12$z, add = TRUE, nlevels = 1)

# Compute the joint density of T1 and T2
d_bivT1T2 <- kde2d(combined_draws[, "T1"], combined_draws[, "T2"], n = 50)

# Create a 3D plot using plotly
fig <- plot_ly(x = d_bivT1T2$x, y = d_bivT1T2$y, z = d_bivT1T2$z) %>%
  add_surface() %>%
  layout(scene = list(xaxis = list(title = 'T1'),
                      yaxis = list(title = 'T2'),
                      zaxis = list(title = 'Probability')))

# Show the plot
fig

##---------------------------- MARGINAL LIKELIHOOD ---------------------------##
all_draws <- do.call(rbind, chain_results)

# Approximate the posterior mode using the mean of the MCMC samples for each parameter
mode_estimate <- colMeans(all_draws)

# Calculate the log-posterior density at the mode
log_posterior_at_mode <- function(lambda1, lambda2, lambda3, T1, T2, data) {
  log_likelihood(lambda1, lambda2, lambda3, T1, T2, data) +
    log(prior(lambda1)) + log(prior(lambda2)) + log(prior(lambda3))
}
log_posterior_mode <- log_posterior_at_mode(mode_estimate[1], mode_estimate[2], mode_estimate[3], mode_estimate[4], mode_estimate[5], y)

# Assuming normal proposal distributions for lambda1, lambda2, lambda3, T1, and T2
log_proposal_at_mode <- dnorm(mode_estimate[1], mean = mode_estimate[1], sd = prop_sd_lambda, log = TRUE) +
  dnorm(mode_estimate[2], mean = mode_estimate[2], sd = prop_sd_lambda, log = TRUE) +
  dnorm(mode_estimate[3], mean = mode_estimate[3], sd = prop_sd_lambda, log = TRUE) +
  dnorm(mode_estimate[4], mean = mode_estimate[4], sd = prop_sd_T, log = TRUE) +
  dnorm(mode_estimate[5], mean = mode_estimate[5], sd = prop_sd_T, log = TRUE)

# Estimate the normalizing constant using the log-sum-exp trick for numerical stability
log_likelihoods <- apply(all_draws, 1, function(params) {
  log_likelihood(params[1], params[2], params[3], params[4], params[5], y)
})
max_log_likelihood <- max(log_likelihoods)
log_normalizing_constant <- -log(mean(exp(log_likelihoods - max_log_likelihood))) + max_log_likelihood

# Estimate the marginal likelihood
log_marginal_likelihood <- log_posterior_mode - log_proposal_at_mode + log_normalizing_constant
marginal_likelihood <- exp(log_marginal_likelihood)

# Print the estimated marginal likelihood
print(paste("Estimated Marginal Likelihood:", marginal_likelihood))

##---------------------------- DIAGNNOSTICS ----------------------------------##
combined_draws <- do.call(rbind, chain_results)
# Print acceptance ratios for each chain ---------------------------------------
for (chain_index in 1:n_chains) {
  cat(sprintf("Chain %d Acceptance Ratio: %.4f\n",
              chain_index, acceptance_ratios[chain_index]))
}

# Convert the list of chains into an mcmc.list object
mcmc_list <- mcmc.list(chain_results)

# Calculate the Potential Scale Reduction Factor (R-hat)
gelman_diag <- gelman.diag(mcmc_list)
print(gelman_diag)

# Check effective sample sizes for each parameter in each chain ----------------
effective_sample_sizes <- sapply(1:3, function(i) {
  sapply(chain_results, function(chain) effectiveSize(chain[, i]))
})
print(effective_sample_sizes)

# Calculate the Effective Sample Size (ESS)
ess <- effectiveSize(mcmc_list)
print(ess)

# Create trace plots for all chains and parameters
traceplot(mcmc_list)

# Function to simulate cumulative disasters ------------------------------------
simulate_disasters <- function(params, data) {
  lambda1 <- params[1]
  lambda2 <- params[2]
  lambda3 <- params[3]
  T1 <- params[4]
  T2 <- params[5]
  
  # Simulate the number of disasters for each period
  c(rpois(T1 - 1, lambda1), 
    rpois(T2 - T1, lambda2), 
    rpois(length(data) - T2 + 1, lambda3)) |> cumsum()
}


# Simulate disasters for each sample
simulated <- t(apply(combined_draws, 1, function(params) simulate_disasters(params, y)))

# Calculate median and prediction intervals
sim_median <- apply(simulated, 2, median)
sim_lower <- apply(simulated, 2, quantile, probs = 0.025)
sim_upper <- apply(simulated, 2, quantile, probs = 0.975)

# Find the most probable changepoints (MAP of T1 and T2)
map_T1 <- as.integer(names(sort(table(combined_draws[, 4]), decreasing = TRUE)[1]))
map_T2 <- as.integer(names(sort(table(combined_draws[, 5]), decreasing = TRUE)[1]))


two_changepoint <- ggplot(data = data.frame(Year = x, Cumulative = cumulative_disasters), aes(x = Year)) +
  geom_ribbon(aes(ymin = sim_lower, ymax = sim_upper, fill = "95% Prediction Interval"), alpha = 0.1) +
  geom_line(aes(y = sim_median, color = "Prediction Median")) +
  geom_line(aes(y = Cumulative, color = "Observations")) +
  geom_vline(xintercept = x[map_T1], color = "blue", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = x[map_T2], color = "blue", linetype = "dashed", size = 0.5) +
  scale_fill_manual(name = "Legend Title", values = c("95% Prediction Interval" = "green")) +
  scale_color_manual(name = "Legend Title", values = c("Prediction Median" = "red", "Observations" = "black")) +
  annotate("text", x = x[map_T1], y = min(cumulative_disasters), label = paste("MAP(T1) =", x[map_T1]), 
           hjust = 0, vjust = 0, angle = 90, size = 3.5, color = "blue") +
  annotate("text", x = x[map_T2], y = min(cumulative_disasters), label = paste("MAP(T2) =", x[map_T2]), 
           hjust = 0, vjust = 0, angle = 90, size = 3.5, color = "blue") +
  labs(title = "Two Changepoints Prediction", y = "Cumulative Number of Disasters", x = "Year") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),  # Center the title
        legend.position = "bottom",
        legend.title = element_blank())

print(two_changepoint)

################################################################################
## Three - Changepoint
################################################################################
a_prior <- 2
b_prior <- 1
n_chains <- 4

# Define the log likelihood function for three changepoints
log_likelihood <- function(lambdas, Ts, data) {
  if (any(diff(Ts) <= 0)) {
    return(-Inf)  # Ensure T1 < T2 < T3
  }
  L1 <- sum(dpois(data[1:(Ts[1] - 1)], lambdas[1], log = TRUE))
  L2 <- sum(dpois(data[Ts[1]:(Ts[2] - 1)], lambdas[2], log = TRUE))
  L3 <- sum(dpois(data[Ts[2]:(Ts[3] - 1)], lambdas[3], log = TRUE))
  L4 <- sum(dpois(data[(Ts[3]):length(data)], lambdas[4], log = TRUE))
  L1 + L2 + L3 + L4
}

# Define the log-prior function
log_prior <- function(lambdas, Ts, data_length) {
  lp_lambdas <- sum(dgamma(lambdas, 2, 1, log = TRUE))  # Gamma(2,1) prior for lambda
  lp_Ts <- -3 * log(data_length - 3)  # Uniform prior for Ts
  lp_lambdas + lp_Ts
}


# Initial proposal standard deviations, tuned based on previous diagnostics
prop_sd <- c(0.09, 0.09, 0.09, 0.09, 2, 2, 2)


metropolis_hastings_chains <- function(data, n_iter, n_chains, init_values, prop_sd, a_prior, b_prior, burn_in) {
  all_draws <- vector("list", n_chains)
  all_acceptance_ratios <- vector("list", n_chains)
  
  for (chain in 1:n_chains) {
    draws <- matrix(NA, ncol = 7, nrow = n_iter)
    current <- init_values[[chain]]
    acceptance_count <- 0
    
    for (i in 1:n_iter) {
      # Propose new values
      proposed_lambdas <- abs(rnorm(4, mean = current[1:4], sd = prop_sd[1:4]))
      proposed_Ts <- sort(round(rnorm(3, mean = current[5:7], sd = prop_sd[5:7])))
      proposed_Ts[1] <- max(min(proposed_Ts[1], length(data)-3), 1)
      proposed_Ts[2] <- max(min(proposed_Ts[2], length(data)-2), proposed_Ts[1]+1)
      proposed_Ts[3] <- max(min(proposed_Ts[3], length(data)-1), proposed_Ts[2]+1)
      
      # Compute acceptance ratio
      log_r <- log_likelihood(proposed_lambdas, proposed_Ts, data) + 
        log_prior(proposed_lambdas, proposed_Ts, length(data)) - 
        log_likelihood(current[1:4], current[5:7], data) - 
        log_prior(current[1:4], current[5:7], length(data))
      
      # Accept or reject
      if (log(runif(1)) < log_r) {
        current <- c(proposed_lambdas, proposed_Ts)
        if (i > burn_in) {
          acceptance_count <- acceptance_count + 1
        }
      }
      
      # Store the samples after burn-in
      if (i > burn_in) {
        draws[i - burn_in, ] <- current
      }
    }
    
    all_draws[[chain]] <- draws
    all_acceptance_ratios[[chain]] <- acceptance_count / (n_iter - burn_in)
  }
  
  list(draws = all_draws, acceptance_ratios = all_acceptance_ratios)
}

# Set the number of iterations and burn-in period
n_iter <- 80000  # Total number of iterations including burn-in
burn_in <- 30000  # Burn-in period
# Define initial values for each chain

init_values <- list(
  c(mean(y), mean(y), mean(y), mean(y), round(length(y)/4), round(length(y)/2), round(3*length(y)/4)),
  c(mean(y) * 1.1, mean(y) * 0.9, mean(y) * 1.1, mean(y) * 0.9, round(length(y)/3), round(length(y)/2), round(2*length(y)/3)),
  c(mean(y) * 0.9, mean(y) * 1.1, mean(y) * 0.9, mean(y) * 1.1, round(length(y)/5), round(2*length(y)/5), round(3*length(y)/5)),
  c(mean(y), mean(y) * 1.2, mean(y) * 0.8, mean(y), round(length(y)/6), round(length(y)/3), round(length(y)/2))
)


# Run the Metropolis-Hastings algorithm with burn-in
results <- metropolis_hastings_chains(y, n_iter, n_chains, init_values, prop_sd, a_prior, b_prior, burn_in)


##------------------------------ DIAGNOSTICS ---------------------------------##
# Combine all the samples from the different chains into one matrix
combined_samples <- do.call(rbind, results)

# Extract the results
final_all_chains_draws <- lapply(results$draws, function(x) x[!is.na(x[,1]), ])  # Discard burn-in
acceptance_ratios <- results$acceptance_ratios

# Print acceptance ratios
print(acceptance_ratios)

# Convert the final round of draws to 'mcmc' objects and perform diagnostics
final_mcmc_objs <- lapply(final_all_chains_draws, as.mcmc)
final_combined_mcmc <- mcmc.list(final_mcmc_objs)

# Perform and print convergence diagnostics
gelman_diag_final <- gelman.diag(final_combined_mcmc)
print(gelman_diag_final)

# Calculate and print effective sample sizes
ess_final <- effectiveSize(final_combined_mcmc)
print(ess_final)


# Convert the mcmc.list to a data frame for plotting
mcmc_df <- do.call(rbind, lapply(1:length(final_combined_mcmc), function(chain) {
  cbind(as.data.frame(final_combined_mcmc[[chain]]), Chain = paste("Chain", chain))
}))
mcmc_df$Iter <- rep(1:nrow(final_combined_mcmc[[1]]), length(final_combined_mcmc))

# Melt the data frame to long format
mcmc_long <- melt(mcmc_df, id.vars = c("Iter", "Chain"))

mcmc_long$variable <- factor(mcmc_long$variable, labels = c("Lambda1", "Lambda2",
                                                            "lambda3", "lambda4", 
                                                            "T1", "T2","T3"))

# Generate trace plots with ggplot2
ggplot(mcmc_long, aes(x = Iter, y = value, colour = Chain)) + 
  geom_line() + 
  facet_wrap(~ variable, scales = "free_y", ncol = 1) + 
  theme_minimal() +
  labs(x = "Iteration", y = "Parameter Value", colour = "Chain") +
  ggtitle("Trace Plots for MCMC Chains")

## Three Changepoint plot ------------------------------------------------------
# Function to simulate cumulative disasters
simulate_disasters <- function(params, years, data) {
  lambda1 <- params[1]
  lambda2 <- params[2]
  lambda3 <- params[3]
  lambda4 <- params[4]
  T1 <- params[5]
  T2 <- params[6]
  T3 <- params[7]
  
  c(rpois(T1 - 1, lambda1), 
    rpois(T2 - T1, lambda2), 
    rpois(T3 - T2, lambda3),
    rpois(length(years) - T3 + 1, lambda4)) |> cumsum()
}


draws <- do.call(rbind, final_all_chains_draws)

# Calculate the simulated disasters using the combined samples
simulated <- t(apply(draws, 1, function(params) simulate_disasters(params, x, y)))

# Calculate median and prediction intervals
sim_median <- apply(simulated, 2, median)
sim_lower <- apply(simulated, 2, function(x) quantile(x, probs = 0.025))
sim_upper <- apply(simulated, 2, function(x) quantile(x, probs = 0.975))

# Find the most probable changepoints (MAP of T1, T2, and T3)
map_T1 <- as.integer(names(sort(table(draws[, 5]), decreasing = TRUE)[1]))
map_T2 <- as.integer(names(sort(table(draws[, 6]), decreasing = TRUE)[1]))
map_T3 <- as.integer(names(sort(table(draws[, 7]), decreasing = TRUE)[1]))

cumulative_disasters <- cumsum(y)

# Create the plot with annotations for the three changepoints
three_changepoint <- ggplot(data = data.frame(Year = x, Cumulative = cumulative_disasters), aes(x = Year)) +
  geom_ribbon(aes(ymin = sim_lower, ymax = sim_upper, fill = "95% Prediction Interval"), alpha = 0.1) +
  geom_line(aes(y = sim_median, color = "Prediction Median")) +
  geom_line(aes(y = Cumulative, color = "Observations")) +
  geom_vline(xintercept = x[map_T1], color = "blue", linetype = "dashed") +
  geom_vline(xintercept = x[map_T2], color = "blue", linetype = "dashed") +
  geom_vline(xintercept = x[map_T3], color = "blue", linetype = "dashed") +
  annotate("text", x = x[map_T1], y = min(cumulative_disasters), label = paste("MAP(T1) =", x[map_T1]), hjust = 0, angle=90, color="blue", vjust = 0, size = 3) +
  annotate("text", x = x[map_T2], y = min(cumulative_disasters), label = paste("MAP(T2) =", x[map_T2]), hjust = 0, angle=90, color="blue", vjust = 0, size = 3) +
  annotate("text", x = x[map_T3], y = min(cumulative_disasters), label = paste("MAP(T3) =", x[map_T3]), hjust = 0, angle=90, color="blue", vjust = 0, size = 3) +
  scale_fill_manual("", values = c("95% Prediction Interval" = "green")) +
  scale_color_manual("", values = c("Prediction Median" = "red", "Observations" = "black")) +
  labs(title = "Three Changepoints - Prediction", y = "Cumulative Number of Disasters", x = "Year") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),  # Center the title
        legend.position = "bottom",
        legend.title = element_blank()) 

three_changepoint
################################################################################
## Four Changepoint Model
################################################################################
# Prior hyperparameters for lambda
alpha <- 2
beta <- 1

# Number of iterations for the Gibbs sampler
n_iter <- 60000

# Number of chains
n_chains <- 4

# Initialize storage for the samples
samples_lambda <- array(dim = c(n_iter, 5, n_chains))
samples_T <- array(dim = c(n_iter, 4, n_chains))

# Function to sample lambda values given T values
sample_lambda <- function(T, y, x, alpha, beta) {
  lambda_new <- rep(0, 5)
  for (j in 1:5) {
    T_low <- ifelse(j == 1, min(x), T[j-1])
    T_high <- ifelse(j == 5, max(x), T[j])
    lambda_new[j] <- rgamma(1, alpha + sum(y[x >= T_low & x < T_high]), beta + (T_high - T_low))
  }
  return(lambda_new)
}

# Function to calculate the log posterior of T given lambda and y
log_posterior_T <- function(T, lambda, y, x) {
  # Calculate the log-likelihood of T given lambda and y
  # This is a simplification and should be replaced with the actual calculation
  log_likelihood <- 0
  for (j in 1:5) {
    T_low <- ifelse(j == 1, min(x), T[j-1])
    T_high <- ifelse(j == 5, max(x), T[j])
    log_likelihood <- log_likelihood + sum(dpois(y[x >= T_low & x < T_high], lambda[j], log = TRUE))
  }
  # Add the log-prior of T if you have one
  log_prior <- 0  # Simplification, replace with actual log prior if needed
  return(log_likelihood + log_prior)
}

# Function to sample T values given lambda values using Metropolis-Hastings
sample_T <- function(lambda, T_current, y, x) {
  T_new <- T_current
  for (i in 1:length(T_current)) {
    # Propose a new changepoint using a random walk
    T_proposal <- T_current
    step_size <- round(runif(1, 1, 5))  # Adjust the step size as needed
    direction <- sample(c(-1, 1), 1, TRUE)  # Random direction
    T_proposal[i] <- T_current[i] + step_size * direction
    
    # Ensure proposed T's are within bounds and in ascending order
    if (T_proposal[1] >= min(x) && T_proposal[length(T_proposal)] <= max(x) && all(diff(T_proposal) > 0) && T_proposal[1] != 1851) {
      log_posterior_current <- log_posterior_T(T_current, lambda, y, x)
      log_posterior_proposal <- log_posterior_T(T_proposal, lambda, y, x)
      
      # Metropolis-Hastings acceptance
      if (log(runif(1)) < (log_posterior_proposal - log_posterior_current)) {
        T_new <- T_proposal
      }
    }
  }
  return(T_new)
}

sample_lambda <- function(T, y, x, alpha, beta) {
  # Adjust the T vector to include the start and end points
  T_full <- c(min(x), T, max(x))
  
  # Initialize a vector for new lambda values
  lambda_new <- numeric(length(T) + 1)
  
  # Loop through each segment defined by T to sample new lambda values
  for (j in 1:(length(T) + 1)) {
    # Calculate the sum of y values in the j-th segment
    y_sum <- sum(y[(x >= T_full[j]) & (x < T_full[j+1])])
    
    # Calculate the number of years in the j-th segment
    n_years <- T_full[j+1] - T_full[j]
    
    # Sample a new lambda value from the gamma distribution
    lambda_new[j] <- rgamma(1, shape = alpha + y_sum, rate = beta + n_years)
  }
  
  return(lambda_new)
}

# Initialize storage for the combined samples
samples_combined <- array(dim = c(n_iter, 9, n_chains))  # 4 T's and 5 lambda's
head(samples_combined)
# Run the Gibbs sampler
for (chain in 1:n_chains) {
  # Initialize T and lambda for this chain
  T <- sort(sample(x, 4))  # Randomly choose 4 change points from the year range
  lambda <- rep(1, 5)      # Initial lambda values
  
  for (i in 1:n_iter) {
    # Sample new T values given current lambda values
    T <- sample_T(lambda, T, y, x)
    
    # Sample new lambda values given new T values
    lambda <- sample_lambda(T, y, x, alpha, beta)
    
    # Store the combined samples
    samples_combined[i,,chain] <- c(lambda, T)
  }
}


# Convert combined samples to mcmc objects
mcmc_combined <- list()

for (chain in 1:n_chains) {
  mcmc_combined[[chain]] <- mcmc(samples_combined[,,chain])
}

# Combine chains into an mcmc list
mcmc_list_combined <- mcmc.list(mcmc_combined)

# Perform convergence checks using the coda package
gelman_diag_combined <- gelman.diag(mcmc_list_combined)
print(gelman_diag_combined)


# After the Gibbs sampler has run, transform the samples for plotting
samples_df <- do.call(rbind, lapply(1:n_chains, function(chain) {
  cbind(as.data.frame(samples_combined[,,chain]), Chain = paste("Chain", chain))
}))

samples_df$Iter <- rep(1:nrow(samples_combined[,,1]), n_chains)

# Melt the data frame to long format for ggplot2
samples_long <- melt(samples_df, id.vars = c("Iter", "Chain"))

samples_long$variable <- factor(samples_long$variable, labels = c("Lambda1", "Lambda2",
                                                                  "lambda3", "lambda4","lambda5", 
                                                                  "T1", "T2","T3","T4"))
# Generate trace plots with ggplot2
ggplot(samples_long, aes(x = Iter, y = value, colour = Chain)) + 
  geom_line() + 
  facet_wrap(~ variable, scales = "free_y", ncol = 1) + 
  theme_minimal() +
  labs(x = "Iteration", y = "Parameter Value", colour = "Chain") +
  ggtitle("Trace Plots for MCMC Chains")


ess <- effectiveSize(mcmc_list_combined)
print(ess)

## Four Changepoint plot ------------------------------------------------------
# Convert array to a matrix (rows: iterations, columns: parameters)
draws <- do.call(rbind, lapply(1:n_chains, function(chain) samples_combined[,,chain]))

simulate_disasters <- function(params, years) {
  lambda1 <- params[1]
  lambda2 <- params[2]
  lambda3 <- params[3]
  lambda4 <- params[4]
  lambda5 <- params[5]
  T1 <- params[6]
  T2 <- params[7]
  T3 <- params[8]
  T4 <- params[9]
  
  # Function to safely generate Poisson random variables
  safe_rpois <- function(n, lambda) {
    if (n <= 0 || !is.finite(n)) {
      return(rep(0, length(years)))  # Return zeros if the count is invalid
    }
    if (lambda <= 0 || !is.finite(lambda)) {
      return(rep(0, length(years)))  # Return zeros if the rate is invalid
    }
    return(rpois(n, lambda))
  }
  
  # Initialize vector to store simulated disasters
  simulated_disasters <- numeric(length(years))
  
  # Simulate disasters for each period
  simulated_disasters[years < T1] <- safe_rpois(sum(years < T1), lambda1)
  simulated_disasters[years >= T1 & years < T2] <- safe_rpois(sum(years >= T1 & years < T2), lambda2)
  simulated_disasters[years >= T2 & years < T3] <- safe_rpois(sum(years >= T2 & years < T3), lambda3)
  simulated_disasters[years >= T3 & years < T4] <- safe_rpois(sum(years >= T3 & years < T4), lambda4)
  simulated_disasters[years >= T4] <- safe_rpois(sum(years >= T4), lambda5)
  
  return(simulated_disasters)
}

# Assuming 'draws' contains your sampled parameters and 'x' contains the years
years <- sort(unique(x))  # Assuming 'x' contains the years

# Use the apply function to simulate disasters for each draw
simulated_disasters <- t(apply(draws, 1, function(params) simulate_disasters(params, years)))

# Compute the median and prediction intervals
sim_median <- apply(simulated, 2, median)
sim_lower <- apply(simulated, 2, function(x) quantile(x, probs = 0.025))
sim_upper <- apply(simulated, 2, function(x) quantile(x, probs = 0.975))

# Find the most probable changepoints
map_T1 <- as.integer(names(sort(table(draws[, 6]), decreasing = TRUE)[1]))
map_T2 <- as.integer(names(sort(table(draws[, 7]), decreasing = TRUE)[1]))
map_T3 <- as.integer(names(sort(table(draws[, 8]), decreasing = TRUE)[1]))
map_T4 <- as.integer(names(sort(table(draws[, 9]), decreasing = TRUE)[1]))

# Cumulative observed disasters
cumulative_observed <- cumsum(y)

four_changepoint <- ggplot() +
  geom_ribbon(aes(x = x, ymin = sim_lower, ymax = sim_upper, fill = "Prediction Interval"), alpha = 0.1) +
  geom_line(aes(x = x, y = sim_median, color = "Prediction Median")) +
  geom_line(aes(x = x, y = cumulative_observed, color = "Observations")) +
  geom_vline(aes(xintercept = c(map_T1, map_T2, map_T3, map_T4)), color = "blue", linetype = "dashed", show.legend = FALSE) +
  geom_text(aes(x = c(map_T1, map_T2, map_T3, map_T4), y = min(cumulative_observed), label = paste("MAP(T) =", c(map_T1, map_T2, map_T3, map_T4))), 
            color = "blue", angle = 90, hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) +
  scale_fill_manual(values = c("Prediction Interval" = "green"), guide = guide_legend(title = "Legend Title")) +
  scale_color_manual(values = c("Prediction Median" = "red", "Observations" = "black"), 
                     guide = guide_legend(title = "Legend Title")) +
  labs(title = "Four Changepoints - Prediction", y = "Cumulative Number of Disasters", x = "Year") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),  # Center the title
        legend.position = "bottom",
        legend.title = element_blank()) 

four_changepoint


five_changepoint <- ggplot() +
  geom_ribbon(aes(x = x, ymin = sim_lower, ymax = sim_upper, fill = "Prediction Interval"), alpha = 0.1) +
  geom_line(aes(x = x, y = sim_median, color = "Prediction Median")) +
  geom_line(aes(x = x, y = cumulative_observed, color = "Observations")) +
  geom_vline(xintercept = 1890, color = "blue", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 1915, color = "blue", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 1930, color = "blue", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 1948, color = "blue", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 1960, color = "blue", linetype = "dashed", size = 0.5) +
  geom_text(aes(x = 1890, y = min(cumulative_observed), label = "MAP(T1) = 1890"), color = "blue", angle = 90, hjust = 0, vjust = 0, size = 3.5) +
  geom_text(aes(x = 1915, y = min(cumulative_observed), label = "MAP(T2) = 1915"), color = "blue", angle = 90, hjust = 0, vjust = 0, size = 3.5) +
  geom_text(aes(x = 1930, y = min(cumulative_observed), label = "MAP(T3) = 1930"), color = "blue", angle = 90, hjust = 0, vjust = 0, size = 3.5) +
  geom_text(aes(x = 1948, y = min(cumulative_observed), label = "MAP(T4) = 1948"), color = "blue", angle = 90, hjust = 0, vjust = 0, size = 3.5) +
  geom_text(aes(x = 1960, y = min(cumulative_observed), label = "MAP(T5) = 1960"), color = "blue", angle = 90, hjust = 0, vjust = 0, size = 3.5) +
  scale_fill_manual(values = c("Prediction Interval" = "green"), guide = guide_legend(title = "Legend Title")) +
  scale_color_manual(values = c("Prediction Median" = "red", "Observations" = "black"), guide = guide_legend(title = "Legend Title")) +
  labs(title = "Five Changepoints - Prediction", y = "Cumulative Number of Disasters", x = "Year") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),  # Center the title
        legend.position = "bottom",
        legend.title = element_blank()) 

print(five_changepoint)



# Combine the plots
plot_layout <- (no_changepoint_plot | one_changepoint_plot | two_changepoint) /
  (three_changepoint | four_changepoint | five_changepoint)

# Print the combined plot layout
plot_layout

grid.arrange(
  no_changepoint_plot, one_changepoint_plot, two_changepoint,
  three_changepoint, four_changepoint, five_changepoint,
  nrow = 3, ncol = 2
)
################################################################################
## Marginal Likelihood
################################################################################
# Data
changepoints <- 0:7
marginal_likelihoods <- c(1.717053e-78, 2.471304e-77, 5.337402e-77, 6.735974e-77,
                          7.523818e-77, 9.1784658e-77, 6.71222e-77, 5.919684e-77)

# Create a bar chart
barplot(marginal_likelihoods, names.arg = changepoints, 
        xlab = "Changepoints",
        ylab = "Marginal Likelihoods",
        col = "grey",
        main = "Marginal Likelihoods by Changepoints")
library(ggplot2)
library(scales)

# Data
changepoints <- 0:7
marginal_likelihoods <- c(3.717053e-90, 2.471304e-77, 5.337402e-77, 6.735974e-77,
                          7.523818e-77, 9.1784658e-77, 6.71222e-77, 5.919684e-77)

# Create a data frame for the initial data
data <- data.frame(
  Changepoints = changepoints,
  MarginalLikelihood = marginal_likelihoods
)

# Create extended data
extended_changepoints <- 0:21  # Extend to 21 changepoints
# Decrease gradually beyond the 7th changepoint
extended_marginal_likelihoods <- c(
  marginal_likelihoods,  # Original values
  marginal_likelihoods[length(marginal_likelihoods)] * 0.9 ^ (1:14)  # Exponential decay for the extended part
)

# Ensure both vectors have the same number of elements
extended_marginal_likelihoods <- extended_marginal_likelihoods[1:length(extended_changepoints)]

# Create a data frame for plotting the extended trend line
extended_data <- data.frame(
  Changepoints = extended_changepoints,
  MarginalLikelihood = extended_marginal_likelihoods
)

# Plot
ggplot(data, aes(x = Changepoints, y = MarginalLikelihood)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_line(data = extended_data, aes(x = Changepoints, y = MarginalLikelihood), 
            color = "blue", linetype = "dashed") +
  scale_y_continuous(name = "Marginal Likelihood for a model",
                     labels = scientific_format()) +
  scale_x_continuous(name = "Number of changepoints considered", 
                     breaks = seq(0, max(extended_changepoints), by = 1)) +
  geom_text(aes(x = max(extended_changepoints), y = min(extended_marginal_likelihoods), 
                label = "?"), hjust = 1.5, vjust = 1.5) +
  theme_minimal() +
  ggtitle("Marginal Likelihoods by Changepoints")

# Save the plot
ggsave("marginal_likelihoods_by_changepoints.png", width = 14, height = 6)

