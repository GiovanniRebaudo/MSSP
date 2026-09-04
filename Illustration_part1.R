# ============================================================
# Figure 2: different random probabilities, same pEPPF
# ============================================================
#
# Models:
# \star: mSSP with iid non-atomic common atoms
# 1: repulsive common atoms
# 2: one fixed spike at zero and M-1 iid slab atoms
# 3: iid atoms + stochastic dependent matching of weights/atoms
#
# Required packages:
# install.packages(c("ggplot2", "patchwork"))

library(ggplot2)
library(patchwork)

set.seed(1234)

# ------------------------------------------------------------
# 1. Hyperparameters
# ------------------------------------------------------------

M       = 5
gamma   = 0.60

# P0 = N(mu0, sd0^2).
mu0     = 0
sd0     = 0.90

# Repulsion parameter in Model 1
delta_A = 0.35

# Matching strength in Model 3 [g(theta) = theta]
kappa   = 50.0

# Gaussian mixture kernel
sigma_kernel = 0.30

# Sample sizes used for the partition entropy
n1 = 50
n2 = 50

# Monte Carlo size 
B = 10000

# Grid used for the random densities
y_grid = seq(-3.5, 3.5, length.out = 300)

# Base measure P0 and matching score g
rP0   = function(n) rnorm(n, mean = mu0, sd = sd0)
g_fun = function(x) x


# ------------------------------------------------------------
# 2. Basic utilities
# ------------------------------------------------------------

rdirichlet1 <- function(alpha) {
  z = rgamma(length(alpha), shape = alpha, rate = 1)
  z / sum(z)
}

# Compute common weights using the representation of model 3
# The complete vectors V_j = (omega_j, (1-omega_j) pi_j) are iid symmetric
# Dirichlet_M vectors.
draw_common_weights <- function(M, gamma) {
  stopifnot(M >= 2L, gamma > 0)

  omega1 = rbeta(1L, shape1 = gamma, shape2 = (M - 1L) * gamma)
  omega2 = rbeta(1L, shape1 = gamma, shape2 = (M - 1L) * gamma)

  pi1 = rdirichlet1(rep(gamma, M - 1L))
  pi2 = rdirichlet1(rep(gamma, M - 1L))

  V1 = c(omega1, (1 - omega1) * pi1)
  V2 = c(omega2, (1 - omega2) * pi2)

  list(
    omega1 = omega1,
    omega2 = omega2,
    pi1 = pi1,
    pi2 = pi2,
    V1 = V1,
    V2 = V2
  )
}

# Compute entropy: 
# sample allocation given the weights and compute entropy
partition_entropy <- function(V1, V2, n1, n2) {
  K  = length(V1)

  z1 = sample.int(K, size = n1, replace = TRUE, prob = V1)
  z2 = sample.int(K, size = n2, replace = TRUE, prob = V2)

  nk = tabulate(c(z1, z2), nbins = K)
  nk = nk[nk > 0]

  p = nk / (n1 + n2)
  -sum(p * log(p))
}

# Evaluate sum_h w_h phi_sigma(y - theta_h) on a grid.
mixture_density <- function(y, w, theta, sigma) {
  K = length(theta)

  kernels = vapply(
    theta,
    function(th) dnorm(y, mean = th, sd = sigma),
    numeric(length(y))
  )

  # vapply returns a vector when K = 1; convert it to a matrix.
  if (K == 1L) {
    kernels = matrix(kernels, ncol = 1L)
  }

  as.vector(kernels %*% w)
}

# Mean of the random density under Gaussian kernels.
mixture_mean <- function(w, theta) {
  sum(w * theta)
}


# ------------------------------------------------------------
# 3. Model 1: exact rejection sampler for repulsive atoms
# ------------------------------------------------------------
#
# Target:
#   prod_h p0(theta_h)
#   prod_{h<l} [1 - exp{-(theta_h-theta_l)^2/(2 delta_A^2)}]
#
# The proposal is iid P0. 

log_repulsion_factor <- function(theta, delta) {
  K = length(theta)
  if (K <= 1L) return(0)

  dif = outer(theta, theta, "-")
  d2  = dif[upper.tri(dif)]^2

  z = -d2 / (2 * delta^2)

  # log{1 - exp(z)}, with z <= 0.
  terms = log1p(-exp(z))

  if (any(!is.finite(terms))) return(-Inf)
  sum(terms)
}

draw_repulsive_atoms <- function(K, delta, rP0, max_tries = 1000000L) {
  if (K == 1L) return(rP0(1L))

  for (iter in seq_len(max_tries)) {
    theta   = rP0(K)
    log_acc = log_repulsion_factor(theta, delta)

    if (log(runif(1)) < log_acc) {
      return(theta)
    }
  }

  stop(
    "Repulsive rejection sampler did not accept. ",
    "Try decreasing delta_A or M."
  )
}


# ------------------------------------------------------------
# 4. Model 2: stochastic matching
# ------------------------------------------------------------

all_permutations <- function(x) {
  if (length(x) == 1L) {
    return(matrix(x, nrow = 1L))
  }

  out = lapply(seq_along(x), function(i) {
    rest = all_permutations(x[-i])
    cbind(x[i], rest)
  })

  do.call(rbind, out)
}

# Cache all permutations because M is small.
perm_cache <- setNames(
  lapply(seq_len(M), function(k) all_permutations(seq_len(k))),
  as.character(seq_len(M))
)

sample_matching_permutation <- function(theta, S, kappa, g_fun, perm_cache) {
  K = length(theta)
  if (K == 1L) return(1L)

  P = perm_cache[[as.character(K)]]

  # Row q contains (theta_{sigma_q(1)}, ..., theta_{sigma_q(K)}).
  idx_rowwise = as.vector(t(P))
  theta_perm  = matrix(
    theta[idx_rowwise],
    nrow = nrow(P),
    ncol = K,
    byrow = TRUE
  )

  gtheta = matrix(
    g_fun(as.vector(theta_perm)),
    nrow = nrow(theta_perm),
    ncol = K
  )

  scores   = as.vector(gtheta %*% S)
  log_prob = kappa * scores
  prob     = exp(log_prob - max(log_prob))
  prob     = prob / sum(prob)

  q = sample.int(nrow(P), size = 1L, prob = prob)
  as.integer(P[q, ])
}


# ------------------------------------------------------------
# 5. Monte Carlo simulation
# ------------------------------------------------------------

model_names = c("Star", "1", "2", "3")
group_names = c("Group 1", "Group 2")

entropy_draws = numeric(B)
omega_draws = matrix(
  NA_real_,
  nrow = B,
  ncol = 2L,
  dimnames = list(NULL, c("omega1", "omega2"))
)

# store densities
# densities[[model]][[group]] is a B x length(y_grid) matrix
densities = setNames(
  lapply(model_names, function(m) {
    list(
      `Group 1` = matrix(NA_real_, nrow = B, ncol = length(y_grid)),
      `Group 2` = matrix(NA_real_, nrow = B, ncol = length(y_grid))
    )
  }),
  model_names
)

# store (mu_1, mu_2) per model
mu_draws = setNames(
  lapply(model_names, function(m) {
    data.frame(mu1 = numeric(B), mu2 = numeric(B))
  }),
  model_names
)

for (b in seq_len(B)) {

  # ----------------------------------------------------------
  # Common complete weights and common joint partition
  # ----------------------------------------------------------
  W = draw_common_weights(M, gamma)

  K  = M
  V1 = W$V1
  V2 = W$V2

  omega_draws[b, ] = c(W$omega1, W$omega2)

  # Since all four models have the same pEPPF, we use the same
  # partition realization for all four models.
  entropy_draws[b] = partition_entropy(V1, V2, n1, n2)


  # ----------------------------------------------------------
  # Star: mSSP
  # ----------------------------------------------------------
  theta_star = rP0(K)

  densities[["Star"]][["Group 1"]][b, ] =
    mixture_density(y_grid, V1, theta_star, sigma_kernel)

  densities[["Star"]][["Group 2"]][b, ] =
    mixture_density(y_grid, V2, theta_star, sigma_kernel)

  mu_draws[["Star"]]$mu1[b] = mixture_mean(V1, theta_star)
  mu_draws[["Star"]]$mu2[b] = mixture_mean(V2, theta_star)


  # ----------------------------------------------------------
  # 1: repulsive common atoms
  # ----------------------------------------------------------
  theta_1 = draw_repulsive_atoms(K, delta_A, rP0)

  densities[["1"]][["Group 1"]][b, ] =
    mixture_density(y_grid, V1, theta_1, sigma_kernel)

  densities[["1"]][["Group 2"]][b, ] =
    mixture_density(y_grid, V2, theta_1, sigma_kernel)

  mu_draws[["1"]]$mu1[b] = mixture_mean(V1, theta_1)
  mu_draws[["1"]]$mu2[b] = mixture_mean(V2, theta_1)


  # ----------------------------------------------------------
  # 2: spike-and-slab
  # ----------------------------------------------------------
  #
  # The first support point is the common fixed spike at zero.

  theta_2 = c(0, rP0(M - 1L))

  densities[["2"]][["Group 1"]][b, ] =
    mixture_density(y_grid, V1, theta_2, sigma_kernel)

  densities[["2"]][["Group 2"]][b, ] =
    mixture_density(y_grid, V2, theta_2, sigma_kernel)

  mu_draws[["2"]]$mu1[b] = mixture_mean(V1, theta_2)
  mu_draws[["2"]]$mu2[b] = mixture_mean(V2, theta_2)


  # ----------------------------------------------------------
  # 3: dependent matching of weights and iid atoms
  # ----------------------------------------------------------
  theta_3 = rP0(K)

  S = V1 + V2

  Sigma = sample_matching_permutation(
    theta = theta_3,
    S = S,
    kappa = kappa,
    g_fun = g_fun,
    perm_cache = perm_cache
  )

  # W_{j, Sigma(h)} = V_{j,h}
  W1_3 = numeric(K)
  W2_3 = numeric(K)

  W1_3[Sigma] = V1
  W2_3[Sigma] = V2

  densities[["3"]][["Group 1"]][b, ] =
    mixture_density(y_grid, W1_3, theta_3, sigma_kernel)

  densities[["3"]][["Group 2"]][b, ] =
    mixture_density(y_grid, W2_3, theta_3, sigma_kernel)

  mu_draws[["3"]]$mu1[b] = mixture_mean(W1_3, theta_3)
  mu_draws[["3"]]$mu2[b] = mixture_mean(W2_3, theta_3)


  if (b %% 250L == 0L) {
    message("Monte Carlo iteration ", b, " / ", B)
  }
}


# ------------------------------------------------------------
# 6. Summaries for the density panels
# ------------------------------------------------------------

summarise_density_matrix <- function(mat, model, group, y_grid) {
  data.frame(
    model = model,
    group = group,
    y     = y_grid,
    mean  = colMeans(mat),
    lower = apply(mat, 2L, quantile, probs = 0.025),
    upper = apply(mat, 2L, quantile, probs = 0.975)
  )
}

density_summary <- do.call(
  rbind,
  lapply(model_names, function(m) {
    do.call(
      rbind,
      lapply(group_names, function(g) {
        summarise_density_matrix(
          densities[[m]][[g]],
          model = m,
          group = g,
          y_grid = y_grid
        )
      })
    )
  })
)

entropy_df <- do.call(
  rbind,
  lapply(model_names, function(m) {
    data.frame(
      model = m,
      entropy = entropy_draws
    )
  })
)

mu_df <- do.call(
  rbind,
  lapply(model_names, function(m) {
    data.frame(
      model = m,
      mu1 = mu_draws[[m]]$mu1,
      mu2 = mu_draws[[m]]$mu2
    )
  })
)

model_titles <- c(
  "Star" = "mSSP",
  "1" = "repulsive atoms",
  "2" = "spike-and-slab",
  "3" = "dependent matching"
)


# ------------------------------------------------------------
# 7. Build the 3 x 4 figure
# ------------------------------------------------------------
theme_illustration <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.title = element_blank()
  )
# Top row: same entropy distribution in every panel.
top_plots <- lapply(model_names, function(m) {
  dat <- entropy_df[entropy_df$model == m, ]
  ggplot(dat, aes(x = entropy)) +
    geom_density(
      fill = "grey80",
      colour = "black",
      linewidth = 0.55,
      adjust = 1
    ) +
    labs(
      title = model_titles[[m]],
      x = "Entropy of pooled cluster frequencies",
      y = "Prior density"
    ) +
    theme_illustration +
    theme(legend.position = "none")
})
# Middle row: random densities
bottom_ymax  = 1.03 * max(density_summary$upper, na.rm = TRUE)
middle_plots = lapply(model_names, function(m) {
  dat <- density_summary[density_summary$model == m, ]
  p <- ggplot(
    dat,
    aes(
      x = y,
      y = mean,
      colour = group,
      fill = group,
      linetype = group
    )
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.18,
      colour = NA
    ) +
    geom_line(linewidth = 0.75) +
    coord_cartesian(
      xlim = range(y_grid),
      ylim = c(0, bottom_ymax)
    ) +
    labs(
      x = expression(y),
      y = "Mixture density"
    ) +
    theme_illustration
  if (m == "Star") {
    p <- p +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.05, 0.95),
        legend.justification = c(0, 1),
        legend.direction = "horizontal",
        legend.background = element_rect(
          fill = "white",
          colour = NA
        )
      ) +
      guides(
        colour = guide_legend(nrow = 1)
      )
  } else {
    p = p +
      theme(legend.position = "none")
  }
  p
})
# Third row: joint distribution of (mu1, mu2)
mu_min  = min(c(mu_df$mu1, mu_df$mu2), na.rm = TRUE)
mu_max  = max(c(mu_df$mu1, mu_df$mu2), na.rm = TRUE)
mu_pad  = 0.05 * (mu_max - mu_min)
mu_lims = c(mu_min - mu_pad, mu_max + mu_pad)
third_plots = lapply(model_names, function(m) {
  dat = mu_df[mu_df$model == m, ]
  ggplot(dat, aes(x = mu1, y = mu2)) +
    geom_point(alpha = 0.6, size = 0.8, colour = "grey80") +
    stat_density_2d(linewidth = 0.45, colour = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    coord_cartesian(xlim = mu_lims, ylim = mu_lims) +
    labs(
      x = expression(mu[1]),
      y = expression(mu[2])
    ) +
    theme_illustration
})
figure_A =
  (
    wrap_plots(top_plots, nrow = 1) /
      wrap_plots(middle_plots, nrow = 1) /
      wrap_plots(third_plots, nrow = 1)
  ) +
  plot_layout(
    heights = c(1, 1.25, 1.15)
  )
print(figure_A)
# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------
ggsave(
  filename = "illustration_A.pdf",
  plot = figure_A,
  width = 12,
  height = 8.0,
  units = "in"
)
ggsave(
  filename = "illustration_A.png",
  plot = figure_A,
  width = 12,
  height = 8.0,
  units = "in",
  dpi = 300
)



# ------------------------------------------------------------
# 9. Optional numerical summaries
# ------------------------------------------------------------

cat("\nMean partition entropy:", mean(entropy_draws), "\n")
cat(
  "95% interval for partition entropy:",
  quantile(entropy_draws, c(0.025, 0.975)),
  "\n"
)

cat("\nEmpirical prior correlations of (mu1, mu2):\n")
for (m in model_names) {
  cat(m, ": ", cor(mu_draws[[m]]$mu1, mu_draws[[m]]$mu2), "\n", sep = "")
}

cat("\nModel parameters:\n")
cat("M =", M, "\n")
cat("gamma =", gamma, "\n")
cat(
  "omega_j iid ~ Beta(",
  gamma,
  ", ",
  (M - 1L) * gamma,
  ")\n",
  sep = ""
)
cat("theoretical E[omega_j] =", 1 / M, "\n")
cat("empirical means of omega_1 and omega_2 =", colMeans(omega_draws), "\n")
cat("delta_A =", delta_A, "\n")
cat("kappa =", kappa, "\n")
cat("P0 = Normal(", mu0, ", ", sd0, "^2)\n", sep = "")
cat("kernel sigma =", sigma_kernel, "\n")
