# ============================================================
# Figure 3: different random probabilities, same latent partition 
# ============================================================
#
# Three-row figure:
#   Row 1: entropy of the final partition Pi induced by ties in
#          the atom locations
#   Row 2: prior mean and 95% pointwise bands of the mixture densities
#          (for Model 6, all categorical covariate configurations are shown)
#   Row 3: joint prior distribution of (mu_1, mu_2)
#
# Models:
# Star: common-atoms finite mSSP
# 4: paired group-specific atoms (full-range borrowing)
# 5: spike-and-slab common atoms (coalescence at zero)
# 6: ANOVA-type covariate-dependent atoms
#
# Required packages:
# install.packages(c("ggplot2", "patchwork"))

library(ggplot2)
library(patchwork)

set.seed(1234)

# ------------------------------------------------------------
# 1. Hyperparameters
# ------------------------------------------------------------

# Common latent mSSP
L     = 5
alpha = 0.60

# Prior sample sizes used to construct the partitions
n1 = 50
n2 = 50
N  = n1 + n2

# Monte Carlo size
B_MC = 10000

# Gaussian mixture kernel
sigma_kernel = 0.30

# Grid for density summaries
y_grid = seq(-3.5, 3.5, length.out = 300)

# ------------------------------------------------------------
# Common marginal atom scale
# ------------------------------------------------------------

mu0 = 0
sd0 = 0.90

rP0 <- function(n) {
  rnorm(n, mean = mu0, sd = sd0)
}

# ------------------------------------------------------------
# 4: paired group-specific atoms
# ------------------------------------------------------------
# (theta_1h, theta_2h) ~ N_2(theta0, Sigma0)

rho_4 = -0.85

theta0_4 = c(mu0, mu0)

Sigma0_4 = sd0^2 * matrix(
  c(
    1, rho_4,
    rho_4, 1
  ),
  nrow = 2,
  byrow = TRUE
)

# ------------------------------------------------------------
# 5: spike-and-slab atoms
# ------------------------------------------------------------

omega_5 = 0.30

# ------------------------------------------------------------
# 6: ANOVA-type covariate-dependent atoms
# ------------------------------------------------------------
#
# x = (v,w), v = 1,...,V, w = 1,...,W
#
# theta_{x,h} = m_h + A_{v,h} + B_{w,h}
#
# All V x W covariate-indexed densities are displayed in the
# Model 6 panel in the second row.

V = 2
W = 2

# Random-effect standard deviations.
sd_m_6 = 0.60
sd_A_6 = 0.45
sd_B_6 = 0.50

# Distribution of observation-specific categorical covariates,
# used both for the final partition and for the covariate-averaged
# density mean used in Row 3.
prob_v = rep(1 / V, V)
prob_w = rep(1 / W, W)

# All Model 6 covariate configurations
configs_6 = expand.grid(
  v = seq_len(V),
  w = seq_len(W)
)

config_labels_6 = paste0(
  "(v=", configs_6$v,
  ", w=", configs_6$w, ")"
)

n_configs_6 = nrow(configs_6)

# Probability of each configuration under the specified covariate law.
config_prob_6 = prob_v[configs_6$v] * prob_w[configs_6$w]
config_prob_6 = config_prob_6 / sum(config_prob_6)


# ------------------------------------------------------------
# 2. Basic utilities
# ------------------------------------------------------------

rdirichlet1 <- function(alpha_vec) {
  z = rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)
  z / sum(z)
}

draw_latent_weights <- function(L, alpha) {
  list(
    pi1 = rdirichlet1(rep(alpha, L)),
    pi2 = rdirichlet1(rep(alpha, L))
  )
}

draw_allocations <- function(pi1, pi2, n1, n2) {
  list(
    z1 = sample.int(
      length(pi1),
      size = n1,
      replace = TRUE,
      prob = pi1
    ),
    z2 = sample.int(
      length(pi2),
      size = n2,
      replace = TRUE,
      prob = pi2
    )
  )
}

# Entropy of a partition represented by arbitrary cluster keys
partition_entropy_from_keys <- function(keys) {
  counts = table(keys)
  p = as.numeric(counts) / length(keys)
  -sum(p * log(p))
}

# Latent partition Pi^0:
# equality of latent component labels across the two groups.
latent_partition_entropy <- function(z1, z2) {
  keys = c(
    paste0("h", z1),
    paste0("h", z2)
  )
  partition_entropy_from_keys(keys)
}

# Evaluate sum_h w_h phi_sigma(y - theta_h)
mixture_density <- function(y, weights, theta, sigma) {

  kernels = vapply(
    theta,
    function(th) {
      dnorm(y, mean = th, sd = sigma)
    },
    numeric(length(y))
  )

  if (length(theta) == 1L) {
    kernels = matrix(kernels, ncol = 1L)
  }

  as.vector(kernels %*% weights)
}

# Mean of a Gaussian location mixture:
# mu = integral y f(y) dy = sum_h w_h theta_h.
mixture_mean <- function(weights, theta) {
  sum(weights * theta)
}

# Multivariate normal sampler
rmvnorm_base <- function(n, mean, Sigma) {
  p = length(mean)

  Z = matrix(
    rnorm(n * p),
    nrow = n,
    ncol = p
  )

  Rchol = chol(Sigma)

  sweep(
    Z %*% Rchol,
    2,
    mean,
    FUN = "+"
  )
}


# ------------------------------------------------------------
# 3. Final partition under each atom-assignment mechanism
# ------------------------------------------------------------

# Star:
# common continuous atoms imply Pi = Pi^0 almost surely
final_entropy_Star <- function(z1, z2) {
  latent_partition_entropy(z1, z2)
}

# 4:
# group-specific continuous atoms imply no cross-group ties.
# Within each group, observations with the same h share the same atom.
final_entropy_4 <- function(z1, z2) {

  keys = c(
    paste0("g1_h", z1),
    paste0("g2_h", z2)
  )

  partition_entropy_from_keys(keys)
}

# 5:
# all latent components whose atom is zero coalesce.
final_entropy_5 <- function(z1, z2, theta) {

  atom_key = ifelse(
    theta == 0,
    "zero",
    paste0("nonzero_h", seq_along(theta))
  )

  keys = c(
    atom_key[z1],
    atom_key[z2]
  )

  partition_entropy_from_keys(keys)
}

# 6:
# theta_{x,h} is shared by observations having the same
# latent component h and the same covariate configuration x=(v,w),
# including across groups.
#
# With continuous random effects, different (x,h) pairs produce
# distinct locations almost surely.
final_entropy_6 <- function(
  z1, z2,
  v1, w1,
  v2, w2
) {

  keys = c(
    paste0(
      "v", v1,
      "_w", w1,
      "_h", z1
    ),
    paste0(
      "v", v2,
      "_w", w2,
      "_h", z2
    )
  )

  partition_entropy_from_keys(keys)
}


# ------------------------------------------------------------
# 4. 6 atom trajectories
# ------------------------------------------------------------

draw_6_effects <- function(
  L, V, W,
  sd_m, sd_A, sd_B
) {

  list(
    m = rnorm(
      L,
      mean = 0,
      sd = sd_m
    ),

    A = matrix(
      rnorm(
        V * L,
        mean = 0,
        sd = sd_A
      ),
      nrow = V,
      ncol = L
    ),

    Beff = matrix(
      rnorm(
        W * L,
        mean = 0,
        sd = sd_B
      ),
      nrow = W,
      ncol = L
    )
  )
}

theta_6 <- function(effects, v, w) {
  effects$m +
    effects$A[v, ] +
    effects$Beff[w, ]
}


# ------------------------------------------------------------
# 5. Monte Carlo simulation
# ------------------------------------------------------------

model_names = c(
  "Star", "4", "5", "6"
)

group_names = c(
  "Group 1", "Group 2"
)

latent_entropy = numeric(B_MC)

# Final tie-induced partition entropy
final_entropy = setNames(
  lapply(
    model_names,
    function(x) numeric(B_MC)
  ),
  model_names
)

# Density draws for Star--5.
# 6 is stored separately for all covariate configurations.
densities = setNames(
  lapply(
    c("Star", "4", "5"),
    function(m) {
      list(
        `Group 1` = matrix(
          NA_real_,
          nrow = B_MC,
          ncol = length(y_grid)
        ),
        `Group 2` = matrix(
          NA_real_,
          nrow = B_MC,
          ncol = length(y_grid)
        )
      )
    }
  ),
  c("Star", "4", "5")
)

# 6:
# array dimensions = Monte Carlo iteration x y-grid x configuration
densities_6_all = list(
  `Group 1` = array(
    NA_real_,
    dim = c(
      B_MC,
      length(y_grid),
      n_configs_6
    )
  ),

  `Group 2` = array(
    NA_real_,
    dim = c(
      B_MC,
      length(y_grid),
      n_configs_6
    )
  )
)

# Row 3:
# joint prior draws of (mu_1, mu_2)
mu_draws = setNames(
  lapply(
    model_names,
    function(m) {
      data.frame(
        mu1 = numeric(B_MC),
        mu2 = numeric(B_MC)
      )
    }
  ),
  model_names
)


for (b in seq_len(B_MC)) {

  # ----------------------------------------------------------
  # Common latent mSSP stage
  # ----------------------------------------------------------

  WW = draw_latent_weights(
    L = L,
    alpha = alpha
  )

  pi1 = WW$pi1
  pi2 = WW$pi2

  ZZ = draw_allocations(
    pi1 = pi1,
    pi2 = pi2,
    n1 = n1,
    n2 = n2
  )

  z1 = ZZ$z1
  z2 = ZZ$z2

  # Same latent partition Pi^0 for Star--6
  latent_entropy[b] =
    latent_partition_entropy(
      z1 = z1,
      z2 = z2
    )


  # ----------------------------------------------------------
  # Star: common non-atomic atoms
  # ----------------------------------------------------------

  theta_Star <- rP0(L)

  densities[["Star"]][["Group 1"]][b, ] =
    mixture_density(
      y = y_grid,
      weights = pi1,
      theta = theta_Star,
      sigma = sigma_kernel
    )

  densities[["Star"]][["Group 2"]][b, ] =
    mixture_density(
      y = y_grid,
      weights = pi2,
      theta = theta_Star,
      sigma = sigma_kernel
    )

  final_entropy[["Star"]][b] =
    final_entropy_Star(
      z1 = z1,
      z2 = z2
    )

  mu_draws[["Star"]]$mu1[b] =
    mixture_mean(
      weights = pi1,
      theta = theta_Star
    )

  mu_draws[["Star"]]$mu2[b] =
    mixture_mean(
      weights = pi2,
      theta = theta_Star
    )


  # ----------------------------------------------------------
  # 4: paired group-specific atoms
  # ----------------------------------------------------------

  theta_pair = rmvnorm_base(
    n = L,
    mean = theta0_4,
    Sigma = Sigma0_4
  )

  theta1_4 = theta_pair[, 1]
  theta2_4 = theta_pair[, 2]

  densities[["4"]][["Group 1"]][b, ] =
    mixture_density(
      y = y_grid,
      weights = pi1,
      theta = theta1_4,
      sigma = sigma_kernel
    )

  densities[["4"]][["Group 2"]][b, ] =
    mixture_density(
      y = y_grid,
      weights = pi2,
      theta = theta2_4,
      sigma = sigma_kernel
    )

  final_entropy[["4"]][b] =
    final_entropy_4(
      z1 = z1,
      z2 = z2
    )

  mu_draws[["4"]]$mu1[b] =
    mixture_mean(
      weights = pi1,
      theta = theta1_4
    )

  mu_draws[["4"]]$mu2[b] =
    mixture_mean(
      weights = pi2,
      theta = theta2_4
    )


  # ----------------------------------------------------------
  # 5: spike-and-slab common atoms
  # ----------------------------------------------------------

  is_zero = rbinom(
    L,
    size = 1,
    prob = omega_5
  ) == 1

  theta_5 = rP0(L)
  theta_5[is_zero] = 0

  densities[["5"]][["Group 1"]][b, ] =
    mixture_density(
      y = y_grid,
      weights = pi1,
      theta = theta_5,
      sigma = sigma_kernel
    )

  densities[["5"]][["Group 2"]][b, ] =
    mixture_density(
      y = y_grid,
      weights = pi2,
      theta = theta_5,
      sigma = sigma_kernel
    )

  final_entropy[["5"]][b] =
    final_entropy_5(
      z1 = z1,
      z2 = z2,
      theta = theta_5
    )

  mu_draws[["5"]]$mu1[b] =
    mixture_mean(
      weights = pi1,
      theta = theta_5
    )

  mu_draws[["5"]]$mu2[b] =
    mixture_mean(
      weights = pi2,
      theta = theta_5
    )


  # ----------------------------------------------------------
  # 6: ANOVA-type covariate-dependent atoms
  # ----------------------------------------------------------

  effects_6 = draw_6_effects(
    L = L,
    V = V,
    W = W,
    sd_m = sd_m_6,
    sd_A = sd_A_6,
    sd_B = sd_B_6
  )

  # Store densities for ALL covariate configurations.
  # Also store their means in this Monte Carlo iteration.
  mu1_6_configs = numeric(n_configs_6)
  mu2_6_configs = numeric(n_configs_6)

  for (cc in seq_len(n_configs_6)) {

    vv = configs_6$v[cc]
    ww = configs_6$w[cc]

    theta_now = theta_6(
      effects = effects_6,
      v = vv,
      w = ww
    )

    densities_6_all[["Group 1"]][b, , cc] =
      mixture_density(
        y = y_grid,
        weights = pi1,
        theta = theta_now,
        sigma = sigma_kernel
      )

    densities_6_all[["Group 2"]][b, , cc] =
      mixture_density(
        y = y_grid,
        weights = pi2,
        theta = theta_now,
        sigma = sigma_kernel
      )

    mu1_6_configs[cc] =
      mixture_mean(
        weights = pi1,
        theta = theta_now
      )

    mu2_6_configs[cc] =
      mixture_mean(
        weights = pi2,
        theta = theta_now
      )
  }

  # For the Model 6 Row-3 summary, use the mean of the density
  # marginalized over the categorical covariate distribution:
  #
  #   mu_j = sum_x Pr(x) * int y f_{j,x}(y) dy.
  #
  # With the default uniform covariate distribution this is simply
  # the average across the four configurations.
  mu_draws[["6"]]$mu1[b] =
    sum(
      config_prob_6 *
        mu1_6_configs
    )

  mu_draws[["6"]]$mu2[b] =
    sum(
      config_prob_6 *
        mu2_6_configs
    )

  # Observation-specific categorical covariates used to construct
  # the final tie-induced partition.
  v1 = sample.int(
    V,
    size = n1,
    replace = TRUE,
    prob = prob_v
  )

  w1 = sample.int(
    W,
    size = n1,
    replace = TRUE,
    prob = prob_w
  )

  v2 = sample.int(
    V,
    size = n2,
    replace = TRUE,
    prob = prob_v
  )

  w2 = sample.int(
    W,
    size = n2,
    replace = TRUE,
    prob = prob_w
  )

  final_entropy[["6"]][b] =
    final_entropy_6(
      z1 = z1,
      z2 = z2,
      v1 = v1,
      w1 = w1,
      v2 = v2,
      w2 = w2
    )


  if (b %% 500L == 0L) {
    message(
      "Monte Carlo iteration ",
      b,
      " / ",
      B_MC
    )
  }
}


# ------------------------------------------------------------
# 6. Density summaries
# ------------------------------------------------------------

summarise_density_matrix <- function(
  mat,
  model,
  group,
  y_grid
) {

  data.frame(
    model = model,
    group = group,
    y = y_grid,
    mean = colMeans(mat),
    lower = apply(
      mat,
      2L,
      quantile,
      probs = 0.025
    ),
    upper = apply(
      mat,
      2L,
      quantile,
      probs = 0.975
    )
  )
}

# Star--4,5
density_summary = do.call(
  rbind,
  lapply(
    c("Star", "4", "5"),
    function(m) {
      do.call(
        rbind,
        lapply(
          group_names,
          function(g) {
            summarise_density_matrix(
              mat = densities[[m]][[g]],
              model = m,
              group = g,
              y_grid = y_grid
            )
          }
        )
      )
    }
  )
)

# 6: all V x W configurations
summarise_density_array <- function(
  arr,
  group,
  y_grid,
  configs,
  config_labels
) {

  out = lapply(
    seq_len(nrow(configs)),
    function(cc) {

      data.frame(
        model = "6",
        group = group,
        config = config_labels[cc],
        v = configs$v[cc],
        w = configs$w[cc],
        y = y_grid,
        mean = colMeans(
          arr[, , cc]
        ),
        lower = apply(
          arr[, , cc],
          2L,
          quantile,
          probs = 0.025
        ),
        upper = apply(
          arr[, , cc],
          2L,
          quantile,
          probs = 0.975
        )
      )
    }
  )

  do.call(rbind, out)
}

density_summary_6_all = rbind(
  summarise_density_array(
    arr = densities_6_all[["Group 1"]],
    group = "Group 1",
    y_grid = y_grid,
    configs = configs_6,
    config_labels = config_labels_6
  ),

  summarise_density_array(
    arr = densities_6_all[["Group 2"]],
    group = "Group 2",
    y_grid = y_grid,
    configs = configs_6,
    config_labels = config_labels_6
  )
)


# ------------------------------------------------------------
# 7. Plotting data frames
# ------------------------------------------------------------

# Row 1: FINAL partition entropy
final_entropy_df = do.call(
  rbind,
  lapply(
    model_names,
    function(m) {
      data.frame(
        model = m,
        entropy = final_entropy[[m]]
      )
    }
  )
)

# Row 3: joint prior means
mu_df = do.call(
  rbind,
  lapply(
    model_names,
    function(m) {
      data.frame(
        model = m,
        mu1 = mu_draws[[m]]$mu1,
        mu2 = mu_draws[[m]]$mu2
      )
    }
  )
)

model_titles = c(
  "Star" = "mSSP",
  "4" = "full-range borrowing",
  "5" = "spike-and-slab",
  "6" = "ANOVA-type"
)


# ------------------------------------------------------------
# 8. Build the 3 x 4 figure
# ------------------------------------------------------------

theme_illustration = theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 10
    ),
    legend.title = element_blank()
  )


# ------------------------------------------------------------
# Row 1:
# entropy of the FINAL partition Pi induced by atom ties
# ------------------------------------------------------------

final_entropy_values = unlist(
  final_entropy,
  use.names = FALSE
)

entropy_min = max(
  0,
  min(
    final_entropy_values,
    na.rm = TRUE
  )
)

entropy_max = max(
  final_entropy_values,
  na.rm = TRUE
)

entropy_pad = 0.04 * max(
  entropy_max - entropy_min,
  1e-8
)

entropy_lims = c(
  max(
    0,
    entropy_min - entropy_pad
  ),
  entropy_max + entropy_pad
)

top_plots = lapply(
  model_names,
  function(m) {

    dat = final_entropy_df[
      final_entropy_df$model == m,
    ]

    ggplot(
      dat,
      aes(x = entropy)
    ) +
      geom_density(
        fill = "grey80",
        colour = "black",
        linewidth = 0.55,
        adjust = 1
      ) +
      coord_cartesian(
        xlim = entropy_lims
      ) +
      labs(
        title = model_titles[[m]],
        x = "Entropy of pooled cluster frequencies",
        y = "Prior density"
      ) +
      theme_illustration +
      theme(
        legend.position = "none"
      )
  }
)


# ------------------------------------------------------------
# Row 2:
# random mixture densities
# ------------------------------------------------------------

density_ymax = 1.03 * max(
  c(
    density_summary$upper,
    density_summary_6_all$upper
  ),
  na.rm = TRUE
)

middle_plots = lapply(
  model_names,
  function(m) {

    # ----------------------------------------
    # Star-4-5
    # ----------------------------------------
    if (m != "6") {

      dat = density_summary[
        density_summary$model == m,
      ]

      p = ggplot(
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
          aes(
            ymin = lower,
            ymax = upper
          ),
          alpha = 0.18,
          colour = NA
        ) +
        geom_line(
          linewidth = 0.75
        ) +
        coord_cartesian(
          xlim = range(y_grid),
          ylim = c(
            0,
            density_ymax
          )
        ) +
        labs(
          x = expression(y),
          y = "Mixture density"
        ) +
        theme_illustration

      # Put the Group 1 / Group 2 legend inside Star only.
      if (m == "Star") {

        p = p +
          theme(
            legend.position = "inside",
            legend.position.inside = c(
              0.05,
              0.95
            ),
            legend.justification = c(
              0,
              1
            ),
            legend.direction = "horizontal",
            legend.background = element_rect(
              fill = "white",
              colour = NA
            )
          ) +
          guides(
            colour = guide_legend(
              nrow = 1
            )
          )

      } else {

        p = p +
          theme(
            legend.position = "none"
          )
      }

      return(p)
    }

    # ----------------------------------------
    # 6:
    # all covariate configurations
    # ----------------------------------------

    dat = density_summary_6_all

    ggplot(
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
        aes(
          ymin = lower,
          ymax = upper
        ),
        alpha = 0.18,
        colour = NA
      ) +
      geom_line(
        linewidth = 0.60
      ) +
      facet_wrap(
        ~ config,
        nrow = V
      ) +
      coord_cartesian(
        xlim = range(y_grid),
        ylim = c(
          0,
          density_ymax
        )
      ) +
      labs(
        x = expression(y),
        y = "Mixture density"
      ) +
      theme_illustration +
      theme(
        legend.position = "none",
        strip.background = element_rect(
          fill = "white"
        ),
        strip.text = element_text(
          size = 7
        ),
        axis.text = element_text(
          size = 7
        ),
        axis.title = element_text(
          size = 8
        ),
        panel.spacing = grid::unit(
          0.12,
          "lines"
        )
      )
  }
)


# ------------------------------------------------------------
# Row 3:
# joint prior distribution of (mu_1, mu_2)
# ------------------------------------------------------------

mu_min = min(
  c(
    mu_df$mu1,
    mu_df$mu2
  ),
  na.rm = TRUE
)

mu_max = max(
  c(
    mu_df$mu1,
    mu_df$mu2
  ),
  na.rm = TRUE
)

mu_pad = 0.05 * max(
  mu_max - mu_min,
  1e-8
)

mu_lims = c(
  mu_min - mu_pad,
  mu_max + mu_pad
)

bottom_plots = lapply(
  model_names,
  function(m) {

    dat = mu_df[
      mu_df$model == m,
    ]

    ggplot(
      dat,
      aes(
        x = mu1,
        y = mu2
      )
    ) +
      geom_point(alpha = 0.6, size = 0.8, colour = "grey80") +
      stat_density_2d(linewidth = 0.45, colour = "black") +
      geom_abline(
        intercept = 0,
        slope = 1,
        linetype = 2
      ) +
      coord_cartesian(
        xlim = mu_lims,
        ylim = mu_lims
      ) +
      labs(
        x = expression(mu[1]),
        y = expression(mu[2])
      ) +
      theme_illustration +
      theme(
        legend.position = "none"
      )
  }
)


# ------------------------------------------------------------
# Combine the 3 x 4 figure
# ------------------------------------------------------------

figure_B =
  (
    wrap_plots(
      top_plots,
      nrow = 1
    ) /
      wrap_plots(
        middle_plots,
        nrow = 1
      ) /
      wrap_plots(
        bottom_plots,
        nrow = 1
      )
  ) +
  plot_layout(
    heights = c(
      1,
      1.40,
      1.15
    )
  )

print(figure_B)


# ------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------

ggsave(
  filename = "illustration_B.pdf",
  plot = figure_B,
  width = 12,
  height = 8,
  units = "in"
)

ggsave(
  filename = "illustration_B.png",
  plot = figure_B,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)


# ------------------------------------------------------------
# 10. Optional numerical summaries / checks
# ------------------------------------------------------------

cat(
  "\nMean entropy of latent partition Pi^0:\n"
)

cat(
  mean(latent_entropy),
  "\n"
)

cat(
  "\nMean entropy of final tie-induced partition Pi:\n"
)

for (m in model_names) {

  cat(
    m,
    ": ",
    mean(final_entropy[[m]]),
    "\n",
    sep = ""
  )
}

cat(
  "\nExpected partition-ordering checks:\n"
)

cat(
  "Star: final entropy equals latent entropy exactly = ",
  isTRUE(
    all.equal(
      final_entropy[["Star"]],
      latent_entropy
    )
  ),
  "\n",
  sep = ""
)

cat(
  "4: proportion H(Pi) >= H(Pi^0) = ",
  mean(
    final_entropy[["4"]] >=
      latent_entropy - 1e-12
  ),
  "\n",
  sep = ""
)

cat(
  "5: proportion H(Pi) <= H(Pi^0) = ",
  mean(
    final_entropy[["5"]] <=
      latent_entropy + 1e-12
  ),
  "\n",
  sep = ""
)

cat(
  "6: proportion H(Pi) >= H(Pi^0) = ",
  mean(
    final_entropy[["6"]] >=
      latent_entropy - 1e-12
  ),
  "\n",
  sep = ""
)

cat(
  "\nEmpirical prior correlations of (mu1, mu2):\n"
)

for (m in model_names) {

  cat(
    m,
    ": ",
    cor(
      mu_draws[[m]]$mu1,
      mu_draws[[m]]$mu2
    ),
    "\n",
    sep = ""
  )
}

cat(
  "\n6 covariate configurations shown in Row 2:\n"
)

print(configs_6)

cat(
  "\n6 Row-3 means are marginalized over these configuration probabilities:\n"
)

print(
  data.frame(
    configs_6,
    probability = config_prob_6
  )
)
