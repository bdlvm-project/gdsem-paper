library(brms)
library(dplyr)
library(here)
library(ggplot2)
library(bayesplot)
library(tikzDevice)
library(patchwork)
library(tarchetypes)
library(targets)
# Analysis based on https://osf.io/bdqe8
fix_chains = function(mod) {
  for (i in seq_along(mod$fit@sim$samples)) {
    pars = grepl("bsp.+neuro", colnames(mod$fit@sim$samples[[i]]))
    samp = mod$fit@sim$samples[[i]][, pars]
    if(mean(samp[,1]) < 0)
      mod$fit@sim$samples[[i]][, pars] = -samp
  }
  mod
}

data_load = \(x){
  d = readRDS(x) %>%
  mutate(
    # add NA for latent variables
    neuro = NA_real_, nega = NA_real_,
    diary_day = as.numeric(diary_day)
  ) %>%
  rename(
    goodmood = good_mood,
    selfesteem = self_esteem
  ) %>%
  mutate(
    goodmood = 4 - goodmood,
    selfesteem = 4 - selfesteem,
  )
  # add censoring indicator
  for (i in seq_len(ncol(d))) {
    d[, i] = unclass(d[, i])
    q = colnames(d)[i]
    if (grepl("_neuro_", q))
      d[, paste0(q, "_cens")] =
        ifelse(d[, i] == 1, "left",
        ifelse(d[, i] == 5, "right",
                            "none"))
    if (q %in% c(
      "stressed", "irritable", "loneliness",
      "goodmood", "selfesteem"))
      d[, paste0(q, "_cens")] =
        ifelse(is.na(d[,i]), "interval",
        ifelse(d[, i] == 0, "left",
        ifelse(d[, i] == 4, "right",
                            "none")))
  }
  d %>%
    filter(!(
      is.na(stressed) |
      is.na(irritable) |
      is.na(loneliness) |
      is.na(goodmood) |
      is.na(selfesteem)))
}

tar_plan(
  # --- load data
  tar_target(
    data_file,
    here::here("data", "diary_gocd2_luna_fazio_05062023.rds"),
    format = "file"
  ),
  data = data_load(data_file) %>%
    mutate(id = 1:n()) %>%
    mutate(
      # add first-measurement indicator
      first_day = diary_day == min(diary_day),
      idx = id[first_day],
      .by = woman
    ),
  # --- formula
  f01 =
    bf(neuro | mi() + subset(first_day) + index(id) ~ 0) +
    bf(bfi_neuro_1 | cens(bfi_neuro_1_cens) + subset(first_day) ~ mi(neuro, idx)) +
    bf(bfi_neuro_2r| cens(bfi_neuro_2r_cens)+ subset(first_day) ~ mi(neuro, idx)) +
    bf(bfi_neuro_3 | cens(bfi_neuro_3_cens) + subset(first_day) ~ mi(neuro, idx)) +
    bf(bfi_neuro_4 | cens(bfi_neuro_4_cens) + subset(first_day) ~ mi(neuro, idx)) +
    bf(bfi_neuro_5r| cens(bfi_neuro_5r_cens)+ subset(first_day) ~ mi(neuro, idx)) +
    bf(bfi_neuro_6r| cens(bfi_neuro_6r_cens)+ subset(first_day) ~ mi(neuro, idx)) +
    bf(bfi_neuro_7 | cens(bfi_neuro_7_cens) + subset(first_day) ~ mi(neuro, idx)) +
    bf(bfi_neuro_8 | cens(bfi_neuro_8_cens) + subset(first_day) ~ mi(neuro, idx)) +
    bf(
      nega | mi() ~ 0 + mi(neuro, idx) + (1 | woman),
      sigma ~ mi(neuro, idx) + (1 | woman)
    ) +
    bf(stressed | cens(stressed_cens) ~ mi(nega)) +
    bf(irritable | cens(irritable_cens) ~ mi(nega)) +
    bf(loneliness | cens(loneliness_cens) ~ mi(nega)) +
    bf(goodmood | cens(goodmood_cens) ~ mi(nega)) +
    bf(selfesteem | cens(selfesteem_cens) ~ mi(nega)) +
    set_rescor(FALSE),
  # --- priors
  p01 = {
    p = get_prior(f01, data)
    p = p[p$resp != "",]
    # If we use ULI, the LV will have the same scale as the items (1-5)
    # so in the most extreme case one item moving 1 -> 2 could go along
    # with another item moving 1 -> 5 which means a loading of 0.25/4
    # depending on which one was used as reference
    p[p$class == "b" & p$coef != "", 1] = "normal(0,2)"
    # In the case of item-specific variance, a theoretically completely
    # unrelated item with a very low mean value would still be able to
    # span the whole scale if given sd = 2
    p[p$class == "sigma", 1] = "gamma(5, 5)"
    # Latent to latent pars
    p[p$resp == "nega" , 1] = "normal(0,2)"
    # The exp-scale scedastic parameters can't be too large
    p[p$dpar == "sigma" , 1] = "normal(0,0.25)"
    # Random intercepts... pending sanity check
    p[p$class == "sd", 1] = "normal(0, 0.25)"
    # set ULI constraint
    p[p$class == "b" & p$resp == "goodmood", 1] = "constant(1)"
    # set UVI constraint
    p[p$class == "sigma" & p$resp == "neuro", 1] = "constant(1)"
    p
  },
  m01_full =
    brm(f01, data, prior = p01,
      backend = "cmdstanr", cores = 4, init = 0),
  m01_final = fix_chains(m01_full),
  m01_plot = {
    draws = as_draws_df(m01_final)[,c(1:44,5768:5770)]
    pars = colnames(draws)
    # item
    rearr = c(1:8,12,9,10,11,13)
    item_intr = which(grepl("b_[a-zA-Z0-9]+_Intercept$", pars))[rearr]
    item_stdv = which(grepl("^sigma_[bsigl].+", pars))[rearr]
    item_load = which(grepl("^bsp_[bsigl][a-zA-Z0-9]+_[a-zA-Z0-9]+$", pars))[c(1:8,13,9:12)]
    # factors
    lvar_pars = setdiff(1:44, c(item_intr, item_stdv, item_load))[c(2,4,1,3,5)]
    # labels
    colnames(draws)[item_intr] = c(paste0("$\\nu_{\\mathrm{",
      c(paste0("Ne},", 1:8), paste0("Em},", 1:5)), "}$"))
    colnames(draws)[item_stdv] = c(paste0("$\\sigma_{\\mathrm{",
      c(paste0("Ne},", 1:8), paste0("Em},", 1:5)), "}$"))
    colnames(draws)[item_load] = c(paste0("$\\lambda_{\\mathrm{",
      c(paste0("Ne},", 1:8), paste0("Em},", 1:5)), "}$"))
    colnames(draws)[lvar_pars] = c(
      "$\\beta_{1 \\mu_\\mathrm{Em}}$",
      "$\\sigma_{\\mu_\\mathrm{Em}}$",
      "$\\beta_{0 \\sigma_\\mathrm{Em}}$",
      "$\\beta_{1 \\sigma_\\mathrm{Em}}$",
      "$\\sigma_{\\sigma_\\mathrm{Em}}$")
    # plots
    mytheme = theme_bw() + theme(
      panel.border = element_blank(), 
      axis.line = element_line(colour = "black")
    )
    lvar_plot = mcmc_intervals(draws[,c(lvar_pars, 45:47)], prob_outer = 0.95) +
      mytheme + scale_x_continuous(breaks = 0.1*c(-3:6), minor_breaks = NULL)
    intr_plot = mcmc_intervals(draws[,c(item_intr, 45:47)], prob_outer = 0.95) + mytheme +
      scale_x_continuous(minor_breaks = NULL) # Not used
    stdv_plot = mcmc_intervals(draws[,c(item_stdv, 45:47)], prob_outer = 0.95) + mytheme +
      scale_x_continuous(minor_breaks = NULL)
    load_plot = mcmc_intervals(draws[,c(item_load, 45:47)], prob_outer = 0.95) + mytheme +
      scale_x_continuous(minor_breaks = NULL)
    # final layout
    (
      (lvar_plot + ggtitle("Latent variable parameters")) /
      ((load_plot + ggtitle("Item parameters")) + stdv_plot)
    ) + plot_layout(heights = c(0.4, 0.6))
  },
  tar_target(
    tikz_plot, {
      path = here("case-study", "case_study.tex")
      tikzDevice::tikz(file = path, width = 5.5, height = 4)
      print(m01_plot)
      dev.off()
      path
    },
    format = "file"
  )
)