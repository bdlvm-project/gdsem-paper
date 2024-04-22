# Preprocessing -----
utils.named_tibble = function(x) {
  objname <- deparse(substitute(x))
  eval(parse(text = paste0(
    "tibble(id = names(x), ", objname, " = unname(x))"
  )))
}

# cross a vector's values with all of the rows in d
utils.cross_vec = function(d, v) {
  objname <- deparse(substitute(v))
  eval(parse(text = paste0(
    "d %>%
    bind_cols(tibble(", objname, " = list(v))) %>%
    unnest(", objname, ")"
  )))
}

# get response variables in brms formula
utils.get_resp = function(x, pattern) {
  if(inherits(x, "bform"))
    r = x$responses
  if(inherits(x, "brmsprior"))
    r = x$resp
  unique(r[grepl(pattern, r)])
}

# generate data with additional handling:
# - set vars matching na_pattern to NA
# - cache the result and return filepath
utils.make_data = \(
  condition, ..., na_pattern, drop = "nonconstant"
) {
  stopifnot(identical(nrow(condition), 1L))
  filepath = with(condition, {
    formula = eval(parse(text = formulas[[1]])) +
      bf(nonconstant ~ 1) + set_rescor(FALSE)
    generator = SBC::SBC_generator_brms(
      formula,
      data_template[[1]][rep(1, num_obs), ],
      prior = priors_gen[[1]],
      stanvars = bdlvm_stanvars, generate_lp = TRUE,
      thin = 50, warmup = 10000, refresh = 2000
    )
    d = SBC::generate_datasets(generator, num_rep)
    d = utils.format_data(d, formula, na_pattern, drop)
    # cache and return path
    filepath = here(
      tar_path_store(), "user", paste0(fullname, ".data")
    )
    qsave(d, filepath)
    filepath # return
  })
  condition %>% mutate(data = filepath)
}

utils.format_data = \(d, formula, na_pattern, drop) {
  data_ok = sapply(d$generated, \(x) all(!is.na(x)))
  # drop simulations with any NA/NaN values
  d$variables = d$variables[data_ok,]
  d$generated = d$generated[data_ok]
  # drop user-specified columns
  d$variables = d$variables[, !grepl(drop, colnames(d$variables))]
  # put missing data where needed
  make_na = utils.get_resp(formula, na_pattern)
  if(!is.null(make_na)) {
    d$na_values = lapply(d$generated,\(x)x[, make_na])
    d$generated = lapply(d$generated,\(x){
      x[, make_na] = NA_real_
      x
    })
  }
  return(d)
}

# Fitting -----
utils.sbc_compute = function(condition, ...) {
  with(condition, {
    formula = eval(parse(text = formulas[[1]])) +
      set_rescor(FALSE)
    prior = priors_fit[[1]]
    prior = prior[prior$source != "ignore",]
    backend = SBC_backend_brms(
      formula,
      template_data = data_template[[1]],
      prior = prior,
      stanvars = bdlvm_stanvars,
      chains = 4, thin = 1, warmup = 1000, iter = 2000
    )
    data = qread(data[[1]])
    SBC::compute_SBC(data, backend, ...)
  })
}

# Postprocessing ----
## Posterior summaries ----
utils.summarize_sbc = function(meta_df, res) {
  data = qs::qread(meta_df$data)
  for (i_sim in seq_along(res$fits)) {
    fit = res$fits[[i_sim]]
    pars = data$variables[i_sim,]
    lvs = data$na_values[[i_sim]]
    # 1. bias + rmse
    res$stats[res$stats$sim_id == i_sim, "bias"] = utils.err_col(fit, pars, 1)
    res$stats[res$stats$sim_id == i_sim, "rmse"] = utils.err_col(fit, pars, 2)
    # 2. chain times
  res$stats[res$stats$sim_id == i_sim, "essb_sec"] = {
    res$stats[res$stats$sim_id == i_sim,
      "ess_bulk"]/sum(rstan::get_elapsed_time(fit$fit))
  }
  res$stats[res$stats$sim_id == i_sim, "esst_sec"] = {
    res$stats[res$stats$sim_id == i_sim,
      "ess_tail"]/sum(rstan::get_elapsed_time(fit$fit))
  }
    # 3. get logliks
    ll_row = utils.ll_row(fit, pars, lvs)
    res$stats = bind_rows(
      res$stats,
      tibble(
        sim_id = i_sim,
        variable = "loglik",
        simulated_value = ll_row$val,
        rank = ll_row$rank,
        mean = ll_row$mean,
        rhat = ll_row$rhat,
        max_rank = unique(res$stats$max_rank)
      )
    )
  }
  meta_df %>%
    bind_cols(res$stats) %>%
    nest(.by = c(fullname, priors_fit), .key = "fit_dx") %>%
    mutate(backend_dx = list(res$backend_diagnostics))
}

### ^Helpers ----
utils.err_col = function(fit, pars, exponent) {
  draws = as.data.frame(as_draws_df(fit))[, colnames(pars)]
  apply(
    apply(draws, 1, `-`, pars),
    1, \(x) mean(x^exponent)
  )^(1/exponent)
}

utils.ll_row = function(fit, pars, lvs) {
    ll_ref = utils.true_ll(fit, pars, lvs)
    ll_sim = utils.post_ll(fit, colnames(lvs))
    list(
      val = ll_ref,
      rank = floor(sum(ll_sim < ll_ref)/10),
      mean = mean(ll_sim),
      rhat = rhat(ll_sim)
    )
  }

utils.true_ll = function(fit, pars, lvs) {
  newdata = fit$data
  # get true latent values
  for(lv in colnames(lvs)) {
    newdata[,lv] = lvs[lv]
  }
  # replace draws for true parameters
  fit$fit@sim$samples = fit$fit@sim$samples[1]
  draws = fit$fit@sim$samples[[1]]
  fit$fit@sim$samples[[1]] = cbind(
    pars,
    draws[1, setdiff(colnames(draws), colnames(pars))]
  )
  # compute true log_lik
  sum(log_lik(
    fit,
    newdata = newdata,
    resp = setdiff(colnames(fit$data), colnames(lvs))
  ))
}

utils.post_ll = function(fit, lvs) {
  matrix(apply(
    log_lik(fit, resp = setdiff(colnames(fit$data), lvs)),
    1, sum
  ), ncol = 4)
}

## Rearrange summaries ----
utils.regroup = function(summ_df) {
  bck_dx = summ_df %>%
    select(fullname, priors_fit, backend_dx) %>%
    unnest(backend_dx) %>%
    mutate(
      id = substr(fullname, 1, 5),
      batch = substr(fullname, 26, 30)
    ) %>%
    mutate(sim_id = as.numeric(as.factor(paste(batch, sim_id)))) %>%
    select(id, priors_fit, sim_id, n_divergent, time = max_chain_time)

  summ_df %>%
    select(fullname, priors_fit, fit_dx) %>%
    unnest(fit_dx) %>%
    mutate(sim_id = as.numeric(as.factor(paste(batch, sim_id)))) %>%
    select(-c(fullname, num_obs, num_rep, batch, priors_gen, data)) %>%
    utils.alt_names %>%
    left_join(bck_dx, by = c("id", "priors_fit", "sim_id")) %>%
    nest(.by = c(id, priors_fit, sim_id, n_divergent, time)) %>%
    nest(.by = c(id, priors_fit)) %>%
    nest(.by = id)
}

### ^Helpers ----
utils.alt_names = function(df) {
  mutate(df,
    param_type = case_when(
      grepl("^bsp_x.+_miy[0-9]+$", variable) ~
        #"Factor loading",
        r"($\lambda$)",
      grepl("^sigma_x[0-9]+i[0-9]+$", variable) ~
        #"Item std. dev.",
        r"($\tau$)",
      grepl("^(sigma_y[0-9]+)|(Intercept_sigma_y[0-9]+)$", variable) ~
        #"Factor std. dev.",
        r"($\beta_{0 \sigma}$)",
      grepl("^bsp_y[0-9]+_.*miy.+$", variable) ~
        #"Slope on mean",
        r"($\beta_\mu$)",
      id == "sim01" & variable == "bsp_sigma_y2_miy1" ~
        #"Slope on std. dev. for sim01",
        r"($\beta_{1 \sigma}$)",
      grepl("^bsp_sigma_y[0-9]+_.*miy.+$", variable) ~
        #"Slope on std. dev.",
        r"($\beta_{\geq 1 \sigma}$)",
      variable == "loglik" ~ "Log-likelihood"
    ) %>%
    factor(levels = c(
      r"($\beta_\mu$)",
      r"($\beta_{0 \sigma}$)",
      r"($\beta_{1 \sigma}$)",
      r"($\beta_{\geq 1 \sigma}$)",
      r"($\lambda$)",
      r"($\tau$)",
      "Log-likelihood"
    ), ordered = TRUE),
  alt_variable = param_type
  )
}

# Plotting ----
utils.rhat_heatmap = function(df) {
  df = df %>%
  mutate(rmean = mean(rhat), .by = sim_id) %>%
  summarise(
    rhat = max(rhat),
    .by = c(sim_id, variable, rmean)) %>%
  mutate(
    ranking = dense_rank(rmean),
    categories = case_when(
      rhat < 1.005 ~ "A", # 1.00
      rhat < 1.015 ~ "B", # 1.01
      rhat < 1.055 ~ "C", # <= 1.05
      TRUE ~ "D"          # >  1.05
  ))
  df$variable = fct_rev(df$variable)
  ggplot(df, aes(x = ranking, y = variable, fill = categories)) +
  geom_tile() +
  scale_y_discrete(expand = c(0,0),
    labels = TeX(unique(as.character(df$variable)))) +
  scale_x_continuous(expand = c(0,0), limits = c(-1e-9,251),
    breaks = (1:5)*50, minor_breaks = 249.9) +
  ggtitle("(b) Convergence") +
  scale_fill_viridis_d(
    name = TeX(r"($\hat{R}$)"),
    labels = c("1.00", "1.01", TeX(r"($\leq 1.05$)"), "> 1.05"),
    direction =-1, option = "mako") +
    theme_bw() + theme(
      axis.title = element_blank(),
      legend.key.size = unit(0.25, "cm"),
      legend.title = element_text(size = 9), 
      legend.text = element_text(size = 7.5),
      legend.margin = margin(t=0,r=0,b=0,l=-8),
      panel.grid = element_blank(),
      panel.grid.minor.x = element_line(color = "grey40", linewidth = 1.1)
    )
}

utils.sbc_ecdf = function(df) {
  grouping = df[, c("variable", "param_type")] %>%
    unique %>%
    summarise(
      param_type = setNames(list(variable), unique(param_type)),
      .by = param_type
    ) %>% pull(param_type)
  SBC::plot_ecdf_diff(
    df, combine_variables = grouping
  ) + facet_grid(
    . ~ (\(x){
      factor(
        TeX(as.character(x), output = "character"),
        levels = c(
          "beta[mu]", "beta[0*sigma]", "beta[1*sigma]",
          "beta[phantom() * {phantom() >= phantom()} * 1*sigma]",
          "lambda", "tau", "'Log-likelihood'"
        )
      )
    })(group),
    labeller = label_parsed
  ) + ggtitle("(a) Calibration") +
  theme_bw() + theme(
    legend.position = "none",
    axis.text.x =
      element_text(angle = 35, size = 8, hjust = 1)
  )
}

utils.recov_box = function(df, metric, title) {
  ggplot(df, aes(x = param_type, y = {{metric}})) +
    geom_boxplot(
      outlier.size = 0.75, outlier.alpha = 0.5, outlier.shape = 1
    ) + theme_bw() + ylab(title) + scale_x_discrete(
      labels = TeX(unique(as.character(df$param_type)))
    ) + theme(
      axis.title.x = element_blank(),
      plot.background = element_blank()
    )
}

utils.make_figures = function(df) {
  fit_dx = df %>% unnest(data) %>%
    filter(priors_fit == "weak") %>%
    unnest(data) %>% unnest(data) %>%
    filter(
      !variable %in% c("loglik", "lprior", "lp__"),
      !is.na(rhat),
      !grepl("^b_.+_Intercept$", variable)
    ) %>%
    mutate(variable = alt_variable)
  
  recov_dx = fit_dx %>%
  summarise(
    rhat = max(rhat),
    bias = mean(bias),
    rmse = mean(rmse),
    .by = c(sim_id, param_type)
  ) %>%
  mutate(max_rhat = max(rhat), .by = sim_id) %>%
  filter(max_rhat < 1.055)

  sbc_dx = df %>% unnest(data) %>%
    filter(priors_fit == "sbc") %>%
    unnest(data) %>% unnest(data) %>% 
    filter(
      !variable %in% c("lprior", "lp__"),
      !is.na(rhat),
      !grepl("^b_.+_Intercept$", variable)
    )

  p = (
    (utils.sbc_ecdf(sbc_dx)) /
    (
      (utils.rhat_heatmap(fit_dx) |
      utils.recov_box(recov_dx, bias, "Bias") +
        ggtitle("(c) Parameter recovery") |
      utils.recov_box(recov_dx, rmse, "RMSE")) +
      plot_layout(widths = c(0.5,0.25,0.25))
    )
  ) + plot_layout(heights = c(0.5, 0.5))
  attr(p, "plot_id") = df$id
  p
}

utils.ess_plot = \(df) {
  df %>% unnest(data) %>%
    filter(priors_fit == "weak") %>%
    unnest(data) %>% unnest(data) %>%
    filter(
      !variable %in% c("loglik", "lprior", "lp__"),
      !is.na(rhat),
      !grepl("^b_.+_Intercept$", variable)
    ) %>%
    mutate(max_rhat = max(rhat),
      .by = c(id, sim_id)) %>%
    filter(max_rhat < 1.055) %>%
    transmute(
      id, sim_id,
      param_type = fct_recode(param_type,
        "$\\beta_{\\geq 1 \\sigma}$" = "$\\beta_{1 \\sigma}$"),
      `Bulk ESS` = ess_bulk, `Tail ESS` = ess_tail,
      `Bulk ESS/sec` = essb_sec, `Tail ESS/sec` = esst_sec,
      `Wall time` = time
    ) %>%
    pivot_longer(
      cols = contains("ess") | contains("time")
    ) %>%
    summarise(
      value = mean(value),
      .by = c(id, sim_id, param_type, name)
    ) %>%
    mutate(
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75),
      .by = c(id, param_type, name)
    ) %>% (\(x) {(
      utils.ess_box(x %>% filter(grepl("ESS$", name))) +
        theme(
          legend.position = "none",
          strip.text.y = element_blank(),
        ) |
      utils.ess_box(x %>% filter(grepl("sec$", name))) +
        theme(legend.position = "none")
      ) / (
      free(guide_area() +
      utils.ess_box(x %>% filter(grepl("^Wall", name)) %>%
        mutate(param_type = "", id = factor(id)) %>% unique
      ) + scale_y_log10(name = "seconds") +
      scale_x_discrete(limits = rev) +
      coord_flip() + theme(
        strip.text.y = element_blank(),
        axis.text.x = element_text(),
        axis.title.x = element_text(),
        axis.ticks.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) + 
      plot_layout(
        guides = "collect",
        widths = c(0.7,0.3)
      ))) + plot_layout(heights = c(0.8, 0.2))
    })
}

utils.ess_box = \(df) {
    ggplot(df, aes(x = id, y = value, fill = id, color = id)) +
    geom_boxplot(
      color = "black", outlier.size = 0.75,
      outlier.alpha = 0.5, outlier.shape = 1
    ) +
    geom_boxplot(
      color = NA, fatten = 0,
      outlier.shape = NA
    ) +
    geom_boxplot(
      color = "#00000055", fill = NA, outlier.shape = NA,
      linewidth = 0.1, fatten = 10
    ) +
    scale_fill_manual(name = "Model",
      labels = c("Two-factor", "Mediation", "Interaction", "Sequential"),
      values = c('#BBCCEE', '#CCEEFF', '#CCDDAA', '#EEEEBB'),
    ) + facet_grid(
      (\(x){
        factor(
          TeX(as.character(x), output = "character"),
          levels = c(
            "beta[mu]", "beta[0*sigma]", "beta[1*sigma]",
            "beta[phantom() * {phantom() >= phantom()} * 1*sigma]",
            "lambda", "tau", "'Log-likelihood'"
          )
        )
      })(param_type) ~ TeX(name, output = "character"),
      labeller = label_parsed
    ) +
    theme_bw() + theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.y = element_text(angle = 0)
    )
}

utils.fig_appendix = function(df) {
  fit_dx = df %>% unnest(data) %>%
    filter(priors_fit == "sbc") %>%
    unnest(data) %>% unnest(data) %>%
    filter(
      !variable %in% c("loglik", "lprior", "lp__"),
      !is.na(rhat),
      !grepl("^b_.+_Intercept$", variable)
    ) %>%
    mutate(variable = alt_variable)
  
  recov_dx = fit_dx %>%
  summarise(
    rhat = max(rhat),
    bias = mean(bias),
    rmse = mean(rmse),
    .by = c(id, sim_id, param_type)
  ) %>%
  mutate(max_rhat = max(rhat), .by = c(id, sim_id)) %>%
  filter(max_rhat < 1.055)

  p = lapply(1:4, \(x) {
    y = c("Two-factor model", "Mediation model",
      "Interaction model", "Sequential model")[x]
    x = c("sim01", "sim02", "sim03", "sim04")[x]
    (utils.rhat_heatmap(fit_dx %>% filter(id == x)) +
    scale_fill_manual(
      name = TeX(r"($\hat{R}$)"),
      labels = c("1.00", "1.01", TeX(r"($\leq 1.05$)"), "> 1.05"),
      values = rev(c("#0B0405FF", "#40498EFF",
        "#38AAACFF", "#DEF5E5FF"))
    ) + ggtitle(y) |
    utils.recov_box(recov_dx %>% filter(id == x), bias, "Bias") |
    utils.recov_box(recov_dx %>% filter(id == x), rmse, "RMSE")) +
    plot_layout(widths = c(0.5,0.25,0.25))
  })

  p[[1]]/p[[2]]/p[[3]]/p[[4]]
}