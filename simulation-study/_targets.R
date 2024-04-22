# Setup ----
library(bdlvm) # available via github: bdlvm-project/bdlvm-pkg
library(brms)
library(cmdstanr)
library(dplyr)
library(forcats)
library(future)
library(ggplot2)
library(here)
library(latex2exp)
library(patchwork)
library(purrr)
library(qs)
library(SBC) # available via github: hyunjimoon/SBC
library(tarchetypes)
library(targets)
library(tidyr)
library(tikzDevice)
plan(multisession)
options(
  SBC.min_chunk_size = 5,
  brms.backend = "cmdstanr"
)
tar_option_set(
  format = "qs", memory = "transient",
  garbage_collection = TRUE
)
tar_source("simulation-study/R")
# Simulation code ----
tar_plan(
  sim_df = simulation_conditions(),
  all_generators = sim_df %>%
    select(
      formulas,
      data_template,
      priors_gen,
      num_obs,
      num_rep,
      fullname
    ),
  tar_target(
    each_generator,
    all_generators,
    pattern = map(all_generators)
  ),
  tar_target(
    datasets,
    utils.make_data(each_generator, na_pattern = "^y[0-9]*$"),
    pattern = map(each_generator)
  ),
  tar_target(
    all_models,
    sim_df %>%
      left_join(datasets) %>%
      unnest(priors_fit) %>%
      select(
        formulas,
        data_template,
        priors_fit,
        data
      )
  ),
  tar_target(
    each_model,
    all_models,
    pattern = map(all_models)
  ),
  tar_target(
    full_results,
    utils.sbc_compute(
      each_model,
      keep_fits = TRUE
    ),
    pattern = map(each_model),
    iteration = "list"
  ),
  metadata_df = sim_df %>%
    left_join(datasets) %>%
    unnest(priors_fit) %>%
    mutate(
      priors_gen = names(priors_gen),
      priors_fit = names(priors_fit),
    ) %>%
    select(
      fullname,
      id,
      num_obs,
      num_rep,
      batch,
      priors_gen,
      priors_fit,
      data
    ),
  tar_target(
    selected_summaries,
    utils.summarize_sbc(metadata_df, full_results),
    pattern = map(metadata_df, full_results)
  ),
  grouped_summaries = utils.regroup(selected_summaries),
  tar_target(
    full_figure,
    utils.make_figures(grouped_summaries),
    pattern = map(grouped_summaries),
    iteration = "list"
  ),
  tar_target(
    figure_pdf, {
      path = here("simulation-study", "plots",
        paste0(attr(full_figure, "plot_id"), ".pdf"))
      ggsave(path, plot = full_figure, width = 7, height = 4,
        device = cairo_pdf)
    },
    pattern = map(full_figure),
    format = "file"
  ),
  ess_figure = utils.ess_plot(grouped_summaries),
  tar_target(
    ess_pdf, {
      path = here("simulation-study", "plots", "ess_figure.pdf")
      ggsave(path, plot = ess_figure, width = 7, height = 7,
        device = cairo_pdf)
    },
    format = "file"
  ),
  tar_target(
    fig_appendix, {
      path = here("simulation-study", "plots", "fig_appendix.pdf")
      ggsave(path, plot = utils.fig_appendix(grouped_summaries),
        width = 7, height = 7, device = cairo_pdf)
    },
    format = "file"
  )
)