simulation_conditions = \() {
  num_obs = c(500)
  num_rep = 5
  num_batch = 50 # multiplies num_rep
  formulas = list(
    sim01 =
      lv(y1 ~ items(x1, 5)) +
      lv(y2 ~ items(x2, 5) + y1, sigma ~ y1),
    sim02 =
      lv(y1 ~ items(x1, 5)) +
      lv(y21 ~ items(x21, 5) + y1, sigma ~ y1) +
      lv(y22 ~ items(x22, 5) + y1, sigma ~ y1) +
      lv(y3 ~ items(x3, 5) + y1 + y21, sigma ~ y1 + y22),
    sim03 =
      lv(y1 ~ items(x1, 5)) +
      lv(y21 ~ items(x21, 5)) +
      lv(y22 ~ items(x22, 5)) +
      lv(y3 ~ items(x3, 5) + y1 + y1:y21, sigma ~ y1 + y1:y22),
    sim04 =
      lv(y1 ~ items(x1, 5)) +
      lv(y2 ~ items(x2, 5) + y1, sigma ~ I(y1^2)) +
      lv(y3 ~ items(x3, 5) + y2, sigma ~ I(y2^2)) +
      lv(y4 ~ items(x4, 5) + y3, sigma ~ I(y3^2)) +
      lv(y5 ~ items(x5, 5) + y4, sigma ~ I(y4^2))
  )
  # Setting up the dataframe ----
  utils.named_tibble(formulas) %>%
    mutate(
      data_template = lapply(formulas, bdlvm_template),
      prior_template = map2(formulas, data_template, bdlvm_model_prior)
    ) %>%
  # Assign generating priors based on id via generator.get_priors()
    mutate(
      priors_gen = map2(id, prior_template, generator.get_priors)
    ) %>%
    unnest(priors_gen) %>%
  # For SBC, use the same priors for fitting (minus nonconstant)
    mutate(
      priors_fit = map2(
        id, prior_template,
        \(id, template) {
        c(
          generator.get_priors(id, template),
          fitting.get_priors(id, template)
        )
      })
    ) %>%
  # Cross over the number of simulations and SBC reps
    utils.cross_vec(num_obs) %>%
    utils.cross_vec(num_rep) %>%
  # Parse the formulas and create calls to make each generator/backend
    mutate(
      # Keep formulas as strings to avoid environment creation
      formulas = map(
        formulas, \(x) {
          paste(
            bdlvm_parse(x, eval = FALSE),
            collapse = "+"
          )
        }
      ),
    ) %>%
  # Batching: repeats reps
    slice(
      rep(1:n(), each = num_batch)
    ) %>%
    mutate(batch = 1:n(), .by = id) %>%
  # Create descriptive names for each condition
    (\(x) {
      fullname <- paste0(
        x$id, "_n", x$num_obs,
        "_GEN.", names(x$priors_gen),
        "_batch.", 1:num_batch
      )
      bind_cols(fullname = fullname, x)
    })
}