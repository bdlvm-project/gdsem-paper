generator.get_priors = function(f, p) {
  # ULI identification constraint
  p[p$source == "fact_incpt",1] <- "constant(0)"
  p[p$source == "first_load",1] <- "constant(1)"
  # We will not estimate item intercepts
  p[p$source == "item_incpt",1] <- "constant(0)"
  if(identical(f, "sim01"))
    list(
    sbc = (\(){
    p[p$source == "fact_sigma",1] <- "gamma(11,11)"
    p[p$source == "fact_dsigma",1] <- "expgamma(11,11)"
    p[p$source == "item_load",1] <- "normal(1,0.3)"
    p[p$source == "item_sigma",1] <- "normal(0.5,0.15)"
    p[p$source == "coef_mean",1] <- "normal(1,0.3)"
    p[p$source == "coef_sigma",1] <- "normal(0.15,0.05)"
    # set lower bounds
    p[p$source == "fact_sigma",8] <- "0.7"
    p[p$source == "item_sigma",8] <- "0.3"
    # item coeffs only take bounds when specified
    # on the generic parameter
    newl <- p[p$source == "item_load",]
    newl[,c(1, 3, 4, 6, 7)] <- ""
    newl[,8] <- "0.7"
    rbind(p, newl)})()
    )
  else if(identical(f, "sim02"))
    list(
    sbc = (\(){
    p[p$source == "fact_sigma",1] <- "gamma(11,11)"
    p[p$source == "fact_dsigma",1] <- "expgamma(11,11)"
    p[p$source == "item_load",1] <- "normal(1,0.3)"
    p[p$source == "item_sigma",1] <- "normal(0.5,0.15)"
    p[p$source == "coef_mean",1] <- "normal(1,0.3)"
    p[p$source == "coef_sigma",1] <- paste0(
    "normal(", c(-1, -1, 1, -1)*0.15, ", 0.05)")
    # set lower bounds
    p[p$source == "fact_sigma",8] <- "0.7"
    p[p$source == "item_sigma",8] <- "0.3"
    newl <- p[p$source == "item_load",]
    newl[,c(1, 3, 4, 6, 7)] <- ""
    newl[,8] <- "0.7"
    rbind(p, newl)})()
    )
  else if(identical(f, "sim03"))
    list(
    sbc = (\(){
    p[p$source == "fact_sigma",1] <- "gamma(11,11)"
    p[p$source == "fact_dsigma",1] <- "expgamma(11,11)"
    p[p$source == "item_load",1] <- "normal(1,0.3)"
    p[p$source == "item_sigma",1] <- "normal(0.5,0.15)"
    p[p$source == "coef_mean",1] <- paste0(
      "normal(", c(1, 0.5), ", 0.3)")
    p[p$source == "coef_sigma",1] <- paste0(
    "normal(", c(0.1, 0.05), ", 0.05)")
    # set lower bounds
    p[p$source == "fact_sigma",8] <- "0.7"
    p[p$source == "item_sigma",8] <- "0.3"
    newl <- p[p$source == "item_load",]
    newl[,c(1, 3, 4, 6, 7)] <- ""
    newl[,8] <- "0.7"
    rbind(p, newl)})()
    )
  else if(identical(f, "sim04"))
    list(
    sbc = (\(){
    p[p$source == "fact_sigma",1] <- "gamma(11,11)"
    p[p$source == "fact_dsigma",1] <- "expgamma(11,11)"
    p[p$source == "item_load",1] <- "normal(1,0.3)"
    p[p$source == "item_sigma",1] <- "normal(0.5,0.15)"
    p[p$source == "coef_mean",1] <- "normal(0,0.2)"
    p[p$source == "coef_sigma",1] <- "normal(0, 0.05)"
    # set lower bounds
    p[p$source == "fact_sigma",8] <- "0.5"
    p[p$source == "item_sigma",8] <- "0.25"
    newl <- p[p$source == "item_load",]
    newl[,c(1, 3, 4, 6, 7)] <- ""
    newl[,8] <- "0.7"
    rbind(p, newl)})()
    )
  else
    stop("No match found")
}