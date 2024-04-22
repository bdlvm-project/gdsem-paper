fitting.get_priors = function(f, p) {
  # ULI identification constraint
  p[p$source == "fact_incpt",1] <- "constant(0)"
  p[p$source == "first_load",1] <- "constant(1)"
  # We will not estimate item intercepts
  p[p$source == "item_incpt",1] <- "constant(0)"
  if(f %in% c("sim01", "sim02", "sim03", "sim04"))
    list(
    # Uncomment to test with blavaan default priors where available
    # blavaan = (\(){
    # p[p$source == "fact_sigma",1] <- "gamma(1,0.5)"
    # p[p$source == "fact_dsigma",1] <- "expgamma(1,0.5)"
    # p[p$source == "item_load",1] <- "normal(0,10)"
    # p[p$source == "item_sigma",1] <- "gamma(1,0.5)"
    # p[p$source == "coef_mean",1] <- "normal(0,10)"
    # p[p$source == "coef_sigma",1] <- "normal(0,2.5)"
    # p})(),
    weak = (\(){
    p[p$source == "fact_sigma",1] <- "gamma(5,5)"
    p[p$source == "fact_dsigma",1] <- "expgamma(5,5)"
    p[p$source == "item_load",1] <- "normal(0,2.5)"
    p[p$source == "item_sigma",1] <- "gamma(2.5,5)"
    p[p$source == "coef_mean",1] <- "normal(0,2.5)"
    p[p$source == "coef_sigma",1] <- "normal(0,0.5)"
    p})()
    )
  else
    stop("No match found")
}