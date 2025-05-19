copula_estimation_settings <- list(
  list(copula_model = "vine",
       vine_fit_settings = list(
         family_set = c("onepar","indep"),
         selcrit="mbicv", psi0=0.9,
         threshold=0.0
       ),
       connection_threshold=NULL
  ),
  list(copula_model = "vine",
       vine_fit_settings = list(
         family_set = c("onepar","indep"),
         selcrit="mbicv", psi0=0.9,
         threshold=0.8
       ),
       connection_threshold=NULL
  ),
  list(copula_model = "vine",
       vine_fit_settings = list(
         family_set = c("onepar","indep"),
         selcrit="mbicv", psi0=0.9
       ),
       connection_threshold=0.05),
  list(copula_model = "vine",
       vine_fit_settings = list(
         family_set = c("onepar","indep"),
         selcrit="mbicv", psi0=0.9,
         threshold=0.8
       ),
       connection_threshold=0.05),
  list(copula_model = "vine",
       vine_fit_settings = list(
         family_set = c("onepar","indep"),
         selcrit="aic", psi0=0.9,
         threshold=0.8
       ),
       connection_threshold=0.05)
)
