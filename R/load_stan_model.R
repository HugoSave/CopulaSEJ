# STAN models (bayesian model for fitting the decoupler marginals) needs to be compiled
# to C++ code before being used. This can take some time. The .stan_cache environment can be used to
# to cache an already compiled model.

.stan_cache <- new.env(parent = emptyenv())

load_stan_model <- function(load_precompiled = TRUE, precompiled_model_path = "beta_hiearchical.rds", save_compiled=TRUE, save_compiled_path = "beta_hiearchical.rds",
                           force_compilation=FALSE) {
  model_name = "beta_hierarchical"
  model_file_names <- paste0(model_name, ".stan")
  # Use cache if model already compiled
  if (!exists(model_name, envir = .stan_cache) || force_compilation) {

    # check if exists saved compiled model
    if (load_precompiled && file.exists(precompiled_model_path) && (!force_compilation)) {
      compiled_model <- readRDS(precompiled_model_path)
    } else {
      stan_file <- system.file("stan", model_file_names, package = "CopulaSEJ")
      if (stan_file == "") {
        stop("Stan file not found: ", model_file_names)
      }
      message("Compiling stan model. This might take a few minutes the first time.")
      compiled_model <- rstan::stan_model(file = stan_file)

    # Optionally save the compiled model
      if (save_compiled) {
        if (is.null(save_compiled_path)) {
          save_compiled_path <- file.path(tempdir(), paste0(model_name, "_compiled.rds"))
        }
        message("Saving compiled model to: ", save_compiled_path)
        saveRDS(compiled_model, save_compiled_path)
      }
    }

    assign(model_name, compiled_model, envir = .stan_cache)
  }

  get(model_name, envir = .stan_cache)
}
