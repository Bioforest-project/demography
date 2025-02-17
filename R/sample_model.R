sample_model <- function(model, mdata) {
  model_path <- file.path("models", paste0(model, ".stan"))
  stan_model <- cmdstan_model(model_path)
  chains_path <- file.path("chains", model)
  unlink(chains_path, recursive = TRUE)
  dir.create(chains_path)
  fit <- stan_model$sample(
    data = mdata,
    chains = 4,
    parallel_chains = 4,
    save_warmup = FALSE,
    output_dir = chains_path
  )
  print(paste(model, "done!"))
}
