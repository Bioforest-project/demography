extract_draws <- function(fit, param) {
  df_draw <- fit$draws(variables = param)
  data.frame(
    variable = rep(paste0(param, "[", seq_len(dim(df_draw)[3]), "]"),
      each = dim(df_draw)[1] * dim(df_draw)[2]
    ),
    value = c(df_draw),
    iter = rep(seq_len(dim(df_draw)[1] * dim(df_draw)[2]), dim(df_draw)[3])
  )
}
