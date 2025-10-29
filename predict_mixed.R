predict_mixed <- function(model, newdata = NULL) {
  
  if(is.null(newdata)) {
    newdata <- model$frame
  }
  
  fixed_formula <- lme4::nobars(model$call$formula)
  random_formula <- lme4::findbars(model$call$formula)[[1]]
  random_vars <- all.vars(random_formula)
  
  grouping <- random_vars[length(random_vars)]
  
  terms_obj <- delete.response(terms(model))
  xlevs <- model.frame(model)[,-1] |>
    select(!(!!grouping))
  xlevs <- lapply(xlevs, levels)
  
  
  beta <- fixef(model)$cond
  # U <- ranef(model)$cond[[1]] |>
  #   as.matrix()
  U <- ranef(model)$cond[[1]] |>
    rownames_to_column(var = "temp") |>
    rename(!!grouping := temp)
  
  z <- newdata |>
    dplyr::mutate("(Intercept)" = 1) |> 
    dplyr::select("(Intercept)", all_of(random_vars)) |>
    as.data.frame()
  
  common_rows_model <- intersect(z[[grouping]], U[[grouping]])
  z <- z |>
    dplyr::filter(.data[[grouping]] %in% common_rows_model)
  U_common <- U |>
    dplyr::filter(.data[[grouping]] %in% common_rows_model)
  
  zu_df <- dplyr::left_join(
    x = z, y = U_common, by = setNames(grouping, grouping)
  ) |>
    dplyr::mutate(random_pred = 0)
  if(nrow(zu_df) > 0) {
    for(i in 1:length(random_vars)) {
      zu_df <- dplyr::mutate(zu_df,
                      random_pred = random_pred + 
                        .data[[names(zu_df)[i]]] * .data[[names(zu_df)[length(random_vars) + 1 + i]]]
      )
    }
  }
  zu_df <- zu_df[,-((length(random_vars) + 2):(2*length(random_vars) + 1))]
  names(zu_df) <- gsub(".x", replacement = "", names(zu_df))
  zu_df <- zu_df |>
    tidyr::unite(join_var, all_of(random_vars), sep = "||&&||", remove = TRUE) |>
    dplyr::select(random_pred, join_var) |>
    unique()
  
  predictions <- model.matrix(terms_obj, data = newdata, xlev = xlevs) %*% beta
  
  random_part <- paste(random_vars[-length(random_vars)], collapse = " + ")
  new_formula <- update(fixed_formula, paste(". ~ . +", random_part))
  
  pred_df <- model.frame(new_formula, data = newdata) |> 
    dplyr::mutate(temp = newdata[[grouping]]) |>
    dplyr::rename(!!grouping := temp) |>
    dplyr::select(all_of(random_vars)) |>
    dplyr::mutate(linear_pred = predictions) |>
    tidyr::unite(join_var, all_of(random_vars), sep = "||&&||", remove = TRUE)
  
  pred_df <- pred_df |>
    dplyr::left_join(zu_df, by = "join_var") |> 
    dplyr::mutate(random_pred = case_when(
      is.na(random_pred) ~ 0,
      .default = random_pred
      )
    ) |>
    dplyr::mutate(
      combined_pred = linear_pred + random_pred,
      pop_prob_prediction = exp(linear_pred)/(1+exp(linear_pred)),
      indiv_prob_prediction = exp(combined_pred)/(1+exp(combined_pred))
    ) |>
    dplyr::select(!join_var)
  
  if(nrow(zu_df) == 0) {
    pred_df <- pred_df |>
      dplyr::select(!indiv_prob_prediction)
  }
  
  pred_list <- list(pred_df |> dplyr::select(last_col()) |> unlist() |> unname(),
                    pred_df)
  
  names(pred_list) <- c(names(pred_df |> dplyr::select(last_col())),
                        "pred_df")
  
  return(pred_list)
}