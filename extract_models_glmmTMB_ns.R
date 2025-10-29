extract_models_glmmTMB_ns <- function(model, LD = NULL, p.values = FALSE,
                                      include_model = TRUE) {
  library(lme4) |> suppressPackageStartupMessages()
  # This function is used to obtain a separate model for each combination of
  #   the categorical variables used to fit the model when there is ONE
  #   quantitative variable.
  # LD is the lethal dose value. If provided, this function will also give
  #   information about the number of observations in each category, the
  #   number below that LD value, and the percentage below that LD value.
  # include_model: whether to output the logit model or not
  
  
  # Check LD
  if (length(LD) > 0) {
    if (LD > 1 & LD <= 100 & length(LD) == 1) {
      LD <- LD / 100
    } else if (LD < 0 | LD > 100 | length(LD) != 1) {
      stop("Invalid LD value. Should be between 0 and 100.")
    }
  }
  
  # Save the variance-covariance matrix if the user wants p-values
  if(p.values == TRUE) {
    if(class(model)[1] == "speedglm") {
      Sigma <- solve(model$XTX) |>
        as.data.frame()
    } else if(class(model)[1] == "glm") {
      Sigma <- solve(t(model$R)%*%model$R) |>
        as.data.frame() # Same as vcov(model)
    } else if(class(model)[1] == "glmerMod") {
      Sigma <- vcov(model) |>
        as.matrix() |> 
        as.data.frame()
    } else if(class(model)[1] == "glmmTMB") {
      Sigma <- vcov(model)$cond |>
        as.data.frame()
    }
  }
  
  # Extract the data from the model
  mod_frame <- model.frame(model)
  mod_vars <- colnames(mod_frame)[-1] # Remove the first element
  if(class(model)[1] == "glmerMod") {
    rand_var <- names(model1@cnms)
    mod_vars <- mod_vars[!(mod_vars %in% rand_var)]
  } else if (class(model)[1] == "glmmTMB") {
    rand_var <- names(model$modelInfo$reTrms$cond$cnms)
    mod_vars <- mod_vars[!(mod_vars %in% rand_var)]
  }
  var_type <- sapply(mod_frame[mod_vars], is.numeric)
  # Extract categorical variables from the model formula
  numeric_var <- mod_vars[var_type == TRUE]
  if(substring(numeric_var[1], 1, 3) == "ns(" |
     substring(numeric_var[1], 1, 12) == "splines::ns(") {
    terms <- names(fixef(model)$cond)
    numeric_place <- grep(paste0("ns(", numeric_var[2]), 
                          terms, fixed = T)
    interaction_place <- which(grepl(":", terms) & 
                                 !grepl("^[^:]*::[^:]*$", terms))
    correct_place <- numeric_place[! (numeric_place %in% 
                    intersect(numeric_place, interaction_place))]
    spline_vars <- terms[correct_place]
    numeric_var <- numeric_var[2]
  }
  
  #categorical_columns <- names(model$xlevels)
  categorical_columns <- mod_vars[var_type == FALSE]
  
  # If there is no quantitative variable, return an informative message
  if(sum(var_type == TRUE) == 0) {
    stop("No quantitative variables found in the model.")
  }
  
  # If there are no categorical variables, return a model for all together
  if (length(categorical_columns) == 0) { ####### UPDATE THIS
    if(class(model)[1] == "glmerMod") {
      intercept <- lme4::fixef(model)[1]
      slope <- lme4::fixef(model)[2]
    } else if (class(model)[1] == "glmmTMB") {
      intercept <- glmmTMB::fixef(model)$cond[1]
      slope <- glmmTMB::fixef(model)$cond[2]
    } else {
      intercept <- coef(model)[1]
      slope <- coef(model)[2]
    }
    if (!is.null(LD)) {
      LD_val <- as.numeric((log(LD / (1 - LD)) - intercept) / slope)
      
      # Get summary information related to the LD value
      n <- nrow(mod_frame)
      x <- length(which(mod_frame[[numeric_var]] <= LD_val))
      n_inds <- n_distinct(mod_frame[[rand_var]])
      x_inds <- n_distinct(mod_frame[[rand_var]] |>
                             _[which(mod_frame[[numeric_var]] <= LD_val)])
      risk <- data.frame(
        "Total Records" = n,
        "Records Below LD" = x,
        "Percent Records Below LD" = round(x / n * 100, 2),
        "Total Individuals" = n_inds,
        "Indivds Below LD" = x_inds,
        "Percent Indivds Below LD" = round(x_inds / n_inds * 100, 2)
      )
      
      models_df <- data.frame(
        Group = "Everyone Combined",
        Model = ifelse(slope >= 0,
                       paste0(
                         round(intercept, 3), " + ", round(slope, 3),
                         "(", numeric_var, ")"
                       ),
                       paste0(
                         round(intercept, 3), " - ", round(abs(slope), 3),
                         "(", numeric_var, ")"
                       )
        ),
        "LD Value" = round(LD_val, 2),
        risk
      )
    } else {
      models_df <- data.frame(
        Group = "Everyone Combined",
        Model = ifelse(slope >= 0,
                       paste0(
                         round(intercept, 3), " + ", round(slope, 3),
                         "(", numeric_var, ")"
                       ),
                       paste0(
                         round(intercept, 3), " - ", round(abs(slope), 3),
                         "(", numeric_var, ")"
                       )
        )
        # Intercept = intercept,
        # Slope = slope
      )
    }
    rownames(models_df) <- NULL
    names(models_df) <- gsub("\\.", " ", names(models_df))
    #stop("No categorical variables found in the model.")
    return(models_df)
  }
  
  # Convert character variables to factors
  for (var in categorical_columns) {
    if (is.character(mod_frame[[var]])) {
      mod_frame[[var]] <- as.factor(mod_frame[[var]])
    }
  }
  
  # Create a data frame to store the models
  models_df <- data.frame()
  
  # Get all combinations of levels of categorical variables
  combinations <- expand.grid(lapply(mod_frame[categorical_columns], levels))
  
  # Extract coefficients from the original model
  if(class(model)[1] == "glmerMod") {
    original_coeffs <- lme4::fixef(model)
  } else if (class(model)[1] == "glmmTMB") {
    original_coeffs <- glmmTMB::fixef(model)$cond
  } else {
    original_coeffs <- coef(model)
  }
  
  # Change the numeric variable from the left or middle to the right
  for(numeric_i in 1:length(spline_vars)) {
    coef_names <- names(original_coeffs)
    left_numeric <- grep(paste0(spline_vars[numeric_i], ":"), 
                         coef_names, fixed = T)
    coef_names[left_numeric] <- paste0(spline_vars[numeric_i], ":") |>
      gsub("", coef_names[left_numeric]) |>
      paste0(":", spline_vars[numeric_i])
    names(original_coeffs) <- coef_names
    if(p.values == TRUE) {
      p_vals <- matrix(rep(0, nrow(combinations) * 4),
                       nrow = nrow(combinations))
    }
  }
  
  for (i in 1:nrow(combinations)) {
    slope <- c()
    for(numeric_i in 1:length(spline_vars)) {
      # Initialize the intercept, slope, and term list
      intercept <- original_coeffs["(Intercept)"]
      slope[numeric_i] <- original_coeffs[spline_vars[numeric_i]]
      term <- list()
      
      # Update the intercept and slope for each categorical variable with no
      #   interactions
      for (j in 1:length(categorical_columns)) {
        variable <- categorical_columns[j]
        level <- as.character(combinations[i, j])
        term[[j]] <- paste(variable, level, sep = "")
        
        intercept_add <- original_coeffs[term[[j]]]
        slope_add <- original_coeffs[paste(term[[j]], 
                                           spline_vars[numeric_i], sep = ":")]
        if (!is.na(intercept_add)) {
          intercept <- intercept + intercept_add
        }
        if (!is.na(slope_add)) {
          slope[numeric_i] <- slope[numeric_i] + slope_add
        }
      }
      
      # Update the intercept and slope for each categorical variable from
      #   the interaction terms
      if (length(term) >= 2) {
        term_ints <- sapply(
          2:length(term),
          function(x) {
            apply(combn(term, x), 2, paste, collapse = ":")
          }
        ) |>
          unlist()
        for (l in 1:length(term_ints)) {
          intercept_add <- original_coeffs[term_ints[l]]
          slope_add <- original_coeffs[paste(term_ints[l],
                                             spline_vars[numeric_i],
                                             sep = ":"
          )]
          if (!is.na(intercept_add)) {
            intercept <- intercept + intercept_add
          }
          if (!is.na(slope_add)) {
            slope[numeric_i] <- slope[numeric_i] + slope_add
          }
        }
      }
      if(p.values == TRUE) {
        p_vals[i,] = extract_pvalues(term, spline_vars[numeric_i], 
                                     original_coeffs, Sigma)
      }
    }
      
      # Add the information to the models_df
      if (!is.null(LD)) {
        # LD_val <- as.numeric((log(LD / (1 - LD)) - intercept) / slope)
        spline_basis <- ns(mod_frame[[numeric_var]], df = length(spline_vars))
        # LD_val <- optimize(
        #   f = function(x) {
        #     abs(sum(slope * predict(spline_basis, newx = x)) - 
        #           log(LD / (1 - LD)) + intercept)
        #   },
        #   interval = c(0, 4)
        # )$minimum
        obj_func <- function(x){
          intercept + sum(slope * predict(spline_basis, newx = x)) - qlogis(LD)
          }
        # LD_val = tryCatch(max(
        #   rootSolve::uniroot.all(
        #     Vectorize(obj_func),
        #     interval = c(0, 4)
        #     )
        #   ),
        #   warning = function(w) {
        #     message("No root, opimizing instead.")
        #     optimize(\(x) abs(obj_func(x)), c(0, 4))$minimum
        #     }
        #   )
        LD_val = max(c(0,
          rootSolve::uniroot.all(
            Vectorize(obj_func),
            interval = c(0, 4)
          ))
        )
        
        
      # Reduce the data frame to the group combo we are working with 
      reduced_data <- mod_frame
      for (k in 1:length(categorical_columns)) {
        reduced_data <- reduced_data[which(
          reduced_data[[categorical_columns[k]]] ==
            as.character(combinations[i, k])
        ), ]
      }
      
      # Get summary information related to the LD value
      n <- nrow(reduced_data)
      x <- length(which(reduced_data[[numeric_var]] <= LD_val))
      n_inds <- n_distinct(reduced_data[[rand_var]])
      x_inds <- n_distinct(reduced_data[[rand_var]] |>
                             _[which(reduced_data[[numeric_var]] <= LD_val)])
      risk <- data.frame(
        "Total Records" = n,
        "Records Below LD" = x,
        "Percent Records Below LD" = round(x / n * 100, 2),
        "Total Individuals" = n_inds,
        "Indivds Below LD" = x_inds,
        "Percent Indivds Below LD" = round(x_inds / n_inds * 100, 2)
      )
      
      model_info <- data.frame(
        Group = paste(t(combinations[i, ]), collapse = " AND "),
        Model = paste0(round(intercept, 3), " + ", 
                       paste0(round(slope,3), "*", 
                              spline_vars, collapse = " + ")),
        "LD Value" = LD_val,
        risk
      )
    } else {
      model_info <- data.frame(
        Group = paste(t(combinations[i, ]), collapse = " AND "),
        # Model = ifelse(slope >= 0,
        #                paste0(
        #                  round(intercept, 3), " + ", round(slope, 3),
        #                  "(", numeric_var, ")"
        #                ),
        #                paste0(
        #                  round(intercept, 3), " - ", round(abs(slope), 3),
        #                  "(", numeric_var, ")"
        #                )
        # )
        Model = paste0(round(intercept, 3), " + ", 
                       paste0(round(slope,3), "*", 
                              spline_vars, collapse = " + "))
        # Intercept = intercept,
        # Slope = slope
      )
    }
    
    models_df <- rbind(models_df, model_info)
  }
  
  
  names(models_df) <- gsub("\\.", " ", names(models_df))
  
  if (!is.null(LD)) {
    models_df$`LD Value` <- ifelse(models_df$`LD Value` <= 0,
                                   "None",
                                   as.character(round(models_df$`LD Value`, 2))
    )
  }
  
  rownames(models_df) <- NULL
  
  if(p.values & is.null(LD)) {
    models_df <- data.frame(models_df, as.data.frame(p_vals))
    names(models_df) <- c("Group", "Model", 
                          "p-value for Intercept (Compare)",
                          "p-value for Slope (Compare)",
                          "p-value for Intercept",
                          "p-value for Slope")
  } else if(p.values & !is.null(LD)) {
    models_df <- data.frame(models_df, as.data.frame(p_vals))
    names(models_df) <- c("Group", "Model", 
                          "LD Value",
                          "Total Records",
                          "Records Below LD",
                          "Percent Records Below LD",
                          "Total Individuals",
                          "Indivds Below LD",
                          "Percent Indivds Below LD",
                          "p-value for Intercept (Compare)",
                          "p-value for Slope (Compare)",
                          "p-value for Intercept",
                          "p-value for Slope")
  }
  
  if(!include_model) {
    models_df <- models_df[,-2]
  }
  
  return(models_df)
}



