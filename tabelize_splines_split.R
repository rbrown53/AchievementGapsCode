tableize_splines_split <- function(model = NULL, model_group = "None", 
                                   mod_variables = NA, ld_value = 25,
                                   dataset = NULL, smaller_table = FALSE) {
  # # resp is the response variable
  # # quant is the quantitative variable (the name should be entered in quotes)
  # # mod_variables is a vector of the names of the x variables wanted in the table. 
  # #  These should be categorical.
  # # model_group is the x variable to split up the tables by
  # # ld_value is the predicted probability of the response variable being true.
  # # latex_tables controls if we want a latex table or an HTML one.
  # # dataset is the data set on which to fit the model. 
  # #   resp, mod_vars, and model_group should be in this data set.
  source('~/My Drive/SUU/Attainment Gaps Fellow/extract_models_glmmTMB_ns.R')
  var_label <- sapply(
    mod_variables,
    function(x) switch(x,
                       "course_college" = "Course College",
                       "student_college" = "Student College",
                       "sex" = "Sex",
                       "gpa_group" = "HS GPA Group",
                       "degree" = "Degree",
                       "student_type" = "Student Type",
                       "course_level" = "Course Level",
                       "hsgpa" = "High School GPA",
                       "sem_gpa" = "Semester GPA",
                       "overall_gpa" = "Overall GPA",
                       "stu_level" = "Student Level (Grad/Undergrad)",
                       "student_class" = "Student Class",
                       "age" = "Age",
                       "marital_status" = "Marital Status",
                       "first_gen" = "First Gen Status",
                       "race" = "Race/Ethnicity",
                       "percent_need_met" = "Percent of Financial Need Met",
                       "term" = "Academic Term",
                       "persist1" = "persistence 1",
                       "persist2" = "persistence 2",
                       "covid" = "COVID",
                       "dfw" = "DFW",
                       "df" = "D or F",
                       "term" = "Term",
                       "course_subject" = "Course Subject",
                       "course_department" = "Course Department",
                       "student_department" = "Student Department"
    )
  )
  for(i in 1:length(var_label)) {
    if(is.null(var_label[[i]])) {
      var_label[[i]] = ""
    }
  }
  model_group_label <- sapply(
    model_group,
    function(x) switch(x,
                       "None" = "None",
                       "course_college" = "Course College",
                       "student_college" = "Student College",
                       "sex" = "Sex",
                       "gpa_group" = "HS GPA Group",
                       "degree" = "Degree",
                       "student_type" = "Student Type",
                       "course_level" = "Course Level",
                       "hsgpa" = "High School GPA",
                       "sem_gpa" = "Semester GPA",
                       "overall_gpa" = "Overall GPA",
                       "stu_level" = "Student Level (Grad/Undergrad)",
                       "student_class" = "Student Class",
                       "age" = "Age",
                       "marital_status" = "Marital Status",
                       "first_gen" = "First Gen Status",
                       "race" = "Race/Ethnicity",
                       "percent_need_met" = "Percent of Financial Need Met",
                       "term" = "Academic Term",
                       "persist1" = "persistence 1",
                       "persist2" = "persistence 2",
                       "covid" = "COVID",
                       "dfw" = "DFW",
                       "df" = "D or F",
                       "term" = "Term",
                       "course_subject" = "Course Subject",
                       "course_department" = "Course Department",
                       "student_department" = "Student Department"
    )
  )
  resp_label = "DFW"
  if(model_group != "None") {
    tables <- list()
    i <- 0
    for(level in levels(factor(dataset[[model_group]]))) {
      i <- i + 1
      data_filtered = dataset[dataset[[model_group]] == level,] |>
        droplevels()
      # formulas_vars <-  all.vars(formula(model$call))
      # formulas_vars <- formulas_vars[! formulas_vars %in% model_group]
      model_formula <- formula(paste0("dfw", " ~ ", mod_variables[1],
                                      " * ", mod_variables[2],
                                      " * splines::ns(hsgpa, df = 4)",
                                      " + (1 + hsgpa | random_id )"))
      data_mod <- try(
        glmmTMB::glmmTMB(
          model_formula,
          data = data_filtered, family = binomial(),
          # sparseX = c(cond = TRUE),
          control = glmmTMB::glmmTMBControl(
            optCtrl = list(iter.max = 2000, eval.max = 2100)
          )
        ),
        silent = TRUE
      )
      
      if(class(data_mod)[1] == "try-error"){
        tables[[i]] <- NULL
      } else {
        if(ld_value > 0) {
          mod_table <- extract_models_glmmTMB_ns(data_mod, ld_value,
                                                 include_model = FALSE)
          mod_table[,1] <- gsub(pattern = "<=", replacement = "&le;", mod_table[,1])
          mod_table[,1] <- gsub(pattern = ">=", replacement = "&ge;", mod_table[,1])
          # if(nrow(unique(data_filtered[mod_vars[!(mod_vars %in% quant)]])) == 1) {
          #   mod_table[1] = unique(data_filtered[mod_vars[!(mod_vars %in% quant)]])
          # }
          table.col.names <- c(
            "Group",
            paste0("LD", ld_value, "\nValue"),
            "Number Records\nin Category",
            paste0("Number Records\nBelow LD", ld_value, " Value"),
            paste0("Percent Records\nBelow LD", ld_value, " Value"),
            "Number Individuals\nin Category",
            "Number Individuals\nBelow LD",
            "Percent Individuals\nBelow LD"
          )
          mod_table <- mod_table %>%
            mutate(
              `LD Value` = kableExtra::cell_spec(
                `LD Value`, color = case_when(
                  `LD Value` > 3.5 & `LD Value` != "None" ~ "darkred",
                  `LD Value` > 3.0 & `LD Value` != "None" ~ "red",
                  `LD Value` < 1.0 | `LD Value` == "None" ~ "blue",
                  .default = "black"
                ),
                bold = case_when(
                  `LD Value` > 3.5 & `LD Value` != "None" ~ T,
                  `LD Value` > 3.0 & `LD Value` != "None" ~ T,
                  `LD Value` < 1.0 | `LD Value` == "None" ~ T,
                  .default = F)
              ),
              `Percent Records Below LD` = kableExtra::cell_spec(
                `Percent Records Below LD`, color = case_when(
                  `Percent Records Below LD` > 50 ~ "darkred",
                  `Percent Records Below LD` > 20 ~ "red",
                  `Percent Records Below LD` < 5 ~ "blue",
                  .default = "black"),
                bold = case_when(`Percent Records Below LD` > 50 ~ T,
                                 `Percent Records Below LD` > 20 ~ T,
                                 `Percent Records Below LD` < 5 ~ T,
                                 .default = F)
              ),
              `Percent Indivds Below LD` = kableExtra::cell_spec(
                `Percent Indivds Below LD`, color = case_when(
                  `Percent Indivds Below LD` > 50 ~ "darkred",
                  `Percent Indivds Below LD` > 20 ~ "red",
                  `Percent Indivds Below LD` < 5 ~ "blue",
                  .default = "black"),
                bold = case_when(`Percent Indivds Below LD` > 50 ~ T,
                                 `Percent Indivds Below LD` > 20 ~ T,
                                 `Percent Indivds Below LD` < 5 ~ T,
                                 .default = F)
              )
            )
        } else {
          mod_table <- extract_models_glmmTMB_ns(data_mod,
                                                 include_model = FALSE)
          mod_table[,1] <- gsub(pattern = "<=", replacement = "&le;", mod_table[,1])
          mod_table[,1] <- gsub(pattern = ">=", replacement = "&ge;", mod_table[,1])
          table.col.names <- c(paste0(var_label, " Group"), 
                               "Fixed Effects\nLogit Model")
        }
        column_width1 <- max(
          strwidth(mod_table[,1], font = 12, units = "in")
        ) - 0.4
        column_width2 <- max(
          strwidth(mod_table[,2], font = 12, units = "in")
        ) - 0.1
        if(column_width1 > 3.6) {
          column_width1 <- column_width1 / 2
        }
        table_caption <- paste0(
          "<center><strong><span style='color:black'>",
          "Table for the \"", 
          level, 
          "\" Category of the \"",
          model_group_label,
          "\" Variable when Predicting ", resp_label,
          "</center></strong></span>"
        )
        
        output_table <- mod_table |>
          knitr::kable(digits = 3,
                       format = "html",
                       caption = table_caption,
                       col.names = table.col.names,
                       align = "lcccccccc",
                       escape = F,
                       linesep = ""
          ) |>
          kableExtra::column_spec(1, width = paste0(column_width1, "in")) |>
          kableExtra::column_spec(2, width = paste0(column_width2, "in")) |>
          kableExtra::kable_styling(
            bootstrap_options = c("striped", "hover", "bordered", "condensed"),
            full_width = T,
            font_size = 12
          )
        tables[[i]] <- output_table
      }
    }
  } else if (model_group == "None") {
    if (class(model)[1] == "try-error") {
      tables <- NULL
    } else {
      if(ld_value > 0) {
        mod_table <- extract_models_glmmTMB_ns(model, ld_value,
                                               include_model = FALSE)
        mod_table[,1] <- gsub(pattern = "<=", replacement = "&le;", mod_table[,1])
        mod_table[,1] <- gsub(pattern = ">=", replacement = "&ge;", mod_table[,1])
        column_width1 <- max(
          strwidth(mod_table[,1], font = 12, units = "in")
        ) - 0.4
        column_width2 <- max(
          strwidth(mod_table[,2], font = 12, units = "in")
        ) 
        if(column_width1 > 3.6) {
          column_width1 <- column_width1 / 2
        }
        table.col.names <- c(
          paste0("Group"),
          paste0("LD", ld_value, "\nValue"),
          "Number Records\nin Category",
          paste0("Number Records\nBelow LD", ld_value, " Value"),
          paste0("Percent Records\nBelow LD", ld_value, " Value"),
          "Number Individuals\nin Category",
          "Number Individuals\nBelow LD",
          "Percent Individuals\nBelow LD"
        )
        mod_table <- mod_table %>%
          mutate(
            `LD Value` = kableExtra::cell_spec(
              `LD Value`, color = case_when(
                `LD Value` > 3.5 & `LD Value` != "None" ~ "darkred",
                `LD Value` > 3.0 & `LD Value` != "None" ~ "red",
                `LD Value` < 1.0 | `LD Value` == "None" ~ "blue",
                .default = "black"
              ),
              bold = case_when(
                `LD Value` > 3.5 & `LD Value` != "None" ~ T,
                `LD Value` > 3.0 & `LD Value` != "None" ~ T,
                `LD Value` < 1.0 | `LD Value` == "None" ~ T,
                .default = F)
            ),
            `Percent Records Below LD` = kableExtra::cell_spec(
              `Percent Records Below LD`, color = case_when(
                `Percent Records Below LD` > 50 ~ "darkred",
                `Percent Records Below LD` > 20 ~ "red",
                `Percent Records Below LD` < 5 ~ "blue",
                .default = "black"),
              bold = case_when(`Percent Records Below LD` > 50 ~ T,
                               `Percent Records Below LD` > 20 ~ T,
                               `Percent Records Below LD` < 5 ~ T,
                               .default = F)
            ),
            `Percent Indivds Below LD` = kableExtra::cell_spec(
              `Percent Indivds Below LD`, color = case_when(
                `Percent Indivds Below LD` > 50 ~ "darkred",
                `Percent Indivds Below LD` > 20 ~ "red",
                `Percent Indivds Below LD` < 5 ~ "blue",
                .default = "black"),
              bold = case_when(`Percent Indivds Below LD` > 50 ~ T,
                               `Percent Indivds Below LD` > 20 ~ T,
                               `Percent Indivds Below LD` < 5 ~ T,
                               .default = F)
            )
          )
      } else {
        mod_table <- extract_models_glmmTMB_ns(model, include_model = FALSE)
        table.col.names <- c("Group")
      }
      table_caption <- paste0(
        "<center><strong><span style='color:black'>",
        "Table for the Entire University when Predicting ", resp_label,
        "</center></strong></span>"
      )
      output_table <- mod_table |>
        knitr::kable(digits = 3,
                     format = "html",
                     caption = table_caption,
                     col.names = table.col.names,
                     align = "lcccccccc",
                     escape = F,
                     linesep = ""
        ) |>
        kableExtra::column_spec(1, width = paste0(column_width1, "in")) |>
        kableExtra::column_spec(2, width = paste0(column_width2, "in")) |>
        kableExtra::kable_styling(
          bootstrap_options = c("striped", "hover", "bordered", "condensed"),
          full_width = F,
          font_size = 12
        )
      tables <- output_table
    }
  }
  
  if (exists("tables")) {
    return(tables)
  } else {
    return(NULL)
  }
}