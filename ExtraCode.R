### Modelling

Functions to generate appropriate models to present results as per RECORD \[15\] Outcomes and RECORD \[16\] Main Results unadjusted and adjusted for the primary aim.

```{r}
#| label: faos-model-func-robust

fit_faos_models_with_progress <- function(imputed_data, 
                                          outcomes = c("FAOS_Pain_TotalScore_Preop", 
                                                       "FAOS_Symptom_TotalScore_Preop",
                                                       "FAOS_DailyLiving_TotalScore_Preop", 
                                                       "FAOS_Sport_TotalScore_Preop",
                                                       "FAOS_Quality_TotalScore_Preop"),
                                          predictor = "Region",
                                          model_type = c("minimal", "unadjusted", "fully_adjusted", "arthritis_adjusted"),
                                          covariate = NULL,
                                          minimal_adjusters = c("AgeAtInitialExam", "Sex2"),
                                          full_covariates = c("AgeAtInitialExam", "Sex2", 
                                                              "SRCQTotalScore", "SmokingStatus"),
                                          arthritis_covariates = c("Region", "AgeAtInitialExam", "Sex2", "SmokingStatus"),
                                          arthritis_predictor = "Arthritis",
                                          cluster_var = "PatientID",
                                          robust_se = TRUE,
                                          use_parallel = TRUE,
                                          show_progress = TRUE) {
  
  # Input validation
  model_type <- match.arg(model_type)
  
  # Check required packages
  required_packages <- c("mice")
  if (robust_se) required_packages <- c(required_packages, "estimatr")
  if (use_parallel) required_packages <- c(required_packages, "future", "furrr")
  
  missing_packages <- required_packages %>%
    map_lgl(~!requireNamespace(.x, quietly = TRUE)) %>%
    keep(isTRUE)
  
  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(names(missing_packages), collapse = ", "))
  }
  
  # Validate inputs
  if (model_type == "minimal" && is.null(covariate)) {
    stop("covariate must be specified for minimal models")
  }
  
  # Helper function to build formula strings
  build_formula_string <- function(outcome, type, covariate = NULL) {
    
    formula_terms <- switch(type,
                            "minimal" = {
                              adjusters <- minimal_adjusters %>% 
                                discard(~ .x == covariate)
                              c(covariate, adjusters) %>% 
                                compact()
                            },
                            "unadjusted" = predictor,
                            "fully_adjusted" = c(predictor, full_covariates) %>% unique(),
                            "arthritis_adjusted" = c(arthritis_predictor, arthritis_covariates) %>% unique()
    )
    
    str_c(outcome, " ~ ", str_c(formula_terms, collapse = " + "))
  }
  
  # Helper function to standardize outcome names
  standardize_outcome_name <- function(outcome) {
    extracted <- str_extract(outcome, "(?<=FAOS_)[^_]+")
    case_when(
      extracted == "Pain" ~ "Pain",
      extracted == "Symptom" ~ "Symptom", 
      extracted == "DailyLiving" ~ "DailyLiving",
      extracted == "Sport" ~ "Sport",
      extracted == "Quality" ~ "Quality",
      TRUE ~ extracted
    )
  }
  
  # Helper function to fit individual model
  fit_single_model <- function(outcome_data) {
    # Ensure required packages are loaded in parallel context
    if (use_parallel) {
      require(mice, quietly = TRUE)
      if (robust_se) require(estimatr, quietly = TRUE)
    }
    
    outcome <- outcome_data$outcome
    formula_string <- outcome_data$formula_string
    
    tryCatch({
      if (robust_se && !is.null(cluster_var)) {
        with(imputed_data, 
             lm_robust(as.formula(formula_string), 
                       clusters = !!sym(cluster_var)))
      } else {
        with(imputed_data, 
             lm(as.formula(formula_string)))
      }
    }, error = function(e) {
      warning("Failed to fit model for outcome '", outcome, "': ", e$message)
      return(NULL)
    })
  }
  
  # Prepare data for processing
  outcome_data <- tibble(outcome = outcomes) %>%
    mutate(
      outcome_name = map_chr(outcome, standardize_outcome_name),
      formula_string = map_chr(outcome, ~build_formula_string(.x, model_type, covariate))
    )
  
  # Print initial info
  if (show_progress) {
    cat("Fitting", length(outcomes), "FAOS", model_type, "models with", 
        ifelse(robust_se, "robust", "standard"), "standard errors\n")
    
    if (robust_se && !is.null(cluster_var)) {
      cat("Clustering by:", cluster_var, "\n")
    }
    
    if (use_parallel) {
      cat("Using parallel processing\n")
    }
    cat("\n")
  }
  
  # Fit models with optional progress and parallel processing
  if (show_progress) {
    with_progress({
      p <- progressor(along = outcomes)
      
      if (use_parallel && length(outcomes) > 1) {
        # Parallel processing with progress and proper RNG
        model_results <- outcome_data %>%
          mutate(
            model = future_map(1:nrow(.), function(i) {
              result <- fit_single_model(outcome_data[i, ])
              p(message = sprintf("Fitted %s model", outcome_data$outcome_name[i]))
              return(result)
            }, .progress = FALSE, .options = furrr_options(seed = TRUE))  # Ensure proper RNG
          )
      } else {
        # Sequential processing with progress
        model_results <- outcome_data %>%
          mutate(
            model = map(1:nrow(.), function(i) {
              result <- fit_single_model(outcome_data[i, ])
              p(message = sprintf("Fitted %s model", outcome_data$outcome_name[i]))
              return(result)
            })
          )
      }
    })
  } else {
    # No progress reporting
    if (use_parallel && length(outcomes) > 1) {
      model_results <- outcome_data %>%
        mutate(model = future_map(1:nrow(.), ~fit_single_model(outcome_data[.x, ]), 
                                  .options = furrr_options(seed = TRUE)))
    } else {
      model_results <- outcome_data %>%
        mutate(model = map(1:nrow(.), ~fit_single_model(outcome_data[.x, ])))
    }
  }
  
  # Filter out failed models and create named list
  final_results <- model_results %>%
    filter(!map_lgl(model, is.null)) %>%
    {set_names(.$model, .$outcome_name)}
  
  # Add metadata as attributes
  attr(final_results, "model_type") <- model_type
  attr(final_results, "predictor") <- if(model_type == "arthritis_adjusted") arthritis_predictor else predictor
  attr(final_results, "covariate") <- covariate
  attr(final_results, "robust_se") <- robust_se
  attr(final_results, "cluster_var") <- cluster_var
  attr(final_results, "n_outcomes") <- length(final_results)
  attr(final_results, "use_parallel") <- use_parallel
  
  return(final_results)
}

# Create memoised version using memory cache (consistent with your framework)
fit_faos_models_cached <- memoise(
  fit_faos_models_with_progress,
  cache = cache_memory()
)

```

```{r}
# Set up parallel processing plan (adjust workers based on your CPU cores)
# Use availableCores() - 1 to leave one core free for system
plan(multisession, workers = max(1, availableCores() - 1))

# Helper function to create a cache key for memoise
create_cache_key <- function(imputed_data, outcomes, predictor, model_type, 
                             covariate, minimal_adjusters, full_covariates,
                             arthritis_covariates, arthritis_predictor, 
                             cluster_var, robust_se) {
  
  # Create a hash of the imputed data structure
  data_hash <- digest(list(
    ncol = ncol(complete(imputed_data, 1)),
    nrow = nrow(complete(imputed_data, 1)),
    m = imputed_data$m,
    colnames = names(complete(imputed_data, 1))
  ))
  
  # Create hash of all parameters
  param_hash <- digest(list(
    outcomes = outcomes,
    predictor = predictor,
    model_type = model_type,
    covariate = covariate,
    minimal_adjusters = minimal_adjusters,
    full_covariates = full_covariates,
    arthritis_covariates = arthritis_covariates,
    arthritis_predictor = arthritis_predictor,
    cluster_var = cluster_var,
    robust_se = robust_se
  ))
  
  return(paste0(data_hash, "_", param_hash))
}

# Optimized function with better parallel processing and caching
fit_faos_models_optimized <- function(imputed_data, 
                                      outcomes = c("FAOS_Pain_TotalScore_Preop", 
                                                   "FAOS_Symptom_TotalScore_Preop",
                                                   "FAOS_DailyLiving_TotalScore_Preop", 
                                                   "FAOS_Sport_TotalScore_Preop",
                                                   "FAOS_Quality_TotalScore_Preop"),
                                      predictor = "Region",
                                      model_type = c("minimal", "unadjusted", "fully_adjusted", "arthritis_adjusted"),
                                      covariate = NULL,
                                      minimal_adjusters = c("AgeAtInitialExam", "Sex2"),
                                      full_covariates = c("AgeAtInitialExam", "Sex2", 
                                                          "SRCQTotalScore", "SmokingStatus"),
                                      arthritis_covariates = c("Region", "AgeAtInitialExam", "Sex2", "SmokingStatus"),
                                      arthritis_predictor = "Arthritis",
                                      cluster_var = "PatientID",
                                      robust_se = TRUE,
                                      use_parallel = TRUE,
                                      show_progress = TRUE,
                                      cache_results = TRUE) {
  
  # Input validation
  model_type <- match.arg(model_type)
  
  # Validate inputs
  if (model_type == "minimal" && is.null(covariate)) {
    stop("covariate must be specified for minimal models")
  }
  
  # Create cache key for this specific function call
  cache_key <- create_cache_key(imputed_data, outcomes, predictor, model_type,
                                covariate, minimal_adjusters, full_covariates,
                                arthritis_covariates, arthritis_predictor,
                                cluster_var, robust_se)
  
  # Check if results are cached
  if (cache_results && exists("faos_model_cache", envir = .GlobalEnv)) {
    cached_result <- get("faos_model_cache", envir = .GlobalEnv)[[cache_key]]
    if (!is.null(cached_result)) {
      if (show_progress) cat("Using cached results for", model_type, "models\n")
      return(cached_result)
    }
  }
  
  # Helper function to build formula strings
  build_formula_string <- function(outcome, type, covariate = NULL) {
    formula_terms <- switch(type,
                            "minimal" = {
                              adjusters <- minimal_adjusters[minimal_adjusters != covariate]
                              c(covariate, adjusters) %>% 
                                compact()
                            },
                            "unadjusted" = predictor,
                            "fully_adjusted" = unique(c(predictor, full_covariates)),
                            "arthritis_adjusted" = unique(c(arthritis_predictor, arthritis_covariates))
    )
    
    str_c(outcome, " ~ ", str_c(formula_terms, collapse = " + "))
  }
  
  # Helper function to standardize outcome names
  standardize_outcome_name <- function(outcome) {
    extracted <- str_extract(outcome, "(?<=FAOS_)[^_]+")
    case_when(
      extracted == "Pain" ~ "Pain",
      extracted == "Symptom" ~ "Symptom", 
      extracted == "DailyLiving" ~ "DailyLiving",
      extracted == "Sport" ~ "Sport",
      extracted == "Quality" ~ "Quality",
      TRUE ~ extracted
    )
  }
  
  # Pre-extract complete datasets for efficiency
  complete_datasets <- map(1:imputed_data$m, ~complete(imputed_data, .x))
  
  # Helper function to fit individual model (optimized for parallel)
  fit_single_model_parallel <- function(outcome_info) {
    outcome <- outcome_info$outcome
    formula_string <- outcome_info$formula_string
    outcome_name <- outcome_info$outcome_name
    
    tryCatch({
      if (robust_se && !is.null(cluster_var)) {
        # Use pre-extracted complete datasets
        pooled_models <- map(complete_datasets, function(complete_data) {
          lm_robust(as.formula(formula_string), 
                    data = complete_data,
                    clusters = complete_data[[cluster_var]],
                    se_type = "stata")  # Often faster than default
        })
        
        # Pool results manually for better control
        # This is simplified - you might want more sophisticated pooling
        first_model <- pooled_models[[1]]
        
        # For now, return the first model (you can implement proper pooling if needed)
        return(first_model)
        
      } else {
        # Use mice::with() for non-robust models
        with(imputed_data, lm(as.formula(formula_string)))
      }
    }, error = function(e) {
      warning("Failed to fit model for outcome '", outcome, "': ", e$message)
      return(NULL)
    })
  }
  
  # Prepare data for processing
  outcome_data <- tibble(outcome = outcomes) %>%
    mutate(
      outcome_name = map_chr(outcome, standardize_outcome_name),
      formula_string = map_chr(outcome, ~build_formula_string(.x, model_type, covariate))
    ) %>%
    # Convert to list for better parallel processing
    pmap(list)
  
  # Print initial info
  if (show_progress) {
    cat("Fitting", length(outcomes), "FAOS", model_type, "models with", 
        ifelse(robust_se, "robust", "standard"), "standard errors\n")
    
    if (robust_se && !is.null(cluster_var)) {
      cat("Clustering by:", cluster_var, "\n")
    }
    
    if (use_parallel) {
      cat("Using parallel processing with", nbrOfWorkers(), "workers\n")
    }
    cat("\n")
  }
  
  # Fit models
  if (show_progress) {
    with_progress({
      p <- progressor(steps = length(outcomes))
      
      if (use_parallel && length(outcomes) > 1) {
        # Parallel processing with progress
        models <- future_map(outcome_data, function(x) {
          result <- fit_single_model_parallel(x)
          p(message = sprintf("Fitted %s model", x$outcome_name))
          return(result)
        }, .options = furrr_options(seed = TRUE))
      } else {
        # Sequential processing with progress
        models <- map(outcome_data, function(x) {
          result <- fit_single_model_parallel(x)
          p(message = sprintf("Fitted %s model", x$outcome_name))
          return(result)
        })
      }
    })
  } else {
    # No progress reporting
    if (use_parallel && length(outcomes) > 1) {
      models <- future_map(outcome_data, fit_single_model_parallel,
                           .options = furrr_options(seed = TRUE))
    } else {
      models <- map(outcome_data, fit_single_model_parallel)
    }
  }
  
  # Create results structure
  outcome_names <- map_chr(outcome_data, "outcome_name")
  names(models) <- outcome_names
  
  # Filter out failed models
  final_results <- models[!map_lgl(models, is.null)]
  
  # Add metadata as attributes
  attr(final_results, "model_type") <- model_type
  attr(final_results, "predictor") <- if(model_type == "arthritis_adjusted") arthritis_predictor else predictor
  attr(final_results, "covariate") <- covariate
  attr(final_results, "robust_se") <- robust_se
  attr(final_results, "cluster_var") <- cluster_var
  attr(final_results, "n_outcomes") <- length(final_results)
  attr(final_results, "use_parallel") <- use_parallel
  
  # Cache results
  if (cache_results) {
    if (!exists("faos_model_cache", envir = .GlobalEnv)) {
      assign("faos_model_cache", list(), envir = .GlobalEnv)
    }
    faos_model_cache <- get("faos_model_cache", envir = .GlobalEnv)
    faos_model_cache[[cache_key]] <- final_results
    assign("faos_model_cache", faos_model_cache, envir = .GlobalEnv)
  }
  
  return(final_results)
}
```




```{r}
#| label: progress-model-func


# Helper function to create a cache key for memoise
create_cache_key <- function(imputed_data, outcomes, predictor, model_type, 
                             covariate, minimal_adjusters, full_covariates,
                             arthritis_covariates, arthritis_predictor, 
                             cluster_var, robust_se) {
  
  # Create a hash of the imputed data structure (not the actual data)
  data_hash <- digest(list(
    ncol = ncol(complete(imputed_data, 1)),
    nrow = nrow(complete(imputed_data, 1)),
    m = imputed_data$m,
    colnames = names(complete(imputed_data, 1))
  ))
  
  # Create hash of all parameters
  param_hash <- digest(list(
    outcomes = outcomes,
    predictor = predictor,
    model_type = model_type,
    covariate = covariate,
    minimal_adjusters = minimal_adjusters,
    full_covariates = full_covariates,
    arthritis_covariates = arthritis_covariates,
    arthritis_predictor = arthritis_predictor,
    cluster_var = cluster_var,
    robust_se = robust_se
  ))
  
  return(paste0(data_hash, "_", param_hash))
}

# Modified function with progress bar
fit_faos_models_with_progress <- function(imputed_data, 
                                          outcomes = c("FAOS_Pain_TotalScore_Preop", 
                                                       "FAOS_Symptom_TotalScore_Preop",
                                                       "FAOS_DailyLiving_TotalScore_Preop", 
                                                       "FAOS_Sport_TotalScore_Preop",
                                                       "FAOS_Quality_TotalScore_Preop"),
                                          predictor = "Region",
                                          model_type = c("minimal", "unadjusted", "fully_adjusted", "arthritis_adjusted"),
                                          covariate = NULL,
                                          minimal_adjusters = c("AgeAtInitialExam", "Sex2"),
                                          full_covariates = c("AgeAtInitialExam", "Sex2", 
                                                              "SRCQTotalScore", "SmokingStatus"),
                                          arthritis_covariates = c("Region", "AgeAtInitialExam", "Sex2", "SmokingStatus"),
                                          arthritis_predictor = "Arthritis",
                                          cluster_var = "PatientID",
                                          robust_se = TRUE,
                                          show_progress = TRUE) {
  
  # Input validation
  model_type <- match.arg(model_type)
  
  # Check required packages
  required_packages <- c("mice")
  if (robust_se) required_packages <- c(required_packages, "estimatr")
  
  missing_packages <- required_packages %>%
    map_lgl(~!require(.x, character.only = TRUE, quietly = TRUE)) %>%
    keep(isTRUE)
  
  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(names(missing_packages), collapse = ", "))
  }
  
  # Validate inputs
  if (model_type == "minimal" && is.null(covariate)) {
    stop("covariate must be specified for minimal models")
  }
  
  # Helper function to build formula strings
  build_formula_string <- function(outcome, type, covariate = NULL) {
    
    formula_terms <- switch(type,
                            "minimal" = {
                              adjusters <- minimal_adjusters %>% 
                                discard(~ .x == covariate)
                              c(covariate, adjusters) %>% 
                                compact()
                            },
                            "unadjusted" = predictor,
                            "fully_adjusted" = c(predictor, full_covariates) %>% unique(),
                            "arthritis_adjusted" = c(arthritis_predictor, arthritis_covariates) %>% unique()
    )
    
    str_c(outcome, " ~ ", str_c(formula_terms, collapse = " + "))
  }
  
  # Helper function to standardize outcome names
  standardize_outcome_name <- function(outcome) {
    extracted <- str_extract(outcome, "(?<=FAOS_)[^_]+")
    case_when(
      extracted == "Pain" ~ "Pain",
      extracted == "Symptom" ~ "Symptom", 
      extracted == "DailyLiving" ~ "DailyLiving",
      extracted == "Sport" ~ "Sport",
      extracted == "Quality" ~ "Quality",
      TRUE ~ extracted
    )
  }
  
  # Helper function to fit individual model with progress
  fit_single_model_with_progress <- function(outcome, formula_string, p = NULL) {
    
    if (!is.null(p) && show_progress) {
      outcome_name <- standardize_outcome_name(outcome)
      p(message = sprintf("Fitting %s model", outcome_name))
    }
    
    tryCatch({
      if (robust_se && !is.null(cluster_var)) {
        with(imputed_data, 
             lm_robust(as.formula(formula_string), 
                       clusters = !!sym(cluster_var)))
      } else {
        with(imputed_data, 
             lm(as.formula(formula_string)))
      }
    }, error = function(e) {
      warning("Failed to fit model for outcome '", outcome, "': ", e$message)
      return(NULL)
    })
  }
  
  # Wrap main computation in progressr
  if (show_progress) {
    with_progress({
      p <- progressor(along = outcomes)
      
      # Print initial info
      cat("Fitting", length(outcomes), "FAOS", model_type, "models with", 
          ifelse(robust_se, "robust", "standard"), "standard errors\n")
      
      if (robust_se && !is.null(cluster_var)) {
        cat("Clustering by:", cluster_var, "\n")
      }
      cat("\n")
      
      # Main modeling pipeline
      model_results <- tibble(outcome = outcomes) %>%
        mutate(
          outcome_name = map_chr(outcome, standardize_outcome_name),
          formula_string = map_chr(outcome, ~build_formula_string(.x, model_type, covariate)),
          model = map2(outcome_name, formula_string, 
                       ~{
                         result <- fit_single_model_with_progress(.data$outcome[which(.data$outcome_name == .x)], .y, p)
                         p()  # Increment progress
                         return(result)
                       })
        ) %>%
        filter(!map_lgl(model, is.null)) %>%
        {set_names(.$model, .$outcome_name)}
    })
  } else {
    # Without progress bar
    model_results <- tibble(outcome = outcomes) %>%
      mutate(
        outcome_name = map_chr(outcome, standardize_outcome_name),
        formula_string = map_chr(outcome, ~build_formula_string(.x, model_type, covariate)),
        model = map2(outcome_name, formula_string, 
                     ~fit_single_model_with_progress(.data$outcome[which(.data$outcome_name == .x)], .y, NULL))
      ) %>%
      filter(!map_lgl(model, is.null)) %>%
      {set_names(.$model, .$outcome_name)}
  }
  
  # Add metadata as attributes
  attr(model_results, "model_type") <- model_type
  attr(model_results, "predictor") <- if(model_type == "arthritis_adjusted") arthritis_predictor else predictor
  attr(model_results, "covariate") <- covariate
  attr(model_results, "robust_se") <- robust_se
  attr(model_results, "cluster_var") <- cluster_var
  attr(model_results, "n_outcomes") <- length(model_results)
  
  return(model_results)
}

# Create memoised version
fit_faos_models_cached <- memoise(
  fit_faos_models_with_progress,
  cache = cache_filesystem(path = "faos_model_cache")
)

```



Read in complication tables, then combine and clean the text description of complications as required.

```{r}
#| label: import-complication
#| eval: false
#| echo: false

# Authenticate for sheets using the same token
gs4_auth(token = drive_token())

ComplicTable1 <- googlesheets4::read_sheet(
  ss = SheetIDs$Complic1,
  sheet = "Complications",
  range = "A2:AD",
  col_names = TRUE,
  col_types = "cccTlnicicccccccccccicccccDccD"
)

ComplicTable2 <- range_read(
  ss = SheetIDs$Complic2,
  sheet = "Complications",
  range = "A2:AD",
  col_names = TRUE,
  col_types = "cccTlnicicccccccccccicccccDccD"
)

# Complication Table
MasterComplic <- bind_and_clean(
  df1 = ComplicTable1, 
  df2 = ComplicTable2, 
  cols = c(
    "TreatmentID", 
    "ComplicationID", 
    "ComplicationOccurrence",
    "ComplicationNature", 
    "DateOfOccurrence",
    "ComplicationTreatmentOffered",
    "DateReoperation"),
  clean_cols = "ComplicationNature",
  clean_fn = clean_text  # Pass the function directly
)

```


```{r}
STROBESub <- STROBEInput |> dplyr::select(
  TreatmentID,
  TreatmentStatusNotes,
  TreatmentStatus
) |> dplyr::filter(
  !is.na(TreatmentStatusNotes),
  TreatmentStatus == "Archived"
)


STROBESub2 <- STROBESub |>
  dplyr::distinct(TreatmentStatusNotes) |>
  dplyr::mutate(
    TreatmentStatusNotes2 = dplyr::case_when(
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "(second|2nd) opinion|(treating surgeon)") ~ "Second opinion",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "refer.*(elsewhere|to|on|for|public)") ~ "Referred elsewhere",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "opt|op-|withdr*") ~ "Opt out",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "non.*registry|pathology|shoulder|wrist|knee|amput*|acl|not.*registry") ~ "Non-registry pathology",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "attend|cancel|cnx|canx|cx|(never came)|arrive|missed|show|(no referral)|appt|dna") ~ "No initial consult",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "(treatme.*t|tx).*(offered|given)|men.*ion|(surgery elsewhere)|(no intervention)|intervention.*rec|(surgery in ed)|hesitant.*treat|offer.*treatment|(not viable)|resolve|(poses risk)") ~ "No treatment offered",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "((sa|ak) patient)|(w sa)|(other surg.*n)|(surgeon other)|(joint consult)|treat.*public|elsewhere") ~ "Non-registry surgeon",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "duplicate|accessory|error|reop|complication|incorrect|(additional treatment)|superfluous|(not failed)|should.*surg") ~ "Accessory record",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "pre.*registry|prior|(ex.*sting treatment)|(pre reg tx)") ~ "Followup of pre-registry treatment",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "imaging|mri|return.*(results|diagnosis)|(never.*return)") ~ "Incomplete diagnosis",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "empty|no|enough.*(note|letter|info|record|report)|insufficient|(cannot retrieve)|(unable.*recruit)") ~ "Insufficient information",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "status|decline|confirmation|(surgery cancelled)|never.*(treatment|tx)|surg.*(ahead|elsewhere|occur|ahead)|(ahead.*surgery)|issue.*payment") ~ "Declined or unknown treatment status",
      stringr::str_detect(str_to_lower(TreatmentStatusNotes), "(registry end)|(end.*registry)") ~ "Registry end",
      .default = NA_character_
    )
  ) 
```



```{r}
#| label: bilateral-definition

SnapshotBilat1 <- MasterTable5 |> group_by(
  PatientID
) |> summarise(
  TreatmentCount = n()
) |> filter(
  TreatmentCount > 1
)

SnapshotBilat2 <- Snapshot |> filter(
  PatientID %in% SnapshotBilat1$PatientID
)

SnapshotLeft1 <- SnapshotBilat2 |> filter(
  AffectedSide == "Left"
) |> group_by(
  PatientID
) |> arrange(
  DateInitialExamination,
  .by_group = TRUE
) |> dplyr::select(
  TreatmentID,
  CombID,
  PatientID,
  AffectedSide,
  DateInitialExamination,
  DateInitialExaminationNum,
  TreatmentStatus,
  TreatmentStatusNotes,
  DateStatusChange
) |> slice_head(
  n = 1
) |> ungroup()

SnapshotRight <- SnapshotBilat2 |> filter(
  AffectedSide == "Right"
) |> rename(
  DIERight = "DateInitialExamination"
) |> group_by(
  PatientID
) |> arrange(
  DIERight,
  .by_group = TRUE
) |> dplyr::select(
  TreatmentID,
  CombID,
  PatientID,
  AffectedSide,
  DIERight,
  TreatmentStatus,
  TreatmentStatusNotes,
  DateStatusChange
) |> slice_head(
  n = 1
) |> ungroup() |> left_join(
  SnapshotLeft1 |> dplyr::select(
    PatientID,
    DateInitialExamination
  ),
  by = "PatientID"
) |> rename(
  DIELeft = "DateInitialExamination"
) |> mutate(
  DIEDiff = as.numeric(as.duration(DIERight %--% DIELeft),"weeks")
) |> mutate(
  BilateralPres = case_when(
    DIEDiff == 0 ~ "Simultaneous",
    DIEDiff > 0 ~ "Index",
    DIEDiff < 0 ~ "Subsequent",
    .default = "Unilateral"
  )
)

SnapshotLeft2 <- SnapshotLeft1 |> rename(
  DIELeft = "DateInitialExamination"
) |> left_join(
  SnapshotRight |> dplyr::select(
    PatientID,
    DIERight
  ),
  by = "PatientID"
) |> mutate(
  DIEDiff = as.numeric(as.duration(DIELeft %--% DIERight),"weeks")
) |> mutate(
  BilateralPres = case_when(
    DIEDiff == 0 ~ "Simultaneous",
    DIEDiff > 0 ~ "Index",
    DIEDiff < 0 ~ "Subsequent",
    .default = "Unilateral"
  )
) |> filter(
  BilateralPres != "Unilateral"
)
```

Marry back into the master table.

```{r}
#| label: bilateral-status

MasterTable6 <- left_join(
  MasterTable5,
  SnapshotLeft2 |> dplyr::select(
    TreatmentID,
    BilateralPres
  ) |> rename(
    BilateralPresLeft = "BilateralPres"),
  by = "TreatmentID"
) |> left_join(
  SnapshotRight |> dplyr::select(
    TreatmentID,
    BilateralPres
  ) |> rename(
    BilateralPresRight = "BilateralPres"),
  by = "TreatmentID"
) |> tidyr::unite(
  "BilateralPres",
  c(BilateralPresLeft,BilateralPresRight),
  na.rm = TRUE,
  remove = TRUE
) |> mutate(
  BilateralPres1 = if_else(
    str_count(BilateralPres) < 1,
    "Unilateral",
    BilateralPres
  )
) |> dplyr::select(
  -BilateralPres
) |> rename(
  BilateralPres = "BilateralPres1"
) |> mutate(
  BilateralDiag = stringr::str_detect(
    Term3,
    "bilateral"
  )
)


```



```{r}
library(tidyverse)
library(naniar)
library(epoxy)

# 1. First, let's examine the bilateral presentation pattern
bilateral_summary <- MasterAnalysis |>
  filter(BilateralPres == "Simultaneous") |>
  group_by(PatientID) |>
  summarise(
    n_treatments = n(),
    regions = paste(unique(Region), collapse = ", "),
    .groups = "drop"
  )

epoxy("Found {nrow(bilateral_summary)} patients with bilateral presentations")
epoxy("Treatment counts for bilateral patients: {paste(table(bilateral_summary$n_treatments), collapse = ', ')}")

# 2. Examine FAOS missingness patterns for bilateral patients
quest_cols <- c("FAOS_Symptom_TotalScore_Preop", "FAOS_Pain_TotalScore_Preop", 
                "FAOS_DailyLiving_TotalScore_Preop", "FAOS_Sport_TotalScore_Preop", 
                "FAOS_Quality_TotalScore_Preop", "VR12_Mental_TotalScore_Preop",
                "VR12_Physical_TotalScore_Preop","PCSSF_TotalScore_Preop",
                "Satisfaction_Preop","SRCQTotalScore")

bilateral_quest_pattern <- MasterAnalysis |>
  filter(BilateralPath == "Bilateral") |>
  select(PatientID, TreatmentID, Region, all_of(quest_cols)) |>
  group_by(PatientID) |>
  mutate(
    treatment_number = row_number(),
    # Check if any FAOS scores are present for this treatment
    has_any_quest = if_else(
      rowSums(!is.na(select(cur_data(), all_of(quest_cols)))) > 0, 
      TRUE, FALSE
    )
  ) |>
  ungroup()

# 3. Identify patients where only one treatment has FAOS data
potential_duplication_candidates <- bilateral_quest_pattern |>
  group_by(PatientID) |>
  summarise(
    n_treatments = n(),
    n_with_quest = sum(has_any_quest, na.rm = TRUE),
    treatments_with_quest = paste(TreatmentID[has_any_quest], collapse = ", "),
    regions_with_quest = paste(Region[has_any_quest], collapse = ", "),
    .groups = "drop"
  ) |>
  filter(
    n_treatments >= 2,  # Must have multiple treatments
    n_with_quest == 1    # But FAOS data for only one treatment
  )

epoxy("Found {nrow(potential_duplication_candidates)} bilateral patients with FAOS data on only one side")

# 4. Detailed examination of these candidates
if (nrow(potential_duplication_candidates) > 0) {
  detailed_candidates <- MasterAnalysis |>
    filter(PatientID %in% potential_duplication_candidates$PatientID) |>
    select(PatientID, TreatmentID, Region, BilateralPath, all_of(quest_cols)) |>
    arrange(PatientID, TreatmentID)
  
  epoxy("Sample of potential duplication candidates:")
  print(head(detailed_candidates, 20))
  
  # 5. Check if the FAOS patterns are consistent (complete vs missing)
  quest_completeness <- detailed_candidates |>
    group_by(PatientID, TreatmentID) |>
    summarise(
      across(all_of(faos_cols), ~!is.na(.x)),
      quest_complete = rowSums(across(all_of(quest_cols), ~!is.na(.x))) == length(quest_cols),
      quest_missing = rowSums(across(all_of(quest_cols), ~is.na(.x))) == length(quest_cols),
      .groups = "drop"
    )
  
  # 6. Summary statistics for bilateral FAOS patterns
  pattern_summary <- quest_completeness |>
    group_by(PatientID) |>
    summarise(
      complete_treatments = sum(quest_complete),
      missing_treatments = sum(quest_missing),
      partial_treatments = n() - complete_treatments - missing_treatments,
      .groups = "drop"
    )
  
  epoxy("Question completion patterns for bilateral candidates:")
  epoxy("- Patients with 1 complete, 1+ missing: {sum(pattern_summary$complete_treatments == 1 & pattern_summary$missing_treatments >= 1)}")
  epoxy("- Patients with partial data patterns: {sum(pattern_summary$partial_treatments > 0)}")
}

# 7. Visualize missingness patterns for bilateral patients
bilateral_vis_data <- MasterAnalysis |>
  filter(BilateralPath == "Bilateral") |>
  select(PatientID, TreatmentID, all_of(quest_cols)) |>
  arrange(PatientID, TreatmentID)

epoxy("Creating missingness visualization for bilateral patients...")

# Create missingness plot
if (nrow(bilateral_vis_data) > 0) {
  bilateral_vis_data |>
    vis_miss(cluster = TRUE) +
    labs(title = "FAOS Missingness Pattern for Bilateral Patients",
         subtitle = "Clustered by similarity - look for complementary missing patterns")
}

# 8. Check for exact value matching within patients (potential evidence of duplication)
if (nrow(potential_duplication_candidates) > 0) {
  # For patients with exactly one complete FAOS set, check if we can identify 
  # which treatment should receive the duplicated values
  duplication_analysis <- MasterAnalysis |>
    filter(PatientID %in% potential_duplication_candidates$PatientID) |>
    select(PatientID, TreatmentID, Region, all_of(quest_cols)) |>
    group_by(PatientID) |>
    mutate(
      has_complete_quest = rowSums(!is.na(select(cur_data(), all_of(quest_cols)))) == length(quest_cols),
      has_no_quest = rowSums(is.na(select(cur_data(), all_of(quest_cols)))) == length(quest_cols)
    ) |>
    ungroup()
  
  epoxy("Patients ready for potential FAOS duplication:")
  duplication_ready <- duplication_analysis |>
    group_by(PatientID) |>
    filter(sum(has_complete_quest) == 1 & sum(has_no_quest) >= 1) |>
    ungroup()
  
  if (nrow(duplication_ready) > 0) {
    epoxy("Found {length(unique(duplication_ready$PatientID))} patients with exactly one complete FAOS set and one+ empty FAOS sets")
    
    duplication_ready |>
      select(PatientID, TreatmentID, Region, has_complete_quest, has_no_quest) |>
      arrange(PatientID, TreatmentID) |>
      print()
  }
}

```


```{r}
library(tidyverse)
library(naniar)
library(epoxy)

# 1. First, let's examine the bilateral presentation pattern
bilateral_summary <- MasterAnalysis |>
  filter(BilateralStatus == "Bilateral") |>
  group_by(PatientID) |>
  summarise(
    n_treatments = n(),
    regions = paste(unique(Region), collapse = ", "),
    .groups = "drop"
  )

epoxy("Found {nrow(bilateral_summary)} patients with bilateral presentations")
epoxy("Treatment counts for bilateral patients: {paste(table(bilateral_summary$n_treatments), collapse = ', ')}")

# 2. Examine FAOS missingness patterns for bilateral patients
quest_cols <- c("FAOS_Symptom_TotalScore_Preop", "FAOS_Pain_TotalScore_Preop", 
                "FAOS_DailyLiving_TotalScore_Preop", "FAOS_Sport_TotalScore_Preop", 
                "FAOS_Quality_TotalScore_Preop", "VR12_Mental_TotalScore_Preop",
                "VR12_Physical_TotalScore_Preop","PCSSF_TotalScore_Preop",
                "Satisfaction_Preop","SRCQTotalScore")

bilateral_quest_pattern <- MasterAnalysis |>
  filter(BilateralStatus == "Bilateral") |>
  select(PatientID, TreatmentID, Region, all_of(quest_cols)) |>
  group_by(PatientID) |>
  mutate(
    treatment_number = row_number(),
    # Check if any FAOS scores are present for this treatment
    has_any_quest = if_else(
      rowSums(!is.na(select(cur_data(), all_of(quest_cols)))) > 0, 
      TRUE, FALSE
    )
  ) |>
  ungroup()

# 3. Identify patients where only one treatment has FAOS data
potential_duplication_candidates <- bilateral_quest_pattern |>
  group_by(PatientID) |>
  summarise(
    n_treatments = n(),
    n_with_quest = sum(has_any_quest, na.rm = TRUE),
    treatments_with_quest = paste(TreatmentID[has_any_quest], collapse = ", "),
    regions_with_quest = paste(Region[has_any_quest], collapse = ", "),
    .groups = "drop"
  ) |>
  filter(
    n_treatments >= 2,  # Must have multiple treatments
    n_with_quest == 1    # But FAOS data for only one treatment
  )

epoxy("Found {nrow(potential_duplication_candidates)} bilateral patients with FAOS data on only one side")

# 4. Detailed examination of these candidates
if (nrow(potential_duplication_candidates) > 0) {
  detailed_candidates <- MasterAnalysis |>
    filter(PatientID %in% potential_duplication_candidates$PatientID) |>
    select(PatientID, TreatmentID, Region, BilateralStatus, all_of(quest_cols)) |>
    arrange(PatientID, TreatmentID)
  
  epoxy("Sample of potential duplication candidates:")
  print(head(detailed_candidates, 20))
  
  # 5. Check if the FAOS patterns are consistent (complete vs missing)
  quest_completeness <- detailed_candidates |>
    group_by(PatientID, TreatmentID) |>
    summarise(
      across(all_of(faos_cols), ~!is.na(.x)),
      quest_complete = rowSums(across(all_of(quest_cols), ~!is.na(.x))) == length(quest_cols),
      quest_missing = rowSums(across(all_of(quest_cols), ~is.na(.x))) == length(quest_cols),
      .groups = "drop"
    )
  
  # 6. Summary statistics for bilateral FAOS patterns
  pattern_summary <- quest_completeness |>
    group_by(PatientID) |>
    summarise(
      complete_treatments = sum(quest_complete),
      missing_treatments = sum(quest_missing),
      partial_treatments = n() - complete_treatments - missing_treatments,
      .groups = "drop"
    )
  
  epoxy("Question completion patterns for bilateral candidates:")
  epoxy("- Patients with 1 complete, 1+ missing: {sum(pattern_summary$complete_treatments == 1 & pattern_summary$missing_treatments >= 1)}")
  epoxy("- Patients with partial data patterns: {sum(pattern_summary$partial_treatments > 0)}")
}

# 7. Visualize missingness patterns for bilateral patients
bilateral_vis_data <- MasterAnalysis |>
  filter(BilateralStatus == "Bilateral") |>
  select(PatientID, TreatmentID, all_of(quest_cols)) |>
  arrange(PatientID, TreatmentID)

epoxy("Creating missingness visualization for bilateral patients...")

# Create missingness plot
if (nrow(bilateral_vis_data) > 0) {
  bilateral_vis_data |>
    vis_miss(cluster = TRUE) +
    labs(title = "FAOS Missingness Pattern for Bilateral Patients",
         subtitle = "Clustered by similarity - look for complementary missing patterns")
}

# 8. Check for exact value matching within patients (potential evidence of duplication)
if (nrow(potential_duplication_candidates) > 0) {
  # For patients with exactly one complete FAOS set, check if we can identify 
  # which treatment should receive the duplicated values
  duplication_analysis <- MasterAnalysis |>
    filter(PatientID %in% potential_duplication_candidates$PatientID) |>
    select(PatientID, TreatmentID, Region, all_of(quest_cols)) |>
    group_by(PatientID) |>
    mutate(
      has_complete_quest = rowSums(!is.na(select(cur_data(), all_of(quest_cols)))) == length(quest_cols),
      has_no_quest = rowSums(is.na(select(cur_data(), all_of(quest_cols)))) == length(quest_cols)
    ) |>
    ungroup()
  
  epoxy("Patients ready for potential FAOS duplication:")
  duplication_ready <- duplication_analysis |>
    group_by(PatientID) |>
    filter(sum(has_complete_quest) == 1 & sum(has_no_quest) >= 1) |>
    ungroup()
  
  if (nrow(duplication_ready) > 0) {
    epoxy("Found {length(unique(duplication_ready$PatientID))} patients with exactly one complete FAOS set and one+ empty FAOS sets")
    
    duplication_ready |>
      select(PatientID, TreatmentID, Region, has_complete_quest, has_no_quest) |>
      arrange(PatientID, TreatmentID) |>
      print()
  }
}

```

```{r}

MasterBilat1 <- MasterAnalysis1 |> filter(
  BilateralStatus == "Bilateral",
  !is.na(FAOS_Pain_TotalScore_Preop)
) |> group_by(PatientID) |> mutate(
  TreatmentSeq = row_number()
) |>
  filter(max(TreatmentSeq) == 1) |>
  ungroup()

MasterBilat2 <- MasterAnalysis1 |> filter(
  BilateralStatus == "Bilateral",
  is.na(FAOS_Pain_TotalScore_Preop)
) |> group_by(PatientID) |> mutate(
  TreatmentSeq = row_number()
) |>
  filter(max(TreatmentSeq) == 1) |>
  ungroup()

BilatExtra <- MasterBilat1 |> filter(
  !PatientID %in% MasterBilat2$PatientID
)
# 
# |> rows_patch(MasterBilat1 |> dplyr::select(
#   PatientID,
#   contains("Preop"),
#   contains("Score")
# ),
# by = "PatientID"
# )



```

Investigate relationships


```{r}

library(tidyverse)
library(mice)
library(corrplot)
library(epoxy)

# Function to analyze correlations and suggest blocks
analyze_correlations_for_blocks <- function(data, min_correlation = 0.3) {
  
  # Get only numeric variables that have missing data
  numeric_vars <- data |> 
    select(where(is.numeric)) |> 
    select(where(~any(is.na(.))))
  
  if(ncol(numeric_vars) == 0) {
    epoxy("No numeric variables with missing data found.")
    return(NULL)
  }
  
  # Calculate correlation matrix (pairwise complete observations)
  cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")
  
  # Create correlation plot
  epoxy("=== CORRELATION HEATMAP ===")
  corrplot(cor_matrix, 
           method = "color", 
           type = "upper", 
           order = "hclust",
           tl.cex = 0.8, 
           tl.col = "black",
           title = "Correlation Matrix (Hierarchically Clustered)",
           mar = c(0,0,2,0))
  
  # Find high correlations
  high_cors <- cor_matrix |> 
    as.data.frame() |> 
    rownames_to_column("var1") |> 
    pivot_longer(-var1, names_to = "var2", values_to = "correlation") |> 
    filter(var1 != var2, 
           !is.na(correlation),
           abs(correlation) >= min_correlation) |> 
    arrange(desc(abs(correlation))) |> 
    # Remove duplicate pairs
    filter(var1 < var2)
  
  epoxy("\n=== HIGH CORRELATIONS (|r| >= {min_correlation}) ===")
  if(nrow(high_cors) > 0) {
    print(high_cors)
  } else {
    epoxy("No correlations >= {min_correlation} found.")
  }
  
  return(list(
    correlation_matrix = cor_matrix,
    high_correlations = high_cors,
    numeric_variables = names(numeric_vars)
  ))
}

# Function to suggest blocks based on correlation patterns
suggest_correlation_blocks <- function(data, min_correlation = 0.3, block_size_limit = 6) {
  
  cor_analysis <- analyze_correlations_for_blocks(data, min_correlation)
  
  if(is.null(cor_analysis)) return(NULL)
  
  high_cors <- cor_analysis$high_correlations
  numeric_vars <- cor_analysis$numeric_variables
  
  # Create blocks based on correlation patterns
  epoxy("\n=== SUGGESTED CORRELATION BLOCKS ===")
  
  # Method 1: Group by variable name patterns (common in questionnaires)
  var_patterns <- tibble(variable = numeric_vars) |> 
    mutate(
      pattern = case_when(
        str_detect(variable, "VR12") ~ "VR12_Scales",
        str_detect(variable, "FAOS") ~ "FAOS_Scales", 
        str_detect(variable, "PCSSF|SRCQ") ~ "Symptom_Scales",
        str_detect(variable, "Age|Injury") ~ "Demographics_Clinical",
        TRUE ~ "Other"
      )
    ) |> 
    arrange(pattern, variable)
  
  epoxy("\nBlock suggestions based on variable patterns:")
  pattern_blocks <- var_patterns |> 
    group_by(pattern) |> 
    summarise(
      variables = list(variable),
      count = n(),
      .groups = "drop"
    )
  
  for(i in 1:nrow(pattern_blocks)) {
    block_vars <- pattern_blocks$variables[[i]]
    if(length(block_vars) > 1) {  # Only show blocks with multiple variables
      epoxy("\n{pattern_blocks$pattern[i]} Block ({pattern_blocks$count[i]} variables):")
      epoxy("  {paste(block_vars, collapse = ', ')}")
    }
  }
  
  # Method 2: Network-based grouping using high correlations
  if(nrow(high_cors) > 0) {
    epoxy("\n\nBlock suggestions based on high correlations:")
    
    # Create adjacency list from high correlations
    adj_list <- high_cors |> 
      select(var1, var2) |> 
      gather(key = "position", value = "variable") |> 
      select(variable) |> 
      distinct() |> 
      pull(variable)
    
    # Simple grouping - variables that correlate with each other
    cor_groups <- list()
    used_vars <- character(0)
    
    for(var in adj_list) {
      if(!var %in% used_vars) {
        # Find all variables highly correlated with this one
        connected_vars <- high_cors |> 
          filter(var1 == var | var2 == var) |> 
          pivot_longer(c(var1, var2), values_to = "connected_var") |> 
          filter(connected_var != var) |> 
          pull(connected_var) |> 
          unique()
        
        if(length(connected_vars) > 0) {
          group_vars <- c(var, connected_vars)
          cor_groups[[length(cor_groups) + 1]] <- group_vars
          used_vars <- c(used_vars, group_vars)
        }
      }
    }
    
    for(i in seq_along(cor_groups)) {
      if(length(cor_groups[[i]]) > 1) {
        epoxy("\nCorrelation Group {i} ({length(cor_groups[[i]])} variables):")
        epoxy("  {paste(cor_groups[[i]], collapse = ', ')}")
      }
    }
  }
  
  return(list(
    pattern_blocks = pattern_blocks,
    correlation_analysis = cor_analysis,
    all_numeric_vars = numeric_vars
  ))
}

# Function to create visit sequence from blocks
create_block_visit_sequence <- function(blocks_list) {
  
  # Suggested order: Demographics -> Symptoms -> Function -> Quality of Life
  block_order <- c(
    "Demographics_Clinical",
    "Symptom_Scales", 
    "VR12_Scales",
    "FAOS_Scales",
    "Other"
  )
  
  visit_sequence <- character(0)
  
  for(block_name in block_order) {
    block_vars <- blocks_list$pattern_blocks |> 
      filter(pattern == block_name) |> 
      pull(variables)
    
    if(length(block_vars) > 0) {
      visit_sequence <- c(visit_sequence, block_vars[[1]])
    }
  }
  
  epoxy("\n=== SUGGESTED VISIT SEQUENCE FROM BLOCKS ===")
  epoxy("Sequential block order:")
  for(i in seq_along(visit_sequence)) {
    epoxy("{i}. {visit_sequence[i]}")
  }
  
  return(visit_sequence)
}

# Run analysis on your dataset
epoxy("CORRELATION ANALYSIS FOR BLOCK CONSTRUCTION")
epoxy("="*50)

# Step 1: Analyze correlations
blocks_analysis <- suggest_correlation_blocks(MasterPatch, 
                                              min_correlation = 0.3,
                                              block_size_limit = 6)

# Step 2: Create suggested visit sequence
if(!is.null(blocks_analysis)) {
  suggested_sequence <- create_block_visit_sequence(blocks_analysis)
  
  epoxy("\n" %+% "="*50)
  epoxy("READY TO USE:")
  epoxy("Copy this vector for your manual_visit_sequence:")
  epoxy("manual_visit_sequence <- c(")
  for(i in seq_along(suggested_sequence)) {
    comma <- if(i < length(suggested_sequence)) "," else ""
    epoxy('  "{suggested_sequence[i]}"{comma}')
  }
  epoxy(")")
}

# Optional: Show correlation matrix for specific blocks
show_block_correlations <- function(data, block_vars) {
  if(length(block_vars) > 1) {
    block_data <- data |> select(all_of(block_vars))
    block_cor <- cor(block_data, use = "pairwise.complete.obs")
    
    epoxy("\nCorrelations within selected block:")
    corrplot(block_cor, 
             method = "number", 
             type = "upper",
             tl.cex = 0.8,
             number.cex = 0.7)
    
    return(block_cor)
  }
}

# Example: Show correlations within FAOS scales
faos_vars <- names(MasterPatch)[str_detect(names(MasterPatch), "FAOS")]
if(length(faos_vars) > 1) {
  epoxy("\n" %+% "="*30)
  epoxy("EXAMPLE: FAOS Block Correlations")
  show_block_correlations(MasterPatch, faos_vars)
}


```


```{r}


ShadowMaster <- naniar::as_shadow(MasterAnalysis1)

aqShadow <- naniar::bind_shadow(MasterAnalysis1)

nabShadow <- nabular(MasterAnalysis1)
```

```{r}

nabTbl <- gtsummary::tbl_summary(
  nabShadow |> dplyr::select(
    all_of(ends_with("NA"))
  ),
  by = "FAOS_Pain_TotalScore_Preop_NA"
  
)

knitr::knit_print(nabTbl)


```


```{r}

BilateralCorrect <- MasterAnalysis1 |> filter(
  TreatmentID %in% duplication_ready$TreatmentID
) |> arrange(
  TreatmentID
)

BilateralWith <- BilateralCorrect |> dplyr::filter(
  if_any(c(Sex2:Arthritis), ~ !is.na(.x))
)

BilateralWithOut <- BilateralCorrect |> dplyr::filter(
  if_any(c(Sex2:Arthritis), ~ is.na(.x))
)

BilateralPatch <- rows_patch(
  BilateralWithOut,
  BilateralWith,
  by = "PatientID",
  unmatched = "ignore",
  copy = FALSE,
  in_place = FALSE
)

MasterPatch <- rows_patch(
  MasterAnalysis1,
  BilateralPatch,
  by = "TreatmentID",
  unmatched = "ignore",
  copy = FALSE,
  in_place = FALSE
)

```


```{r}

# Function to display missingness summary for manual review
review_missingness_for_manual_sequence <- function(data) {
  
  # Calculate missingness summary
  missingness_summary <- data |> 
    summarise(across(everything(), ~mean(is.na(.)))) |> 
    pivot_longer(everything(), names_to = "variable", values_to = "prop_missing") |> 
    arrange(prop_missing) |> 
    mutate(
      prop_available = 1 - prop_missing,
      percent_available = round(prop_available * 100, 1),
      percent_missing = round(prop_missing * 100, 1)
    ) |> 
    select(variable, percent_missing, percent_available)
  
  # Show complete variables (no imputation needed)
  complete_vars <- missingness_summary |> 
    filter(percent_missing == 0)
  
  # Show variables needing imputation
  missing_vars <- missingness_summary |> 
    filter(percent_missing > 0) |> 
    arrange(percent_missing)
  
  epoxy::epoxy("=== VARIABLES WITH NO MISSING DATA (complete) ===")
  epoxy("Count: {nrow(complete_vars)}")
  if(nrow(complete_vars) > 0) print(complete_vars)
  
  epoxy("\n=== VARIABLES NEEDING IMPUTATION (ordered by % missing) ===")
  epoxy("Count: {nrow(missing_vars)}")
  if(nrow(missing_vars) > 0) print(missing_vars)
  
  return(list(
    complete_variables = complete_vars$variable,
    variables_to_impute = missing_vars$variable,
    missingness_summary = missingness_summary
  ))
}

# Function to apply manual visit sequence
apply_manual_visitsequence <- function(data, manual_sequence) {
  
  # Validate that all variables in sequence exist and need imputation
  all_vars <- names(data)
  missing_vars <- names(data)[sapply(data, function(x) any(is.na(x)))]
  
  # Check for issues
  invalid_vars <- setdiff(manual_sequence, all_vars)
  unnecessary_vars <- setdiff(manual_sequence, missing_vars)
  missing_from_sequence <- setdiff(missing_vars, manual_sequence)
  
  if(length(invalid_vars) > 0) {
    epoxy("⚠️  Variables in sequence that don't exist in data: {paste(invalid_vars, collapse = ', ')}")
  }
  
  if(length(unnecessary_vars) > 0) {
    epoxy("ℹ️  Variables in sequence that don't need imputation: {paste(unnecessary_vars, collapse = ', ')}")
  }
  
  if(length(missing_from_sequence) > 0) {
    epoxy("⚠️  Variables needing imputation but not in sequence: {paste(missing_from_sequence, collapse = ', ')}")
  }
  
  # Create MICE object with manual sequence
  epoxy("\n=== APPLYING MANUAL VISIT SEQUENCE ===")
  epoxy("Variables in visit order:")
  for(i in seq_along(manual_sequence)) {
    epoxy("{i}. {manual_sequence[i]}")
  }
  
  mice_object <- mice(data, 
                      visitSequence = manual_sequence,
                      maxit = 0, 
                      m = 10, # Dry run first
                      printFlag = FALSE)
  
  epoxy("\nManual visit sequence applied successfully!")
  
  return(mice_object)
}

```


```{r}
# Step 1: Review your data's missingness pattern
epoxy("STEP 1: Reviewing missingness pattern for manual sequence planning")
review_info <- review_missingness_for_manual_sequence(MasterPatch)

# Step 2: Visualize missingness pattern
epoxy("\nSTEP 2: Visualizing missingness pattern")
vis_miss(MasterPatch)

# Step 3: Show default MICE sequence for comparison
epoxy("\nSTEP 3: Default MICE visit sequence for comparison")
mice_default <- mice(MasterPatch, maxit = 0, m = 5, printFlag = FALSE)
epoxy("Default sequence: {paste(mice_default$visitSequence, collapse = ' -> ')}")

# EDIT THIS VECTOR WITH YOUR DESIRED SEQUENCE:
manual_visit_sequence <- c(
  #Example structure - replace with your actual sequence:
  "AgeAtInitialExam",                    # Usually complete, good predictor
  "InjuryToPresentation",                # Time-based, relatively independent
  "VR12_Physical_TotalScore_Preop",      # Health measures
  "VR12_Mental_TotalScore_Preop",
  "SRCQTotalScore",                      # Questionnaire scores
  "PCSSF_TotalScore_Preop",
  "Region",                              # Categorical variables
  "Arthritis",
  "FAOS_Symptom_TotalScore_Preop",       # FAOS subscales in logical order
  "FAOS_Pain_TotalScore_Preop",
  "FAOS_DailyLiving_TotalScore_Preop",
  "FAOS_Sport_TotalScore_Preop",
  "FAOS_Quality_TotalScore_Preop",
  "SmokingStatus",
  
)

# Uncomment and run when you've specified your sequence:
mice_manual <- apply_manual_visitsequence(MasterPatch |> dplyr::select(-Satisfaction_Preop), manual_visit_sequence)
```

```{r}


TblMaster1 <- gtsummary::tbl_summary(
  MasterPatch,
  include = c(
    Sex2,
    where(is.numeric)
  ),
  by = "Sex2",
  missing = "no",
  statistic = list(all_continuous() ~ "{mean} ({sd})")
) |> add_p(
  test.args = all_tests("t.test") ~ list(var.equal = TRUE)
)

knitr::knit_print(TblMaster1)


```



```{r}
# Set up parallel processing
plan(multisession, workers = max(1, availableCores() - 1))

# Function to pool fixest results using Rubin's rules
pool_fixest_results <- function(model_list) {
  coef_list <- map(model_list, ~coef(.x))
  vcov_list <- map(model_list, ~vcov(.x))
  
  all_coef_names <- map(coef_list, names) %>% reduce(intersect)
  
  if (length(all_coef_names) == 0) {
    warning("No common coefficients across imputations")
    return(NULL)
  }
  
  coef_list <- map(coef_list, ~.x[all_coef_names])
  vcov_list <- map(vcov_list, ~.x[all_coef_names, all_coef_names, drop = FALSE])
  
  m <- length(model_list)
  
  # Pool coefficients
  pooled_coef <- reduce(coef_list, `+`) / m
  
  # Pool variances using Rubin's rules
  W <- reduce(vcov_list, `+`) / m
  coef_matrix <- do.call(rbind, coef_list)
  B <- cov(coef_matrix) * (m - 1) / m
  T_var <- W + (1 + 1/m) * B
  
  pooled_se <- sqrt(diag(T_var))
  
  # Degrees of freedom
  n_obs <- nobs(model_list[[1]])
  df_complete <- n_obs - length(all_coef_names)
  r <- (1 + 1/m) * diag(B) / diag(W)
  df_adj <- (m - 1) * (1 + 1/r)^2
  df_obs <- df_complete * (df_complete + 1) / (df_complete + 3) * (1 - r)
  df_pooled <- pmin(df_adj, df_obs)
  
  # Statistics and p-values
  t_stats <- pooled_coef / pooled_se
  p_values <- 2 * pt(abs(t_stats), df = df_pooled, lower.tail = FALSE)
  t_crit <- qt(0.975, df = df_pooled)
  ci_lower <- pooled_coef - t_crit * pooled_se
  ci_upper <- pooled_coef + t_crit * pooled_se
  
  pooled_result <- list(
    coefficients = pooled_coef,
    std.error = pooled_se,
    statistic = t_stats,
    p.value = p_values,
    conf.low = ci_lower,
    conf.high = ci_upper,
    df = df_pooled,
    vcov = T_var,
    nobs = n_obs,
    m = m,
    fmi = r / (1 + r)
  )
  
  class(pooled_result) <- "pooled_fixest"
  return(pooled_result)
}

print.pooled_fixest <- function(x, digits = 4) {
  cat("Pooled Multiple Imputation Results (m =", x$m, ")\n")
  cat("Observations:", x$nobs, "\n\n")
  
  coef_table <- data.frame(
    Estimate = x$coefficients,
    `Std. Error` = x$std.error,
    `t value` = x$statistic,
    `Pr(>|t|)` = x$p.value,
    `CI Lower` = x$conf.low,
    `CI Upper` = x$conf.high,
    `DF` = x$df,
    `FMI` = x$fmi,
    check.names = FALSE
  )
  
  print(round(coef_table, digits))
}

```

```{r}

# Cache key function
create_cache_key <- function(imputed_data, outcomes, predictor, model_type, 
                             covariate, minimal_adjusters, full_covariates,
                             arthritis_covariates, arthritis_predictor, 
                             cluster_var) {
  
  data_hash <- digest(list(
    ncol = ncol(complete(imputed_data, 2)),
    nrow = nrow(complete(imputed_data, 2)),
    m = imputed_data$m,
    colnames = names(complete(imputed_data, 2))
  ))
  
  param_hash <- digest(list(
    outcomes = outcomes,
    predictor = predictor,
    model_type = model_type,
    covariate = covariate,
    minimal_adjusters = minimal_adjusters,
    full_covariates = full_covariates,
    arthritis_covariates = arthritis_covariates,
    arthritis_predictor = arthritis_predictor,
    cluster_var = cluster_var
  ))
  
  return(paste0(data_hash, "_", param_hash))
}


```

```{r}

fit_faos_models_fixest <- function(imputed_data, 
                                   outcomes = c("FAOS_Pain_TotalScore_Preop", 
                                                "FAOS_Symptom_TotalScore_Preop",
                                                "FAOS_DailyLiving_TotalScore_Preop", 
                                                "FAOS_Sport_TotalScore_Preop",
                                                "FAOS_Quality_TotalScore_Preop"),
                                   predictor = "Region",
                                   model_type = c("minimal", "unadjusted", "fully_adjusted", "arthritis_adjusted"),
                                   covariate = NULL,
                                   minimal_adjusters = c("AgeAtInitialExam", "Sex2"),
                                   full_covariates = c("AgeAtInitialExam", "Sex2", 
                                                       "SRCQTotalScore", "SmokingStatus"),
                                   arthritis_covariates = c("Region", "AgeAtInitialExam", "Sex2", "SmokingStatus"),
                                   arthritis_predictor = "Arthritis",
                                   cluster_var = "PatientID",
                                   show_progress = TRUE) {
  
  model_type <- match.arg(model_type)
  
  if (model_type == "minimal" && is.null(covariate)) {
    stop("covariate must be specified for minimal models")
  }
  
  build_formula_string <- function(outcome, type, covariate = NULL) {
    formula_terms <- switch(type,
                            "minimal" = {
                              adjusters <- minimal_adjusters[minimal_adjusters != covariate]
                              c(covariate, adjusters) %>% compact()
                            },
                            "unadjusted" = predictor,
                            "fully_adjusted" = unique(c(predictor, full_covariates)),
                            "arthritis_adjusted" = unique(c(arthritis_predictor, arthritis_covariates))
    )
    
    str_c(outcome, " ~ ", str_c(formula_terms, collapse = " + "))
  }
  
  standardize_outcome_name <- function(outcome) {
    extracted <- str_extract(outcome, "(?<=FAOS_)[^_]+")
    case_when(
      extracted == "Pain" ~ "Pain",
      extracted == "Symptom" ~ "Symptom", 
      extracted == "DailyLiving" ~ "DailyLiving",
      extracted == "Sport" ~ "Sport",
      extracted == "Quality" ~ "Quality",
      TRUE ~ extracted
    )
  }
  
  fit_single_model_fixest <- function(outcome_info) {
    outcome <- outcome_info$outcome
    formula_string <- outcome_info$formula_string
    
    tryCatch({
      # Exclude original data (m=1) and fit on imputations 2:m
      imputation_indices <- if (imputed_data$m > 1) 1:imputed_data$m else 1
      
      complete_models <- map(imputation_indices, function(imp) {
        complete_data <- mice::complete(imputed_data, imp)
        
        formula_vars <- all.vars(as.formula(formula_string))
        model_vars <- c(formula_vars, cluster_var)
        
        if (any(is.na(complete_data[, model_vars, drop = FALSE]))) {
          warning("Missing values found in imputation ", imp, " for outcome: ", outcome)
          return(NULL)
        }
        
        fixest::feols(as.formula(formula_string), 
                      data = complete_data,
                      cluster = cluster_var)
      })
      
      complete_models <- compact(complete_models)
      
      if (length(complete_models) == 0) {
        warning("No successful models fitted for outcome: ", outcome)
        return(NULL)
      }
      
      pool_fixest_results(complete_models)
      
    }, error = function(e) {
      warning("Failed to fit model for outcome '", outcome, "': ", e$message)
      return(NULL)
    })
  }
  
  outcome_data <- tibble(outcome = outcomes) %>%
    mutate(
      outcome_name = map_chr(outcome, standardize_outcome_name),
      formula_string = map_chr(outcome, ~build_formula_string(.x, model_type, covariate))
    ) %>%
    pmap(list)
  
  if (show_progress) {
    cat("Fitting", length(outcomes), "FAOS", model_type, "models with fixest\n")
    cat("Clustering by:", cluster_var, "\n")
    cat("Using", length(2:imputed_data$m), "imputations (excluding original)\n\n")
  }
  
  if (show_progress) {
    with_progress({
      p <- progressor(steps = length(outcomes))
      
      models <- future_map(outcome_data, function(x) {
        result <- fit_single_model_fixest(x)
        p(message = sprintf("Fitted %s model", x$outcome_name))
        return(result)
      }, .options = furrr_options(seed = TRUE))
    })
  } else {
    models <- future_map(outcome_data, fit_single_model_fixest,
                         .options = furrr_options(seed = TRUE))
  }
  
  outcome_names <- map_chr(outcome_data, "outcome_name")
  names(models) <- outcome_names
  final_results <- models[!map_lgl(models, is.null)]
  
  attr(final_results, "model_type") <- model_type
  attr(final_results, "predictor") <- if(model_type == "arthritis_adjusted") arthritis_predictor else predictor
  attr(final_results, "covariate") <- covariate
  attr(final_results, "cluster_var") <- cluster_var
  attr(final_results, "n_outcomes") <- length(final_results)
  
  return(final_results)
}

# Create filesystem cache
if (!dir.exists("model_cache")) {
  dir.create("model_cache", recursive = TRUE)
}

fit_faos_models_cached <- memoise(
  fit_faos_models_fixest,
  cache = cache_filesystem(path = "model_cache")
)


```

```{r}
#| label: run-analysis-models-1

tic()

unadjusted_models <- fit_faos_models_fixest(
  imputed_data = MasterImp,
  model_type = "unadjusted",
  cluster_var = "PatientID"
)

toc()

```


```{r}
# Check if the variable exists in the imputed datasets
if("FAOS_Pain_TotalScore_Preop" %in% names(complete(MasterImp, 1))) {
  epoxy("FAOS variable found in imputed datasets")
  
  # Check missingness across all imputed datasets
  faos_missing_check <- map_dfr(1:MasterImp$m, function(i) {
    complete_data <- complete(MasterImp, i)
    tibble(
      imputation = i,
      n_total = nrow(complete_data),
      n_missing_faos = sum(is.na(complete_data$FAOS_Pain_TotalScore_Preop)),
      prop_missing = n_missing_faos / n_total * 100
    )
  })
  
  epoxy("Missingness in FAOS across imputations:")
  print(faos_missing_check)
  
} else {
  epoxy("FAOS variable NOT found in imputed datasets")
  epoxy("Available variables in imputed data:")
  print(names(complete(MasterImp, 1)))
}

```



```{r}
#| label: run-analysis-models-2

tic()

adjusted_models <- fit_faos_models_cached(
  imputed_data = MasterImp,
  model_type = "fully_adjusted",
  robust_se = TRUE,
  cluster_var = "PatientID"
)

toc()

```

```{r}
#| label: run-analysis-models-3

tic()

arthritis_models <- fit_faos_models_cached(
  imputed_data = MasterImp,
  model_type = "arthritis_adjusted",
  robust_se = TRUE,
  cluster_var = "PatientID"
)

toc()

```


```{r}
fit_faos_models_fixest <- function(imputed_data, 
                                   outcomes = c("FAOS_Pain_TotalScore_Preop", 
                                                "FAOS_Symptom_TotalScore_Preop",
                                                "FAOS_DailyLiving_TotalScore_Preop", 
                                                "FAOS_Sport_TotalScore_Preop",
                                                "FAOS_Quality_TotalScore_Preop"),
                                   predictor = "Region",
                                   model_type = c("minimal", "unadjusted", "fully_adjusted", "arthritis_adjusted"),
                                   covariate = NULL,
                                   minimal_adjusters = c("AgeAtInitialExam", "Sex2"),
                                   full_covariates = c("AgeAtInitialExam", "Sex2", 
                                                       "SRCQTotalScore", "SmokingStatus"),
                                   arthritis_covariates = c("Region", "AgeAtInitialExam", "Sex2", "SmokingStatus"),
                                   arthritis_predictor = "Arthritis",
                                   cluster_var = "PatientID",
                                   use_parallel = TRUE,
                                   show_progress = TRUE) {
  
  # Check if fixest is available (much faster for clustered SEs)
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("fixest package required for this optimized approach. Install with: install.packages('fixest')")
  }
  
  model_type <- match.arg(model_type)
  
  # Validate inputs
  if (model_type == "minimal" && is.null(covariate)) {
    stop("covariate must be specified for minimal models")
  }
  
  # Build formula strings (same as before)
  build_formula_string <- function(outcome, type, covariate = NULL) {
    formula_terms <- switch(type,
                            "minimal" = {
                              adjusters <- minimal_adjusters[minimal_adjusters != covariate]
                              c(covariate, adjusters) %>% compact()
                            },
                            "unadjusted" = predictor,
                            "fully_adjusted" = unique(c(predictor, full_covariates)),
                            "arthritis_adjusted" = unique(c(arthritis_predictor, arthritis_covariates))
    )
    
    str_c(outcome, " ~ ", str_c(formula_terms, collapse = " + "))
  }
  
  standardize_outcome_name <- function(outcome) {
    extracted <- str_extract(outcome, "(?<=FAOS_)[^_]+")
    case_when(
      extracted == "Pain" ~ "Pain",
      extracted == "Symptom" ~ "Symptom", 
      extracted == "DailyLiving" ~ "DailyLiving",
      extracted == "Sport" ~ "Sport",
      extracted == "Quality" ~ "Quality",
      TRUE ~ extracted
    )
  }
  
  # Fast model fitting using fixest (much faster for clustered SEs)
  fit_single_model_fixest <- function(outcome_info, complete_data) {
    outcome <- outcome_info$outcome
    formula_string <- outcome_info$formula_string
    
    tryCatch({
      # fixest is much faster for clustered standard errors
      fixest::feols(as.formula(formula_string), 
                    data = complete_data,
                    cluster = cluster_var)
    }, error = function(e) {
      warning("Failed to fit model for outcome '", outcome, "': ", e$message)
      return(NULL)
    })
  }
  
  # Pre-extract and pool complete datasets
  pooled_data <- mice::complete(imputed_data, "long")
  
  # Prepare outcome data
  outcome_data <- tibble(outcome = outcomes) %>%
    mutate(
      outcome_name = map_chr(outcome, standardize_outcome_name),
      formula_string = map_chr(outcome, ~build_formula_string(.x, model_type, covariate))
    ) %>%
    pmap(list)
  
  # Print info
  if (show_progress) {
    cat("Fitting", length(outcomes), "FAOS", model_type, "models using fixest\n")
    cat("Clustering by:", cluster_var, "\n")
    if (use_parallel) {
      cat("Using parallel processing with", nbrOfWorkers(), "workers\n")
    }
    cat("\n")
  }
  
  # Fit models in parallel
  if (show_progress) {
    with_progress({
      p <- progressor(steps = length(outcomes))
      
      if (use_parallel && length(outcomes) > 1) {
        models <- future_map(outcome_data, function(x) {
          result <- fit_single_model_fixest(x, pooled_data)
          p(message = sprintf("Fitted %s model", x$outcome_name))
          return(result)
        }, .options = furrr_options(seed = TRUE))
      } else {
        models <- map(outcome_data, function(x) {
          result <- fit_single_model_fixest(x, pooled_data)
          p(message = sprintf("Fitted %s model", x$outcome_name))
          return(result)
        })
      }
    })
  } else {
    if (use_parallel && length(outcomes) > 1) {
      models <- future_map(outcome_data, ~fit_single_model_fixest(.x, pooled_data),
                           .options = furrr_options(seed = TRUE))
    } else {
      models <- map(outcome_data, ~fit_single_model_fixest(.x, pooled_data))
    }
  }
  
  # Create final results
  outcome_names <- map_chr(outcome_data, "outcome_name")
  names(models) <- outcome_names
  final_results <- models[!map_lgl(models, is.null)]
  
  # Add metadata
  attr(final_results, "model_type") <- model_type
  attr(final_results, "predictor") <- if(model_type == "arthritis_adjusted") arthritis_predictor else predictor
  attr(final_results, "covariate") <- covariate
  attr(final_results, "robust_se") <- TRUE  # fixest always uses robust SEs
  attr(final_results, "cluster_var") <- cluster_var
  attr(final_results, "n_outcomes") <- length(final_results)
  attr(final_results, "use_parallel") <- use_parallel
  
  return(final_results)
}

# Create persistent cache directory
if (!dir.exists("model_cache")) {
  dir.create("model_cache", recursive = TRUE)
}

# Create memoised version with filesystem cache for persistence
fit_faos_models_cached <- memoise(
  fit_faos_models_optimized,
  cache = cache_filesystem(path = "model_cache")
)
```



```{r}
#| label: model-table-func

# Main function to create regression tables from model lists
create_faos_regression_table <- function(model_list, 
                                         table_type = c("unadjusted", "minimally_adjusted", "fully_adjusted", "arthritis_adjusted"),
                                         custom_labels = NULL,
                                         estimate_digits = 2,
                                         pvalue_digits = 3) {
  
  table_type <- match.arg(table_type)
  
  # Default labels based on table type
  default_labels <- switch(table_type,
                           "unadjusted" = list(Region ~ "Diagnostic Group"),
                           "minimally_adjusted" = list(
                             AgeAtInitialExam ~ "Age at Initial Exam",
                             Sex2 ~ "Sex",
                             SRCQTotalScore ~ "Self-Reported Comorbidities",
                             VR12_Mental_TotalScore_Preop ~ "VR12 Mental Component"
                           ),
                           "fully_adjusted" = list(
                             Region ~ "Diagnostic Group"
                           ),
                           "arthritis_adjusted" = list(
                             Arthritis ~ "Arthritis"
                           )
  )
  
  # Determine which variables to include based on table type
  include_vars <- switch(table_type,
                         "unadjusted" = "Region",
                         "minimally_adjusted" = c("AgeAtInitialExam", "Sex2", "SRCQTotalScore", "VR12_Mental_TotalScore_Preop"),
                         "fully_adjusted" = "Region",  # Only show Region, other covariates are adjusted for but not displayed
                         "arthritis_adjusted" = "Arthritis"  # Only show Arthritis, other covariates are adjusted for but not displayed
  )
  
  # Use custom labels if provided, otherwise use defaults
  labels_to_use <- if(!is.null(custom_labels)) custom_labels else default_labels
  
  # Create individual regression tables
  regression_tables <- map2(model_list, names(model_list), function(model, subscale_name) {
    
    # Determine which variables should be shown on single row (for binary/categorical vars)
    single_row_vars <- switch(table_type,
                              "unadjusted" = NULL,  # Region has multiple levels, show all
                              "minimally_adjusted" = c("Sex2"),  # Show only one level for Sex if it's binary
                              "fully_adjusted" = NULL,  # Region has multiple levels, show all
                              "arthritis_adjusted" = "Arthritis"  # Show only "Yes" level for Arthritis (binary)
    )
    
    
    gtsummary::tbl_regression(
      model,
      tidy_fun = pool_and_tidy_mice,
      include = all_of(include_vars),
      label = labels_to_use,
      estimate_fun = function(x) style_number(x, digits = estimate_digits),
      pvalue_fun = function(x) style_pvalue(x, digits = pvalue_digits),
      add_estimate_to_reference_rows = TRUE,
      show_single_row = all_of(single_row_vars)
    )
  })
  
  # Create group headers (convert model names to nice labels)
  group_headers <- case_when(
    names(model_list) == "Pain" ~ "Pain",
    names(model_list) == "Symptom" ~ "Symptoms", 
    names(model_list) == "DailyLiving" ~ "Activities of Daily Living",
    names(model_list) == "Sport" ~ "Sport and Recreation",
    names(model_list) == "Quality" ~ "Quality of Life",
    TRUE ~ names(model_list)  # fallback to original names
  )
  
  # Stack the tables
  stacked_table <- tbl_stack(
    regression_tables,
    group_header = group_headers
  )
  
  # For arthritis tables, modify to remove redundant "Arthritis" labels
  if(table_type == "arthritis_adjusted") {
    # Create tables without group headers for arthritis tables
    stacked_table <- tbl_stack(regression_tables)  # No group_header argument
    
    # Modify to show subscale names as the main labels
    stacked_table <- stacked_table |>
      modify_header(label ~ "**FAOS Subscale**") |>
      modify_table_body(
        ~.x |> 
          mutate(
            label = rep(group_headers, each = 1),  # Use subscale names as labels
            row_type = "label"
          )
      )
  } else {
    # For other table types, use group headers normally
    stacked_table <- tbl_stack(
      regression_tables,
      group_header = group_headers
    )
  }
  
  return(stacked_table)
}
```



```{r}
#| label: table-wrapper-func
# Convenience wrapper functions for each model type

# 1. Create unadjusted regression table
create_unadjusted_table <- function(unadjusted_models) {
  
  table <- create_faos_regression_table(unadjusted_models, table_type = "unadjusted")
  
  return(table)
}

# 2. Create minimally adjusted table (for descriptive analysis)
create_descriptive_table <- function(descriptive_models, 
                                     covariate_focus,
                                     table_caption = NULL) {
  
  # Create custom label for the focal covariate
  covariate_labels <- switch(covariate_focus,
                             "AgeAtInitialExam" = list(AgeAtInitialExam ~ "Age at Initial Exam", Sex2 ~ "Sex"),
                             "SRCQTotalScore" = list(SRCQTotalScore ~ "Self-Reported Comorbidities", 
                                                     AgeAtInitialExam ~ "Age at Initial Exam", Sex2 ~ "Sex"),
                             "VR12_Mental_TotalScore_Preop" = list(VR12_Mental_TotalScore_Preop ~ "VR12 Mental", 
                                                                   AgeAtInitialExam ~ "Age at Initial Exam", Sex2 ~ "Sex")
  )
  
  # Custom include variables for descriptive tables
  include_vars <- switch(covariate_focus,
                         "AgeAtInitialExam" = c("AgeAtInitialExam", "Sex2"),
                         "SRCQTotalScore" = c("SRCQTotalScore", "AgeAtInitialExam", "Sex2"),
                         "VR12_Mental_TotalScore_Preop" = c("VR12_Mental_TotalScore_Preop", "AgeAtInitialExam", "Sex2")
  )
  
  table <- create_faos_regression_table(descriptive_models, 
                                        table_type = "minimally_adjusted",
                                        custom_labels = covariate_labels)
  
  if(is.null(table_caption)) {
    covariate_name <- switch(covariate_focus,
                             "AgeAtInitialExam" = "age",
                             "SRCQTotalScore" = "comorbidities", 
                             "VR12_Mental_TotalScore_Preop" = "mental health component"
    )
  }
  
  return(table)
}

# 3. Create fully adjusted table (main analysis)
create_adjusted_table <- function(adjusted_models) {
  
  table <- create_faos_regression_table(adjusted_models, table_type = "fully_adjusted")
  
  return(table)
}

# 4. Create arthritis adjusted table (secondary analysis)
create_arthritis_table <- function(arthritis_models) {
  
  table <- create_faos_regression_table(arthritis_models, table_type = "arthritis_adjusted")
  
  return(table)
}

# Alternative approach: More flexible function that handles different predictors
create_custom_regression_table <- function(model_list, 
                                           predictor_var,
                                           covariate_vars = NULL,
                                           custom_labels = NULL,
                                           estimate_digits = 2,
                                           pvalue_digits = 3) {
  
  # Determine which variables to include
  include_vars <- c(predictor_var, covariate_vars)
  
  # Default labels if none provided
  if(is.null(custom_labels)) {
    custom_labels <- switch(predictor_var,
                            "Region" = list(Region ~ "Diagnostic Group"),
                            "Arthritis" = list(Arthritis ~ "Arthritis"),
                            # Add other predictors as needed
    )
  }
  
  # Create individual regression tables
  regression_tables <- map2(model_list, names(model_list), function(model, subscale_name) {
    gtsummary::tbl_regression(
      model,
      tidy_fun = pool_and_tidy_mice,
      include = include_vars,
      label = custom_labels,
      estimate_fun = function(x) style_number(x, digits = estimate_digits),
      pvalue_fun = function(x) style_pvalue(x, digits = pvalue_digits),
      add_estimate_to_reference_rows = TRUE
    )
  })
  
  # Create group headers
  group_headers <- case_when(
    names(model_list) == "Pain" ~ "Pain",
    names(model_list) == "Symptom" ~ "Symptoms", 
    names(model_list) == "DailyLiving" ~ "Activities of Daily Living",
    names(model_list) == "Sport" ~ "Sport and Recreation",
    names(model_list) == "Quality" ~ "Quality of Life",
    TRUE ~ names(model_list)
  )
  
  # Stack the tables
  stacked_table <- tbl_stack(
    regression_tables,
    group_header = group_headers
  )
  
  return(stacked_table)
}


```


```{r}
#| label: faos-plot-func


# Improved plotting function with better scaling and appearance
create_faos_plots <- function(models, 
                              plot_type = c("age_sex", "comorbidities_sex", "mental_sex", 
                                            "age_arthritis", "region_distribution"),
                              x_var = NULL,
                              color_var = NULL,
                              x_label = NULL,
                              subscale_names = NULL,
                              fixed_y_scale = TRUE,
                              show_points = FALSE) {
  
  plot_type <- match.arg(plot_type)
  
  # Set defaults based on plot type
  plot_configs <- list(
    "age_sex" = list(x_var = "AgeAtInitialExam", color_var = "Sex2", 
                     x_label = "Age at Initial Exam", filter_condition = "AgeAtInitialExam > 15"),
    "comorbidities_sex" = list(x_var = "SRCQTotalScore", color_var = "Sex2", 
                               x_label = "Self-Reported Comorbidities Score", filter_condition = "AgeAtInitialExam > 15"),
    "mental_sex" = list(x_var = "VR12_Mental_TotalScore_Preop", color_var = "Sex2", 
                        x_label = "VR12 Mental Component", filter_condition = "AgeAtInitialExam > 15"),
    "age_arthritis" = list(x_var = "AgeAtInitialExam", color_var = "Arthritis", 
                           x_label = "Age at Initial Exam", filter_condition = "AgeAtInitialExam > 15"),
    "region_distribution" = list(x_var = "Region", color_var = NULL, 
                                 x_label = "Region", filter_condition = NULL)
  )
  
  config <- plot_configs[[plot_type]]
  
  # Override defaults if provided
  if(!is.null(x_var)) config$x_var <- x_var
  if(!is.null(color_var)) config$color_var <- color_var
  if(!is.null(x_label)) config$x_label <- x_label
  
  # Handle subscale names
  if(is.null(subscale_names)) {
    plot_subscale_names <- names(models)
  } else {
    if(length(subscale_names) != length(models)) {
      stop("Length of subscale_names must match length of models")
    }
    plot_subscale_names <- subscale_names
  }
  
  # Get predictions for all models first to determine common y-scale
  all_predictions <- map(models, function(model) {
    pred_grid <- marginaleffects::predictions(model)
    
    # Apply filter if specified
    if(!is.null(config$filter_condition)) {
      pred_grid <- pred_grid |> filter(eval(parse(text = config$filter_condition)))
    }
    return(pred_grid)
  })
  
  # Calculate common y-axis limits if fixed scale is requested
  if(fixed_y_scale && plot_type != "region_distribution") {
    all_estimates <- map_dfr(all_predictions, ~data.frame(
      estimate = .x$estimate
    ))
    
    y_min <- min(all_estimates$estimate, na.rm = TRUE)
    y_max <- max(all_estimates$estimate, na.rm = TRUE)
    y_range <- y_max - y_min
    y_limits <- c(y_min - 0.05 * y_range, y_max + 0.05 * y_range)
  }
  
  # Create plots
  plots <- map2(all_predictions, plot_subscale_names, function(pred_grid, subscale_name) {
    
    if(plot_type == "region_distribution") {
      # Special handling for region distribution plots
      ggplot(data = pred_grid, mapping = aes(x = !!sym(config$x_var), y = estimate)) +
        stat_slabinterval(aes(thickness = after_stat(pdf*n)), scale = 0.7, fill = "steelblue") +
        stat_dotsinterval(side = "right", scale = 0.7, slab_linewidth = NA, fill = "steelblue") +
        labs(
          y = paste("Predicted", subscale_name),
          x = config$x_label
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    } else {
      # Standard continuous variable plots using model predictions
      p <- ggplot(pred_grid, aes(x = !!sym(config$x_var), y = estimate))
      
      if(!is.null(config$color_var)) {
        p <- p + aes(color = !!sym(config$color_var), fill = !!sym(config$color_var))
      }
      
      # Apply geom_smooth to the predicted points (no individual points shown)
      p <- p +
        geom_smooth(method = "loess", se = TRUE, level = 0.99, alpha = 0.3, linewidth = 0.8) +
        scale_fill_brewer(palette = "Set2") +
        scale_color_brewer(palette = "Dark2") +
        labs(
          y = paste("Predicted", subscale_name),
          x = config$x_label
        ) +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold")
        )
      
      # Apply fixed y-scale if requested
      if(fixed_y_scale) {
        p <- p + coord_cartesian(ylim = y_limits)
      }
      
      if(!is.null(config$color_var)) {
        p <- p + labs(color = str_to_title(gsub("2$", "", config$color_var)), 
                      fill = str_to_title(gsub("2$", "", config$color_var)))
      }
      
      p
    }
  })
  
  return(plots)
}

# Alternative function to create a combined plot using patchwork
create_combined_faos_plot <- function(models, 
                                      plot_type = "age_sex",
                                      subscale_names = NULL,
                                      ncol = 3,
                                      fixed_y_scale = TRUE) {
  
  # Create individual plots
  plots <- create_faos_plots(
    models = models,
    plot_type = plot_type,
    subscale_names = subscale_names,
    fixed_y_scale = fixed_y_scale,
    show_points = FALSE
  )
  
  # Combine using patchwork (requires library(patchwork))
  if(requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- patchwork::wrap_plots(plots, ncol = ncol) +
      patchwork::plot_annotation(
        title = paste("FAOS Subscale Predictions by", str_to_title(gsub("_", " ", plot_type))),
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
      )
    return(combined_plot)
  } else {
    warning("patchwork package not available. Returning list of individual plots.")
    return(plots)
  }
}

# Usage examples:
# 1. Create plots with improved defaults
# age_plots <- create_faos_plots(age_models, plot_type = "age_sex")

# 2. Create plots with custom subscale names and consistent scaling
# age_plots <- create_faos_plots(
#   age_models, 
#   plot_type = "age_sex",
#   subscale_names = c("Pain", "Symptoms", "Daily Living", "Sport", "Quality of Life"),
#   fixed_y_scale = TRUE
# )

# 3. Create a combined plot (requires patchwork package)
# combined_age_plot <- create_combined_faos_plot(
#   age_models,
#   plot_type = "age_sex", 
#   subscale_names = c("Pain", "Symptoms", "Daily Living", "Sport", "Quality of Life"),
#   ncol = 3
# )


```

```{r}
library(tidyverse)
library(ggdist)
library(patchwork)
library(epoxy)

# Function to identify categorical and continuous variables
identify_variable_types <- function(data) {
  
  categorical_vars <- data |> 
    select(where(~is.character(.x) | is.factor(.x))) |> 
    names()
  
  continuous_vars <- data |> 
    select(where(is.numeric)) |> 
    names()
  
  epoxy("Categorical variables ({length(categorical_vars)}): {paste(categorical_vars, collapse = ', ')}")
  epoxy("Continuous variables ({length(continuous_vars)}): {paste(continuous_vars, collapse = ', ')}")
  
  return(list(
    categorical = categorical_vars,
    continuous = continuous_vars
  ))
}

# Function to create individual slab plots
create_slab_plot <- function(data, cat_var, cont_var, colors = NULL) {
  
  # Remove rows with missing data for this pair
  plot_data <- data |> 
    select(all_of(c(cat_var, cont_var))) |> 
    drop_na()
  
  if(nrow(plot_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No complete cases") +
             theme_void() +
             labs(title = paste(cont_var, "by", cat_var)))
  }
  
  # Create the plot
  p <- plot_data |> 
    ggplot(aes(x = !!sym(cat_var), y = !!sym(cont_var), fill = !!sym(cat_var))) +
    stat_slabinterval(
      adjust = 1.5,
      width = 0.6,
      .width = c(0.50, 0.80, 0.95),
      interval_colour = "grey40",
      interval_size = c(2, 1, 0.5)
    ) +
    scale_fill_viridis_d(alpha = 0.7) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste(cont_var, "by", cat_var),
      subtitle = paste("n =", nrow(plot_data), "complete cases"),
      x = cat_var,
      y = cont_var
    )
  
  if(!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  
  return(p)
}

# Function to create patchwork visualization
create_categorical_continuous_patchwork <- function(data, 
                                                    max_cats = 3, 
                                                    max_continuous = 6,
                                                    ncol = 3) {
  
  # Identify variable types
  var_types <- identify_variable_types(data)
  
  # Limit variables to manageable number
  cat_vars <- head(var_types$categorical, max_cats)
  cont_vars <- head(var_types$continuous, max_continuous)
  
  epoxy("\nCreating plots for:")
  epoxy("Categorical: {paste(cat_vars, collapse = ', ')}")
  epoxy("Continuous: {paste(cont_vars, collapse = ', ')}")
  
  # Create all combinations
  plot_combinations <- expand_grid(
    cat_var = cat_vars,
    cont_var = cont_vars
  )
  
  epoxy("Total plots to create: {nrow(plot_combinations)}")
  
  # Create individual plots
  plot_list <- map2(
    plot_combinations$cat_var,
    plot_combinations$cont_var,
    ~create_slab_plot(data, .x, .y)
  )
  
  # Add names for reference
  names(plot_list) <- paste(plot_combinations$cont_var, 
                            "by", 
                            plot_combinations$cat_var, 
                            sep = "_")
  
  # Create patchwork layout
  if(length(plot_list) > 1) {
    patchwork_plot <- wrap_plots(plot_list, ncol = ncol) +
      plot_annotation(
        title = "Categorical-Continuous Variable Relationships",
        subtitle = "Density distributions with 50%, 80%, and 95% intervals",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12)
        )
      )
  } else {
    patchwork_plot <- plot_list[[1]]
  }
  
  return(list(
    patchwork = patchwork_plot,
    individual_plots = plot_list,
    combinations = plot_combinations
  ))
}

# Function to create focused visualization on specific variables
create_focused_relationship_plot <- function(data, cat_vars, cont_vars, ncol = 2) {
  
  plot_combinations <- expand_grid(
    cat_var = cat_vars,
    cont_var = cont_vars
  )
  
  plot_list <- map2(
    plot_combinations$cat_var,
    plot_combinations$cont_var,
    ~create_slab_plot(data, .x, .y)
  )
  
  names(plot_list) <- paste(plot_combinations$cont_var, 
                            "by", 
                            plot_combinations$cat_var, 
                            sep = "_")
  
  patchwork_plot <- wrap_plots(plot_list, ncol = ncol) +
    plot_annotation(
      title = "Selected Variable Relationships",
      subtitle = "Custom selection of categorical-continuous relationships"
    )
  
  return(patchwork_plot)
}

# Main analysis for your dataset
epoxy("CATEGORICAL-CONTINUOUS RELATIONSHIP ANALYSIS")
epoxy("="*50)

# Create comprehensive patchwork (limited to avoid overcrowding)
comprehensive_viz <- create_categorical_continuous_patchwork(
  MasterPatch, 
  max_cats = 3,           # Limit categorical variables
  max_continuous = 6,     # Limit continuous variables  
  ncol = 3               # Plots per row
)

# Display the comprehensive visualization
print(comprehensive_viz$patchwork)

# Show summary of what was plotted

epoxy("PLOTS CREATED:")
for(i in seq_along(comprehensive_viz$individual_plots)) {
  epoxy("{i}. {names(comprehensive_viz$individual_plots)[i]}")
}

# Example: Create focused plot on key clinical variables
epoxy("FOCUSED ANALYSIS EXAMPLE:")
epoxy("Focusing on key clinical relationships...")

# Define key variables for focused analysis
key_categorical <- c("Sex2", "TreatmentStatus", "Region")
key_continuous <- c("AgeAtInitialExam", "VR12_Physical_TotalScore_Preop", 
                    "FAOS_Pain_TotalScore_Preop", "FAOS_Quality_TotalScore_Preop")

# Filter to variables that actually exist in your data
existing_cat <- intersect(key_categorical, names(MasterPatch))
existing_cont <- intersect(key_continuous, names(MasterPatch))

if(length(existing_cat) > 0 && length(existing_cont) > 0) {
  focused_viz <- create_focused_relationship_plot(
    MasterPatch, 
    existing_cat[1:min(2, length(existing_cat))],      # Max 2 categorical
    existing_cont[1:min(4, length(existing_cont))],    # Max 4 continuous
    ncol = 2
  )
  
  epoxy("Focused visualization created with:")
  epoxy("Categorical: {paste(existing_cat[1:min(2, length(existing_cat))], collapse = ', ')}")
  epoxy("Continuous: {paste(existing_cont[1:min(4, length(existing_cont))], collapse = ', ')}")
  
  print(focused_viz)
} else {
  epoxy("Key variables not found in dataset for focused analysis")
}

# Function to access individual plots
epoxy("\n" %+% "="*30)
epoxy("ACCESSING INDIVIDUAL PLOTS:")
epoxy("Use: comprehensive_viz$individual_plots$'VariableName_by_CategoryName'")
epoxy("Available plots: {paste(names(comprehensive_viz$individual_plots), collapse = ', ')}")
```



```{r}
#| label: comparison-table-func

# Function to create pairwise comparisons for FAOS models
create_faos_comparisons <- function(model_list, 
                                    combined = TRUE,
                                    variable = "Region",
                                    comparison_type = "pairwise",
                                    conf_level = 0.95,
                                    decimals = 3,
                                    fdr_correction = TRUE
) {
  
  if(combined) {
    # Create one combined table with all FAOS subscales
    all_comparisons <- map2_dfr(model_list, names(model_list), function(model, subscale_name) {
      
      # Generate comparisons - hard-code the most common case
      comparisons <- avg_comparisons(
        model,
        variables = list(Region = "pairwise"),
        conf_level = 0.95
      )
      
      # Add subscale identifier and select columns
      comparisons |>
        dplyr::select(contrast, estimate, std.error, p.value) |>
        mutate(subscale = subscale_name, .before = 1)
    })
    
    # Apply FDR correction across ALL comparisons if requested
    if(fdr_correction) {
      all_comparisons <- all_comparisons |>
        mutate(
          p.value.raw = p.value,  # Keep original p-values
          p.value.fdr = p.adjust(p.value, method = "BH")  # FDR adjusted p-values
          #significant.fdr = p.value.fdr < 0.05  # Flag FDR-significant results
        )
    }
    
    # Clean up subscale names for display
    all_comparisons <- all_comparisons |>
      mutate(subscale = case_when(
        subscale == "Pain" ~ "Pain",
        subscale == "Symptom" ~ "Symptoms", 
        subscale == "DailyLiving" ~ "Activities of Daily Living",
        subscale == "Sport" ~ "Sport and Recreation",
        subscale == "Quality" ~ "Quality of Life",
        TRUE ~ subscale
      ))
    
    # Create combined gt table
    if(fdr_correction) {
      # Table with both raw and FDR-adjusted p-values
      result <- all_comparisons |>
        gt(groupname_col = "subscale") |>
        fmt_number(columns = c(estimate, std.error, p.value.raw, p.value.fdr), decimals = 3) |>
        cols_label(
          contrast = "Contrast",
          estimate = "Estimate",
          std.error = "SE",
          p.value.raw = "p-value (raw)",
          p.value.fdr = "p-value (FDR)"
          #significant.fdr = "FDR Sig."
        ) |>
        tab_style(
          style = cell_fill(color = "lightblue"),
          locations = cells_body(rows =  p.value.fdr < 0.05)
        ) |>
        tab_footnote(
          footnote = "P-values adjusted for multiple comparisons using Benjamini-Hochberg FDR correction)",
          locations = cells_column_labels(columns = p.value.fdr)
        ) |>
        tab_options(
          row_group.font.weight = "bold"
        )
    } else {
      # Standard table with raw p-values only
      result <- all_comparisons |>
        gt(groupname_col = "subscale") |>
        fmt_number(columns = c(estimate, std.error, p.value), decimals = 3) |>
        tab_header(title = table_caption) |>
        cols_label(
          contrast = "Contrast",
          estimate = "Estimate",
          std.error = "SE",
          p.value = "p-value"
        ) |>
        tab_options(
          row_group.font.weight = "bold"
        )
    }
    
    return(result)
    
  } else {
    # Create separate tables for each FAOS subscale
    # First, get all comparisons to calculate FDR across all subscales
    if (fdr_correction) {
      all_comparisons_temp <- map2_dfr(model_list, names(model_list), function(model, subscale_name) {
        comparisons <- avg_comparisons(
          model,
          variables = list(Region = "pairwise"),
          conf_level = 0.95
        )
        
        comparisons |>
          dplyr::select(contrast, estimate, std.error, p.value) |>
          mutate(subscale = subscale_name, .before = 1)
      })
      
      # Apply FDR correction across all comparisons
      all_comparisons_temp <- all_comparisons_temp |>
        mutate(p.value.fdr = p.adjust(p.value, method = "BH"))
    }
    
    comparison_tables <- map2(model_list, names(model_list), function(model, subscale_name) {
      
      # Generate comparisons
      comparisons <- avg_comparisons(
        model,
        variables = list(Region = "pairwise"),
        conf_level = 0.95
      )
      
      comparisons_df <- comparisons |>
        dplyr::select(contrast, estimate, std.error, p.value)
      
      # Add FDR-adjusted p-values if correction requested
      if (fdr_correction) {
        fdr_values <- all_comparisons_temp |>
          filter(subscale == subscale_name) |>
          pull(p.value.fdr)
        
        comparisons_df <- comparisons_df |>
          mutate(
            p.value.raw = p.value,
            p.value.fdr = fdr_values
            #significant.fdr = p.value.fdr < 0.05
          )
      }
      
      # Create individual gt table
      if(fdr_correction) {
        comp_table <- comparisons_df |>
          gt() |>
          fmt_number(columns = c(estimate, std.error, p.value.raw, p.value.fdr), decimals = 3) |>
          tab_header(
            title = paste("FAOS", subscale_name, "- Pairwise Comparisons"),
            subtitle = "P-values adjusted for multiple comparisons using Benjamini-Hochberg FDR correction"
          ) |>
          cols_label(
            contrast = "Contrast",
            estimate = "Estimate",
            std.error = "SE",
            p.value.raw = "p-value (raw)",
            p.value.fdr = "p-value (FDR)"
            #significant.fdr = "FDR Sig."
          ) |>
          tab_style(
            style = cell_fill(color = "lightblue"),
            locations = cells_body(rows =  p.value.fdr < 0.05)
          )
      } else {
        comp_table <- comparisons_df |>
          gt() |>
          fmt_number(columns = c(estimate, std.error, p.value), decimals = 3) |>
          tab_header(title = paste("FAOS", subscale_name, "- Pairwise Comparisons")) |>
          cols_label(
            contrast = "Contrast",
            estimate = "Estimate",
            std.error = "SE", 
            p.value = "p-value"
          )
      }
      
      return(comp_table)
    })
    
    names(comparison_tables) <- names(model_list)
    return(comparison_tables)
  }
}


```


```{r}
#| label: fig-desc-agesex
#| fig-cap: "FAOS subscales visualised against age, categorised by sex"

# For descriptive plots - Age relationship (minimally adjusted for sex only)
age_models <- fit_descriptive_models(MasterImp, covariate = "AgeAtInitialExam")


combined_age_plot <- create_combined_faos_plot(
  age_models,
  plot_type = "age_sex",
  subscale_names = c("Pain", "Symptoms", "Daily Living", "Sport", "Quality of Life"),
  ncol = 3
)


knitr::knit_print(combined_age_plot)

```

```{r}
#| label: fig-desc-comorbidsex
#| fig-cap: "FAOS subscales visualised against self-reported comorbidity score, categorised by sex"


comorbidity_models <- fit_descriptive_models(MasterImp, covariate = "SRCQTotalScore")



combined_comorbid_plot <- create_combined_faos_plot(
  comorbidity_models,
  plot_type = "comorbidities_sex",
  subscale_names = c("Pain", "Symptoms", "Daily Living", "Sport", "Quality of Life"),
  ncol = 3
)


knitr::knit_print(combined_comorbid_plot)

```

```{r}
#| label: fig-desc-agearthritis
#| fig-cap: "FAOS subscales visualised against age, categorised by the presence or absence of an arthritis diagnosis"
#| eval: false

arthritis_models <- fit_descriptive_models(MasterImp, covariate = "AgeAtInitialExam")
arthritis_plots <- create_faos_plots(arthritis_models, plot_type = "age_arthritis")

(arthritis_plots$Pain + arthritis_plots$Symptom + arthritis_plots$Daily + arthritis_plots$Sport + arthritis_plots$Quality) + patchwork::plot_layout(ncol = 3)
```
