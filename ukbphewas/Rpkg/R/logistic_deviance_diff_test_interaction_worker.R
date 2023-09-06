# logistic_deviance_diff_test_worker for interaction test

library(fastglm)

args <- commandArgs(trailingOnly = TRUE)

worker_id <- args[1]
dealer_addr <- args[2]
monitor_addr <- args[3]


fit_fastglm <- function(data_matrix) {
  fastglm(
    # We add the intercept to the x matrix.
    x = cbind(1, data_matrix[, -1]),
    y = data_matrix[, 1],
    family = binomial(),
    method = 2
  )
}


fastglm_anova <- function(base, augmented) {
  d_dev <- base$deviance - augmented$deviance
  d_df <- base$df.residual - augmented$df.residual
  p <- 1 - pchisq(d_dev, d_df)

  # Results
  list(
    resid_deviance_base = base$deviance,
    resid_deviance_augmented = augmented$deviance,
    deviance = d_dev,
    p = p
  )
}


logistic_deviance_diff_test_worker <- function(worker_id, ...) {
  # Read the analysis configuration object.
  conf <- fromJSON(file = "analysis_configuration.json")

  # Read the covariables.
  covars <- get_xs(conf)

  # Read the XPCs.
  xpcs <- read.csv(
    conf$binary_conf$xpcs_path,
    colClasses = read_csv_filter_columns(
      conf$binary_conf$xpcs_path,
      c("sample_id", conf$binary_conf$augmented_variables)
    )
  )
  covars <- merge(covars, xpcs, by = "sample_id")

  # For now, sample_id needs to be a string. Eventually, we will infer this
  # better.
  covars$sample_id <- as.character(covars$sample_id)

  # Create a header file.
  cat(paste0("variable_id,analysis_type,sex_subset,",
             "n_cases,n_controls,",
             "n_excluded_from_controls,",
             "deviance_base,n_params_base,",
             "deviance_augmented,n_params_augmented,",
             "deviance_diff,p,",
             "deviance_inter,n_params_inter,",
             "deviance_diff_inter,p_inter\n"),
      file = "header_summary.csv")

  cat("R: Opened output file.\n")
  output_file <- file(
    paste0("results_worker_", worker_id, "_summary.csv"), open = "wt"
  )

  model_file <- file(
    paste0("results_worker_", worker_id, "_model.json"), open = "wt"
  )

  # Check if sex is in the independant variables.
  # This is important later when the variability in the sex variable is
  # tested after filtering.
  base_cols <- char_to_terms(conf$model_rhs)
  aug_cols <- conf$binary_conf$augmented_variables
  aug_int_cols <- conf$binary_conf$augmented_interaction_variables

  # The final columns (for the two last tests)
  base_aug_cols <- c(base_cols, aug_cols)
  base_aug_int_cols <- c(base_aug_cols, aug_int_cols)

  unique_inter_cols <- unique(
    purrr::flatten(stringr::str_split(aug_int_cols, "[*]"))
  )

  # Callback for when data is to be processed from the queue.
  do.work <- function(metadata, data) {
    data <- deserialize(data)

    cases <- as.character(data[data[, "y"] == 1, "sample_id"])
    to_exclude <- as.character(data[is.na(data[, "y"]), "sample_id"])

    covars$y <- as.numeric(covars$sample_id %in% cases)
    covars[covars$sample_id %in% to_exclude, "y"] <- NA

    n_excluded <- sum(is.na(covars[, "y"]))

    aug_model_rhs <- paste(base_aug_cols, collapse = " + ")
    aug_formula <- as.formula(paste0(
      "y ~ ", aug_model_rhs
    ))
    int_model_rhs <- paste(base_aug_int_cols, collapse = " + ")
    int_formula <- as.formula(paste0(
      "y ~ ", int_model_rhs
    ))
     # Augmented
    aug_data_matrix <- model.frame(aug_formula, data = covars)
    keep <- complete.cases(aug_data_matrix)

    aug_data_matrix <- aug_data_matrix[keep, ]
    col_drop_li <- drop_columns_with_no_variance(aug_data_matrix)
    aug_data_matrix <- col_drop_li$mat
    dropped_cols <- col_drop_li$dropped_cols

    if (length(dropped_cols) != 0) {
      # Did we drop a required column for interaction?
      if (length(intersect(dropped_cols, unique_inter_cols)) > 0) {
        return()
      }
      base_cols <- base_cols[!(base_cols %in% dropped_cols)]
      base_aug_cols <- base_aug_cols[!(base_aug_cols %in% dropped_cols)]
    }

    n_cases <- sum(aug_data_matrix[, "y"] == 1)
    n_controls <- sum(aug_data_matrix[, "y"] == 0)

    if (n_cases < conf$binary_conf$min_num_cases) {
      return()
    }

     # Interaction
    int_data_matrix <- model.frame(int_formula, data = covars)
    keep <- complete.cases(int_data_matrix)

    int_data_matrix <- int_data_matrix[keep, ]
    col_drop_li <- drop_columns_with_no_variance(int_data_matrix)
    int_data_matrix <- col_drop_li$mat
    dropped_cols <- col_drop_li$dropped_cols

    if (length(dropped_cols) != 0) {
      aug_cols <- aug_cols[!(aug_cols %in% dropped_cols)]
      aug_int_cols <- aug_int_cols[!(aug_int_cols %in% dropped_cols)]
    }

    if (n_cases != sum(int_data_matrix[, "y"] == 1)) {
      print("N cases not equal")
      return()
    }
    if (n_controls != sum(int_data_matrix[, "y"] == 0)) {
      print("N controls not equal")
      return()
    }

    # Creating the matrices for the base and augmented models
    aug_data_matrix <- as.matrix(aug_data_matrix)
    fit_base <- fit_fastglm(aug_data_matrix[, c("y", base_cols)])
    fit_aug <- fit_fastglm(aug_data_matrix)

    # Creating the matrix for the interaction model (adding the interaction
    # columns)
    for (expression in aug_int_cols) {
      int_data_matrix[gsub("[*]", ":", expression)] <- eval(
        parse(text = expression),
        envir = int_data_matrix
      )
    }
    int_data_matrix <- as.matrix(int_data_matrix)
    fit_int <- fit_fastglm(int_data_matrix)

    # Anova
    lrt_aug <- fastglm_anova(fit_base, fit_aug)
    lrt_inter <- fastglm_anova(fit_aug, fit_int)

    line <- paste(
      metadata$variable_id,
      metadata$analysis_type,
      ifelse(is.null(metadata$sex_subset), "BOTH", metadata$sex_subset),
      n_cases,
      n_controls,
      n_excluded,
      lrt_aug$resid_deviance_base,
      length(base_cols) + 1,          # Add 1 for intercept.
      lrt_aug$resid_deviance_augmented,
      length(base_aug_cols) + 1,
      lrt_aug$deviance,
      lrt_aug$p,
      lrt_inter$resid_deviance_augmented,
      length(base_aug_int_cols) + 1,  # Add 1 for intercept
      lrt_inter$deviance,
      lrt_inter$p,
      sep = ","
    )

    writeLines(line, con = output_file)

    # We will write the model coefficients in case they are needed.
    infer_df_aug <- data.frame(
      variable = names(fit_aug$coefficients),
      beta = fit_aug$coefficients,
      se = fit_aug$se,
      z = fit_aug$coefficients / fit_aug$se
    )

    infer_df_int <- data.frame(
      variable = names(fit_int$coefficients),
      beta = fit_int$coefficients,
      se = fit_int$se,
      z = fit_int$coefficients / fit_int$se
    )

    infer_df_aug$p <- 2 * pnorm(-abs(infer_df_aug$z))
    infer_df_aug$nlog10p <- -2 * pnorm(-abs(infer_df_aug$z), log.p = TRUE) / log(10)
    infer_df_aug[1, "variable"] <- "intercept"

    infer_df_int$p <- 2 * pnorm(-abs(infer_df_int$z))
    infer_df_int$nlog10p <- -2 * pnorm(-abs(infer_df_int$z), log.p = TRUE) / log(10)
    infer_df_int[1, "variable"] <- "intercept"

    cat(toJSON(
      list(
        variable_id = metadata$variable_id,
        analysis_type = metadata$analysis_type,
        model_aug_fit = infer_df_aug,
        model_inter_fit = infer_df_int
      )
    ), file = model_file)
    cat("\n", file = model_file)

  }

  # Start the main loop.
  Worker(worker_id, ..., callback = do.work)

  cat("R: closing output file.\n")
  close(output_file)

}

logistic_deviance_diff_test_worker(worker_id, dealer_addr, monitor_addr)
