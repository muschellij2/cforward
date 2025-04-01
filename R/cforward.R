

#' Estimate Out-of-Sample Concordance
#'
#' @param train A data set to perform model training.
#' @param test A data set to estimate concordance, from fit model with `train`.
#' Set to `train` if estimating on the same data
#' @param event_time Character vector of length 1 with event times, passed to
#' \code{\link[survival]{Surv}}
#' @param event_status Character vector of length 1 with event status, passed to
#' \code{\link[survival]{Surv}}
#' @param weight_column Character vector of length 1 with weights for
#' model.  If no weights are available, set to `NULL`
#' @param all_variables Character vector of variables to put in the
#' model.  All must be in `data`.
#' @param cfit_args Arguments passed to \code{\link[survival]{concordancefit}}.  If
#' `strata` is to be passed, set `strata_column` in this list.
#' @param ... Additional arguments to pass to \code{\link[survival]{coxph}}
#'
#' @return A list of concordance and the model fit with the training data
#' @export
estimate_concordance = function(
    train,
    test = train,
    event_time = "event_time_years",
    event_status = "mortstat",
    weight_column = "WTMEC4YR_norm",
    all_variables = NULL,
    cfit_args = list(),
    ...
) {
    form   <- paste(all_variables, collapse=" + ")
    y = paste0("survival::Surv(", event_time, ", ", event_status, ") ~ ")
    args = list(...)
    if ("data" %in% names(args)) {
        warning("data set in ... for estimate_concordance, will be overridden")
    }
    args$formula = as.formula(paste0(y, form))
    args$data = train
    if ("weights" %in% names(args)) {
        warning(paste0("Weights were specified in argument, ",
                       "overridden by weight_column"))
    }
    args$weights = train[[weight_column]]

    fit_k <- do.call(survival::coxph, args = args)
    eta_k <- predict(fit_k, type='lp', newdata=test)

    # allows for "1" for all_variables
    check_vars = c(all_variables, event_status, event_time)
    check_vars = intersect(check_vars, colnames(test))
    if (anyNA(test[check_vars])) {
        warning(paste0("Missing elements in the test data! ",
                       "Please remove prior to running"))
    }
    cfit_args$y = survival::Surv(test[[event_time]],
                                 test[[event_status]])
    cfit_args$x = eta_k
    cfit_args$weights = test[[weight_column]]
    cfit_args$reverse = TRUE
    if ("strata_column" %in% names(cfit_args) &&
        !"strata" %in% names(cfit_args)) {
        cfit_args$strata = test[[cfit_args$strata_column]]
        cfit_args$strata_column = NULL
    }
    checker = function(x, n) {
        if (is.null(x)) {
            x = rep(TRUE, n)
        } else {
            x = !is.na(x)
        }
        x
    }
    n = length(cfit_args$x)
    keep = checker(cfit_args$x, n ) &
        checker(cfit_args$y, n) &
        checker(cfit_args$weights, n) &
        checker(cfit_args$strata, n)
    cfit_args$y = cfit_args$y[keep]
    cfit_args$x = cfit_args$x[keep]
    cfit_args$weights = cfit_args$weights[keep]
    cfit_args$strata = cfit_args$strata[keep]

    cfit = do.call(survival::concordancefit, args = cfit_args)
    stopifnot(!is.null(cfit))
    fit_k$residuals = fit_k$weights = fit_k$linear.predictors = NULL
    fit_k$call$data = fit_k$y = NULL
    fit_k$call$weights = fit_k$call$strata = NULL
    fit_k$call[[1]] = as.name("survival::coxph")
    list(concordance = cfit$concordance,
         model = fit_k)
}




#' Forward Selection Based on C-Index/Concordance
#'
#' @param data A data set to perform model selection and cross-validation.
#' @param event_time Character vector of length 1 with event times, passed to
#' \code{\link[survival]{Surv}}
#' @param event_status Character vector of length 1 with event status, passed to
#' \code{\link[survival]{Surv}}
#' @param weight_column Character vector of length 1 with weights for
#' model.  If no weights are available, set to `NULL`
#' @param variables Character vector of variables to perform selection.
#' Must be in `data`.
#' @param included_variables Character vector of variables
#' forced to have in the model.  Must be in `data`
#' @param n_folds Number of folds for Cross-validation.  If you want to run on
#' the full data, set to 1
#' @param seed Seed set before folds are created.
#' @param max_model_size maximum number of variables in the model.  Selection
#' will stop if reached. Note, this does not correspond to the number
#' of coefficients, due to categorical variables.
#' @param verbose print diagnostic messages
#' @param cfit_args Arguments passed to \code{\link[survival]{concordancefit}}.  If
#' `strata` is to be passed, set `strata_column` in this list.
#' @param ... Additional arguments to pass to \code{\link[survival]{coxph}}
#' @param save_memory save only a minimal amount of information, discard
#' the fitted models
#' @param c_threshold threshold for concordance.  If the difference in the best
#' concordance and this one does not reach a certain threshold, break.
#'
#' @return A list of lists, with elements of:
#' \describe{
#' \item{full_concordance}{Concordance when fit on the full data}
#' \item{models}{Cox model from full data set fit, stripped of large memory
#' elements}
#' \item{cv_concordance}{Cross-validated Concordance}
#' \item{included_variables}{Variables included in the model, other than
#' those being selection upon}
#' }
#' @export
#' @importFrom stats as.formula predict
#'
#' @examples
#' variables = c("gender",
#'               "age_years_interview", "education_adult")
#'
#' res = cforward(nhanes_example,
#'                event_time = "event_time_years",
#'                event_status = "mortstat",
#'                weight_column = "WTMEC4YR_norm",
#'                variables = variables,
#'                included_variables = NULL,
#'                n_folds = 5,
#'                c_threshold = 0.02,
#'                seed = 1989,
#'                max_model_size = 50,
#'                verbose = TRUE)
#' conc = sapply(res, `[[`, "best_concordance")
#'
#'
#'
#' res = cforward(nhanes_example,
#'                event_time = "event_time_years",
#'                event_status = "mortstat",
#'                weight_column = "WTMEC4YR_norm",
#'                variables = variables,
#'                included_variables = NULL,
#'                n_folds = 5,
#'                seed = 1989,
#'                max_model_size = 50,
#'                verbose = TRUE)
#' conc = sapply(res, `[[`, "best_concordance")
#' threshold = 0.01
#' included_variables = names(conc)[c(1, diff(conc)) > threshold]
#'
#' new_variables = c("diabetes", "stroke")
#' second_level = cforward(nhanes_example,
#'                event_time = "event_time_years",
#'                event_status = "mortstat",
#'                weight_column = "WTMEC4YR_norm",
#'                variables = new_variables,
#'                included_variables = included_variables,
#'                n_folds = 5,
#'                seed = 1989,
#'                max_model_size = 50,
#'                verbose = TRUE)
#' second_conc = sapply(second_level, `[[`, "best_concordance")
#' result = second_level[[which.max(second_conc)]]
#' final_model = result$models[[which.max(result$cv_concordance)]]
cforward = function(
    data,
    event_time = "event_time_years",
    event_status = "mortstat",
    weight_column = "WTMEC4YR_norm",
    variables = NULL,
    included_variables = NULL,
    n_folds = 10,
    seed = 1989,
    max_model_size = 50,
    c_threshold = NULL,
    verbose = TRUE,
    cfit_args = list(),
    save_memory = FALSE,
    ...) {

    all_variables = c(variables, included_variables)
    check_data(
        data,
        event_time,
        event_status,
        weight_column,
        all_variables
    )
    uindvars = unique(unlist(all_variables))

    ######################################
    # setting seed outside of folds in case it affects
    # coxph or something else
    ######################################
    set.seed(seed)
    data = make_folds(data, event_status, n_folds, verbose = verbose)

    ######################################
    # check folds
    ######################################
    u_folds = sort(unique(data$fold_number))
    if (n_folds != length(u_folds)) {
        warning("n_folds not equal to number of folds")
    }
    n_folds = length(u_folds)

    ######################################
    # Another data check- no NA
    ######################################
    all_vars = c(uindvars, event_status, event_time, weight_column,
                 "fold_number")
    if (anyNA(data[all_vars])) {
        warning("NA elements to data")
    }

    ### Vector to multiply vector of k-fold cross validated weights by to
    ## get an estimate of cross-valided AUC.
    ## This accounts for potentially unequal sizes in the test datasets
    # k_id <- vapply(inx_ls, length, numeric(1))/nrow(df_analysis)
    # this is fine because there are no NAs
    # so they shouldn't be dropped from model

    all_lists = NULL
    current_concordance = 0
    while(length(variables) > 0) {
        L = cforward_one(data,
                         event_time = event_time,
                         event_status = event_status,
                         weight_column = weight_column,
                         variables = variables,
                         included_variables = included_variables,
                         verbose = verbose,
                         cfit_args = cfit_args,
                         save_memory = save_memory,
                         ...)
        next_var = choose_next_variable(L$cv_concordance)
        L$best_concordance = next_var
        if (save_memory) {
            L$model = NULL
        }
        next_var = names(next_var)
        all_lists = c(all_lists, list(L))
        included_variables = c(included_variables, next_var)
        variables = setdiff(variables, included_variables)
        if (length(included_variables) >= max_model_size) {
            break
        }
        if (verbose) {
            message(paste0("next variable: ", next_var))
            message(paste0("number of variables: ", length(included_variables),
                           "\n"))
        }
        if (!is.null(c_threshold) &&
            (L$best_concordance - current_concordance) < c_threshold) {
            break
        }
        current_concordance = L$best_concordance

    }
    # names(all_lists) = paste0("level_", seq_along(all_lists))
    return(all_lists)
}




#' @rdname cforward
#' @export
cforward_one = function(
    data,
    event_time = "event_time_years",
    event_status = "mortstat",
    weight_column = "WTMEC4YR_norm",
    variables,
    included_variables = NULL,
    verbose = TRUE,
    cfit_args = list(),
    save_memory = FALSE,
    ...) {

    fold_number = n = NULL
    rm(list=c("fold_number", "n"))
    variables = unique(variables)
    # don't add it if already there
    isect = intersect(variables, included_variables)
    stopifnot(length(isect) == 0)

    u_folds = sort(unique(data$fold_number))
    n_folds = length(u_folds)


    ## create empty matrices/vectors to store in-sample and cross validated C-index results
    ##   C_IS_ij Store the in-sample C-index
    concordance_full  <- rep(NA, length(variables))
    names(concordance_full) = variables
    ##   Matrix to store the cross-validated C-index at each of the k-folds
    concordance_cross_mat <- matrix(NA, nrow = length(variables), ncol = length(u_folds))
    colnames(concordance_cross_mat) = u_folds
    rownames(concordance_cross_mat) = variables

    k_id = data %>%
        dplyr::count(fold_number) %>%
        dplyr::arrange(fold_number) %>%
        dplyr::pull(n)
    k_id = k_id / nrow(data)
    names(k_id) = u_folds
    stopifnot(length(k_id) == length(u_folds))

    all_models = vector(mode = "list",
                        length = length(variables))
    names(all_models) = variables
    if (save_memory) {
        all_models = NULL
    }
    i = 1
    ## loop over the number independent variables
    for (i in seq_along(variables)) {
        ivar = variables[[i]]

        run_vars = unique(c(ivar, included_variables))
        full_conc = estimate_concordance(
            data,
            data,
            event_time = event_time,
            event_status = event_status,
            weight_column = weight_column,
            all_variables = run_vars,
            cfit_args = cfit_args,
            ...)
        concordance_full[[ivar]] = full_conc$concordance
        if (!save_memory) {
            all_models[[ivar]] = full_conc$model
        }

        for (k in u_folds){
            xdf_train <- data %>%
                dplyr::filter(fold_number != k)
            xdf_test <- data %>%
                dplyr::filter(fold_number == k)
            # if you're doing it on the full data
            if (n_folds == 1 && nrow(xdf_train) == 0) {
                xdf_train = xdf_test
            }

            conc = estimate_concordance(
                xdf_train,
                xdf_test,
                event_time = event_time,
                event_status = event_status,
                weight_column = weight_column,
                all_variables = run_vars,
                cfit_args = cfit_args,
                ...)

            concordance_cross_mat[i, as.character(k)]  <- conc$concordance
        }
        if (verbose) {
            cc = concordance_cross_mat[i, ] %*% k_id
            cc = c(cc)
            message(paste0("Concordance for ",
                           ivar, ": ", round(cc, 3)))
        }
    }

    ## average C-index over the k-folds
    concordance_cross  <- concordance_cross_mat %*% k_id
    names(concordance_cross) = variables

    L = list(
        full_concordance = concordance_full)
    L$models = all_models
    L$cv_concordance = concordance_cross
    L$included_variables = included_variables
    return(L)
}

