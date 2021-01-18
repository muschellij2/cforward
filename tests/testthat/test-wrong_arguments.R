testthat::test_that("Passing in weights and data arguments", {
  variables = c("gender",
                "age_years_interview", "education_adult")


  res = cforward(nhanes_example,
                 event_time = "event_time_years",
                 event_status = "mortstat",
                 weight_column = "WTMEC4YR_norm",
                 variables = variables,
                 included_variables = NULL,
                 n_folds = 1, # when fold is 1
                 seed = 1989,
                 max_model_size = 2,
                 verbose = FALSE)
  # should be the same - no test/train split
  testthat::expect_equal(
    res[[1]]$full_concordance,
    c(res[[1]]$cv_concordance)
  )

  testthat::expect_warning({
    res = cforward(nhanes_example,
                   event_time = "event_time_years",
                   event_status = "mortstat",
                   weight_column = "WTMEC4YR_norm",
                   variables = variables,
                   included_variables = NULL,
                   n_folds = 2,
                   weights = 5, # weights are specific
                   seed = 1989,
                   max_model_size = 50,
                   verbose = FALSE)
  }, regexp = "were specified")


  testthat::expect_error({
    res = cforward(nhanes_example,
                   event_time = "event_time_years",
                   event_status = "mortstat",
                   weight_column = "WTMEC4YR_norm",
                   variables = "blahblahbald",
                   included_variables = NULL,
                   n_folds = 2,
                   weights = 5,
                   seed = 1989,
                   max_model_size = 50,
                   verbose = FALSE)
  }, regexp = "Independent variables")


  testthat::expect_warning({
    estimate_concordance(nhanes_example, nhanes_example,
                         data = nhanes_example,
                         all_variables = "1")
  }, regexp = "data set in.*overridden.*")


  x = nhanes_example
  x$strata = rbinom(nrow(x), size = 1, prob = 0.2)
  res = cforward(x,
                 event_time = "event_time_years",
                 event_status = "mortstat",
                 weight_column = "WTMEC4YR_norm",
                 variables = variables,
                 included_variables = NULL,
                 n_folds = 2,
                 cfit_args = list(strata_column = "strata"),
                 seed = 1989,
                 max_model_size = 50,
                 verbose = FALSE)


  x$gender[1] = NA
  testthat::expect_warning({
    res = cforward(x,
                   event_time = "event_time_years",
                   event_status = "mortstat",
                   weight_column = "WTMEC4YR_norm",
                   variables = variables,
                   included_variables = NULL,
                   n_folds = 2,
                   seed = 1989,
                   max_model_size = 50,
                   verbose = FALSE)
  }, regexp = "elements in the test")

})
