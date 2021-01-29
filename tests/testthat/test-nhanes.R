testthat::context("NHANES Example")
testthat::test_that("Example for cforward", {
  variables = c("gender",
                "age_years_interview", "education_adult")

  res = cforward(nhanes_example,
                 event_time = "event_time_years",
                 event_status = "mortstat",
                 weight_column = "WTMEC4YR_norm",
                 variables = variables,
                 included_variables = NULL,
                 n_folds = 5,
                 seed = 1989,
                 max_model_size = 50,
                 verbose = TRUE)

  testthat::expect_equal(
    res[[1]]$full_concordance,
    c(gender = 0.551197833594443,
      age_years_interview = 0.626403395590797,
      education_adult = 0.5687432356904))

  conc = sapply(res, `[[`, "best_concordance")
  testthat::expect_equal(
    conc,
    c(age_years_interview = 0.62884268049712, gender = 0.642862998045428,
      education_adult = 0.651845223650031)
  )
  threshold = 0.01
  included_variables = names(conc)[c(1, diff(conc)) > threshold]
  testthat::expect_equal(
    included_variables,
    c("age_years_interview", "gender")
  )

  new_variables = c("diabetes", "stroke")
  second_level = cforward(nhanes_example,
                          event_time = "event_time_years",
                          event_status = "mortstat",
                          weight_column = "WTMEC4YR_norm",
                          variables = new_variables,
                          included_variables = included_variables,
                          n_folds = 5,
                          seed = 1989,
                          max_model_size = 50,
                          verbose = TRUE)
  second_conc = sapply(second_level, `[[`, "best_concordance")
  testthat::expect_equal(
    second_conc,
    c(diabetes = 0.688798272969669, stroke = 0.695749500392469)
  )
  result = second_level[[which.max(second_conc)]]
  final_model = result$models[[which.max(result$cv_concordance)]]

  testthat::expect_equal(
    coef(final_model),
    c(strokeYes = 0.755580151540889, age_years_interview = 0.0527018992699052,
      genderFemale = -0.455832074189688, diabetesYes = 0.991726669674368
    )
  )

})


testthat::test_that("Variables and included have overlap", {
  variables = c("gender",
                "age_years_interview", "education_adult")

  testthat::expect_error({
    res = cforward(nhanes_example,
                   event_time = "event_time_years",
                   event_status = "mortstat",
                   weight_column = "WTMEC4YR_norm",
                   variables = variables,
                   included_variables = variables,
                   n_folds = 5,
                   seed = 1989,
                   max_model_size = 50,
                   verbose = FALSE)
  }, regexp = "isect")

  x = nhanes_example
  x = make_folds(x)
  testthat::expect_warning({
    testthat::expect_message({
      res = cforward(x,
                     event_time = "event_time_years",
                     event_status = "mortstat",
                     weight_column = "WTMEC4YR_norm",
                     variables = variables,
                     included_variables = NULL,
                     n_folds = 5,
                     seed = 1989,
                     max_model_size = 50,
                     verbose = TRUE)
    }, regexp = "in data")
  }, regexp = "n_folds not ")

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
  }, regexp = "not equal")
})


testthat::test_that("Restrict maximum size", {
  variables = c("gender",
                "age_years_interview", "education_adult", "diabetes", "stroke")
  small_model = cforward(nhanes_example,
                         event_time = "event_time_years",
                         event_status = "mortstat",
                         weight_column = "WTMEC4YR_norm",
                         variables = variables,
                         included_variables = NULL,
                         n_folds = 5,
                         seed = 1989,
                         max_model_size = 3,
                         verbose = TRUE)
  testthat::expect_length(small_model, 3L)
})



testthat::test_that("Restrict c_threshold", {
  variables = c("gender",
                "age_years_interview", "education_adult")
  res = cforward(nhanes_example,
                 event_time = "event_time_years",
                 event_status = "mortstat",
                 weight_column = "WTMEC4YR_norm",
                 variables = variables,
                 included_variables = NULL,
                 n_folds = 5,
                 c_threshold = 0.02,
                 seed = 1989,
                 save_memory = TRUE,
                 max_model_size = 50,
                 verbose = TRUE)
  testthat::expect_length(res, 2L)
})
