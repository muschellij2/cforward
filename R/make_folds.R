#' @rdname cforward
#' @export
make_folds = function(data,
                      event_status = "mortstat",
                      n_folds = 10,
                      verbose = TRUE){
  fold_number = NULL
  rm(list="fold_number")
  if (!"fold_number" %in% colnames(data)) {
    if (verbose) {
      message(paste0("Creating ", n_folds, " Folds"))
    }
    data = as.data.frame(data)
    ss = split(data, data[[event_status]])
    ss = lapply(ss, function(x) {
      folds = rep(1:n_folds, length.out = ceiling(nrow(x)/n_folds)*n_folds)
      folds = sample(folds)
      x$fold_number = folds[seq(nrow(x))]
      x
    })
    data = unsplit(ss, data[[event_status]])
    data = tibble::as_tibble(data)
  } else {
    msg = "fold_number in data, using these folds"
    if (verbose) {
      message(msg)
    }
  }
  stopifnot(!anyNA(data$fold_number))
  if (verbose) {
    data %>%
      dplyr::count(fold_number) %>%
      print()
  }
  data
}
