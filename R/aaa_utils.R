check_data = function(data,
                      event_time = "event_time_years",
                      event_status = "mortstat",
                      weight_col = "WTMEC4YR_norm",
                      ind_vars = NULL) {
  stopifnot(event_status %in% colnames(data),
            event_time %in% colnames(data),
            weight_col %in% colnames(data),
            !is.null(ind_vars))
  # !!
  # stopifnot(is.list(ind_vars))
  uindvars = unique(unlist(ind_vars))
  sd = setdiff(uindvars, colnames(data))
  if (length(sd) > 0) {
    msg = paste0("Independent variables not in data: ",
                 paste(sd, collapse = ", "))
    stop(msg)
  }
}


choose_next_variable = function(cv_concordance) {
  ind = sort(cv_concordance, decreasing = TRUE)
  stopifnot(length(which(cv_concordance == ind[1])) <= 1)
  best_var  = ind[1]
  best_var
}
