# FUNCTION TO GET 'BASIS' MATRIX OF TIME ACCRUED IN EACH ARM ------------

#' Generate Piecewise Constant Cumulative Hazard Matrix
#'
#' This helper function takes in a vector of event times, and
#'   generates a matrix dividing each time into a vector giving the time accrued within each
#'   interval between consecutive elements of a vector of knots.
#'
#' @param y vector of event times.
#' @param knots increasing vector of cutpoints
#'
#' @return a numeric matrix, with rows corresponding to elements of y.
#' @export
pw_cum_mat <- function(y, knots){
  # browser()

  if(knots[1] != 0){knots <- c(0,knots)}
  knots_diff <- diff(knots)

  #count of the number of intervals in each list
  num_int <- c(length(knots))
  n <- length(y)

  #matrix with each row being an observation, and each column being an interval, capturing how much time is accrued in each interval
  knots_mat <- matrix(c(knots_diff,0),nrow=n,ncol=num_int,byrow = TRUE)

  #an n length vector of integers indicating which interval the observation is in
  cut_cats <- findInterval(x = y, vec = knots)

  #an n length vector capturing the residual time accrued in the final interval of each observation
  y_last <- y-knots[cut_cats]

  #a loop through each observation to finalize the matrix of time intervals
  for(i in 1:n){
    knots_mat[i,cut_cats[i]] <- y_last[i]
    knots_mat[i,-c(1:cut_cats[i])] <- 0
  }
  return(knots_mat)
}
