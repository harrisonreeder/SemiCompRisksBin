#' Function to return quadrature points and weights
#'
#' This function stores quadrature points used for numerical integration by
#' Gauss-Hermite rule. The values are pre-computed and stored for
#' number of nodes equal to 7, 11, 15, 21, 31, 41, 51, or 61.
#'
#' @param n_quad number of quadrature points
#' @return A list with two named elements (\code{points} and \code{weights}), each
#'   of which is a numeric vector with length equal to the number of
#'   quadrature nodes
get_ghquad_pointsweights <- function(n_quad) {
  if (!is.numeric(n_quad) || (length(n_quad) > 1L)) {
    stop("'n_quad' should be a numeric vector of length 1.")
  }
  if (n_quad == 15) {
    list( points = c(
      -4.499990707309391553664, -3.669950373404452534729, -2.967166927905603248489,
      -2.325732486173857745454, -1.719992575186488932416, -1.136115585210920666319,
      -0.565069583255575748526, 0, 0.565069583255575748526,
      1.136115585210920666319, 1.719992575186488932416, 2.32573248617385774545,
      2.967166927905603248489, 3.669950373404452534729, 4.499990707309391553664),
      weights = c(
        1.522475804253517020161E-9, 1.059115547711066635775E-6,
        1.00004441232499868127E-4, 0.002778068842912775896079,
        0.03078003387254608222868, 0.1584889157959357468838,
        0.4120286874988986270259, 0.5641003087264175328526,
        0.4120286874988986270259, 0.1584889157959357468838,
        0.03078003387254608222868, 0.00277806884291277589608,
        1.00004441232499868127E-4, 1.059115547711066635775E-6,
        1.52247580425351702016E-9) )
  } else{"only 15 gauss-hermite nodes available for now."}
}
