# Generates a D-Vine Matrix of dimension d. An order can be specified with the elements argument.
# Author: Daniel Kraus
gen_dvine_mat <- function(d = length(elements), elements = seq_len(d)){
    if (!setequal(seq_len(d), elements))
        stop("'elements' must consist of the numbers from 1 to d.")
    mat <- diag(elements)
    for (i in seq_len(d - 1) + 1) {
        mat[d + 2 - i, 1:(d + 1 - i)] <- elements[i:d]
    }

    # reverse for vinecopulib notation
    mat[d:1, ]
}
