# Generates a D-Vine Matrix of dimension d. An order can be specified with the elements argument.
# Author: Daniel Kraus

DVineMatGen <- function(d=length(elements),elements=1:d){
    if (!setequal(1:d,elements))
        stop("'elements' has to consist of the numbers from 1 to d.")
    Mat <- diag(elements)
    for (i in 2:d) {
        Mat[d + 2 - i, 1:(d + 1 - i)] <- elements[i:d]
    }

    Mat
}
