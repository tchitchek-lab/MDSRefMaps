#' @title Setting of the algorithm constants
#' 
#' @description
#' This function is used to set the different constants of the MDSRefMaps algorithm.
#'
#' @details
#' The 'ACC_THRESHOLD' parameter can be used to contrain the object accelerations.
#'
#' @param K a numeric value indicating the spring strength between two particles
#' @param F a numeric value indicating the friction of the springs
#' @param DELTA_T a numeric value indicating the time between two steps of the algorithm
#' @param MASS a numeric value indicating the mass of each particle
#' @param ACC_THRESHOLD a numeric value indicating the acceleration threshold (a value of 0 disactivate this acceleration threshold)
#' 
#' @export
setConstants = function(K=1, F=0.1, DELTA_T=0.001, MASS=10, ACC_THRESHOLD=0) {
    if (K < 0)             {stop("the numeric parameter K must be >= 0")}
    if (F < 0)             {stop("the numeric parameter F must be >= 0")}
    if (DELTA_T < 0)       {stop("the numeric parameter DELTA_T must be >= 0")}
    if (MASS < 0)          {stop("the numeric parameter MASS must be >= 0")}
    if (ACC_THRESHOLD < 0) {stop("the numeric parameter ACC_THRESHOLD must be >= 0")}
	invisible(.Call('MDSRefMaps_setCSpace', PACKAGE = 'MDSRefMaps', K, F, MASS, DELTA_T, ACC_THRESHOLD))
}


#' @title Computation of a distance matrix using the Euclidean metric
#' 
#' @description 
#' This function computes the distance matrix from a numeric matrix using the Euclidean metric. 
#'
#' Each row of the input matrix must corresponds to an object and each column must corresponds to an attribute. 
#' 
#' @details
#' This function has been implemented in C++ for fast computations and can be executed from R using this wrapper.
#'
#' @param data a numeric matrix. Rows must correspond to the particles and columns must correspond to the attributes
#' 
#' @return a numeric matrix containing the Euclidean distances betweens all particles
#' 
#' @export
distEuclidean = function(data) {
    if (isTRUE(nrow(data) < 1)) {stop("data must be a numeric matrix")}
    if (isTRUE(ncol(data) < 1)) {stop("data must be a numeric matrix")}
    dist = matrix(rep(0,nrow(data)*nrow(data)),ncol=nrow(data))
    invisible(.Call('MDSRefMaps_distC_euclidean', PACKAGE = 'MDSRefMaps', data, dist))
    return(dist)
}


#' @title Computation of a distance matrix using the Manhattan metric
#' 
#' @description
#' This function computes the distance matrix from a numeric matrix using the Manhattan metric. 
#'
#' Each row of the input matrix must corresponds to an object and each column must corresponds to an attribute. 
#' 
#' @details
#' This function has been implemented in C++ for fast computations and can be executed from R using this wrapper.
#'
#' @param data a numeric matrix. Rows must correspond to the objects and columns must correspond to the attributes
#' 
#' @return a numeric matrix containing the Manhattan distances betweens all objects
#' 
#' @export
distManhattan = function(data) {
    if (isTRUE(nrow(data) < 1)) {stop("data must be a numeric matrix")}
    if (isTRUE(ncol(data) < 1)) {stop("data must be a numeric matrix")}
    dist = matrix(rep(0,nrow(data)*nrow(data)),ncol=nrow(data))
    invisible(.Call('MDSRefMaps_distC_manhattan', PACKAGE = 'MDSRefMaps', data, dist))
    return(dist)
}


#' @title Computation of the Entourage Score
#' 
#' @description
#' This function is used to compute the Entourage Score between two distance matrices, which quantifies the local quality of a MDS representation. 
#' 
#' The Entourage Score corresponds to the normalised number of identical nearest neighbours for each object in the two distance matrices.
#' 
#' @details
#' This function has been implemented in C++ for fast computations and can be executed from R using this wrapper.
#'
#' @param dist1 a numeric matrix of the first distance matrix
#' @param dist2 a numeric matrix of the second distance matrix
#' @param k a numeric indicating the number of nearest neighbours to compare
#' 
#' @return a numeric value of the Entourage Score
#' 
#' @export
computeEntourageScore = function(dist1, dist2, k=3) {
    if (isTRUE(nrow(dist1) != ncol(dist1))) {stop("dist1 must be a numeric matrix of distance")}
    if (isTRUE(nrow(dist2) != ncol(dist2))) {stop("dist1 must be a numeric matrix of distance")}
    if (nrow(dist1) != nrow(dist2)) {stop("dist1 must be the same length as dist2")}
    size              = nrow(dist1)
    names_values      = c(1:size)
    rownames(dist1)   = names_values
    colnames(dist1)   = names_values
    rownames(dist2)   = names_values
    colnames(dist2)   = names_values
    nref = apply(dist1, 1, function(x){names(x[order(x)])[1:k]})
    nfin = apply(dist2, 1, function(x){names(x[order(x)])[1:k]})
    g    = unlist(lapply(c(1:size), function(x){length(intersect(nref[x],nfin[x]))}))
    ent  = sum(g)/(size*k)
    return(ent)
}


#' @title Computation of the Kruskal Stress
#' 
#' @description
#' This function is used to compute the Kruskal Stress between two distance matrices, which quantifies the global quality of a MDS representation.
#'  
#' The Kruskal Stress corresponds to the quantify of information lost during the dimentionality reduction process.
#'
#' @details
#' This function has been implemented in C++ for fast computations and can be executed from R using this wrapper.
#'
#' @param dist1 a numeric matrix of the first distance matrix
#' @param dist2 a numeric matrix of the second distance matrix
#' 
#' @return a numeric value the Kruskal Stress
#' 
#' @export
computeKruskalStress = function(dist1, dist2) {
    if (isTRUE(nrow(dist1) != ncol(dist1))) {stop("dist1 must be a numeric matrix of distance")}
    if (isTRUE(nrow(dist2) != ncol(dist2))) {stop("dist1 must be a numeric matrix of distance")}
    if (nrow(dist1) != nrow(dist2)) {stop("dist1 must be the same length as dist2")}
    stress = invisible(.Call('MDSRefMaps_computeC_stress', PACKAGE = 'MDSRefMaps', dist1, dist2))
    return(stress)
}
