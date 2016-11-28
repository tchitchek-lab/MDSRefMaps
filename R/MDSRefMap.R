#' @title Construction of a MDS Reference Map
#' 
#' @description
#' This function computes a MDS Reference Map based on a distance matrix. 
#'
#' A MDS Reference Map corresponds to a regular MDS representation on which additional objects can be projected.
#' MDS Reference Maps can be computed based on the Euclidean or Manhattan metrics using the 'metric' parameter.
#' The initialization space of object positions can be specified using the 'init' parameter.
#'
#' @details
#' The RefMaps algorithm implements SVD-MDS algorithm which is based on a molecular dynamic approach (Becavin et al.).
#' This metric performs a dimensionality reduction of the original space by modeling objects by particles and pairwise distances between them by repulsion and attraction forces.
#' SVD-MDS metric use Verlet algorithm (Loup Verlet in 1967) to compute the MDS representation.
#' Algorithm constants can be specified via the 'setConts()' function.
#'
#' This implementation allows to used incomplete distance matrices (distance matrices with missing values modeled by NA). Furthermore, distance matrices can be computed based on the Euclidean or Manhattan metrics. 
#'
#' This implementation has been implemented in C++ to handle large sets of high-dimensional objects. init. KS+entourage.
#'
#' @param dist a numeric matrix with all pairwise distances between objects of the representation
#' @param k a numeric value specifying the desired number of dimensions in the resulting Reference Map representation
#' @param init a character value or a numeric matrix specifying how the objects are positioned in the initial configuration. Possible character values are 'rand', 'center', 'svd' (please refer to the details section for more details). Object positions in the initial configuration can be explicitly specified using a numeric matrix where the rows correspond to the objects and where the columns correspond to the MDS dimensions (in k  dimensions).
#' @param metric a character indicating the distance metric to use ("euclidean" or "manhattan")
#' @param max_it a numeric defining the maximal number of steps the algorithm can perform
#' @param stress_sd_th a numeric defining the threshold for the standard deviation of Kruskal Stress
#' @param stack_length a numeric defining the length of the Kruskal Stress stack (used to compute the standard deviation of the Kruskal Stress)
#' @param verbose a boolean enabling the display of debug information at each step of the algorithm
#' 
#' @return a list of 3 elements containing the position of the objects ('points' element), the Kruskal Stress ('stress' element), and the Entourage Score ('entourage' element)
#' 
#' @export
MDSReferenceMap = function(dist, k = 2, init = "svd", metric = "euclidean", max_it = 6*10^6, stress_sd_th = 10^-4, stack_length = 500, verbose = TRUE) {
	
	if (class(dist) != "dist" && is.null(nrow(dist))) {
        stop("dist must be a dist object or a distance matrix")
		}
    if (!is.matrix(init) && !(init %in% c("rand","svd"))) {
        stop("init must be a matrix, a character rand or a character svd")
		}
    if (!(metric %in% c("euclidean","manhattan"))) {
        stop("metric must be a character euclidean or a character manhattan")
		}
    if (k < 1) {
		stop("k must be >= 1")
		}
    if (max_it < 1) {
		stop("max_it must be >= 1")
		}
    if (stress_sd_th < 0) {
		stop("stress_sd_th must be >= 0")
		}
    if (stack_length < 2) {
		stop("stack_length must be >= 2")
		}
    if (!(verbose %in% c(TRUE,FALSE))) {
		stop("verbose must be a logical")
		}
    
    dist = as.matrix(dist)
    if (nrow(dist) != ncol(dist)) {
		stop("dist must be a dist object or a distance matrix")
		}
		
    manhattan = if (metric == "euclidean") {0} else {1}
    if (is.matrix(init)) {
        if (verbose) {print("Initialization with the given init matrix")}
        positions = init    
    } else if (init == "rand") {
        if (verbose) {print("Initialization with random positions")}
        np        = dim(dist)[1]
        positions = matrix(stats::runif(np * k, min = -1, max = 1), nrow = np, ncol = k)
    } else if (init == "svd") {
        if (verbose) {print("Initialization with a SVD on the dist matrix, may takes times")}
        svd       = svd(dist)
        positions = t(svd$u)[,c(1:k)]
    }
    
    stress_best    = Inf
    ref_point      = 0
    positions_best = matrix(positions, ncol = ncol(positions))
    if (verbose) {print("Computation of the Reference Map")}
    if (any(is.na(dist))) {
        print("NA value(s) detected in the input distance matrix, computations will be slower")
        invisible(.Call('MDSRefMaps_NA_core', PACKAGE = 'MDSRefMaps', positions_best, stress_best, positions, dist, manhattan, ref_point, max_it, stack_length, stress_sd_th, verbose))
    }else {
        invisible(.Call('MDSRefMaps_core', PACKAGE = 'MDSRefMaps', positions_best, stress_best, positions, dist, manhattan, ref_point, max_it, stack_length, stress_sd_th, verbose))
    }
	
    res           = c()
    res$points    = positions_best
    res$stress    = stress_best
    dist_res      = if (manhattan) {distManhattan(positions_best)} else {distEuclidean(positions_best)}
    res$entourage = computeEntourageScore(dist, dist_res)
    
    return(res)
}


#' @title Construction of a MDS Projection
#' 
#' @description
#' This function computes a MDS Projection based on MDS Reference Map and a distance matrix.
#'
#' A MDS Projection consists on a MDS Reference Map on which additional objects have been overlayed.
#' MDS Projections can be computed based on the Euclidean or Manhattan metrics using the 'metric' parameter.
#' The initialization space of object positions can be specified using the 'init' parameter.
#' 
#' @details
#' Use RefMaps algorithm to compute the projection of additional objects.
#' Efficient to compare multiple projection, for example made on subsets of large biologicals datasets.
#'
#' This implementation allows to used incomplete distance matrices (distance matrices with missing values modeled by NA).
#' Furthermore, distance matrices can be computed based on the Euclidean or Manhattan metrics. 
#'
#' This implementation has been implemented in C++ to handle large sets of high-dimensional objects.
#'
#' @param refmap a result from RefMap function
#' @param dist a numeric matrix with all pairwise distances between objects of the representation
#' @param k a numeric value specifying the desired number of dimensions in the resulting MDS Projection representation
#' @param init a character or a numeric matrix specifying how the objects are positioned in the initial configuration. Possible character values are 'rand', 'center', 'svd' (please refer to the details section for more details). Object positions in the initial configuration can be explicitly specified using a numeric matrix where the rows correspond to the objects and where the columns correspond to the MDS dimensions.
#' @param metric a character indicating the distance metric to use ("euclidean" or "manhattan")
#' @param max_it a numeric defining the maximal number of steps the algorithm can perform
#' @param stress_sd_th a numeric defining the threshold for the standard deviation of Kruskal Stress
#' @param stack_length a numeric defining the length of the Kruskal Stress stack (used to compute the standard deviation of the Kruskal Stress)
#' @param verbose a boolean enabling the display of debug information at each step of the algorithm
#' 
#' @return a list of 3 elements containing the position of the objects ('points' element), the Kruskal Stress ('stress' element), and the Entourage Score ('entourage' element)
#' 
#' @export
MDSProjection = function(refmap, dist, k = 2, init = "svd", metric = "euclidean", max_it = 6*10^6, stress_sd_th = 10^-4, stack_length = 500, verbose = TRUE) {
	
	if (is.null(refmap$points)) {
        stop("refmap$points must be a matrix of a MDS representation")
		}
    if (class(dist) != "dist" && is.null(nrow(dist))) {
        stop("dist must be a dist object or a distance matrix")
		}
    if (!is.matrix(init) && !(init %in% c("rand","svd")))  {
        stop("init must be a matrix, a character rand or a character svd")
		}
    if (!(metric %in% c("euclidean","manhattan"))) {
        stop("metric must be a character euclidean or a character manhattan")}
    if (k < 1) {
		stop("k must be >= 1")
		}
    if (max_it < 1) {
		stop("max_it must be >= 1")
		}
    if (stress_sd_th < 0) {
		stop("stress_sd_th must be >= 0")
		}
    if (stack_length < 2) {
		stop("stack_length must be >= 2")
		}
    if (!(verbose %in% c(TRUE,FALSE))) {
		stop("verbose must be a logical")
		}
	
	dist = as.matrix(dist)
    if (nrow(dist) != ncol(dist)) {
		stop("dist must be a dist object or a distance matrix")
		}
    manhattan = if (metric == "euclidean") {0} else {1}
    
    ref_point = nrow(refmap$points)
    if (is.matrix(init)) {
        positions = init    
    } else if (init == "rand") {
        if (verbose) {print("Initialization with random positions")}
        np        = dim(dist)[1]
        positions = matrix(stats::runif(np * k, min = -5, max = 5), nrow = np, ncol = k)
    } else if (init == "svd") {
        if (verbose) {print("Initialization with a SVD on the dist matrix, may takes times")}
        svd       = svd(dist)
        positions = t(svd$u)[,c(1:k)]
    }
    
    positions[c(1:nrow(refmap$points)),c(1:ncol(refmap$points))] = refmap$points
    positions_best = matrix(positions, ncol = ncol(positions))
    stress_best    = Inf
	if (verbose) {print("Computation of the Projection")}
    if (any(is.na(dist))) {
        print("NA value(s) detected in the input distance matrix, computations will be slower")
        invisible(.Call('MDSRefMaps_NA_core', PACKAGE = 'MDSRefMaps', positions_best, stress_best, positions, dist, manhattan, ref_point, max_it, stack_length, stress_sd_th, verbose))
    }else {
        invisible(.Call('MDSRefMaps_core', PACKAGE = 'MDSRefMaps', positions_best, stress_best, positions, dist, manhattan, ref_point, max_it, stack_length, stress_sd_th, verbose))
    }
    
    res           = c()
    res$points    = positions_best
    res$stress    = stress_best
    dist_res      = if(manhattan) {distManhattan(positions_best)} else {distEuclidean(positions_best)}
    res$entourage = computeEntourageScore(dist, dist_res)
    
    return(res)
    
}


#' @title Construction of a MDS Reference Map with a MDS Projection
#' 
#' @description
#' This function computes a MDS Reference Map with an associated MDS Projection based on a distance matrix.
#'
#' A MDS Projection consists on a MDS Reference Map on which additional objects have been overlayed.
#' MDS Projections can be computed based on the Euclidean or Manhattan metrics using the 'metric' parameter.
#' The initialization space of object positions can be specified using the 'init' parameter.
#'
#' @details
#' A numeric vector is used to specify the objects to use as reference and the objects to add on the Reference Map.
#'
#' @param dist a numeric matrix with all pairwise distances between objects of the representation
#' @param ref_values a numeric vector indicating how the objects must be considered. Objects assigned with 1 are treated as reference objects, objects assigned with 0 are treated as projection objects and objects assigned with NA are not used
#' @param k a numeric value specifying the desired number of dimensions in the resulting MDS Projection representation
#' @param init a character or a numeric matrix specifying how the objects are positioned in the initial configuration. Possible character values are 'rand', 'center', 'svd' (please refer to the details section for more details). Object positions in the initial configuration can be explicitly specified using a numeric matrix where the rows correspond to the objects and where the columns correspond to the MDS dimensions.
#' @param metric a character indicating the distance metric to use ("euclidean" or "manhattan")
#' @param max_it a numeric indicating the maximal number of steps the algorithm can perform
#' @param stress_sd_th a numeric defining the threshold for the standard deviation of Kruskal Stress
#' @param stack_length a numeric defining the length of the Kruskal Stress stack (used to compute the standard deviation of the Kruskal Stress)
#' @param verbose a boolean enabling the display of debug information at each step of the algorithm
#' 
#' @return a list of 3 elements containing the position of the objects ('points' element), the Kruskal Stress ('stress' element), and the Entourage Score ('entourage' element)
#' 
#' @export
MDSReferenceMapwithProjection = function(dist, ref_values, k = 2,  init = "svd", metric = "euclidean", max_it = 6*10^6, stress_sd_th = 10^-4, stack_length = 500, verbose = TRUE) {
	
	if (class(dist) != "dist" && is.null(nrow(dist))) {
        stop("dist must be a dist object or a distance matrix")
		}
    if (is.null(ref_values) || !(ref_values %in% c(0,1,NA))) {
        stop("ref_values must be a numeric vector with 1 for ref points and 0 for projected points, NA for the object to skip")
		}
    if (!is.matrix(init) && !(init %in% c("rand","svd"))) {
        stop("init must be a matrix, a character rand or a character svd")}
    if (!(metric %in% c("euclidean","manhattan"))) {
        stop("metric must be a character euclidean or a character manhattan")}
    if (k < 1) {
		stop("k must be >= 1")
		}
    if (max_it < 1) {
		stop("max_it must be >= 1")
		}
    if (stress_sd_th < 0) {
		stop("stress_sd_th must be >= 0")
		}
    if (stack_length < 2) {
		stop("stack_length must be >= 2")
		}
    if (!(verbose %in% c(TRUE,FALSE))) {
		stop("verbose must be a logical")
		}
	
	dist = as.matrix(dist)
    if (nrow(dist) != ncol(dist)) {
		stop("dist must be a dist object or a distance matrix")
		}
    if (length(ref_values) != nrow(dist)) {
		stop("ref_values must be the same length as the distances objects")
		}
		
    dist       = dist[!is.na(ref_values), !is.na(ref_values)]
    ref_values = stats::na.omit(ref_values)
    dist_ref   = dist[ref_values == 1, ref_values == 1]
    res        = MDSReferenceMap(dist_ref, k, init, metric, max_it, stress_sd_th, stack_length, verbose)
    
    if (0 %in% ref_values) {
        rownames(dist)   = c(1:nrow(dist))
        dist             = dist[order(-ref_values),order(-ref_values)]
        disordered_names = rownames(dist)
        res              = MDSProjection(res, dist, k, init, metric, max_it, stress_sd_th, stack_length, verbose)
        res$points       = res$points[order(as.numeric(disordered_names)),]
    }

    return(res)
}