#include <Rcpp.h>
#include <queue>
using namespace Rcpp;

unsigned int ref_point = 0;
double DIST_2 = 0;

double K = 1;
double F = 0.1;
double MASS = 10;
double DELTA_T = 0.01;
double TWO_TIMES_DELTA_T = 2 * DELTA_T;
double DELTA_T_2 = pow(DELTA_T, 2.0);
double DELTA_T_2_DIV_MASS = DELTA_T_2 / MASS;
double ACC_THRESHOLD = 0;

NumericMatrix positions;
NumericMatrix positions_prev;
NumericMatrix accelerations;
NumericMatrix velocities;
NumericMatrix dist_initial;
NumericMatrix dist_current;

std::queue<double> stress_stack_queue;
std::queue<double> stress_trace_queue;

void setCSpace(double K_value, double F_value, double MASS_value, double DELTA_T_value, double ACC_THRESHOLD_value);
RcppExport SEXP MDSRefMaps_setCSpace(SEXP K_valueSEXP, SEXP F_valueSEXP, SEXP MASS_valueSEXP, SEXP DELTA_T_valueSEXP, SEXP ACC_THRESHOLD_valueSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type K_value(K_valueSEXP);
    Rcpp::traits::input_parameter< double >::type F_value(F_valueSEXP);
    Rcpp::traits::input_parameter< double >::type MASS_value(MASS_valueSEXP);
    Rcpp::traits::input_parameter< double >::type DELTA_T_value(DELTA_T_valueSEXP);
    Rcpp::traits::input_parameter< double >::type ACC_THRESHOLD_value(ACC_THRESHOLD_valueSEXP);
    setCSpace(K_value, F_value, MASS_value, DELTA_T_value, ACC_THRESHOLD_value);
    return R_NilValue;
END_RCPP
}

void distC_euclidean(NumericMatrix x, NumericMatrix dist);
RcppExport SEXP MDSRefMaps_distC_euclidean(SEXP xSEXP, SEXP distSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dist(distSEXP);
    distC_euclidean(x, dist);
    return R_NilValue;
END_RCPP
}

void distC_manhattan(NumericMatrix x, NumericMatrix dist);
RcppExport SEXP MDSRefMaps_distC_manhattan(SEXP xSEXP, SEXP distSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dist(distSEXP);
    distC_manhattan(x, dist);
    return R_NilValue;
END_RCPP
}

double computeC_stress(NumericMatrix dist1, NumericMatrix dist2);
RcppExport SEXP MDSRefMaps_computeC_stress(SEXP dist1SEXP, SEXP dist2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type dist1(dist1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dist2(dist2SEXP);
    __result = Rcpp::wrap(computeC_stress(dist1, dist2));
    return __result;
END_RCPP
}

void core(NumericMatrix positions_best, NumericVector stress_best, NumericMatrix positionsV, NumericMatrix distanceV, bool manhattan, unsigned int ref_point_value, unsigned int max_iterations, unsigned int stress_stack_length, double stress_stack_sd_th, bool verbose);
RcppExport SEXP MDSRefMaps_core(SEXP positions_bestSEXP, SEXP stress_bestSEXP, SEXP positionsVSEXP, SEXP distanceVSEXP, SEXP manhattanSEXP, SEXP ref_point_valueSEXP, SEXP max_iterationsSEXP, SEXP stress_stack_lengthSEXP, SEXP stress_stack_sd_thSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type positions_best(positions_bestSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stress_best(stress_bestSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type positionsV(positionsVSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type distanceV(distanceVSEXP);
    Rcpp::traits::input_parameter< bool >::type manhattan(manhattanSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ref_point_value(ref_point_valueSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type stress_stack_length(stress_stack_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type stress_stack_sd_th(stress_stack_sd_thSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    core(positions_best, stress_best, positionsV, distanceV, manhattan, ref_point_value, max_iterations, stress_stack_length, stress_stack_sd_th, verbose);
    return R_NilValue;
END_RCPP
}

void NA_core(NumericMatrix positions_best, NumericVector stress_best, NumericMatrix positionsV, NumericMatrix distanceV, bool manhattan, unsigned int ref_point_value, unsigned int max_iterations, unsigned int stress_stack_length, double stress_stack_sd_th, bool verbose);
RcppExport SEXP MDSRefMaps_NA_core(SEXP positions_bestSEXP, SEXP stress_bestSEXP, SEXP positionsVSEXP, SEXP distanceVSEXP, SEXP manhattanSEXP, SEXP ref_point_valueSEXP, SEXP max_iterationsSEXP, SEXP stress_stack_lengthSEXP, SEXP stress_stack_sd_thSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type positions_best(positions_bestSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stress_best(stress_bestSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type positionsV(positionsVSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type distanceV(distanceVSEXP);
    Rcpp::traits::input_parameter< bool >::type manhattan(manhattanSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ref_point_value(ref_point_valueSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type stress_stack_length(stress_stack_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type stress_stack_sd_th(stress_stack_sd_thSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    NA_core(positions_best, stress_best, positionsV, distanceV, manhattan, ref_point_value, max_iterations, stress_stack_length, stress_stack_sd_th, verbose);
    return R_NilValue;
END_RCPP
}

void setCSpace(double K_value, double F_value, double MASS_value, double DELTA_T_value, double ACC_THRESHOLD_value){
    K = K_value;
    F = F_value;
    MASS = MASS_value;
    DELTA_T = DELTA_T_value;
    TWO_TIMES_DELTA_T = 2 * DELTA_T;
    DELTA_T_2 = pow(DELTA_T, 2.0);
    DELTA_T_2_DIV_MASS = DELTA_T_2 / MASS;
    ACC_THRESHOLD = ACC_THRESHOLD_value;
    return;
}

// Compute the euclidean distance of a matrix
// EXPORT
void distC_euclidean(NumericMatrix x, NumericMatrix dist) {
    
    unsigned int outrows = x.nrow();
    for (unsigned int i = 0; i < outrows - 1; i++) {
        NumericVector row_i = x.row(i);
        for (unsigned int j = i + 1; j < outrows; j ++) {
            double d = sqrt(sum(pow(row_i-x.row(j), 2.0)));
            dist(j,i) = d;
            dist(i,j) = d;
        }
    }
    return;
}

// Compute the manhattan distance of a matrix
// EXPORT
void distC_manhattan(NumericMatrix x, NumericMatrix dist) {
    
    unsigned int outrows = x.nrow();
    for (unsigned int i = 0; i < outrows - 1; i++) {
        NumericVector row_i = x.row(i);
        for (unsigned int j = i + 1; j < outrows; j ++) {
            double d = sum(sqrt(pow(row_i-x.row(j), 2.0)));
            dist(j,i) = d;
            dist(i,j) = d;
        }
    }

    return;
}

// Compute the Kurskal stress based on the difference between the distance. Skip NA Values.
// EXPORT
double computeC_stress(NumericMatrix dist1, NumericMatrix dist2) {

    double sum_diff_sq = 0.0;
	double sum_i_sq    = 0.0;

    unsigned int nbrow = dist1.nrow();

    for (unsigned int i = 0; i < nbrow - 1; i++) {
        for (unsigned int j = i + 1; j < nbrow; j++) {
            if (dist1(i,j) != dist1(i,j)){continue;}
            if (dist2(i,j) != dist2(i,j)){continue;}
            sum_diff_sq = sum_diff_sq + pow((dist2(i,j)-dist1(i,j)),2);
			sum_i_sq    = sum_i_sq + pow(dist1(i,j),2);
        }
    }

    return (sqrt(sum_diff_sq/sum_i_sq)*100);
}

// Compute the standard_deviation of a queue
double standard_deviation(){
    double sum = 0.0;
    std::queue<double> tmp_queue;
    tmp_queue = std::queue<double>(stress_stack_queue);
    
    while (tmp_queue.size() > 0) {
        sum += tmp_queue.front();
        tmp_queue.pop();
    }
    
    double mean = sum/stress_stack_queue.size();
    double temp = 0.0;
    
    tmp_queue = std::queue<double>(stress_stack_queue);
    while (tmp_queue.size() > 0) {
        temp += (mean-tmp_queue.front())*(mean-tmp_queue.front());
        tmp_queue.pop();
    }
    
    return sqrt(temp/stress_stack_queue.size());
}

// Set the values of a NumericVector from a queue 
void nvFromQueue(NumericVector stress_stack, std::queue<double> data){
    unsigned int i = 0;
    while (data.size() > 0) {
        stress_stack(i) = data.front();
        data.pop();
        i=i+1;
    }
}

// Set the value of a NumericMatrix from another NumericMatrix
void setNMValues(NumericMatrix x, const NumericMatrix y){
    
    unsigned int nbrow = x.nrow();
    unsigned int nbcol = x.ncol();
    for (unsigned int i = 0; i < nbrow; i++) {
        for (unsigned int j = 0; j < nbcol; j++) {
            x(i,j) = y(i,j);
        }
    }
    
    return;
}

// Set global the sum of the distance to the square
void setDist2() {

    DIST_2 = 0;

    unsigned int nbrow = dist_initial.nrow();

    for (unsigned int i = 0; i < nbrow - 1; i++) {
        for (unsigned int j = i + 1; j < nbrow; j++) {
            DIST_2 = DIST_2 + pow(dist_initial(i,j),2);
        }
    }
    
    return;
}

// Set global the sum of the distance to the square
void NA_setDist2() {

    DIST_2 = 0;

    unsigned int nbrow = dist_initial.nrow();

    for (unsigned int i = 0; i < nbrow - 1; i++) {
        for (unsigned int j = i + 1; j < nbrow; j++) {
            if (dist_initial(i,j) != dist_initial(i,j)){continue;}
            DIST_2 = DIST_2 + pow(dist_initial(i,j),2);
        }
    }
    
    return;
}

// Compute the positions and velocities matrix.
void compute_p_v_ACC_THRESHOLD() {
    
    double positions_tmp_value = 0;
    double min_acc = -ACC_THRESHOLD; // acceleration negative
    double max_acc =  ACC_THRESHOLD;
    
    unsigned int nbrow = positions.nrow();
    unsigned int nbcol = positions.ncol();
    for (unsigned int i = ref_point; i < nbrow; i++) {
        for (unsigned int j = 0; j < nbcol; j++) {
            if(accelerations(i,j) > max_acc) {
                accelerations(i,j) = max_acc;
            }else if (accelerations(i,j) < min_acc) {
                accelerations(i,j) = min_acc;
            }
            positions_tmp_value = positions(i,j);
            positions(i,j) = 2*positions(i,j)-positions_prev(i,j)+accelerations(i,j)*DELTA_T_2_DIV_MASS;
            velocities(i,j) = (positions(i,j)-positions_prev(i,j))/TWO_TIMES_DELTA_T;
            positions_prev(i,j) = positions_tmp_value;
        }
    }
    
    return;
}

// Compute the positions and velocities matrix.
void compute_p_v() {
    
    double positions_tmp_value = 0;
    
    unsigned int nbrow = positions.nrow();
    unsigned int nbcol = positions.ncol();
    for (unsigned int i = ref_point; i < nbrow; i++) {
        for (unsigned int j = 0; j < nbcol; j++) {
            positions_tmp_value = positions(i,j);
            positions(i,j) = 2*positions(i,j)-positions_prev(i,j)+accelerations(i,j)*DELTA_T_2_DIV_MASS;
            velocities(i,j) = (positions(i,j)-positions_prev(i,j))/TWO_TIMES_DELTA_T;
            positions_prev(i,j) = positions_tmp_value;
        }
    }
    
    return;
}

// Update the accelerations matrix based on the difference between the distance.
double compute_stress_accel() {

    double sum_diff_sq = 0;
    accelerations = NumericMatrix(positions.nrow(),positions.ncol());

    unsigned int nbrow = positions.nrow();
    unsigned int nbcol = positions.ncol();    

    for (unsigned int i = 0; i < nbrow - 1; i++) {
        unsigned int j = i;
        while (j < ref_point) {
            double diff = dist_current(i,j)-dist_initial(i,j);
            sum_diff_sq = sum_diff_sq + pow(diff,2);
            j++;
        }
        while (j < nbrow) {
            double diff = dist_current(i,j)-dist_initial(i,j);
            sum_diff_sq = sum_diff_sq + pow(diff,2);
            for (unsigned int k = 0; k < nbcol; k++) {
                double force       = K*diff;
                accelerations(i,k) = accelerations(i,k) + -force*(positions(i,k)-positions(j,k)) - F*velocities(i,k);
                accelerations(j,k) = accelerations(j,k) + force*(positions(i,k)-positions(j,k)) - F*velocities(j,k);
            }
            j++;
        }
    }
    
    return (sqrt(sum_diff_sq/DIST_2)*100);
}

// Skip NA Values. Update the accelerations matrix based on the difference between the distance.
double NA_compute_stress_accel() {

    double sum_diff_sq = 0;
    accelerations = NumericMatrix(positions.nrow(),positions.ncol());

    unsigned int nbrow = positions.nrow();
    unsigned int nbcol = positions.ncol();    

    for (unsigned int i = 0; i < nbrow - 1; i++) {
        unsigned int j = i;
        while (j < ref_point) {
            if (dist_initial(i,j) != dist_initial(i,j)){
                j++;
                continue;
            }
            double diff = dist_current(i,j)-dist_initial(i,j);
            sum_diff_sq = sum_diff_sq + pow(diff,2);
            j++;
        }
        while (j < nbrow) {
            if (dist_initial(i,j) != dist_initial(i,j)){
                j++;
                continue;
            }
            double diff = dist_current(i,j)-dist_initial(i,j);
            sum_diff_sq = sum_diff_sq + pow(diff,2);
            for (unsigned int k = 0; k < nbcol; k++) {
                double force       = K*diff;
                accelerations(i,k) = accelerations(i,k) + -force*(positions(i,k)-positions(j,k)) - F*velocities(i,k);
                accelerations(j,k) = accelerations(j,k) + force*(positions(i,k)-positions(j,k)) - F*velocities(j,k);
            }
            j++;
        }
    }
    
    return (sqrt(sum_diff_sq/DIST_2)*100);
}

// Compute a MDS using the verlet algorithm. Can ignore the reference points (ref_point).
// EXPORT
void core(NumericMatrix positions_best, NumericVector stress_best, NumericMatrix positionsV, NumericMatrix distanceV, bool manhattan, unsigned int ref_point_value, unsigned int max_iterations, unsigned int stress_stack_length, double stress_stack_sd_th, bool verbose) {
    
    // set local positions, dist_initial, dist_current
    ref_point = ref_point_value;
    positions = positionsV;
    dist_initial = distanceV;
    setDist2();
    dist_current = NumericMatrix(positions.nrow(),positions.nrow());
    if (manhattan) {distC_manhattan(positions, dist_current);}else {distC_euclidean(positions, dist_current);}

    // initialize positions_prev
    positions_prev = NumericMatrix(positions.nrow(),positions.ncol());
    setNMValues(positions_prev, positions);
    
    // initialize accelerations and velocities
    accelerations = NumericMatrix(positions.nrow(),positions.ncol());
    velocities = NumericMatrix(positions.nrow(),positions.ncol());
    
    // initialize stress
    stress_stack_queue = std::queue<double>(); 
    double stress = compute_stress_accel();
    stress_best(0) = stress;
    stress_stack_queue.push(stress);
    double sd_stress_stack = 1;
    if(verbose){
        Rcout << "Starting to compute - stress : " << stress << std::endl;
    }
    unsigned int t = 1;
    while ((t < max_iterations)  &&  (sd_stress_stack > stress_stack_sd_th) ) {
    
        if(verbose){
            Rcout << "Step : " << t << " - stress : " << stress << " - sd : " << sd_stress_stack << std::endl;
        }

        if (ACC_THRESHOLD == 0) {compute_p_v();} else {compute_p_v_ACC_THRESHOLD();}
        
        if (manhattan) {distC_manhattan(positions, dist_current);}else {distC_euclidean(positions, dist_current);}
        
        stress = compute_stress_accel();
        stress_stack_queue.push(stress);
        
        if (t >= stress_stack_length) {
            stress_stack_queue.pop();
            sd_stress_stack = standard_deviation();
        }
        
        if (stress < stress_best(0)) {
            stress_best(0) = stress;
            setNMValues(positions_best, positions);
        }
        
        t = t + 1;
    }
    if(verbose){
        Rcout << "Last step : " << t << " - stress : " << stress << " - sd : " << sd_stress_stack << std::endl;
        Rcout << "Finished to compute - best stress : " << stress_best(0) << std::endl;
        Rcout << "For the parameters : K=" << K << " F=" << F << " MASS=" << MASS << " DELTA_T=" << DELTA_T << " DELTA_T_2=" << DELTA_T_2 << " ACC_THRESHOLD=" << ACC_THRESHOLD << " iMax=" << max_iterations << " sdMin=" << stress_stack_sd_th << std::endl;
    }

    return;
}

// Skip NA Values. Compute a MDS using the verlet algorithm. Can ignore the reference points (ref_point).
// EXPORT
void NA_core(NumericMatrix positions_best, NumericVector stress_best, NumericMatrix positionsV, NumericMatrix distanceV, bool manhattan, unsigned int ref_point_value, unsigned int max_iterations, unsigned int stress_stack_length, double stress_stack_sd_th, bool verbose) {
    
    //Set local positions, dist_initial, dist_current
    ref_point = ref_point_value;
    positions = positionsV;
    dist_initial = distanceV;
    NA_setDist2();
    dist_current = NumericMatrix(positions.nrow(),positions.nrow());
    if (manhattan) {distC_manhattan(positions, dist_current);}else {distC_euclidean(positions, dist_current);}

    // initialize positions_prev
    positions_prev = NumericMatrix(positions.nrow(),positions.ncol());
    setNMValues(positions_prev, positions);
    
    // initialize accelerations and velocities
    accelerations = NumericMatrix(positions.nrow(),positions.ncol());
    velocities = NumericMatrix(positions.nrow(),positions.ncol());
    
    // initialize stress
    stress_stack_queue = std::queue<double>(); 
    double stress = NA_compute_stress_accel();
    stress_best(0) = stress;
    stress_stack_queue.push(stress);
    double sd_stress_stack = 1;
    if(verbose){
        Rcout << "Starting to compute - stress : " << stress << std::endl;
    }
    unsigned int t = 1;
    while ((t < max_iterations)  &&  (sd_stress_stack > stress_stack_sd_th) ) {
    
        if(verbose){
            Rcout << "Step : " << t << " - stress : " << stress << " - sd : " << sd_stress_stack << std::endl;
        }

        if (ACC_THRESHOLD == 0) {compute_p_v();} else {compute_p_v_ACC_THRESHOLD();}
        
        if (manhattan) {distC_manhattan(positions, dist_current);}else {distC_euclidean(positions, dist_current);}
        
        stress = NA_compute_stress_accel();
        stress_stack_queue.push(stress);
        
        if (t >= stress_stack_length) {
            stress_stack_queue.pop();
            sd_stress_stack = standard_deviation();
        }
        
        if (stress < stress_best(0)) {
            stress_best(0) = stress;
            setNMValues(positions_best, positions);
        }
        
        t = t + 1;
    }
    if(verbose){
        Rcout << "Last step : " << t << " - stress : " << stress << " - sd : " << sd_stress_stack << std::endl;
        Rcout << "Finished to compute - best stress : " << stress_best(0) << std::endl;
        Rcout << "For the parameters : K=" << K << " F=" << F << " MASS=" << MASS << " DELTA_T=" << DELTA_T << " DELTA_T_2=" << DELTA_T_2 << " ACC_THRESHOLD=" << ACC_THRESHOLD << " iMax=" << max_iterations << " sdMin=" << stress_stack_sd_th << std::endl;
    }

    return;
}