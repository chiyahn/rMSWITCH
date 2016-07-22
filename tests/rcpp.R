# ==========================================
# Examples to get used to rcpp
# ==========================================
#install.packages("Rcpp")
#install.packages("testthat")
#install.packages("inline")
library(Rcpp)
library(testthat)
library(inline)
# 1. scalar input & scalar output
cppFunction("int add(int x, int y, int z) {  return (x + y + z); }")
expect_equal(add(1,2,3), 1+2+3)

# 2. no input & scalar output
cppFunction("int gimmeOne() { return 1; }")
gimmeOne()

cppFunction("int gimmeOneIfGreaterThanOne(int x) { 
            if (x > 1)
              return 1; 
            return 0;
            }")
expect_equal(gimmeOneIfGreaterThanOne(2), 1)
expect_equal(gimmeOneIfGreaterThanOne(0), 0)

# 3. vector input & scalar output
cppFunction("double sumOfSequences(NumericVector x) {
            int n = x.size();
            double sum = 0;
            for (int i = 0; i < n; i++)
              sum += x[i];
            return sum;
            }")
expect_equal(sumOfSequences(seq(1:5)), sum(seq(1:5)))

# 4. vector input & vector output
cppFunction("NumericVector elementwiseSubtract(double x, NumericVector ys) {
              int n = ys.size();
              NumericVector value = NumericVector(n);
              for (int i = 0; i < n; i++)
                value[i] = ys[i] - x;
              return value;
            }")
cppFunction("NumericVector elementwiseSubtractRLike(double x, NumericVector ys) {
              return (ys-x);
            }")
expect_equal(elementwiseSubtract(1, seq(1:5)), (seq(1:5)) - 1)
expect_equal(elementwiseSubtractRLike(1, seq(1:5)), (seq(1:5)) - 1)


# 5. matrix input & vector output
cppFunction("NumericVector rowSum(NumericMatrix mat) {
            int nrow = mat.nrow();
            int ncol = mat.ncol();
            NumericVector vec = NumericVector(nrow);
            for (int i = 0; i < nrow; i++)
              for (int j = 0; j < ncol; j++)
                vec[i] += mat(i,j);
            return vec;
            }")
mat <- matrix(seq(1:12), nrow = 4)
expect_equal(rowSum(mat), apply(mat, 1, sum))

# 6. RcppArmadillo: eigenvector
GetStationaryDistribution <- cxxfunction(signature(matrix="matrix"),
                                         plugin="RcppArmadillo", body='
arma::mat transition_probs = Rcpp::as<arma::mat>(matrix);
arma::cx_vec eigval;
arma::cx_mat eigvec;
arma::eig_gen(eigval, eigvec, transition_probs.t()); // left eigenv, so take t.
int M = transition_probs.n_cols;
arma::vec stationary_dist(M);

// Need to find which eigenvector has a eigenval. of one and extract real parts
// find index
int stationary_dist_index = 0;
for (int i = 0; i < M; i++) 
  if (eigval(i).real() == 1)
    stationary_dist_index = i;
// extract real
for (int i = 0; i < M; i++) 
  stationary_dist(i) = eigvec(i,stationary_dist_index).real(); 

stationary_dist = abs(stationary_dist) / sum(abs(stationary_dist));
return wrap(stationary_dist);')
P <- matrix(c(1/4,1/5,3/4,4/5), ncol = 2)
mu <- GetStationaryDistribution(P)
expect_equal(t(mu) %*% P, t(mu)) # check if it is indeed stationary
expect_equal(sum(mu), 1) # must sum up to 1
P <- matrix(c(1/4,1/5,1/6,2/4,3/5,2/6,1/4,1/5,3/6), ncol = 3)
mu <- GetStationaryDistribution(P)
expect_equal(t(mu) %*% P, t(mu)) # check if it is indeed stationary
expect_equal(sum(mu), 1) # must sum up to 1

