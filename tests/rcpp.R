# ==========================================
# Examples to get used to rcpp
# ==========================================

#install.packages("Rcpp")
#install.packages("testthat")
library(Rcpp)
library(testthat)
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

sourceCpp("C:\\Users\\chiyahn\\Dropbox\\Work\\June15\\rmrs\\tests\\rowSumC.cpp")
expect_equal(rowSumC(mat), apply(mat, 1, sum))

# 6. sample sigma
