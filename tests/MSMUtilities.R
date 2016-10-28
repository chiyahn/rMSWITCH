library(testthat)



test_that("OrderTransitionMatrix", {
  mu1.order <- order(c(1,-1))   # 2, 1
  mu2.order <- order(c(-1,1))   # 1, 2
  mu3.order <- order(c(3,-2,1)) # 2, 3, 1
  mu4.order <- order(c(-1,3,2)) # 1, 3, 2
  mu5.order <- order(c(1,2,3))  # 1, 2, 3
  transition.matrix1 <- matrix(c(0.4,0.7,0.6,0.3), ncol = 2)
  transition.matrix2 <- matrix(c(0.4,0.5,0.2,0.6,0.3,0.3,0,0.2,0.5), ncol = 3)
  expect_equal(OrderTransitionMatrix(transition.matrix1, mu1.order),
               matrix(c(0.3,0.6,0.7,0.4), ncol = 2))
  expect_equal(OrderTransitionMatrix(transition.matrix1, mu2.order),
               transition.matrix1)
  expect_equal(OrderTransitionMatrix(transition.matrix2, mu3.order),
               matrix(c(0.3,0.3,0.6,0.2,0.5,0,0.5,0.2,0.4), ncol = 3))
  expect_equal(OrderTransitionMatrix(transition.matrix2, mu4.order),
               matrix(c(0.4,0.2,0.5,0,0.5,0.2,0.6,0.3,0.3), ncol = 3))
  expect_equal(OrderTransitionMatrix(transition.matrix2, mu5.order),
               transition.matrix2)
})