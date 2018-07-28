context("Vectorization operators")
library(gmvarkit)

a1 <- 1
a2 <- 1:4
a3 <- 1:7^2
A1 <- matrix(a1, nrow=1)
A2 <- matrix(a2, nrow=2, byrow=FALSE)
A3 <- matrix(a3, nrow=7, byrow=FALSE)


test_that("vec and unvec works correctly", {
  expect_equal(vec(A1), a1)
  expect_equal(vec(A2), a2)
  expect_equal(vec(A3), a3)
  expect_equal(unvec(d=1, a=a1), A1)
  expect_equal(unvec(d=2, a=a2), A2)
  expect_equal(unvec(d=7, a=a3), A3)
})

B1 <- matrix(1)
B2 <- matrix(c(1, 0.5, 0.5, 2), nrow=2, byrow=FALSE)
B3 <- matrix(c(1, 0.3, 0.2, 0.3, 2, 0.4, 0.2, 0.4, 3), nrow=3, byrow=FALSE)
b1 <- 1
b2 <- c(1, 0.5, 2)
b3 <- c(1, 0.3, 0.2, 2, 0.4, 3)

test_that("vech and unvech works correctly", {
  expect_equal(vech(B1), b1)
  expect_equal(vech(B2), b2)
  expect_equal(vech(B3), b3)
  expect_equal(unvech(d=1, a=b1), B1)
  expect_equal(unvech(d=2, a=b2), B2)
  expect_equal(unvech(d=3, a=b3), B3)
})
