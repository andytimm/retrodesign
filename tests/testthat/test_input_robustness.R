library(testthat)
library(retrodesign)

context("Make sure functions respond correctly to positive, negative,
        and impossible inputs")

test_that("all functions reject negative standard error", {
  expect_error(retro_design_closed_form(.5,-1))
  expect_error(retrodesign(.5,-1))
  expect_error(retro_design_closed_form(list(.5,1,2),-1))
  expect_error(retrodesign(list(.5,1,2),-1))
  expect_error(type_s(.5,-1))
  expect_error(type_m(.5,-1))
  expect_error(type_s(list(.5,1,2),-1))
  expect_error(type_m(list(.5,1,2),-1))

})

test_that("type M error is always a positive ratio", {
  # This is assumed to be true by Lu et al. (2018) and Gelman & Carlin (2014),
  # but making it work explicitly in the package required some modifications.
  expect_gte(retrodesign(-.5,1)$type_m,0)
  expect_gte(retro_design_closed_form(-.5,1)$type_m,0)
  expect_gte(unlist(retro_design_closed_form(list(-.5,-1,-2),1)$type_m[2]),0)
  expect_gte(unlist(retrodesign(list(-.5,-1,-2),1)$type_m[2]),0)
  expect_gte(unlist(type_m(-.5,1)),0)
  expect_gte(unlist(type_m(list(-.5,-1,-2),1)$type_m[2]),0)

})

test_that("type M error is always a positive ratio", {
  # This is assumed to be true by Lu et al. (2018) and Gelman & Carlin (2014),
  # but making it work explicitly in the package required some modifications.
  expect_gte(retrodesign(-.5,1)$type_m,0)
  expect_gte(retro_design_closed_form(-.5,1)$type_m,0)
  expect_gte(unlist(retro_design_closed_form(list(-.5,-1,-2),1)$type_m[2]),0)
  expect_gte(unlist(retrodesign(list(-.5,-1,-2),1)$type_m[2]),0)
  expect_gte(unlist(type_m(-.5,1)),0)
  expect_gte(unlist(type_m(list(-.5,-1,-2),1)$type_m[2]),0)

})

test_that("type S error works for negative numbers", {
  # Type S should be P(negative estimate|positive true effect), and vice
  # versa for negative true effects.
  # we can easily test this by symmetry
  expect_equal(retrodesign(-.5,1)$typeS,
                   retrodesign(.5,1)$typeS)
  expect_equal(retro_design_closed_form(-.5,1)$typeS,
                retro_design_closed_form(.5,1)$typeS)
  expect_equal(unlist(retro_design_closed_form(list(-.5,-1,-2),1)$type_s[2]),
                   unlist(retro_design_closed_form(list(.5,1,2),1)$type_s[2]))
  expect_equal(unlist(retrodesign(list(-.5,-1,-2),1)$type_s[2]),
                   unlist(retrodesign(list(.5,1,2),1)$type_s[2]))
  expect_equal(unlist(type_s(-.5,1)),
                   unlist(type_s(.5,1)))
  expect_equal(unlist(type_s(list(-.5,-1,-2),1)$type_s[2]),
                   unlist(type_s(list(.5,1,2),1)$type_s[2]))

})

test_that("All functions are robust to vector input", {
  skip_on_cran()
  # This comes by design from the update the new non-central t-dist
  # code, but was not true in the original paper code. That was a known
  # limitation of the paper code and it was never used incorrectly in the paper
  # analysis, but a package should anticipate common mistakes.
  # Note that to make test runtime more reasonable, we're using a fairly lax tol,
  # and skipping on CRAN to be respectful of their compute.
  expect_equal(retrodesign::retrodesign(c(10,rep(0.1,100)), 1, alpha = 0.05,
                                        df = Inf, n.sims = 10^5)$type_m[1],

               retrodesign::retrodesign(c(10), 1, alpha = 0.05,
                                        df = Inf, n.sims = 10^5)$type_m,

               tolerance = .001
  )

  expect_equal(retrodesign::retro_design_closed_form(c(.5,rep(0.1,100)), 1,
                                                     alpha = 0.05)$type_m[1],

               retrodesign::retro_design_closed_form(c(.5), 1,
                                                     alpha = 0.05)$type_m,

               tolerance = .001
  )

  # This is slightly lower tolerance just t
  expect_equal(retrodesign::type_m(c(10,rep(0.1,100)), 1,
                                                     alpha = 0.05)$type_m[1],

               retrodesign::type_m(c(10), 1,
                                                     alpha = 0.05)$type_m,

               tolerance = .001
  )


})
