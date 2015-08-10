library(phylodyn)
context("Uncategorized tests")

test_that("unif_traj produces expected output", {
  expect_equal(unif_traj(-3:3, level = 5), rep(5, 7))
})

test_that("coalsim produces larger coalescent times when effective population is lower", {
  gen1 <- coalsim(s_times = 0, n_sampled = 10, traj = unif_traj, upper = 1, level = 1);
  gen2 <- coalsim(s_times = 0, n_sampled = 10, traj = unif_traj, upper = 100, level = 100);
  expect_true(gen1$coal_times[9] > gen2$coal_times[9])
})

test_that("gen_INLA_args functions agree on isochronous data", {
  data <- coalsim(s_times = 0, n_sampled = 10, traj = unif_traj, upper = 1, level = 1)
  args_old <- gen_INLA_args_old(coal_times = data$coal_times, s_times = data$s_times, n_sampled = data$n_sampled)
  args_new <- gen_INLA_args(s_times = data$s_times, n_sampled = data$n_sampled, coal_times = data$coal_times)
  expect_identical(args_old, args_new)
})

test_that("gen_INLA_args functions agree on serial data", {
  data <- coalsim(s_times = 0:9, n_sampled = rep(1, 10), traj = unif_traj, upper = 1, level = 1)
  args_old <- gen_INLA_args_old(coal_times = data$coal_times, s_times = data$s_times, n_sampled = data$n_sampled)
  args_new <- gen_INLA_args(s_times = data$s_times, n_sampled = data$n_sampled, coal_times = data$coal_times)
  expect_identical(args_old, args_new)
})