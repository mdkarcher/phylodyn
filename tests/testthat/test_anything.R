library(phylodyn)
context("Uncategorized tests")

test_that("unif_traj produces expected output", {
  expect_equal(unif_traj(-3:3, level = 5), rep(5, 7))
})

test_that("coalsim produces larger coalescent times when effective population is higher", {
  gen1 <- coalsim(samp_times = 0, n_sampled = 10, traj = unif_traj, level = 1);
  gen2 <- coalsim(samp_times = 0, n_sampled = 10, traj = unif_traj, level = 100);
  expect_true(gen1$coal_times[9] < gen2$coal_times[9])
})

test_that("BNPR produces correct INLA arguments", {
  samp_times = 0
  n_sampled  = 4
  coal_times = c(0.25, 0.75, 1.0)
  grid = c(0.0, 0.5, 1.0)
  answer = data.frame(rbind(c(0.25, 1, 2.25, log(2.25)), c(0.75, 2, 1, log(1))))
  names(answer) = c("time", "event", "E", "E_log")
  
  data1 = coal_stats(grid = grid, samp_times = samp_times, coal_times = coal_times, n_sampled = n_sampled)
  data1 = with(data1, condense_stats(time = time, event = event, E = E))
  expect_equal(data1, answer)
  
  if (requireNamespace("INLA", quietly = TRUE)) {
    data2 = BNPR(data = list(samp_times = samp_times, n_sampled = n_sampled, coal_times= coal_times), lengthout = 2)$data
    expect_equal(data2$y, answer$event)
    expect_equal(data2$time, answer$time)
    expect_equal(data2$E_log, answer$E_log)
  }
})