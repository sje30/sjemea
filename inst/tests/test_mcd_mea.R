## Test the MCD reader.

context("Check MCD txt file without pre and post voltages")
data = mcd.data.to.array(
    system.file("testdata", "mcd_exported_spike.txt", package="sjemea"))

expect_equal(data$spikes[[1]], c(66.108, 537.129))
expect_equal(data$spikes[[2]], c(102.857, 132.813, 218.690))
expect_equal(data$spikes[[3]], c(402.257))
expect_equal(data$channels, c(31, 41, 51))

context("Check MCD txt file wit pre and post voltages")
data = mcd.data.to.array(
  system.file("testdata", "mcd_exported_pre_post_spike.txt", package="sjemea"),
  pre.spike.ms = 1,
  post.spike.ms = 2)

expect_equal(data$spikes[[1]], c(66.108, 537.129))
expect_equal(data$spikes[[2]], c(102.857, 132.813, 218.690))
expect_equal(data$spikes[[3]], c(402.257))
expect_equal(data$channels, c(31, 41, 51))
