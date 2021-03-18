test_that("check_binaries works", {
  expect_true(check_binaries(binaries = "zip"))
  expect_false(check_binaries(binaries = "@#O#OASiw_oiwe238"))
})
