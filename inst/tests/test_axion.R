context("Guessing numbers of wells")
expect_equal(axion.guess.well.number("D3_33"), 48)
expect_equal(axion.guess.well.number("A3_37"), 12)
##expect_equal(axion.guess.well.number("D3_39"), 12) #does this fail?
expect_that(axion.guess.well.number("A3_39"), throws_error())
