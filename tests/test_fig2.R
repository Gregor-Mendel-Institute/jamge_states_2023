

test_that("Fig2 Mains are same as used for revision manuscript", {
  expect_doppelganger("Fig2A", fig2a)
  expect_doppelganger("Fig2B", fig2b)
  expect_doppelganger("Fig2C", fig2c)
  expect_doppelganger("Fig2D", fig2d)
  expect_doppelganger("Fig2E", fig2e)
  expect_doppelganger("Fig2F", fig2f)
})

test_that("Suppl. 2 are same as used for revision manuscript", {  
  # B missing
  expect_doppelganger("FigS2C", figS2c)
  expect_doppelganger("FigS2D", figS2d)
  expect_doppelganger("FigS2E", figS2e)
 # expect_doppelganger("FigS2F", figS2f) # pie
  expect_doppelganger("FigS2G", figS2g) 
  expect_doppelganger("FigS2H", figS2h) 
  })
