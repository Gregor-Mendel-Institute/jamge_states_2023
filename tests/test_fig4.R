

test_that("Fig3 Mains are same as used for revision manuscript", {
  expect_doppelganger("Fig4B1", fig4b1)
  expect_doppelganger("Fig4B2", fig4b2)
  expect_doppelganger("Fig4C", fig4c)
})

test_that("Suppl.3 are same as used for revision manuscript", {  
  expect_doppelganger("FigS4C", figS4c)
  expect_doppelganger("FigS4D1", figS4d1)
  expect_doppelganger("FigS4D2", figS4d2)
  expect_doppelganger("FigS4D3", figS4d3)
  expect_doppelganger("FigS4E", figS4e)
  })
