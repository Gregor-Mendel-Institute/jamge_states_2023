

test_that("Fig3 Mains are same as used for revision manuscript", {
  expect_doppelganger("Fig3A", fig3a)
  expect_doppelganger("Fig3B", fig3b)
  expect_doppelganger("Fig3C", fig3c)
  expect_doppelganger("Fig3D", fig3d)
  expect_doppelganger("Fig3E", fig3e)
  expect_doppelganger("Fig3F", fig3f)
})

test_that("Suppl.3 are same as used for revision manuscript", {  
  expect_doppelganger("FigS3A", figS3a)
  expect_doppelganger("FigS3B", figS3b)
  expect_doppelganger("FigS3C", figS3c)
  expect_doppelganger("FigS3D", figS3d)
  expect_doppelganger("FigS3E1", figS3e1)
  expect_doppelganger("FigS3E2", figS3e2)
  expect_doppelganger("FigS3E3", figS3e3)
  expect_doppelganger("FigS3E4", figS3e4)
  expect_doppelganger("FigS3F", figS3f)
  })
