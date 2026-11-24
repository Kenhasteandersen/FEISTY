#
# Global run on a small fraction of the cells from COBALT (fixtures/global_cells.csv).
#
# Environment variables:
#   FEISTY_TEST_NCELLS : number of cells. Default is all cells in the fixture.
#   FEISTY_TEST_GLOBAL : "true" also runs calcGlobalCarbonSequestration() on a small area. Needs the
#                        repository with "LU factored TM.Rdata"; its path is given in FEISTY_ROOT.
#

cells <- readCells()
nCells <- as.integer(Sys.getenv("FEISTY_TEST_NCELLS", nrow(cells)))
cells <- cells[unique(round(seq(1, nrow(cells), length.out = min(nCells, nrow(cells))))), ]
fished <- c(1, 3, 5)

test_that("global: biomass and yield in the cells agree with the reference", {
  ref <- utils::read.csv(test_path("fixtures", "reference_global.csv"))
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    for (Fmax in c(0, 0.3)) {
      sim <- simulateCell(cell, Fmax, if (Fmax > 0) fished)
      label <- sprintf("cell %s F=%g", cell$id, Fmax)
      expect_true(all(is.finite(sim$u)), label = label)
      expect_true(min(sim$u) > -1e-6, label = label)
      s <- simSummary(sim)
      r <- ref[ref$id == cell$id & ref$Fmax == Fmax, ]
      expect_equal(s$biomass, r$biomass, tolerance = 2e-3, info = label)
      expect_equal(s$yield,   r$yield,   tolerance = 2e-3, info = label)
    }
  }
})

test_that("global: R and Fortran derivatives are identical in all cells", {
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    p <- setupVertical2(szprod = cell$szprod, lzprod = cell$lzprod, dfbot = cell$dfbot, depth = cell$depth,
                        photic = cell$photic, Tp = cell$Tp, Tm = cell$Tm, Tb = cell$Tb)
    expect_same_derivs(setFishing(p, 0.3), paste("cell", cell$id))
  }
})

test_that("global: the nearest cell is used", {
  cell <- cells[2, ]
  p <- getParametersPosition(cell$lat + 0.1, cell$lon - 0.1, glob = cells)
  expect_equal(p$szprod, cell$szprod)
  expect_equal(p$depth, cell$depth)
  # negative longitudes
  w <- cells[which(cells$lon > 180)[1], ]
  expect_equal(getParametersPosition(w$lat, w$lon - 360, glob = cells)$szprod, w$szprod)
})

test_that("global: fishing is corrected for temperature and depth", {
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    sim <- simulateCell(cell, 0.3, fished, tEnd = 1)
    for (g in fished) {
      Temp <- c(cell$Tp, cell$Tm, cell$Tp, cell$Tm, cell$Tb)[g]
      expected <- 0.3*1.8^((Temp - 15)/10)*(if (cell$depth > 1000) 0.1 else 1)
      expect_equal(sim$p$Fmax[g], expected)
    }
  }
})

test_that("global: cells without Tm can have fishing on all groups", {
  noTm <- cells[is.na(cells$Tm), ]
  skip_if(nrow(noTm) == 0)
  sim <- try(suppressWarnings(simulateCell(noTm[1, ], 0.3, 1:5, tEnd = 5)), silent = TRUE)
  known_issue(is.list(sim) && all(is.finite(sim$u)),
              "simulatePosition() gives NaN when groups 2 or 4 are fished in cells where Tm is NA")
})

test_that("global: calcGlobalCarbonSequestration on a small area", {
  skip_if_not(identical(Sys.getenv("FEISTY_TEST_GLOBAL"), "true"), "FEISTY_TEST_GLOBAL is not true")
  root <- Sys.getenv("FEISTY_ROOT")
  skip_if_not(file.exists(file.path(root, "LU factored TM.Rdata")), "no LU factored TM.Rdata in FEISTY_ROOT")
  old <- setwd(root)
  on.exit(setwd(old))
  TM <- loadTransportMatrix()
  run <- function(...) {
    utils::capture.output(res <- calcGlobalCarbonSequestration(TM, lon = c(330, 336), lat = c(40, 46),
                                                               bPrintStatus = FALSE, nCores = 2, tEnd = 50, ...))
    res
  }
  free <- run()
  fish <- run(Fmax = 0.3, ixGroups = fished)
  expect_true(all(is.finite(free$Biomass)))
  expect_gt(sum(free$Biomass), 0)
  expect_gt(as.numeric(free$TotInject), 0)
  expect_gt(as.numeric(free$TotSeq), 0)
  expect_equal(sum(free$Yield), 0)
  expect_gt(sum(fish$Yield), 0)
})
