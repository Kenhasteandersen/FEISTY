#
# Carbon fluxes, injection, projection onto the grid of the transport matrix, and sequestration
#

rho <- 9   # gWW per gC
pathways <- c("Fecal", "Carcass", "Repro", "Respiration")
lastTimes <- function(sim) sim$t >= 0.6*max(sim$t)
byGroup <- function(x, p) as.vector(rowsum(x, rep(seq_along(p$ix), lengths(p$ix))))

cells <- readCells()
nCells <- as.integer(Sys.getenv("FEISTY_TEST_NCELLS", nrow(cells)))
cells <- cells[unique(round(seq(1, nrow(cells), length.out = min(nCells, nrow(cells))))), ]

# Simulations used by several tests: the cells without and with fishing.
sims <- list()
for (i in seq_len(nrow(cells))) for (Fmax in c(0, 0.3))
  sims[[sprintf("%s F=%g", cells$id[i], Fmax)]] <-
    calcCarbonFluxes(simulateCell(cells[i, ], Fmax, if (Fmax > 0) c(1, 3, 5), tEnd = 300))
injections <- lapply(sims, calcCarbonInjection)

TM <- local({
  e <- new.env()
  utils::data("CTL", package = "FEISTY", envir = e)
  e$TM
})

# ---- Fluxes ------------------------------------------------------------------------------------

test_that("carbon: the fluxes equal the losses in the model", {
  for (n in names(sims)) {
    sim <- sims[[n]]
    flux <- byGroup(sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration, sim$p)
    expect_close(flux, colMeans(sim$totLoss[lastTimes(sim), ]), 1e-6, paste(n, "losses"))
    carcass <- colMeans(sim$B[lastTimes(sim), ]) * sim$p$mort0[sim$p$ixFish]
    expect_close(sim$fluxCarcass, carcass, 1e-12, paste(n, "carcasses"))
    for (x in paste0("flux", pathways))
      expect_true(all(sim[[x]] >= -1e-8), label = paste(n, x))
  }
})

test_that("carbon: the components of respiration add up", {
  for (n in names(sims)) {
    sim <- sims[[n]]
    expect_close(sim$fluxRespirationBasal + sim$fluxRespirationSDA + sim$fluxRespirationRepro,
                 sim$fluxRespiration, 1e-12, n)
    # Half of the unassimilated food is feces and half is SDA
    expect_close(sim$fluxRespirationSDA, sim$fluxFecal, 1e-12, paste(n, "SDA"))
  }
})

test_that("carbon: what is eaten from the resources is lost, died, or fished", {
  # plus the change of fish biomass (from the derivatives), as the simulations are not exactly in steady state.
  for (n in names(sims)) {
    sim <- sims[[n]]
    ix <- which(lastTimes(sim))
    eaten  <- sum(colMeans((sim$mortpred[, sim$p$ixR]*sim$R)[ix, ]))
    lost   <- sum(sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration + sim$fluxCarcass)
    fished <- sum(colMeans(sim$yield[ix, ]))
    growth <- mean(sapply(ix, function(i) sum(derivativesFEISTYR(0, sim$u[i, ], sim$p)$deriv[sim$p$ixFish])))
    expect_close(lost + fished + growth, eaten, 1e-6, paste(n, "budget"))
  }
})

test_that("carbon: the losses with effective temperature (bET) equal the losses in the model", {
  sim <- calcCarbonFluxes(simulateFEISTY(setupBasic2(), tEnd = 100))
  flux <- sum(sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration)
  known_issue(abs(flux/sum(colMeans(sim$totLoss[lastTimes(sim), ])) - 1) < 1e-6,
              "calcCarbonFluxes() does not use the effective temperature of large demersals (bET, depth < 200 m)")
})

# ---- Injection ---------------------------------------------------------------------------------

test_that("carbon: the injection conserves each flux, in total and for each stage", {
  for (n in names(sims)) {
    inject <- injections[[n]]
    for (x in pathways) {
      flux <- sims[[n]][[paste0("flux", x)]]
      expect_close(sum(inject[[x]])*rho, sum(flux), 1e-10, paste(n, x))
      expect_close(colSums(inject[[paste0(x, "_stage")]])*rho, flux, 1e-10, paste(n, x, "stages"))
      expect_close(rowSums(inject[[paste0(x, "_stage")]]), inject[[x]], 1e-10, paste(n, x, "depths"))
    }
    expect_close(inject$total, inject$Fecal + inject$Carcass + inject$Repro + inject$Respiration,
                 1e-12, paste(n, "total"))
    expect_close(inject$RespirationBasal_stage + inject$RespirationSDA_stage + inject$RespirationRepro_stage,
                 inject$Respiration_stage, 1e-12, paste(n, "respiration"))
  }
})

test_that("carbon: the injection profile is well-formed", {
  for (n in names(sims)) {
    inject <- injections[[n]]
    expect_identical(as.numeric(inject$z), as.numeric(0:sims[[n]]$p$bottom))
    expect_length(inject$total, length(inject$z))
    expect_true(all(is.finite(inject$total)), label = n)
    expect_true(all(inject$total >= -1e-12), label = n)
  }
})

# A simulation with only one size class at depth z0 that only has one flux. Same temperature everywhere.
pointSource <- function(flux, z0, Temp = 10, depth = 3000) {
  p <- setupVertical2(depth = depth, Tp = Temp, Tm = Temp, Tb = Temp)
  p$depthDay[] <- 0
  p$depthDay[z0 + 1, ] <- 1
  p$depthNight <- p$depthDay
  zero <- rep(0, length(p$ixFish))
  sim <- list(p = p)
  for (x in c(pathways, "RespirationBasal", "RespirationSDA", "RespirationRepro"))
    sim[[paste0("flux", x)]] <- zero
  sim[[paste0("flux", flux)]][1] <- 1
  if (flux == "Respiration") sim$fluxRespirationBasal[1] <- 1
  calcCarbonInjection(sim)
}

test_that("carbon: sinking detritus is remineralized as in the analytical solution", {
  z0 <- 400; depth <- 3000
  speed <- c(Fecal = 1000, Carcass = 2000, Repro = 1000)   # m/day
  for (x in names(speed)) for (Temp in c(0, 10, 20)) {
    inject <- pointSource(x, z0, Temp, depth)
    profile <- inject[[x]]*rho
    k <- 0.25*2^((Temp - 10)/10) / speed[[x]]   # per meter
    label <- sprintf("%s at %g degrees", x, Temp)
    expect_true(all(profile[inject$z <= z0] == 0), label = paste(label, "nothing above the source"))
    water <- inject$z > z0 & inject$z < depth
    expect_close(profile[water], k*exp(-k*(inject$z[water] - z0)), 2e-3, paste(label, "profile"))
    expect_close(profile[length(profile)], exp(-k*(depth - z0)), 2e-3, paste(label, "reaching the bottom"))
  }
})

test_that("carbon: respiration is injected where the fish are", {
  inject <- pointSource("Respiration", 400)
  expect_equal(inject$Respiration[401]*rho, 1)
  expect_equal(sum(inject$Respiration[-401]), 0)
})

test_that("carbon: carcasses and cold water give deeper injection", {
  meanDepth <- function(...) { x <- pointSource(...); sum(x$z*x$total)/sum(x$total) }
  expect_gt(meanDepth("Carcass", 100), meanDepth("Fecal", 100))
  expect_gt(meanDepth("Fecal", 100, Temp = 2), meanDepth("Fecal", 100, Temp = 20))
  expect_gt(meanDepth("Fecal", 600), meanDepth("Fecal", 100))
})

test_that("carbon: injection with setupVertical", {
  sim <- calcCarbonFluxes(simulateFEISTY(setupVertical(), tEnd = 10))
  known_issue(!inherits(try(calcCarbonInjection(sim), silent = TRUE), "try-error"),
              "calcCarbonInjection() fails for setupVertical, which has no Tp, Tm, and Tb")
})

test_that("carbon: injections agree with the reference", {
  ref <- utils::read.csv(test_path("fixtures", "reference_carbon.csv"))
  for (n in names(sims)) {
    s <- carbonSummary(sims[[n]], injections[[n]])
    r <- ref[ref$case == n, names(s)]
    expect_equal(nrow(r), 1)
    expect_equal(unlist(s), unlist(r), tolerance = 2e-3, info = n)
  }
})

# ---- Projection onto the grid of the transport matrix ------------------------------------------

column <- function(cell) TM$M3d[which.min((TM$grid$yt - cell$lat)^2), which.min((TM$grid$xt - cell$lon)^2), ]

test_that("carbon: the projection to the transport matrix conserves the injection", {
  for (n in names(sims)) {
    cell <- cells[cells$id == sub(" F=.*", "", n), ]
    inject <- injections[[n]]
    res <- FEISTY:::project_injection_to_TM(inject, cell$lat, cell$lon, TM$grid)
    expect_close(sum(res$inject), sum(inject$total), 1e-10, n)
    # In each layer:
    for (k in seq_along(TM$grid$zt)) {
      top <- TM$grid$zw[k]
      ix <- inject$z <= top + TM$grid$dzt[k] & (if (k == 1) inject$z >= top else inject$z > top)
      expect_equal(res$inject[k], sum(inject$total[ix]), tolerance = 1e-8, info = paste(n, "layer", k))
    }
    # Dry layers are moved to the deepest wet layer. Nothing if it is all dry.
    wet <- column(cell)
    res <- FEISTY:::project_injection_to_TM(inject, cell$lat, cell$lon, TM$grid, M3d_col = wet)
    expect_true(all(res$inject[wet == 0] == 0), label = n)
    expect_close(sum(res$inject) + 1, (if (any(wet == 1)) sum(inject$total) else 0) + 1, 1e-10, paste(n, "wet"))
  }
})

test_that("carbon: injection below the sea floor of the transport matrix goes to the deepest wet layer", {
  grid <- TM$grid
  inject <- list(z = 0:1000, total = numeric(1001))
  inject$total[c(10, 100, 1000) + 1] <- c(1, 2, 4)
  wet <- as.numeric(grid$zt < 80)
  res <- FEISTY:::project_injection_to_TM(inject, 0, 180, grid, M3d_col = wet)
  expect_equal(res$inject[1], 1)
  expect_equal(res$inject[max(which(wet == 1))], 6)
  expect_equal(sum(res$inject), 7)
})

test_that("carbon: a depth on the boundary between two layers belongs to the upper layer", {
  grid <- list(zt = c(5, 15, 25), zw = c(0, 10, 20), dzt = c(10, 10, 10), xt = 0, yt = 0)
  inject <- list(z = 0:30, total = replace(numeric(31), c(0, 10, 11, 20, 30) + 1, c(1, 2, 4, 8, 16)))
  expect_equal(FEISTY:::project_injection_to_TM(inject, 0, 0, grid)$inject, c(1 + 2, 4 + 8, 16))
})

# ---- Sequestration -----------------------------------------------------------------------------

# A small ocean of 2 x 2 columns with 3 layers of the same thickness. Columns do not interact. The water 
# mixes with the layer above with the rates k12 and k23, and the surface is a strong sink.
# In steady state: c2 = (q2 + q3)/k12, and c3 = c2 + q3/k23
smallOcean <- function(k12, k23, dz = 100, area = 1e12) {
  d <- c(2, 2, 3)
  id <- array(seq_len(prod(d)), d)
  A <- Matrix::Matrix(0, prod(d), prod(d), sparse = TRUE)
  for (i in 1:2) for (j in 1:2) {
    b <- id[i, j, ]
    A[b[1], b[1]] <- -k12 - 1e10; A[b[1], b[2]] <- k12
    A[b[2], b[1]] <- k12;         A[b[2], b[2]] <- -k12 - k23; A[b[2], b[3]] <- k23
    A[b[3], b[2]] <- k23;         A[b[3], b[3]] <- -k23
  }
  list(A = Matrix::lu(A), M3d = array(1, d), msk = list(pkeep = seq_len(prod(d)), hkeep = 1:4),
       grid = list(xt = c(1, 3), yt = c(-1, 1), zt = dz*(1:3 - 0.5), zw = dz*(0:2), dzt = rep(dz, 3),
                   Areat = matrix(area, 2, 2), DXT3d = array(sqrt(area), d), DYT3d = array(sqrt(area), d),
                   DZT3d = array(dz, d)))
}

test_that("carbon: sequestration agrees with the analytical solution", {
  k12 <- 0.05; k23 <- 0.002; dz <- 100; area <- 1e12
  ocean <- smallOcean(k12, k23, dz, area)
  inject <- array(0, c(2, 2, 3))   # gC/m2/yr
  inject[1, 1, ] <- c(5, 2, 1)
  inject[2, 1, ] <- c(0, 0, 3)
  inject[1, 2, ] <- c(7, 0, 0)
  res <- calcCarbonSequestration(ocean, inject, photic = dz)

  q <- inject/dz
  c2 <- (q[, , 2] + q[, , 3])/k12
  c3 <- c2 + q[, , 3]/k23
  expect_equal(res$Cseq[, , 2], c2, tolerance = 1e-6)
  expect_equal(res$Cseq[, , 3], c3, tolerance = 1e-6)
  expect_lt(max(abs(res$Cseq[, , 1])), 1e-6)    # the surface does not sequester
  expect_equal(res$Cseq_per_area, (c2 + c3)*dz, tolerance = 1e-6)
  expect_equal(as.numeric(res$TotInject), sum(inject)*area/1e15, tolerance = 1e-10)
  expect_equal(as.numeric(res$TotInject_below_euphotic), sum(inject[, , 2:3])*area/1e15, tolerance = 1e-10)
  expect_equal(as.numeric(res$TotSeq), sum(c2 + c3)*dz*area/1e15, tolerance = 1e-6)
  expect_equal(as.numeric(res$TotSeqTime), sum(c2 + c3)*dz/sum(inject[, , 2:3]), tolerance = 1e-6)
})

test_that("carbon: sequestration is linear and needs injection", {
  ocean <- smallOcean(0.05, 0.002)
  set.seed(1)
  a <- array(runif(12), c(2, 2, 3)); b <- array(runif(12), c(2, 2, 3))
  seq <- function(x) as.numeric(calcCarbonSequestration(ocean, x)$TotSeq)
  expect_equal(seq(a + 2*b), seq(a) + 2*seq(b), tolerance = 1e-8)
  expect_equal(seq(0*a), 0)
})

test_that("carbon: the euphotic depth can vary between columns", {
  ocean <- smallOcean(0.05, 0.002)
  inject <- array(1, c(2, 2, 3))
  photic <- matrix(c(10, 100, 200, 300), 2, 2)   # 3, 2, 1, and 0 layers have centers below
  res <- calcCarbonSequestration(ocean, inject, photic = photic)
  expect_equal(res$inject_below_euphotic, matrix(c(3, 2, 1, 0), 2, 2))
  expect_equal(as.numeric(calcCarbonSequestration(ocean, inject, photic = 100)$TotInject_below_euphotic),
               8*1e12/1e15)
})
