#
# Update the reference values used by the tests. Only do this after an intended change of the model;
# check the differences with "git diff" before committing. Run from the root of the repository, with the 
# new version of FEISTY installed:
#
#   Rscript tests/update_reference.R
#
library(FEISTY)
library(testthat)
source(test_path("helper-FEISTY.R"))

cases <- referenceCases()
ref <- do.call(rbind, lapply(names(cases), function(k) cbind(case = k, simSummary(cases[[k]]()))))
cells <- readCells()
global <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  do.call(rbind, lapply(c(0, 0.3), function(Fmax)
    cbind(id = cells$id[i], Fmax = Fmax, simSummary(simulateCell(cells[i, ], Fmax, if (Fmax > 0) c(1, 3, 5))))))))

carbon <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  do.call(rbind, lapply(c(0, 0.3), function(Fmax) {
    sim <- calcCarbonFluxes(simulateCell(cells[i, ], Fmax, if (Fmax > 0) c(1, 3, 5), tEnd = 300))
    cbind(case = sprintf("%s F=%g", cells$id[i], Fmax), carbonSummary(sim, calcCarbonInjection(sim)))
  }))))

write <- function(x, file) {
  ix <- !(names(x) %in% c("case", "id", "Fmax", "group"))
  x[ix] <- signif(x[ix], 8)
  write.csv(x, test_path("fixtures", file), row.names = FALSE)
}
write(ref, "reference.csv")
write(global, "reference_global.csv")
write(carbon, "reference_carbon.csv")
