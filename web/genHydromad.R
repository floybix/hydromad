## generate the web site

## setwd("X:/Packages/hydromad/web")
## Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs8.63/bin/gswin32c.exe")
## source("genHydromad.R", echo = TRUE)

library(hydromad)

spec <- list()

spec[["modelling framework"]] <-
    list(
         list("hydromad"),
         list("methods...", helpname = "hydromad.object", do.example = FALSE,
              codefile = "update.hydromad.R"),
         list("simulate", helpname = "simulate.hydromad"),
         list("runlist", do.example = FALSE),
         list("hydromad.options", do.example = FALSE, codefile = "options.R")
         )

spec[["assessment"]] <-
    list(
         list("summary", helpname = "summary.hydromad", do.example = FALSE),
         list("eventseq"),
         list("xyplot", helpname = "xyplot.hydromad"),
         list("xyplot.runlist"),
         list("errormasscurve", helpname = "xyplot.hydromad"),
         list("qqmath", helpname = "xyplot.hydromad"),
         list("fitStat", do.example = FALSE)
         )

spec[["calibration"]] <-
    list(
         list("objFunVal", do.example = FALSE),
         list("fitBySampling", do.example = FALSE),
         list("fitByOptim", do.example = FALSE),
         list("fitBySCE", do.example = FALSE),
         list("fitByDE", do.example = FALSE),
         list("fitByDream", do.example = FALSE)
         )

spec[["soil moisture accounting"]] <-
    list(
         list("scalar"),
         list("cwi", helpname = "IHACRES.CWI.model", codefile = "cwi.R"),
         list("cmd", helpname = "IHACRES.CMD.model", codefile = "cmd.R"),
         list("bucket"),
         list("awbm"),
         list("sacramento", height = 700),
         list("snow", height = 500),
         list("runoffratio"),
         list("dbm")
         )

spec[["routing"]] <-
    list(
         list("armax", codefile = "armax.sim.R", do.example = FALSE),
         list("expuh", codefile = "expuh.sim.R", do.example = FALSE),
         list("lambda", codefile = "lambda.sim.R", do.example = FALSE),
         list("powuh", codefile = "powuh.sim.R", do.example = FALSE)
         )

spec[["routing fitting"]] <-
    list(
         list("armax.ls.fit", do.example = FALSE),
         list("armax.sriv.fit", do.example = FALSE),
         list("armax.inverse.fit", do.example = FALSE),
         #"lambda.inverse.fit"
         list("tryModelOrders", do.example = FALSE),
         list("estimateDelay", do.example = FALSE)
         )

spec[["utilities"]] <-
    list(
         list("convertFlow", do.example = FALSE),
         list("rollccf"),
         list("parameterSets", do.example = FALSE),
         list("observed", do.example = FALSE),
         list("SCEoptim", do.example = FALSE)
         )

spec[["datasets"]] <-
    list(
         list("BinghamTrib", codefile = "../data/BinghamTrib.R"),
         list("Canning", codefile = NA),
         list("Cotter", codefile = NA),
         #list("Molonglo"),
         list("Murrindindi", codefile = NA),
         list("Queanbeyan", codefile = NA),
         list("SalmonBrook", codefile = NA),
         list("Wye", codefile = NA),
         list("HydroTestData", codefile = "../data/HydroTestData.R")
         )

#spec[["wetlands"]] <-
#    list(
#         list("swimp"),
#         list("poweroid")
#         )

source("http://latticeextra.r-forge.r-project.org/generate.R")
source("../../latticeextra/www/generate.R")

## stop on errors
lattice.options(panel.error = NULL)

generateWebsite("hydromad", spec = spec, 
                man.src.dir = "../man/",
                imageSrcBase = "",
                codeSrcSpec = "http://github.com/floybix/hydromad/tree/master/R/%s",
                themeNames = "custom_theme_2",
                do.examples = TRUE)
