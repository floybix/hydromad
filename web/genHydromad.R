## generate the web site

## setwd("X:/Packages/hydromad/web")
## Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs8.63/bin/gswin32c.exe")
## source("genHydromad.R", echo = TRUE)

library(hydromad)

spec <- list()

spec[["modelling framework"]] <-
    list(
         list("hydromad", do.example = TRUE),
         list("methods...", helpname = "hydromad.object",
              codefile = "update.hydromad.R"),
         list("simulate", helpname = "simulate.hydromad", do.example = TRUE),
         list("runlist"),
         list("hydromad.options", codefile = "options.R")
         )

spec[["assessment"]] <-
    list(
         list("summary", helpname = "summary.hydromad"),
         list("eventseq", do.example = TRUE),
         list("xyplot", helpname = "xyplot.hydromad", do.example = TRUE),
         list("xyplot.runlist", do.example = TRUE),
         list("errormasscurve", -1, helpname = "xyplot.hydromad", do.example = TRUE),
         list("qqmath", -2, helpname = "xyplot.hydromad", do.example = TRUE),
         list("fitStat")
         )

spec[["calibration"]] <-
    list(
         list("objFunVal"),
         list("fitBySampling"),
         list("fitByOptim"),
         list("fitBySCE"),
         list("fitByDE"),
         list("fitByDream")
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
         list("armax", codefile = "armax.sim.R"),
         list("expuh", codefile = "expuh.sim.R"),
         list("lambda", codefile = "lambda.sim.R"),
         list("powuh", codefile = "powuh.sim.R")
         )

spec[["routing fitting"]] <-
    list(
         list("armax.ls.fit"),
         list("armax.sriv.fit"),
         list("armax.inverse.fit"),
         #"lambda.inverse.fit"
         list("tryModelOrders"),
         list("estimateDelay")
         )

spec[["utilities"]] <-
    list(
         list("convertFlow"),
         list("rollccf", do.example = TRUE),
         list("parameterSets"),
         list("observed"),
         list("SCEoptim")
         )

spec[["datasets"]] <-
    list(
         list("BinghamTrib", codefile = "../data/BinghamTrib.R", do.example = TRUE),
         list("Canning", codefile = NA, do.example = TRUE),
         list("Cotter", codefile = NA, do.example = TRUE),
         #list("Molonglo"),
         list("Murrindindi", codefile = NA, do.example = TRUE),
         list("Queanbeyan", codefile = NA, do.example = TRUE),
         list("SalmonBrook", codefile = NA, do.example = TRUE),
         list("Wye", codefile = NA, do.example = TRUE),
         list("HydroTestData", codefile = "../data/HydroTestData.R", do.example = TRUE)
         )

#spec[["wetlands"]] <-
#    list(
#         list("swimp"),
#         list("poweroid")
#         )

#source("http://latticeextra.r-forge.r-project.org/generate.R")
source("../../latticeextra/www/generate.R")

## stop on errors
lattice.options(panel.error = NULL)

generateWebsite("hydromad", spec = spec, 
                man.src.dir = "../man/",
                imageSrcBase = "",
                codeSrcSpec = "http://github.com/floybix/hydromad/tree/master/R/%s",
                themes = list(custom_theme_2 = custom.theme.2()),
                do.examples = FALSE)
