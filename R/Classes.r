setClass("mpmData",
         representation(t1Files = "character",
                        mtFiles = "character",
                        pdFiles = "character",
                        TR = "numeric",
                        TE = "numeric",
                        FA = "numeric",
                        smoothPar = "numeric",
                        maskFile = "character",
                        mask = "array",
                        model = "integer",
                        sdim = "integer",
                        nFiles = "integer"),
         contains = c("array")
)

setClass("ESTATICSModel",
         representation(isConv = "array",
                        invCov = "array",
                        bi = "array",
                        TEScale = "numeric",
                        dataScale = "numeric",
                        smoothPar = "numeric",
                        mask = "array",
                        model = "integer",
                        sdim = "integer",
                        nFiles = "integer"),
         contains = c("array")
)
