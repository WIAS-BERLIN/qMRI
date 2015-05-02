setClass("mpmData",
         representation(t1Files = "character",
                        mtFiles = "character",
                        pdFiles = "character",
                        TR = "numeric",
                        TE = "numeric",
                        FA = "numeric",
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
                        TEScale = "numeric",
                        dataScale = "numeric",
                        mask = "array",
                        model = "integer",
                        sdim = "integer",
                        nFiles = "integer"),
         contains = c("array")
)
