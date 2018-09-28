# Test CreateOmics wrapper
# source Test_ResponseConversion.R first

CreateOmics <- function(assayData_df,
                         pathwayCollection_ls,
                         response = NULL,
                         respType = c("none", "survival", "regression", "categorical")){

  # browser()

  respType <- match.arg(respType)

  if(respType != "none" && is.null(response)){
    stop(paste0("Response must be specified for type = ", respType))
  }
  if(respType == "none" && !is.null(response)){
    stop("Response type required when a response is given.")
  }

  if(!is.null(response)){
    respClean <- convertResponse(response, type = respType)
  }

  switch (respType,
    none = {

      message("Creating object of class OmicsPathway.")
      CreateOmicsPath(assayData_df = assayData_df,
                       pathwayCollection_ls = pathwayCollection_ls)

    },
    survival = {

      message("Creating object of class OmicsSurv.")
      CreateOmicsSurv(assayData_df = assayData_df,
                       pathwayCollection_ls = pathwayCollection_ls,
                       eventTime_num = respClean$time,
                       eventObserved_lgl = respClean$dead)

    },
    regression = {

      message("Creating object of class OmicsReg.")
      CreateOmicsReg(assayData_df = assayData_df,
                      pathwayCollection_ls = pathwayCollection_ls,
                      response_num = respClean)

    },
    categorical = {

      message("Creating object of class OmicsCateg.")
      CreateOmicsCateg(assayData_df = assayData_df,
                        pathwayCollection_ls = pathwayCollection_ls,
                        response_fact = respClean)

    }
  )

}


# Test
CreateOmics()
CreateOmics(response = "12")
CreateOmics(respType = "s")
# Pathway
CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
             pathwayCollection_ls = gene_set_ls)
# Survival
CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
             pathwayCollection_ls = gene_set_ls,
             response = joinedExperiment_df[, 2:3],
             respType = "s")
# Regression
CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
             pathwayCollection_ls = gene_set_ls,
             response = joinedExperiment_df[, 2],
             respType = "r")
# Categorical
CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
             pathwayCollection_ls = gene_set_ls,
             response = joinedExperiment_df[, 3],
             respType = "c")
