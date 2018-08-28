# Test create_Omics wrapper
# source Test_ResponseConversion.R first

create_Omics <- function(assayData_df,
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
      create_OmicsPath(assayData_df = assayData_df,
                       pathwayCollection_ls = pathwayCollection_ls)

    },
    survival = {

      message("Creating object of class OmicsSurv.")
      create_OmicsSurv(assayData_df = assayData_df,
                       pathwayCollection_ls = pathwayCollection_ls,
                       eventTime_num = respClean$time,
                       eventObserved_lgl = respClean$dead)

    },
    regression = {

      message("Creating object of class OmicsReg.")
      create_OmicsReg(assayData_df = assayData_df,
                      pathwayCollection_ls = pathwayCollection_ls,
                      response_num = respClean)

    },
    categorical = {

      message("Creating object of class OmicsCateg.")
      create_OmicsCateg(assayData_df = assayData_df,
                        pathwayCollection_ls = pathwayCollection_ls,
                        response_fact = respClean)

    }
  )

}


# Test
create_Omics()
create_Omics(response = "12")
create_Omics(respType = "s")
# Pathway
create_Omics(assayData_df = joinedExperiment_df[, -c(1:3)],
             pathwayCollection_ls = gene_set_ls)
# Survival
create_Omics(assayData_df = joinedExperiment_df[, -c(1:3)],
             pathwayCollection_ls = gene_set_ls,
             response = joinedExperiment_df[, 2:3],
             respType = "s")
# Regression
create_Omics(assayData_df = joinedExperiment_df[, -c(1:3)],
             pathwayCollection_ls = gene_set_ls,
             response = joinedExperiment_df[, 2],
             respType = "r")
# Categorical
create_Omics(assayData_df = joinedExperiment_df[, -c(1:3)],
             pathwayCollection_ls = gene_set_ls,
             response = joinedExperiment_df[, 3],
             respType = "c")
