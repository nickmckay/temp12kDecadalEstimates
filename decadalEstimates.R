library(lipdR)
library(geoChronR)
library(tidyverse)
if(!file.exists("loaded.RData")){
  D <- readLipd("~/GitHub/Temperature12k/ScientificDataAnalysis/lipdFilesWithEnsembles/")
  save(list = c("D"),file = "loaded.RData")
}else{
  load("loaded.RData")
}

#extract timeseries
TS <- extractTs(D,whichtables = "meas") #the ensembles have already
ts <- ts2tibble(TS)



addIn <- function(L){
  ts <- extractTs(L,whichtables = "meas") %>% ts2tibble()
  #add in interpretations
  ensVars <- which(grepl(pattern = "_ensemble",ts$paleoData_variableName))
  if(length(ensVars) == 0){
    return(L)
  }
  if(! "temp12kTempEnsembleSource" %in% names(ts)){
    ts$temp12kTempEnsembleSource <- NA
  }
  if(! "paleoData_fromTSid" %in% names(ts)){
    ts$paleoData_fromTSid <- NA
  }

  for(i in 1:length(ensVars)){
    vname <- ts$paleoData_variableName[ensVars[i]]
    tsidi <- which(endsWith(ts$paleoData_TSid,str_remove(vname,"_ensemble")))
    if(length(tsidi) == 1){
      ts$paleoData_inCompilation[ensVars[i]] <- "Temp12kEnsemble"
      #copy over metadata
      iii <- which(grepl(pattern = "interpretation[0-9]_",names(ts)))
      ts[ensVars[i],iii] <- ts[tsidi,iii]

      ts$temp12kTempEnsembleSource[ensVars[i]] <- "fromDatabase"
      ts$paleoData_fromTSid[ensVars[i]] <-  ts$paleoData_TSid[tsidi]
    }else{
      stop("This is bad")
    }
  }

return(as.lipd(ts))
}
D2 <- vector(mode = "list",length = length(D))
for(i in 1:length(D)){
D2[[i]] <- addIn(D[[i]])
}

TS <- extractTs(D2,whichtables = "meas") #the ensembles have already
ts <- ts2tibble(TS)

dv <- seq(0,12000,by = 10)

tsg <- dplyr::filter(ts,
                     tolower(paleoData_inCompilation) == "temp12kensemble",
                     tolower(paleoData_units) == "degc")

,
                     tolower(interpretation1_seasonalityGeneral) %in% c("annual","summeronly","winteronly"))



tsNe <- dplyr::filter(ts,
                     tolower(paleoData_inCompilation) == "temp12k",
                     tolower(paleoData_units) == "degc",
                     tolower(interpretation1_seasonalityGeneral) %in% c("annual","summeronly","winteronly"))


test <- which(!tsNe$dataSetName %in% tsg$dataSetName)

unique(tsNe$dataSetName[test])

# #extract timeseries
# TS <- extractTs(D)
#
#
# sg <- pullTsVariable(TS,variable = "interpretation1_seasonalityGeneral")
# ic <- pullTsVariable(TS,"paleoData_inCompilation")
# u <- pullTsVariable(TS,"paleoData_units")
#
#
# #filter by compilation and seasonality
# tu <- which(tolower(ic) == "temp12kensemble" & (tolower(sg) == "annual" | tolower(sg) == "summeronly" | tolower(sg) == "winteronly") & tolower(u) == "degc")
#
# fTS <- TS[tu]
#


decadalEstimate <- function(paleoData_TSid,ageEnsemble,paleoData_values,...){
  dv <- seq(0,12000,by = 10)
  out <- quantile2d(x = ageEnsemble,y = paleoData_values,x.bin = dv,probs = pnorm(-1:1))
  rmse <- (out$quants[,3]-out$quants[,1])/2 #average 1 sigma range
  mse <- rmse^2
  return(list(tsid = paleoData_TSid,year = out$x.bin,median = out$quants[,2], mse = mse))
}
sameLength <- which(purrr::map_dbl(tsg$paleoData_values,nrow) == purrr::map_dbl(tsg$ageEnsemble,nrow))
diffLength <- which(purrr::map_dbl(tsg$paleoData_values,nrow) != purrr::map_dbl(tsg$ageEnsemble,nrow))

allOut <- purrr::pmap(tsg[sameLength,],decadalEstimate)
#loopy
dv <- seq(0,12000,by = 10)

msemat <- medians <- matrix(nrow = length(dv),ncol = nrow(tsg))
probs = pnorm(-1:1)
for(i in 1:nrow(tsg)){
  ageEnsemble <- tsg$ageEnsemble[[i]]
  paleoData_values <- tsg$paleoData_values[[i]]
  #out <- try(quantile2d(x = ageEnsemble,y = paleoData_values,x.bin = dv,probs = pnorm(-1:1)))
  #out <- try(bin2d(x = ageEnsemble,y = paleoData_values,x.bin = dv,interpolate = FALSE,))
  if(!is(out,"try-error")){
    rmse <- (out$quants[,3]-out$quants[,1])/2 #average 1 sigma range
    mse <- rmse^2

    medians[,i] <- out$quants[,2]
    msemat[,i] <- mse
  }
  print(i)
}

medTib <- tibble::as.tibble(medians) %>% setNames(tsg$paleoData_TSid)
mseTib <- tibble::as.tibble(msemat) %>% setNames(tsg$paleoData_TSid)

medTib <- dplyr::bind_cols(data.frame(year = seq(0,12000,by = 10)),medTib)
mseTib <- dplyr::bind_cols(data.frame(year = seq(0,12000,by = 10)),mseTib)

readr::write_csv(x = medTib,file = "medianMatrix.csv")
readr::write_csv(x = mseTib,file = "mseMatrix.csv")

write_sheet_retry(medTib,ss = "1tir8Amwr5QPxrafNq5WY17IP2xi9elgTdZwO6TfbkfE",sheet = "medians",timeout = 60)
write_sheet_retry(msemat,ss = "1tir8Amwr5QPxrafNq5WY17IP2xi9elgTdZwO6TfbkfE",sheet = "mse")

tsg$


