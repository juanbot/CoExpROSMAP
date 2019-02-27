
#' initDb - Initialization of the package so the ROSMAP networks can be used with
#' CoExpNets
#'
#'
#' @param mandatory If this parameter is `TRUE` then the networks will be added no matter whether they were already there.
#' 
#'
#' @return No value 
#' @export
#'
#' @examples
initDb = function(mandatory=F){
  the.dir = system.file("", "extdata", package = "CoExpROSMAP")
  tissues = c("notad","probad","ad","allsamples")
  nets = c("netnotad.8.it.50.rds",
           "netprobad.11.it.50.rds",
           "netad.8.it.50.rds",
           "netallsamples.6.it.50.rds")
  residuals = c("fpkm.casectrl.qc.qn.combat.covs.svas2.res.cogdx.notad.rds",
                "fpkm.casectrl.qc.qn.combat.covs.svas2.res.cogdx.probad.rds",
                "fpkm.casectrl.qc.qn.combat.covs.svas2.res.cogdx.ad.rds",
                "fpkm.casectrl.qc.qn.combat.covs.svas2.res.rds")
  for(tissue in tissues){
    net = nets[which(tissues == tissue)]
    net = paste0(the.dir,"/",net)
    residual = residuals[which(tissues == tissue)]
    CoExpNets::addNet(which.one="CoExpROSMAP",
           tissue=tissue,
           netfile=net,
           ctfile=paste0(net,".celltype.csv"),
           gofile=paste0(net,"_gprof.csv"),
           exprdatafile=paste0(the.dir,"/",
                              residual),
           overwrite=mandatory)
  }
}

#' Title
#'
#' @param tissue 
#' @param which.one 
#'
#' @return
#' @export
#'
#' @examples
getCovariates = function(tissue=tissue,which.one="CoExpROSMAP"){
 
  expr.data = CoExpNets::getExprDataFromTissue(tissue=tissue,which.one=which.one)
  the.dir = system.file("", "extdata", package = "CoExpROSMAP")
  key = utils::read.csv(paste0(the.dir,"/ROSMAP_IDkey.csv"))
  covs = utils::read.csv(paste0(the.dir,"/ROSMAP_clinical.csv"))
  ids = rownames(expr.data)
  mask = key$projid[match(ids,key$mrna_id)]
  nonmatchingids = ids[is.na(mask)]
  goodids = NULL
  for(id in nonmatchingids){
    subids = stringr::str_split(id,"_")
    recid = paste0(subids[[1]][1],"_",subids[[1]][2])
    goodids = c(goodids,recid)
  }
  mask[is.na(mask)] = key$projid[match(goodids,key$mrna_id)]
  samples = mask
  #samples = rosmap.fromRNAseqID2ProjectID(rownames(expr.data))
  gender = as.factor(covs$msex[match(samples,covs$projid)])
  pmi = as.numeric(covs$pmi[match(samples,covs$projid)])
  braaksc = as.factor(covs$braaksc[match(samples,covs$projid)])
  cogdx = as.factor(covs$cogdx[match(samples,covs$projid)])
  educ = as.numeric(covs$educ[match(samples,covs$projid)])
  ceradsc = as.factor(covs$ceradsc[match(samples,covs$projid)])
  age = as.character(covs$age_death[match(samples,covs$projid)])
  
  #Impute
  pmi[is.na(pmi)] = mean(pmi[!is.na(pmi)])
  age[grep("90\\+",age)] = "90"
  age = as.numeric(age)
  race = as.factor(covs$race[match(samples,covs$projid)])
  
  batch = stringr::str_split(rownames(expr.data),"_")
  batch = as.factor(unlist(lapply(batch,function(x){return(x[[3]])})))
  
  toreturn = data.frame(batch,gender,pmi,age,race,braaksc,cogdx,educ,ceradsc)
  
  rownames(toreturn) = rownames(expr.data)
  return(toreturn)
}
