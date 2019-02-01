
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
    residual = residuals[which(tissues == tissue)]
    CoExpNets::addNet(which.one="CoExpROSMAP",
           tissue=tissue,
           netfile=net,
           ctfile=paste0(the.dir,"/",net,"_celltype.csv"),
           gofile=paste0(the.dir,"/",net,"_gprof.csv"),
           exprdatafile=paste0(the.dir,"/",
                              residual),
           overwrite=mandatory)
  }
}
