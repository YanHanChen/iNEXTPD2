#' Interpolation (rarefaction) and extrapolation of Chao et al.’s (2010) phylogenetic diversity and mean phylogenetic diversity (phylogenetic Hill numbers)
#'
#' Function \code{iNEXTPD} computes phylogenetic diversity estimates for rarefied samples and extrapolated samples
#' along with confidence intervals and related coverage estimates based on Chao et al.’s (2010) phylogenetic
#' diversity (PD) and mean phylogenetic diversity (meanPD). See Chao et al. (2010, 2015, 2016) and
#' Hsieh and Chao (2017) for pertinent background and methodologies.
#' @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).\cr
#' Abundance data: a species-by-site matrix/data.frame of species abundances. The row (species) names of
#' data must match the species names in the phylogenetic tree and thus cannot be missing.\cr
#' Incidence raw data: species-by-site raw incidence matrix/data.frame. When there are N assemblages
#' and thus N matrices, users must first merge the N matrices by species identity to obtain a large
#' merged incidence matrix, where the rows of the matrix refer to all species presented in the merged
#' data. The row (species) names of data must match the species names in the phylogenetic tree and
#' thus cannot be missing.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
#' @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
#' then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#' the pooled assemblage. Default is \code{NULL}.
#' @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
#' @param endpoint a positive integer specifying the endpoint for the rarefaction and extrapolation range.
#' If \code{NULL}, then \code{endpoint} = double of the reference sample size in each assemblage. It is ignored if \code{size} is given.
#' @param knots a positive integer specifying the number of equally-spaced knots between 1 and the \code{endpoint}. Default is 40.
#' @param size a sequence of positive integers specifying the sample sizes for which PD or meanPD estimates will be calculated.
#' If \code{NULL}, then estimates will be calculated for those sample sizes determined by the specified/default \code{endpoint}
#' and \code{knots}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
#' Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @import ape
#' @import ggplot2
#' @import dplyr
#' @import tidytree
#' @importFrom stats rbinom
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom phyclust get.rooted.tree.height
#' @return Returns a list containing two tables:
#' \itemize{
#'  \item{\code{$size_based}: size-based PD or meanPD estimates along with their confidence intervals
#'  (if \code{nboot > 0}) and relevant statistics information.} \cr
#'  \item{\code{$coverage_based}: coverage-based diversity estimates along with confidence intervals
#'  (if \code{nboot > 0}) and relevant statistics information.}
#'  }
#' @examples
#' \donttest{
#' # Datatype: abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- iNEXTPD(data = data, tree = tree, datatype = "abundance", q = c(0, 1, 2), nboot = 30)
#' out
#'
#' # Datatype: incidence_raw data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out <- iNEXTPD(data = data, nT = nT, datatype = "incidence_raw", tree = tree, q = c(0, 1, 2), nboot = 30)
#' out
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
#' Chao, A., Chiu C.-H. and Jost L. (2016). Phylogenetic diversity measures and their decomposition: a framework based on Hill numbers. pp. 141-172 in Pellens R. and Grandcolas P. (eds)
#' \emph{Biodiversity Conservation and Phylogenetic Systematics: Preserving our Evolutionary Heritage in an Extinction Crisis}, Springer. \cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
#' @export
iNEXTPD <- function(data,nT,datatype = "abundance",tree,q = c(0,1,2),reftime=NULL,type = 'PD',endpoint = NULL,knots = 40,size = NULL,nboot = 50,conf = 0.95) {
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  # qq <- 0:2
  # if(is.na(pmatch(q, qq) == T) == T) stop("invalid order of q, we only compute q = 0, 1 or 2", call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of key in sampling units", call. = FALSE)
    n <- 0
    for(i in 1:length(nT)){
      mydata[[i]] <- data[,(n+1):(n+nT[i])]
      n <- n+nT[i]
    }
    if(is.null(names(nT))) {
      names(mydata) <- paste0("assemblage",1:length(nT))
    }else{
      names(mydata) = names(nT)
    }
  }else{
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  
  ###########output1
  #atime <- Sys.time()
  # if(class(mydata) == "list"){
  #   infos <- lapply(mydata, function(x){
  #     datainf(data = x, datatype, phylotr = mytree,reft = reft) %>% mutate(Reference.time = reft)
  #     }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(names(mydata),each = length(reft))) %>%
  #     select(Assemblage,n,S.obs,PD.obs,`f1*`,`f2*`,g1,g2,Reference.time)
  # }else{
  #   return(NULL)
  # }
  #btime <- Sys.time()
  #print(paste0('Info time:',btime-atime))
  ############output2
  if(length(knots)!=length(mydata)) knots <- rep(knots,length(mydata))
  if(is.null(size)){
    if(is.null(endpoint)){
      if(datatype == "abundance") {
        endpoint <- sapply(mydata, function(x) 2*sum(x))
      }else if(datatype == "incidence_raw"){
        endpoint <- sapply(mydata, function(x) 2*ncol(x))
      }
    }else{
      if(length(endpoint)!=length(mydata)){
        endpoint <- rep(endpoint,length(mydata))
      }
    }
    size <- lapply(1:length(mydata),function(i){
      if(datatype == "abundance") {
        ni <- sum(mydata[[i]])
      }else if(datatype == "incidence_raw"){
        ni <- ncol(mydata[[i]])
      }
      
      if(endpoint[i] <= ni){
        mi <- floor(seq(1,endpoint[i],length.out = knots[i]))
      }else{
        mi <- floor(c(seq(1,ni,length.out = floor(knots[i]/2)),
                      seq(ni+1,endpoint[i],length.out = knots[i]-floor(knots[i]/2))))
      }
      unique(mi)
    })
  }else{
    if(class(size)=="numeric"|class(size)=="integer"|class(size)=="double"){
      size <- list(size = size)
    }
    if(length(size)!=length(mydata)) size <- lapply(1:length(mydata), function(x) size[[1]])
    size <- lapply(1:length(mydata),function(i){
      if(datatype == "abundance") {
        ni <- sum(mydata[[i]])
      }else if(datatype == "incidence_raw"){
        ni <- ncol(mydata[[i]])
      }
      if( sum(size[[i]] == ni) == 0 ) mi <- sort(c(ni,size[[i]]))
      else mi <- size[[i]]
      unique(mi)
    })
  }
  
  FUN <- function(e){
    if(class(mydata)=="list"){
      # temp = inextPD(datalist=mydata,phylotr=mytree,datatype,q,nboot,conf,m = size,reftime,cal = type)
      inextPD(datalist = mydata,datatype = datatype,phylotr = mytree,q = q,reft = reftime,m=size,
              cal = type,nboot=nboot,conf = conf,unconditional_var = TRUE)
    }else{
      return(NULL)
    }
  }
  #atime <- Sys.time()
  #RE.table <- FUN(3)
  ans <- tryCatch(FUN(e), error = function(e){return()})
  #btime <- Sys.time()
  #print(paste0('R/E time:',btime-atime))
  ###############output3
  # TYPE <- 1:3
  # FUN2 <- function(e){
  #   if(is.na(sum(pmatch(plot.type, TYPE))) == F){
  #     temp2 <- lapply(plot.type, function(j) RE_plot(RE.table, datatype, j))
  #     allname <- c("RE.plot.size", "RE.plot.C", "RE.plot.sizeC")
  #     names(temp2) <- allname[plot.type]
  #     temp2
  #   }else{
  #     return("invalid plot type", call. = FALSE)
  #   }
  # }
  # #atime <- Sys.time()
  # RE.plot <- tryCatch(FUN2(e), error = function(e){return()})
  #btime <- Sys.time()
  #print(paste0('plot time:',btime-atime))
  # ans <- list(summary = infos, inext = RE.table, figure = RE.plot)
  #class(ans) <- c("iNEXTPD")
  ans
  
}

#' Computes asymptotic estimates for phylogenetic diversity and mean phylogenetic diversity (phylogenetic Hill numbers)
#'
#' Function \code{PhdAsy} computes asymptotic phylogenetic diversity estimates with respect to specified/default
#' diversity order q and reference time to infer true phylogenetic diversity (PD) or phylogenetic Hill numbers (meanPD). See Chao et al. (2015) and Hsieh and Chao (2017) for the statistical estimation detail.
#' @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).
#' See the function \code{\link{iNEXTPD}} for details.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
#' @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{seq(0, 2, by = 0.25)}.
#' @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
#' then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#' the pooled assemblage. Default is \code{NULL}.
#' @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
#' Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @import ape
#' @import ggplot2
#' @import dplyr
#' @import tidytree
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom phyclust get.rooted.tree.height
#' @return Returns a table of estimated asymptotic phylogenetic diversity estimates (\code{type = "PD"}) or
#' phylogenetic Hill numbers (\code{type = "meanPD"}) with respect to specified/default order \code{q} and
#' reference time specified in the argument \code{reftime}.
#' @examples
#' # Datatype: abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdAsy(data = data, datatype = "abundance", tree = tree,
#' q = seq(0, 2, by = 0.25), nboot = 30)
#' out
#'
#' # Datatype: incidence_raw data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out <- PhdAsy(data = data, nT = nT, datatype = "incidence_raw",
#' tree = tree, q = seq(0, 2, by = 0.25))
#' out
#' 
#' @references
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
#' @export
PhdAsy <- function(data,nT,datatype = "abundance",tree,q = seq(0,2,by = 0.25),reftime = NULL,type = 'PD',nboot = 50,conf = 0.95){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  #if (length(q) == 1) stop("length of q should be greater than one", call. = FALSE)
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.", call. = FALSE)
  if (sum(q<0)>=1) stop("q must be a positive number", call. = FALSE)
  if ((conf < 0) | (conf < 0) | (is.numeric(conf)==F)) stop('conf"(confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of key(nT) in sampling units", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(nT)){
      mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
      ntmp <- ntmp+nT[i]
    }
    if(is.null(names(nT))) {
      names(mydata) <- paste0("assemblage",1:length(nT))
    }else{
      names(mydata) = names(nT)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  
  # if(class(mydata) == "list"){
  #   infos <- sapply(mydata, function(x){
  #     datainf(data = x, datatype, phylotr = mytree,reft = reft)})
  # }else{
  #   return(NULL)
  # }
  
  FUN = function(e){
    AsyPD(datalist = mydata,datatype = datatype,phylotr = mytree,q = q,reft = reftime,cal = type,nboot,conf)# mytree is pooled tree of class phylo
  }
  #out <- FUN(3)
  ans <- tryCatch(FUN(e), error = function(e){return()})
  return(ans)
}

#' Computes observed phylogenetic diversity and phylogenetic Hill numbers
#'
#' Function \code{PhdObs} computes empirical or observed phylogenetic diversity (PD) and phylogenetic Hill
#' numbers (meanPD, mean phylogenetic diversity) for specified/default order \code{q} and reference
#' time specified in the argument \code{reftime}. See Chao et al. (2010) for details of PD and meanPD.
#' @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).
#' See the function \code{\link{iNEXTPD}} for details.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
#' @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{seq(0, 2, by = 0.25)}.
#' @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
#' then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#' the pooled assemblage. Default is \code{NULL}.
#' @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
#' Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @import ape
#' @import ggplot2
#' @import dplyr
#' @import tidytree
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom phyclust get.rooted.tree.height
#' @return Returns a table of empirical (observed) phylogenetic diversity (\code{type = "PD"}) or phylogenetic Hill number (\code{type= "meanPD"})
#' for specified/default order q and reference time.
#' @examples
#' # Datatype: abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdObs(data = data, datatype = "abundance", tree = tree,
#' q = seq(0, 2, by = 0.25))
#' out
#'
#' # Datatype: incidence_raw data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out <- PhdObs(data = data, nT = nT, datatype = "incidence_raw",
#' tree = tree, q = seq(0, 2, by = 0.25))
#' out
#' 
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
#' @export
PhdObs <- function(data,nT,datatype = "abundance",tree,q = seq(0, 2, by = 0.25),reftime = NULL,type = "PD",
                   nboot = 50,conf = 0.95){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.", call. = FALSE)
  if (sum(q<0)>0) stop("q must be a positive number", call. = FALSE)
  # if ((profile != "q") & (profile != "time")) stop("invalid profile", call. = FALSE)
  # if (length(tprofile_times) == 1 & is.null(tprofile_times)==F) stop("length of time should be greater than one", call. = FALSE)
  # if (sum(tprofile_times<0)>=1 & is.null(tprofile_times)==F) stop("time must be a positive number", call. = FALSE)
  # if (is.null(knots) ==F) {
  #   if ((knots < 0) | (is.numeric(knots)==F) | (knots%%1>0)) {
  #     stop('knot must be a nonnegative integer, We use "knots" = 50 to calculate!', call. = FALSE)
  #   }
  # }
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of key(nT) in sampling units", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(nT)){
      mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
      ntmp <- ntmp+nT[i]
    }
    if(is.null(names(nT))) {
      names(mydata) <- paste0("assemblage",1:length(nT))
    }else{
      names(mydata) = names(nT)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  #=====old version=====
  # FUN = function(e){
  #   ###########data information
  #   if(class(mydata) == "list"){
  #     infos <- sapply(mydata, function(x){
  #       datainf(data = x, datatype, phylotr = mytree,reft = reftime)})
  #   }else{
  #     return(NULL)
  #   }
  #
  #   if(profile == "q") {
  #
  #     if(is.null(reft)){
  #       if (datatype=="incidence_raw") {
  #         da <- lapply(mydata, rowSums) %>% do.call(cbind, .) %>% rowSums()
  #       }else if (datatype=="abundance") {
  #         da <- do.call(cbind, mydata) %>% rowSums()
  #       }
  #       aL <- phyBranchAL_Abu(phylo = mytree,data = da,"abundance",
  #                             refT = reftime)$treeNabu %>%
  #         select(branch.abun,branch.length,tgroup)
  #       PD2 <- PD.qprofile(aL,q = 2, cal =  "PD",nt = sum(da))
  #       Q <- reftime-(reftime^2)/PD2
  #       reft = sort(c('Q'= Q, 'reftime' = reftime))
  #     }else{
  #       names(reft) <- NULL
  #       reft <- sort(reft)
  #     }
  #
  #     temp <- Phdqtable(datalist = mydata, phylotr = mytree, q, cal = type, datatype, nboot, conf, reft)
  #     ans <- list(summary = infos, forq_table = temp, forq_figure = Plotq(temp, type))
  #     #class(ans) <- c("PhdObs")
  #     return(ans)
  #   }
  #   if(profile == "time") {
  #     if (is.null(tprofile_times)) {
  #       tprofile_times <- seq(0.01, reftime, length.out = 15) %>% unique() %>% sort
  #     } else {
  #       tprofile_times <- c(tprofile_times, 0.01, reftime) %>% unique() %>% sort
  #     }
  #     temp <- Phdttable(datalist = mydata, phylotr = mytree, times = tprofile_times,cal = type, datatype, nboot, conf)
  #     if (knots==0) {
  #       ans <- list(summary = infos, fortime_table = temp[[1]], fortime_figure = Plott(temp[[1]], type, temp[[2]]))
  #     } else {
  #       AUC <- AUC_one_table(datalist = mydata,phylotr = mytree,knot = knots,cal = type,datatype = datatype,nboot, conf,reft_max = max(tprofile_times))
  #       ans <- list(summary = infos, fortime_table = temp[[1]], fortime_figure = Plott(temp[[1]], type, temp[[2]]), AUC_table = AUC)
  #     }
  #     #class(ans) <- c("PhdObs")
  #     return(ans)
  #   }
  # }
  #=====new version=====
  FUN <- function(e){
    EmpPD(datalist = mydata,datatype = datatype,phylotr = mytree,q = q,reft = reftime,cal = type,nboot,conf)
  }
  #temp <- FUN(3)
  ans <- tryCatch(FUN(e), error = function(e){return()})
  return(ans)
}


#' Computes phylogenetic diversity for specified values of sample coverage
#'
#' Function \code{estimatePD} computes Chao et al.’s (2010, 2016) phylogenetic diversity (PD, effective total branch lengths,
#' for diversity order q = 0, 1 and 2) and mean phylogenetic diversity (meanPD, phylogenetic Hill
#' numbers or the effective number of lineages, q = 0, 1 and 2) at specified values of sample coverage. See Chao et al. (2010, 2015) and Hsieh and Chao (2017) for formulas and interpretations.
#' Use the function \code{iNEXTPD} to compute PD or meanPD for specified sample sizes.
#' @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data). See the function \code{\link{iNEXTPD}} for details.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
#' @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
#' then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#' the pooled assemblage. Default is \code{NULL}.
#' @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
#' @param level a positive value < 1 or sequence specifying particular values of sample coverage.
#' If \code{NULL}, then \code{level} is set to be the minimum coverage value among all samples extrapolated up to double the reference sample sizes. Default is \code{NULL}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
#' Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @import ape
#' @import dplyr
#' @import tidytree
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats optimize
#' @importFrom phyclust get.rooted.tree.height
#' @return Returns a table of the computed phylogenetic diversity (PD or meanPD) for specified/default diversity orders \code{q} and reference times
#' for the user-specified values of sample coverage. The corresponding sample sizes and sample coverage values are also provided.
#' @examples
#' # Datatype: abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- estimatePD(data = data, datatype = "abundance", tree = tree)
#' out
#'
#' # Datatype: incidence_raw data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out <- estimatePD(data = data, nT = nT, datatype = "incidence_raw", tree = tree)
#' out
#' 
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
#' Chao, A., Chiu C.-H. and Jost L. (2016). Phylogenetic diversity measures and their decomposition: a framework based on Hill numbers. pp. 141-172 in Pellens R. and Grandcolas P. (eds)
#' \emph{Biodiversity Conservation and Phylogenetic Systematics: Preserving our Evolutionary Heritage in an Extinction Crisis}, Springer. \cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
#' @export
estimatePD <- function(data,nT,datatype = "abundance",tree,q = c(0,1,2),reftime=NULL,type = 'PD',level = NULL,nboot = 50,conf = 0.95){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.", call. = FALSE)
  if (sum(q<0)>=1) stop("q must be a positive number", call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf"(confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  #if (length(level)>1) stop('Currently, we only accept one fixed level of coverage.')
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(nT)){
      mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
      ntmp <- ntmp+nT[i]
    }
    if(is.null(names(nT))) {
      names(mydata) <- paste0("assemblage",1:length(nT))
    }else{
      names(mydata) = names(nT)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[pool.name,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  
  if(is.null(level)){
    if(datatype=='abundance'){
      level <- sapply(mydata,function(x){
        ni <- sum(x)
        Coverage(data = x,datatype = datatype,m = 2*ni,nt = ni)
      })
      
    }else if(datatype=='incidence_raw'){
      level <- sapply(mydata,function(x){
        ni <- ncol(x)
        Coverage(data = x,datatype = datatype,m = 2*ni,nt = ni)
      })
    }
    level <- min(level)
  }
  
  out <- invChatPD(datalist = mydata, datatype = datatype,phylotr = mytree, q = q,
                   reft = reftime, cal = type,level = level, nboot, conf)
  return(out)
}


#' Summarizes phylogenetic data information.
#'
#' Function \code{PDInfo} summarizes phylogenetic data statistics for specified/default reference time.
#' @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).
#' See the function \code{\link{iNEXTPD}} for details.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
#' @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
#' then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#' the pooled assemblage. Default is \code{NULL}.
#' @return Returns a table of phylogenetic data information, including reference sample size (\code{n}), number of sampling units (\code{nT}),
#' number of observed species (\code{S.obs}), observed total branch length, i.e., Faith’s PD (\code{PD.obs}), the first two rare species
#' frequency counts and their branch length sums, at specified/default reference time specified in the argument \code{reftime}.
#' See Chao et al. (2010, 2015) and Hsieh and Chao (2017) for formulas and interpretations.
#' @examples
#' # Datatype: abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PDInfo(data = data, datatype = "abundance", tree = tree)
#' out
#'
#' # Datatype: incidence_raw data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out <- PDInfo(data = data, nT = nT, datatype = "incidence_raw", tree = tree)
#' out
#' 
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
#' @export
PDInfo <- function(data,nT,datatype = "abundance", tree,reftime=NULL){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of key(nT) in sampling units", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(nT)){
      mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
      ntmp <- ntmp+nT[i]
    }
    if(is.null(names(nT))) {
      names(mydata) <- paste0("assemblage",1:length(nT))
    }else{
      names(mydata) = names(nT)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  
  if(datatype=='abundance'){
    infos <- lapply(mydata, function(x){
      datainf(data = x, datatype, phylotr = mytree,reft = reftime) %>% mutate(Reference.time = reftime)
    }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(names(mydata),each = length(reftime))) %>%
      select(Assemblage,n,S.obs,PD.obs,`f1*`,`f2*`,g1,g2,Reference.time)
  }else if (datatype=='incidence_raw'){
    infos <- lapply(mydata, function(x){
      datainf(data = x, datatype, phylotr = mytree,reft = reftime) %>% mutate(Reference.time = reftime)
    }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(names(mydata),each = length(reftime))) %>%
      select(Assemblage,`nT`,S.obs,PD.obs,`Q1*`,`Q2*`,R1,R2,Reference.time)
  }
  
  return(infos)
  
}


#' Plots the outcome of \code{iNEXTPD} using the \code{ggplot2} package.
#'
#' Function \code{ggiNEXTPD} plots the outcome of \code{iNEXTPD} using the \code{ggplot2} package.
#' @param outcome the outcome of the function \code{iNEXTPD}
#' @param plot.type a positive integer or sequence specifying types of curves. There are three
#' types of plots: sample-size-based rarefaction and extrapolation curve (\code{$RE.plot.size}, \code{plot.type = 1});
#' sample coverage curve as a function of sample size (\code{$RE.plot.sizeC}, \code{plot.type = 2}); coverage-based
#' rarefaction and extrapolation curve (\code{$RE.plot.C}, \code{plot.type = 3}). Default is \code{c(1,2,3)}.
#' @return Returns different plots of the estimated diversity curves based on the \code{ggplot2} package. Three choices of \code{plot.type}:
#' \itemize{
#'  \item{\code{$RE.plot.size}: size-based rarefaction and extrapolation curve (\code{plot.type = 1}),
#'  sampling curve depicting phylogenetic diversity estimates as a function of sample size.} \cr
#'  \item{\code{$RE.plot.sizeC}: sample completeness curve (\code{plot.type = 2}), the curve of sample coverage as a function of sample size.} \cr
#'  \item{\code{$RE.plot.C}: coverage-based rarefaction and extrapolation curve (\code{plot.type = 3}), sampling curve depicting phylogenetic diversity estimates as a function of
#'  sample coverage.}
#'  }
#' @examples
#' \donttest{
#' # Datatype: abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- iNEXTPD(data = data, tree = tree, datatype = "abundance", q = c(0, 1, 2), nboot = 30)
#' ggiNEXTPD(out)
#'
#' # Datatype: incidence_raw data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out <- iNEXTPD(data = data, nT = nT, datatype = "incidence_raw", tree = tree, q = c(0, 1, 2), nboot = 30)
#' ggiNEXTPD(out)
#' }
#' 
#' @export
ggiNEXTPD <- function(outcome,plot.type = 1:3){
  TYPE <- c(1,2,3)
  if(is.na(sum(pmatch(plot.type, TYPE))) == F){
    temp2 <- lapply(plot.type, function(j) RE_plot(outcome, j))
    allname <- c('RE.plot.size', 'RE.plot.sizeC','RE.plot.C')
    names(temp2) <- allname[plot.type]
    temp2
  }
}

#' Plots time-profile and q-profile based on the outcome of \code{PhdObs} or \code{PhdAsy} using the \code{ggplot2} package.
#'
#' Function \code{ggtqplotPD} plots time-profile (depicting phylogenetic diversity as a function of time) and q-profile
#' (depicting phylogenetic diversity as a function of order q) based on the outcome of the functions
#' \code{PhdObs} or \code{PhdAsy} using the \code{ggplot2} package.
#' @param outcome the outcome of the functions \code{PhdObs} or \code{PhdAsy}.
#' @param profile specifying the type of profile: \code{profile = "q"} for order q profile and \code{profile = "time"} for time profile.
#' Default is \code{"q"}.
#' @return plot of the PD or meanPD empirical or estimated asymptotic curves using the \code{ggplot2} package.
#' @examples
#' # Observed q-profile plot for abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdObs(data = data, datatype = "abundance", tree = tree, q = seq(0, 2, by = 0.25))
#' ggtqplotPD(out, profile = "q")
#'
#' # Observed time-profile plot for abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdObs(data = data, datatype = "abundance", tree = tree, q = c(0, 1, 2),
#' reftime = seq(0.1, 325, length.out = 40))
#' ggtqplotPD(out, profile = "time")
#'
#' # Asymptotic q-profile plot for incidence_raw data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out <- PhdAsy(data = data, datatype = "incidence_raw", nT = nT, tree = tree,
#' q = seq(0, 2, by = 0.25))
#' ggtqplotPD(out, profile = "q")
#'
#' # Asymptotic time-profile plot for incidence_raw data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out <- PhdAsy(data = data, datatype = "incidence_raw", nT = nT, tree = tree,
#' q = c(0, 1, 2), reftime = seq(0.1, 82.8575, length.out = 40))
#' ggtqplotPD(out, profile = "time")
#' 
#' @export
ggtqplotPD <- function(outcome,profile = 'q'){
  if(profile=='q'){
    Plotq(out = outcome)
  }else if (profile=='time'){
    Plott(out = outcome)
  }else{
    stop("Invalid profile type, please use either q or time.", call. = FALSE)
  }
}


#' @useDynLib iNEXTPD2
#' @importFrom Rcpp sourceCpp
NULL
