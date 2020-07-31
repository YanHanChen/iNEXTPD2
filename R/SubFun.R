#===============PhDObs==================
datainf <- function(data, datatype, phylotr,reft){
  if(datatype == "abundance"){
    new <- phyBranchAL_Abu(phylotr,data,datatype,reft)
    #new$treeNabu$branch.length <- new$BLbyT[,1]
    data <- data[data>0]
    ai <- new$treeNabu$branch.abun
    Lis <- new$BLbyT
    I1 <- which(ai==1);I2 <- which(ai==2)
    f1 <- length(I1);f2 <- length(I2)
    a1 <- sapply(1:ncol(Lis),function(i){
      Li = Lis[,i]
      PD_obs <- sum(Li)
      g1 <- sum(Li[I1])
      g2 <- sum(Li[I2])
      c(PD_obs,g1,g2)
    }) %>% matrix(nrow = 3) %>% t()
    a1 <- tibble(n = sum(data),S.obs = length(data),PD.obs = a1[,1],
                 'f1*' = f1,'f2*' = f2, g1 = a1[,2],g2 = a1[,3])
  }else if(datatype == 'incidence_raw'){
    new <- phyBranchAL_Inc(phylotr,data,datatype,reft)
    #new$treeNabu$branch.length <- new$BLbyT[,1]
    data <- data[rowSums(data)>0,colSums(data)>0,drop=F]
    ai <- new$treeNabu$branch.abun
    Lis <- new$BLbyT
    I1 <- which(ai==1);I2 <- which(ai==2)
    f1 <- length(I1);f2 <- length(I2)
    a1 <- sapply(1:ncol(Lis),function(i){
      Li = Lis[,i]
      PD_obs <- sum(Li)
      g1 <- sum(Li[I1])
      g2 <- sum(Li[I2])
      c('PD.obs' = PD_obs, 'g1' = g1, "g2" = g2)
    }) %>% matrix(nrow = 3) %>% t()
    a1 <- tibble('T' = ncol(data),S.obs = nrow(data),PD.obs = a1[,1],
                 'Q1*' = f1,'Q2*' = f2, R1 = a1[,2],R2 = a1[,3])
    names(a1) <- c("T", "S.obs", "PD.obs", "Q1*", "Q2*", "R1", "R2")
  }
  return(a1)
  
}
PD.qprofile=function(aL, q, cal, nt) {
  #aL is a table of 3 columns: abundance, branch lengths and characters specifying the group of each node.
  Abun <- unlist(aL[,1])
  bL <- unlist(aL[,2])
  t_bar <- sum(Abun / nt * bL)
  
  Sub = function(q) {
    PD=ifelse(q==1, exp(-sum(bL*Abun/t_bar/nt*log(Abun/t_bar/nt)) ),
              sum(bL*(Abun/t_bar/nt)^q)^(1/(1-q)))
    ifelse(cal=="PD",PD,PD/t_bar)
  }
  sapply(q, Sub)
}
Phdqtable <- function(datalist, phylotr, q, cal, datatype, nboot, conf, reft){
  qtile <- qnorm(1-(1-conf)/2)
  # all assemblages.
  nms <- names(datalist)
  if(datatype=="abundance"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]]
      x <- x[x>0]
      n <- sum(x)
      emp <- c(t(PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, q=q,cal = cal,nt = n)))
      if(nboot>0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft, BLs = aL$BLbyT )
        Li_b <- sapply(1:length(reft), function(reft_){
          Li_b <- Boots$Li[,reft_]
          Li_b[Li_b>reft[reft_]] <- reft[reft_]
          Li_b
        })
        Li_b <- matrix(Li_b,ncol=ncol(aL$BLbyT))
        ses <- sapply(1:nboot,function(B){
          c(t(PD.Tprofile(ai = Boots$boot_data[,B],Lis = Li_b, q=q,cal = cal,nt = n)))
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if(datatype=="incidence_raw"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]]
      x <- x[rowSums(x)>0,colSums(x)>0]
      n <- ncol(x)
      # For incidence type, we use occurence frequencies instead of raw data since we already have aL table.
      emp <- c(t(PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, q=q,cal = cal,nt = n)))
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft, BLs = aL$BLbyT,splunits = n)
        Li_b <- sapply(1:length(reft), function(reft_){
          Li_b <- Boots$Li[,reft_]
          Li_b[Li_b>reft[reft_]] <- reft[reft_]
          Li_b
        })
        Li_b <- matrix(Li_b,ncol=ncol(aL$BLbyT))
        ses <- sapply(1:nboot,function(B){
          c(t(PD.Tprofile(ai = Boots$boot_data[,B],Lis = Li_b, q=q,cal = cal,nt = n)))
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }
  
  odr <- rep(q,length(reft)*length(nms))
  reftime <- rep(rep(reft,each = length(q)),length(nms))
  nms_tmp <- rep(nms,each = length(q)*length(reft))
  Outputforq <- tibble(Order.q = odr, Empirical = out[,1],LCL = out[,2], UCL = out[,3],
                       reftime = reftime,Assemblage = nms_tmp)
  Outputforq <- Outputforq %>% mutate (method = ifelse(cal=="PD", "PD", "meanPD"))
  return(Outputforq)
}
color_nogreen <- function(n) {
  all <- c("red", "blue", "orange", "purple", "pink", "cyan", "brown", "yellow")
  all[1:n]
}
PD.Tprofile=function(ai,Lis, q, cal, nt) {
  isn0 <- ai>0
  ai <- ai[isn0]
  Lis <- Lis[isn0,,drop=F]
  #q can be a vector
  t_bars <- as.numeric(ai %*% Lis/nt)
  
  
  bL <- sapply(1:length(t_bars), function(i) Lis[,i]/t_bars[i])
  pAbun <- ai/nt
  
  out <- sapply(q, function(j){
    if(j==1) as.numeric(-(pAbun*log(pAbun)) %*% bL) %>% exp()
    else as.numeric((pAbun)^j %*% bL) %>% .^(1/(1-j))
  }) %>% matrix(.,ncol = length(q))
  if(cal=="PD"){
    out <- sapply(1:length(t_bars), function(i){
      out[i,]*t_bars[i]
    }) %>% matrix(.,ncol = length(t_bars)) %>% t()
  }
  out
}
Phdttable <- function(datalist, phylotr, times, cal, datatype, nboot, conf){
  # Note 200117: currently, the reference time is automatically fixed at tree height of pooled species.
  qtile <- qnorm(1-(1-conf)/2)
  # all assemblages.
  if (datatype=="incidence_raw") {
    H_max <- get.rooted.tree.height(phylotr)
    da <- lapply(datalist, rowSums) %>% do.call(cbind, .) %>% rowSums()
  }
  if (datatype=="abundance") {
    H_max <- get.rooted.tree.height(phylotr)
    da <- do.call(cbind, datalist) %>% rowSums()
  }
  # if(abs(H_max - reft)<=1e-4) H_max <- reft
  # Note 200204: currently, to compute Q, we treat incidence data as abundance and absolutely pool
  aL <- phyBranchAL_Abu(phylo = phylotr,data = da,"abundance",refT = H_max)$treeNabu %>%
    select(branch.abun,branch.length,tgroup)
  PD2 <- PD.qprofile(aL,q = 2, cal =  "PD",nt = sum(da))
  Q <- H_max-(H_max^2)/PD2
  nms <- names(datalist)
  
  q_int <- c(0, 1, 2)
  if(datatype=="abundance"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = times)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, q=q_int,cal = cal,nt = n) %>%
        c()
      
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = times, BLs = aL$BLbyT )
        Lis_b <- Boots$Li
        Lis_b <- sapply(1:length(times),function(l){
          tmp <- Lis_b[,l]
          tmp[tmp>times[l]] <- times[l]
          tmp
        })
        f0 <- Boots$f0
        #tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Lis_b)-length(x)-f0))
        ses <- sapply(1:nboot, function(B){
          x_b <- Boots$boot_data[,B]
          isn0 <- as.vector(x_b>0)
          Lis_b_tmp <- Lis_b[isn0,]
          #tgroup_B_tmp <- tgroup_B[isn0]
          x_b <- x_b[isn0]
          out_b <- PD.Tprofile(ai = x_b,Lis=Lis_b_tmp,q=q_int,cal = cal, nt = n) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if(datatype=="incidence_raw"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = times)
      x <- datalist[[i]] %>% .[rowSums(.)>0,]
      n <- ncol(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT,q=q_int,cal = cal,nt = n) %>%
        c()
      
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = times,
                           BLs = aL$BLbyT,splunits = n)
        Lis_b <- Boots$Li
        Lis_b <- sapply(1:length(times),function(l){
          tmp <- Lis_b[,l]
          tmp[tmp>times[l]] <- times[l]
          tmp
        })
        f0 <- Boots$f0
        ses <- sapply(1:nboot, function(B){
          x_b <- Boots$boot_data[,B]
          isn0 <- as.vector(x_b>0)
          Lis_b_tmp <- Lis_b[isn0,]
          x_b <- x_b[isn0]
          out_b <- PD.Tprofile(ai = x_b,Lis=Lis_b_tmp,q=q_int,cal = cal,nt = n) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }
  
  Outputfort <- tibble(time = rep(times,length(q_int)*length(datalist)),
                       Empirical = out[,1],LCL = out[,2], UCL = out[,3],
                       Order.q = rep(rep(q_int, each=length(times)),length(datalist)),
                       Assemblage = rep(nms, each=length(times)*length(q_int)),
                       method=ifelse(cal=="PD", "PD", "meanPD"))
  out = list(fort = Outputfort, Q_Height=c(Q, H_max))
  return(out)
}
EmpPD <- function(datalist,datatype, phylotr, q, reft, cal, nboot, conf){
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  if(datatype=="abundance"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, q=q,cal = cal,nt = n) %>%
        c()
      
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft, BLs = aL$BLbyT )
        Li_b <- Boots$Li
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        f0 <- Boots$f0
        ses <- sapply(1:nboot, function(B){
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,cal = cal, nt = n) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if(datatype=="incidence_raw"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[rowSums(.)>0,]
      n <- ncol(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT,q=q,cal = cal,nt = n) %>%
        c()
      
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft,
                           BLs = aL$BLbyT,splunits = n)
        Li_b <- Boots$Li
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        f0 <- Boots$f0
        ses <- sapply(1:nboot, function(B){
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,cal = cal,nt = n) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }
  Output <- tibble(Assemblage = rep(nms, each=length(reft)*length(q)),
                   Order.q = rep(rep(q, each=length(reft)),length(datalist)),
                   qPD = out[,1],qPD.LCL = out[,2], qPD.UCL = out[,3],
                   Reference.time = rep(reft,length(q)*length(datalist)),
                   Method='Empirical',
                   Type=ifelse(cal=="PD", "PD", "meanPD")) %>%
    arrange(Reference.time)
  return(Output)
}
AUC_one_table <- function(datalist, phylotr, knot, cal, datatype, nboot, conf, reft_max) {
  qtile <- qnorm(1-(1-conf)/2)
  times_AUC <- seq(0.01, reft_max, length.out = knot)
  nms <- names(datalist)
  q_int <- c(0, 1, 2)
  if(datatype=="abundance"){
    AUC <- lapply(1:length(datalist),function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = times_AUC)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT,
                         q=q_int,cal = cal, nt = n)
      # print(paste0("emp:",dim(emp)," AUC diff:",length(c(diff(times_AUC),0))))
      LA <-  c(diff(times_AUC),0) %*% emp %>% as.numeric()
      RA <-   c(0,diff(times_AUC)) %*% emp %>% as.numeric()
      auc <- colMeans(rbind(LA,RA))
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = times_AUC, BLs = aL$BLbyT )
        Lis_b <- Boots$Li
        Lis_b <- sapply(1:length(times_AUC),function(l){
          tmp <- Lis_b[,l]
          tmp[tmp>times_AUC[l]] <- times_AUC[l]
          tmp
        })
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Lis_b)-length(x)-f0))
        
        ses <- sapply(1:nboot, function(B){
          x_b <- Boots$boot_data[,B]
          isn0 <- as.vector(x_b>0)
          Lis_b_tmp <- Lis_b[isn0,]
          tgroup_B_tmp <- tgroup_B[isn0]
          x_b <- x_b[isn0]
          out_b <- PD.Tprofile(ai = x_b,Lis=Lis_b_tmp,q=q_int,cal = cal,nt = n)
          LA_b <-  c(diff(times_AUC),0) %*% out_b %>% as.numeric()
          RA_b <-   c(0,diff(times_AUC)) %*% out_b %>% as.numeric()
          auc_b <- colMeans(rbind(LA_b,RA_b))
          auc_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(auc))
      }
      output <- cbind(auc,auc-qtile*ses,auc+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if (datatype=="incidence_raw"){
    AUC <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = times_AUC)
      x <- datalist[[i]] %>% .[rowSums(.)>0,]
      n <- ncol(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT,
                         q=q_int,cal = cal,nt = n)
      LA <-  c(diff(times_AUC),0) %*% emp %>% as.numeric()
      RA <-   c(0,diff(times_AUC)) %*% emp %>% as.numeric()
      auc <- colMeans(rbind(LA,RA))
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = times_AUC,
                           BLs = aL$BLbyT,splunits = n)
        Lis_b <- Boots$Li
        Lis_b <- sapply(1:length(times_AUC),function(l){
          tmp <- Lis_b[,l]
          tmp[tmp>times_AUC[l]] <- times_AUC[l]
          tmp
        })
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Lis_b)-nrow(x)-f0))
        ses <- sapply(1:nboot, function(B){
          x_b <- Boots$boot_data[,B]
          isn0 <- as.vector(x_b>0)
          Lis_b_tmp <- Lis_b[isn0,]
          tgroup_B_tmp <- tgroup_B[isn0]
          x_b <- x_b[isn0]
          out_b <- PD.Tprofile(ai = x_b,Lis=Lis_b_tmp,q=q_int,
                               cal = cal,nt = n)
          LA_b <-  c(diff(times_AUC),0) %*% out_b %>% as.numeric()
          RA_b <-   c(0,diff(times_AUC)) %*% out_b %>% as.numeric()
          auc_b <- colMeans(rbind(LA_b,RA_b))
          auc_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(auc))
      }
      output <- cbind(auc,auc-qtile*ses,auc+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }
  
  
  AUC <- tibble(Order.q = rep(q_int,length(nms)), Empirical = AUC[,1],LCL = AUC[,2], UCL = AUC[,3],
                Assemblage = rep(nms,each = length(q_int)))
  AUC
}
#=====old version=====
# Plott <- function(out, cal, Q_Height){
#   fort <- out
#   fort$Order.q = paste0("q = ", fort$Order.q)
#   fort$Order.q <- factor(fort$Order.q)
#   Assemblage <- unique(fort$Assemblage)
#   Q = Q_Height[1]
#   root = Q_Height[2]
#   # print(fort)
#   # print(Q)
#   if (cal=="PD") {
#     if(length(Assemblage)==1){
#       p2 <- ggplot(fort, aes(x=time, y=Empirical)) + geom_line(size=1.5,aes(color=Order.q))+
#         geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Order.q),linetype = 0,alpha=0.3)
#       p2 <-  p2 +xlab("time")+ylab("Phylogenetic Diversity")+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
#         geom_point(size=3, data=subset(fort, time%in%c(Q, root)), aes(x=time, y=Empirical, color=Order.q))+
#         annotate('text',x=Q, y=0.1,label="Q" ,parse = TRUE,size=5, color = "gray") +
#         annotate('text',x=root, y=0.1, label="root",parse = TRUE,size=5, color = "gray") +
#         geom_vline(xintercept = c(Q,root), linetype = "longdash",size=0.5, color = "gray")
#     }else{
#       p2 <- ggplot(fort, aes(x=time, y=Empirical, color=Assemblage, linetype=Assemblage)) + geom_line(size=1.5)  +
#         geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Assemblage),linetype = 0,alpha=0.3)+
#         scale_color_manual(values = color_nogreen(length(unique(fort$Assemblage))))+
#         scale_fill_manual(values = color_nogreen(length(unique(fort$Assemblage))))+
#         theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
#         geom_point(size=3, data=subset(fort, time%in%c(Q, root)), aes(x=time, y=Empirical, color=Assemblage))+
#         annotate('text',x=Q, y=0.1,label="Q" ,parse = TRUE,size=5, color = "gray") +
#         annotate('text',x=root, y=0.1, label="root",parse = TRUE,size=5, color = "gray") +
#         geom_vline(xintercept = c(Q,root), linetype = "longdash",size=0.5, color = "gray") +
#         facet_wrap(~Order.q, scales = "free")
#       p2 <-  p2 +xlab("time")+ylab("Phylogenetic Diversity")
#     }
#   } else {
#     if(length(Assemblage)==1){
#       p2 <- ggplot(fort, aes(x=time, y=Empirical)) + geom_line(size=1.5,aes(color=Order.q))+
#         geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Order.q),linetype = 0,alpha=0.3)
#       p2 <-  p2 +xlab("time")+ylab("Phylogenetic Hill numbers")+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
#         geom_point(size=3, data=subset(fort, time%in%c(0.01, Q, root)), aes(x=time, y=Empirical, color=Order.q))+
#         annotate('text',x=Q, y=0.1,label="Q" ,parse = TRUE,size=5, color = "gray") +
#         annotate('text',x=0.01, y=0.1,label=0.01 ,parse = TRUE,size=5, color = "gray") +
#         annotate('text',x=root, y=0.1, label="root",parse = TRUE,size=5, color = "gray") +
#         geom_vline(xintercept = c(0.01, Q,root), linetype = "longdash",size=0.5, color = "gray")
#     }else{
#       p2 <- ggplot(fort, aes(x=time, y=Empirical, color=Assemblage, linetype=Assemblage)) + geom_line(size=1.5)  +
#         geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Assemblage),linetype = 0,alpha=0.3)+
#         scale_color_manual(values = color_nogreen(length(unique(fort$Assemblage))))+
#         scale_fill_manual(values = color_nogreen(length(unique(fort$Assemblage))))+
#         theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
#         geom_point(size=3, data=subset(fort, time%in%c(0.01, Q, root)), aes(x=time, y=Empirical, color=Assemblage))+
#         annotate('text',x=Q, y=0.1,label="Q" ,parse = TRUE,size=5, color = "gray") +
#         annotate('text',x=0.01, y=0.1,label=0.01 ,parse = TRUE,size=5, color = "gray") +
#         annotate('text',x=root, y=0.1, label="root",parse = TRUE,size=5, color = "gray") +
#         geom_vline(xintercept = c(0.01, Q,root), linetype = "longdash",size=0.5, color = "gray") +
#         facet_wrap(~Order.q, scales = "free")
#       p2 <-  p2 +xlab("time")+ylab("Phylogenetic Hill numbers")
#     }
#   }
#   return(p2)
# }
#=====new version=====
Plott <- function(out){
  fort <- out
  fort$Order.q <- factor(paste0("q = ", fort$Order.q),levels = paste0("q = ", unique(fort$Order.q)))
  Assemblage <- unique(fort$Assemblage)
  ylab_ <- paste0(unique(fort$Method)," ",unique(fort$Type))
  if(length(Assemblage)==1){
    p2 <- ggplot(fort, aes(x=Reference.time, y=qPD)) + theme_bw() +geom_line(size=1.5,aes(color=Order.q))+
      geom_ribbon(aes(ymin=qPD.LCL,ymax=qPD.UCL,fill=Order.q),linetype = 0,alpha=0.2)+
      xlab("Reference time")+ylab(ylab_)+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))
  }else{
    p2 <- ggplot(fort, aes(x=Reference.time, y=qPD, color=Assemblage, linetype=Assemblage)) + theme_bw() + geom_line(size=1.5)  +
      geom_ribbon(aes(ymin=qPD.LCL,ymax=qPD.UCL,fill=Assemblage,alpha=0.2),linetype = 0,alpha=0.2)+
      scale_color_manual(values = color_nogreen(length(unique(fort$Assemblage))))+
      scale_fill_manual(values = color_nogreen(length(unique(fort$Assemblage))))+
      theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
      facet_wrap(~Order.q, scales = "free")+xlab("Reference time")+ylab(ylab_)
  }
  return(p2)
}
Plotq <- function(out){
  forq <- out
  forq$Reference.time <- factor(paste0('Ref.time = ',as.character(round(forq$Reference.time,4))),
                                levels = unique(paste0('Ref.time = ',as.character(round(forq$Reference.time,4)))))
  Assemblage <- unique(forq$Assemblage)
  ylab_ <- paste0(unique(forq$Method)," ",unique(forq$Type))
  q1 <- unique(forq$Order.q[(forq$Order.q %% 1)==0])
  if(length(Assemblage)==1){
    p1 <- ggplot(forq, aes(x=Order.q, y=qPD, color=Reference.time)) + theme_bw() + geom_line(size=1.5)+
      geom_ribbon(aes(ymin=qPD.LCL,ymax=qPD.UCL,fill=Reference.time),linetype = 0,alpha=0.2)
    #lai 1006
    p1 <-  p1 +xlab("Order q")+ylab(ylab_) + theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
      geom_point(size=3, data=subset(forq, Order.q%in%q1), aes(x=Order.q, y=qPD, color=Reference.time))
  }else{
    p1 <- ggplot(forq, aes(x=Order.q, y=qPD, color=Assemblage, linetype=Assemblage)) + theme_bw() + geom_line(size=1.5)  +
      geom_ribbon(aes(ymin=qPD.LCL,ymax=qPD.UCL,fill=Assemblage),linetype = 0,alpha=0.2)+
      scale_color_manual(values = color_nogreen(length(unique(forq$Assemblage))))+
      scale_fill_manual(values = color_nogreen(length(unique(forq$Assemblage))))+
      theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
      geom_point(size=3, data=subset(forq, Order.q%in%q1), aes(x=Order.q, y=qPD, color=Assemblage))+
      facet_wrap(~Reference.time, scales = "free")
    p1 <-  p1 +xlab("Order q")+ylab(ylab_)
  }
  return(p1)
}
#===============PhDAsy==================
AsyPD <- function(datalist, datatype, phylotr, q,reft, cal,nboot, conf){#change final list name
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  tau_l <- length(reft)
  if(datatype=="abundance"){
    Estoutput <- lapply(datalist,function(x){
      #atime <- Sys.time()
      x <- x[x>0]
      n <- sum(x)
      aL <- phyBranchAL_Abu(phylo = phylotr,data = x,datatype,refT = reft)
      #aL$treeNabu$branch.length <- aL$BLbyT[,1]
      #aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      est <- PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,q = q,nt = n,cal = cal)
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft, BLs = aL$BLbyT )
        Li_b <- Boots$Li
        f0 <- Boots$f0
        # tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Li_b)-length(x)-f0))
        # aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          outb <- PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q = q,nt = n,cal = cal)
          return(outb)
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(est))
      }
      est <- tibble(Order.q = rep(q,tau_l), qPD = est,
                    qPD.LCL = est - qtile*ses, qPD.UCL = est + qtile*ses,
                    Reference.time = rep(reft,each = length(q)))
      est
    })
  }else if(datatype=="incidence_raw"){
    Estoutput <- lapply(datalist,function(x){
      #atime <- Sys.time()
      x <- x[rowSums(x)>0,colSums(x)>0]
      n <- ncol(x)
      aL <- phyBranchAL_Inc(phylo = phylotr,data = x,datatype,refT = reft)
      # aL$treeNabu$branch.length <- aL$BLbyT[,1]
      # aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      est <- PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,q = q,nt = n,cal = cal)
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype = datatype,nboot = nboot,
                           splunits = n,reft = reft, BLs = aL$BLbyT )
        Li_b <- Boots$Li
        f0 <- Boots$f0
        # tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Li_b)-nrow(x)-f0))
        # aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          outb <- PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q = q,nt = n,cal = cal)
          return(outb)
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(est))
      }
      est <- tibble(Order.q = rep(q,tau_l), qPD = est,
                    qPD.LCL = est - qtile*ses, qPD.UCL = est + qtile*ses,
                    Reference.time = rep(reft,each = length(q)))
      est
    })
  }
  Estoutput <- do.call(rbind,Estoutput) %>%
    mutate(Assemblage = rep(names(datalist),each = length(q)*tau_l),Method = 'Asymptotic',
           Type=ifelse(cal=="PD", "PD", "meanPD")) %>%
    select(Assemblage,Order.q,qPD,qPD.LCL, qPD.UCL, Reference.time,Method, Type) %>%
    arrange(Reference.time)
  Estoutput$qPD.LCL[Estoutput$qPD.LCL<0] = 0
  return(Estoutput)
}
Asy_plot = function(output, type, method=NULL){##add title
  if(is.null(method) == T){
    ylab_ <- "Phylogenetic diversity"
  }else{
    if( substr(method,1,1) != "1" ){
      ylab_ <- paste("Phylogenetic", method ,"diversity")
    }else{
      ylab_ <- method
    }
  }
  if(is.null(method) == T & type == 1) {
    title_ <- "Phylogenetic diversity (asymptotic)"
  } else if(is.null(method) == T & type == 2) {
    title_ <- "Phylogenetic diversity (observed)"
  } else if( substr(method,1,1) != "1"  & type == 1) {
    title_ <- paste("Phylogenetic", method ,"diversity", "(asymptotic)")
  } else if( substr(method,1,1) != "1"  & type == 2) {
    title_ <- paste("Phylogenetic", method ,"diversity", "(observed)")
  } else if( substr(method,1,1) == "1"  & type == 1) {
    title_ <- paste(method, "(asymptotic)")
  } else if( substr(method,1,1) == "1"  & type == 2) {
    title_ <- paste(method, "(observed)")
  } else if(type == 3) {
    title_ <- paste(method, "(asymptotic)")
  } else {
    title_ <- paste(method, "(observed)")
  }
  if(ncol(output) %in% c(2, 5)) output = cbind(output, "Beta")
  if(ncol(output) == 6){
    colnames(output) = c("x", "y", "se", "LCL", "UCL","Assemblage")
    Assemblage <- unique(output[,6]) %>% unlist
  }else{
    colnames(output) = c("x", "y","Assemblage")
    Assemblage <- unique(output[,3]) %>% unlist
  }
  q <- unlist(output$x)
  q1<-q[(round(q) - q) ==0]
  if(length(Assemblage) == 1){
    p <- ggplot(output,aes(x=x,y=y))+geom_line(size=1.5,color="#F8766D")+xlab("Order q")+
      geom_point(size=3, data=subset(output, x%in%q1),color="#F8766D")+
      ylab(ylab_)+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))
    if(ncol(output) == 6) p <- p + geom_ribbon(aes(ymin=LCL,ymax=UCL),alpha=0.2,fill="#F8766D")
  }else{
    p <- ggplot(output,aes(x=x,y=y,color=Assemblage,linetype=Assemblage))+geom_line(size=1.5)+xlab("Order q")+
      scale_color_manual(values = color_nogreen(length(unique(output$Assemblage))))+
      geom_point(size=3, data=subset(output, x%in%q1))+
      ylab(ylab_)+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))
    if(ncol(output) == 6) p <- p + geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Assemblage),alpha=0.2, colour=NA)+scale_fill_manual(values = color_nogreen(length(unique(output$Assemblage))))
  }
  p <- p+ggtitle(title_)
  return(p)
}
#=====old version=====
# PhD.q.est = function(aL, q, datatype, nt){
#   t_bar <-  sum(aL[,1] / nt * aL[,2])
#   aL <- aL %>% select(branch.abun, branch.length)
#   PD_obs <- sum(aL$branch.length)
#   fg1 <- aL %>% filter(branch.abun==1)
#   fg2 <- aL %>% filter(branch.abun==2)
#   f1 <- nrow(fg1);g1 <- sum(fg1$branch.length)
#   f2 <- nrow(fg2);g2 <- sum(fg2$branch.length)
#   if(f2 > 0){
#     A = 2*f2/((nt-1)*f1+2*f2)
#   }else if(f2 == 0 & f1 > 0){
#     A = 2/((nt-1)*(f1-1)+2)
#   }else{
#     A = 1
#   }
#   tmpaL <- aL %>% group_by(branch.abun, branch.length) %>% summarise(n_node = n()) %>% as.matrix()
#
#   if(sum(abs(q-round(q))!=0)>0 | max(q)>2) {
#     deltas <- sapply(0:(nt-1), function(k){
#       del_tmp <- tmpaL[tmpaL[,1]<=(nt-k),,drop=FALSE]
#       delta(del_tmpaL = del_tmp,k,nt)
#     })
#   }
#   Sub <- function(q){
#     if(q==0){
#       ans <- PD_obs+Dq0(nt,f1,f2,g1,g2)
#     }else if(q==1){
#       h2 <- Dq1_2(nt,g1,A)
#       h1 <- aL %>% filter(branch.abun<=(nt-1)) %>%
#         mutate(diga = digamma(nt)-digamma(branch.abun)) %>% apply(., 1, prod) %>% sum(.)/nt
#       #print(paste0("A:",A," g1:",g1," h1:", h1, " h2:",h2))
#       h <- h1+h2
#       ans <- t_bar*exp(h/t_bar)
#     }else if(q==2){
#       ans <- Dq2(as.matrix(tmpaL),nt,t_bar)
#     }else{
#       # timea <- Sys.time()
#       k <- 0:(nt-1)
#       a <- (choose(q-1,k)*(-1)^k*deltas) %>% sum
#       b <- ifelse(g1==0|A==1,0,(g1*((1-A)^(1-nt))/nt)*(A^(q-1)-sum(choose(q-1,k)*(A-1)^k)))
#       ans <- ((a+b)/(t_bar^q))^(1/(1-q))
#       # timeb <- Sys.time()
#       # print(timeb-timea)
#     }
#     return(ans)
#   }
#   est <- sapply(q, Sub)
#   # if(datatype=="abundance")info <- c("n" = nt, "S.obs" = length(Abun), "PD.obs" = PD_obs, "f1*" = f1,
#   #                                    "f2*" = f2, "g1" = g1, "g2" = g2)
#   # else if (datatype == "incidence_raw") info <- c("T" = nt, "S.obs" = length(inci_freq), "PD.obs" = PD_obs, "Q1*" = f1,
#   #                                                 "Q2*" = f2, "R1" = g1, "R2" = g2)
#
#   est
# }
#=====new version=====
PhD.q.est = function(ai,Lis, q, nt, cal){
  t_bars <- as.numeric(t(ai) %*% Lis/nt)
  I1 <- which(ai==1);I2 <- which(ai==2)
  f1 <- length(I1);f2 <- length(I2)
  if(f2 > 0){
    A = 2*f2/((nt-1)*f1+2*f2)
  }else if(f2 == 0 & f1 > 0){
    A = 2/((nt-1)*(f1-1)+2)
  }else{
    A = 1
  }
  S <- length(ai)
  if(1 %in% q){
    ai_h1_I <- ai<=(nt-1)
    h1_pt2 <- rep(0,S)
    ai_h1 <- ai[ai_h1_I]
    h1_pt2[ai_h1_I] <- tibble(ai = ai) %>% .[ai_h1_I,] %>% mutate(diga = digamma(nt)-digamma(ai)) %>%
      apply(., 1, prod)/nt
  }
  if(2 %in% q){
    q2_pt2 <- unlist(ai*(ai-1)/nt/(nt-1))
  }
  if(sum(abs(q-round(q))!=0)>0 | max(q)>2) {
    deltas_pt2 <- sapply(0:(nt-1), function(k){
      ai_delt_I <- ai<=(nt-k)
      deltas_pt2 <- rep(0,S)
      deltas_pt2[ai_delt_I] <- delta_part2(ai = ai[ai_delt_I],k = k,n = nt)
      deltas_pt2
    }) %>% t() # n x S matrix of delta (2nd part)
  }
  Sub <- function(q,g1,g2,PD_obs,t_bar,Li){
    if(q==0){
      ans <- PD_obs+Dq0(nt,f1,f2,g1,g2)
    }else if(q==1){
      h2 <- Dq1_2(nt,g1,A)
      h1 <- sum(Li*h1_pt2)
      h <- h1+h2
      ans <- t_bar*exp(h/t_bar)
    }else if(q==2){
      #ans <- Dq2(as.matrix(tmpaL),nt,t_bar)
      ans <- t_bar^2/sum(Li*q2_pt2)
    }else{
      # timea <- Sys.time()
      k <- 0:(nt-1)
      deltas <- as.numeric(deltas_pt2 %*% Li)
      a <- (choose(q-1,k)*(-1)^k*deltas) %>% sum
      b <- ifelse(g1==0|A==1,0,(g1*((1-A)^(1-nt))/nt)*(A^(q-1)-sum(choose(q-1,k)*(A-1)^k)))
      ans <- ((a+b)/(t_bar^q))^(1/(1-q))
      # timeb <- Sys.time()
      # print(timeb-timea)
    }
    return(ans)
  }
  est <- sapply(1:ncol(Lis),function(i){
    Li = Lis[,i]
    t_bar <- t_bars[i]
    PD_obs <- sum(Li)
    g1 <- sum(Li[I1])
    g2 <- sum(Li[I2])
    est <- sapply(q, function(q_) Sub(q = q_,g1 = g1,g2 = g2,PD_obs = PD_obs,t_bar = t_bar,Li = Li))
  })
  if(cal=='PD'){
    est <- as.numeric(est)
  }else if (cal=='meanPD'){
    est <- as.numeric(sapply(1:length(t_bars), function(i){
      est[,i]/t_bars[i]
    }))
  }
  return(est)
}
Boots.one = function(phylo, aL, datatype, nboot,reft, BLs, splunits = NULL){
  if(datatype=='abundance'){
    data <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    names(data) <- rownames(BLs)[1:length(data)]
    n <- sum(data)
    f1 <- sum(data==1)
    f2 <- sum(data==2)
    f0 <- ceiling( ifelse( f2>0 , ((n-1) / n) * (((f1)^2) / (2*f2) ), ((n-1) / n) * (f1)*(f1-1) / 2 ) )
    c <- ifelse(f2>0, 1 - (f1/n)*((n-1)*f1/((n-1)*f1+2*f2)),
                1 - (f1/n)*((n-1)*(f1-1)/((n-1)*(f1-1)+2)))
    lambda <- (1-c) / sum((data/n)*(1- (data/n) )^n)
    
    p_hat <- (data/n) * (1-lambda*(1- (data/n) )^n)
    p_hat0 <- rep( (1-c) / f0 , f0 );names(p_hat0) <- paste0("notob",1:length(p_hat0))
    g0_hat <- sapply(1:length(reft), function(i){
      g1 <- BLs[,i][aL$branch.abun==1] %>% sum
      g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
      g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
      g0_hat
    })
    # g1 <- aL$branch.length[aL$branch.abun==1] %>% sum
    # g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
    # g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
    ###Notice that the species order of pL_b doesn't change even that of data changes. (property of phyBranchAL_Abu)
    pL_b <- phyBranchAL_Abu(phylo, p_hat, datatype,refT = reft[1])
    pL_b$treeNabu$branch.length <- pL_b$BLbyT[,1]
    pL_b <- pL_b$treeNabu
    pL_b[length(p_hat)+1,"branch.abun"] <- 1
    Li <- BLs
    Li <- rbind(Li[1:length(data),,drop=F],matrix(g0_hat/f0,nrow = f0,ncol = ncol(Li),byrow = T), Li[-(1:length(data)),,drop=F])
    p_hat <- c(p_hat,p_hat0,unlist(pL_b[-(1:length(data)),"branch.abun"]))
    boot_data <- sapply(p_hat,function(p) rbinom(n = nboot,size = n,prob = p)) %>% t()
  }else if(datatype=='incidence_raw'){
    n <- splunits
    data <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    names(data) <- rownames(BLs)[1:length(data)]
    u <- sum(data)
    f1 <- sum(aL$branch.abun==1)
    f2 <- sum(aL$branch.abun==2)
    f0 <- ceiling( ifelse( f2>0 , ((n-1) / n) * (((f1)^2) / (2*f2) ), ((n-1) / n) * (f1)*(f1-1) / 2 ) )
    c <- ifelse(f2>0, 1 - (f1/u)*((n-1)*f1/((n-1)*f1+2*f2)),
                1 - (f1/u)*((n-1)*(f1-1)/((n-1)*(f1-1)+2)))
    lambda <- u/n*(1-c) / sum((data/n)*(1- (data/n) )^n)
    p_hat <- (data/n) * (1-lambda*(1- (data/n) )^n)
    p_hat0 <- rep( (u/n) * (1-c) / f0 , f0 );names(p_hat0) <- paste0("notob",1:length(p_hat0))
    g0_hat <- sapply(1:length(reft), function(i){
      g1 <- BLs[,i][aL$branch.abun==1] %>% sum
      g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
      g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
      g0_hat
    })
    # g1 <- aL$branch.length[aL$branch.abun==1] %>% sum
    # g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
    # g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
    pL_b <- phy_BranchAL_IncBootP(phylo = phylo, pdata = p_hat, refT = reft[1])
    pL_b <- pL_b$treeNincBP
    pL_b[length(p_hat)+1,"branch.incBP"] <- 1
    #pL_b$treeNincBP$branch.length <- pL_b$BLbyT[,1] # delete for now since length is a list instead of a matrix.
    #data_iB <- unlist(aL$branch.abun)
    #pL_b <- (data_iB/n) * (1-lambda*(1- (data_iB/n) )^n)
    Li <- BLs
    Li <- rbind(Li[1:length(data),,drop=F],matrix(g0_hat/f0,nrow = f0,ncol = ncol(Li),byrow = T),
                Li[-(1:length(data)),,drop=F])
    p_hat <- c(p_hat,p_hat0,unlist(pL_b[-(1:length(data)),"branch.incBP"]))
    boot_data <- sapply(p_hat,function(p) rbinom(n = nboot,size = n,prob = p)) %>% t()
  }
  return(list(boot_data=boot_data,Li = Li, f0 = f0))
}
#===============inextPD==================
inextPD = function(datalist, datatype, phylotr, q,reft, m, cal, nboot, conf=0.95, unconditional_var=TRUE){
  # m is a list
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  if(datatype=="abundance"){
    Estoutput <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      #====conditional on m====
      qPDm <- PhD.m.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,m = m[[i]],
                        q = q,nt = n,cal = cal) %>% as.numeric()
      covm = Coverage(x, datatype, m[[i]],n)
      #====unconditional====
      if(unconditional_var){
        goalSC <- unique(covm)
        qPD_unc <- unique(invChatPD_abu(x = x,ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                                        q = q,Cs = goalSC,n = n,cal = cal))
      }
      if(nboot>1){
        Boots <- Boots.one(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT)
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Li_b)-length(x)-f0))
        #aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        if(unconditional_var){
          ses <- sapply(1:nboot, function(B){
            # atime <- Sys.time()
            ai_B <- Boots$boot_data[,B]
            isn0 <- ai_B>0
            qPDm_b <-  PhD.m.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                 m=m[[i]],q=q,nt = n,cal = cal) %>% as.numeric()
            covm_b <- Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]],n)
            qPD_unc_b <- unique(invChatPD_abu(x = ai_B[isn0&tgroup_B=="Tip"],
                                              ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                              q = q,Cs = goalSC,n = n,cal = cal))$qPD
            # btime <- Sys.time()
            # print(paste0("Est boot sample",B,": ",btime-atime))
            return(c(qPDm_b,covm_b,qPD_unc_b))
          }) %>% apply(., 1, sd)
        }else{
          ses <- sapply(1:nboot, function(B){
            # atime <- Sys.time()
            ai_B <- Boots$boot_data[,B]
            isn0 <- ai_B>0
            qPDm_b <-  PhD.m.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                 m=m[[i]],q=q,nt = n,cal = cal) %>% as.numeric()
            covm_b <- Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]],n)
            # btime <- Sys.time()
            # print(paste0("Est boot sample",B,": ",btime-atime))
            return(c(qPDm_b,covm_b))
          }) %>% apply(., 1, sd)
        }
      }else{
        if(unconditional_var){
          ses <- rep(NA,length(c(qPDm,covm,qPD_unc$qPD)))
        }else{
          ses <- rep(NA,length(c(qPDm,covm)))
        }
      }
      
      ses_pd <- ses[1:length(qPDm)]
      ses_cov <- ses[(length(qPDm)+1):(length(qPDm)+length(covm))]
      m_ <- rep(m[[i]],each = length(q)*length(reft))
      method <- ifelse(m[[i]]>n,'Extrapolation',ifelse(m[[i]]<n,'Rarefaction','Observed'))
      method <- rep(method,each = length(q)*length(reft))
      orderq <- rep(q,length(reft)*length(m[[i]]))
      SC_ <- rep(covm,each = length(q)*length(reft))
      SC.LCL_ <- rep(covm-qtile*ses_cov,each = length(q)*length(reft))
      SC.UCL_ <- rep(covm+qtile*ses_cov,each = length(q)*length(reft))
      reft_ <- rep(rep(reft,each = length(q)),length(m[[i]]))
      out_m <- tibble(Assemblage = nms[i], m=m_,Method=method,Order.q=orderq,
                      qPD=qPDm,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
                      SC=SC_,SC.LCL=SC.LCL_,SC.UCL=SC.UCL_,
                      Reference.time = reft_,
                      Type=ifelse(cal=="PD", "PD", "meanPD")) %>%
        arrange(Reference.time,Order.q,m)
      out_m$qPD.LCL[out_m$qPD.LCL<0] <- 0;out_m$SC.LCL[out_m$SC.LCL<0] <- 0
      out_m$SC.UCL[out_m$SC.UCL>1] <- 1
      if(unconditional_var){
        ses_pd_unc <- ses[-(1:(length(qPDm)+length(covm)))]
        out_C <- qPD_unc %>% mutate(qPD.LCL = qPD-qtile*ses_pd_unc,qPD.UCL = qPD+qtile*ses_pd_unc,
                                    Type=ifelse(cal=="PD", "PD", "meanPD"),
                                    Assemblage = nms[i])
        id_C <- match(c('Assemblage','goalSC','SC','m', 'Method', 'Order.q', 'qPD', 'qPD.LCL','qPD.UCL','Reference.time',
                        'Type'), names(out_C), nomatch = 0)
        out_C <- out_C[, id_C] %>% arrange(Reference.time,Order.q,m)
        out_C$qPD.LCL[out_C$qPD.LCL<0] <- 0
      }else{
        out_C <- NULL
      }
      return(list(size_based = out_m, coverage_based = out_C))
    })
  }else if(datatype=="incidence_raw"){
    Estoutput <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[rowSums(.)>0,colSums(.)>0]
      n <- ncol(x)
      #====conditional on m====
      qPDm <- PhD.m.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,m = m[[i]],
                        q = q,nt = n,cal = cal) %>% as.numeric()
      covm = Coverage(x, datatype, m[[i]], n)
      #====unconditional====
      if(unconditional_var){
        goalSC <- unique(covm)
        qPD_unc <- unique(invChatPD_inc(x = rowSums(x),ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                                        q = q,Cs = goalSC,n = n,cal = cal))
        #colnames(qPD_unc)[which(colnames(qPD_unc)=='t')] <- 'm'
      }
      if(nboot>1){
        Boots <- Boots.one(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT,n)
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Li_b)-nrow(x)-f0))
        if(unconditional_var){
          ses <- sapply(1:nboot, function(B){
            # atime <- Sys.time()
            ai_B <- Boots$boot_data[,B]
            isn0 <- ai_B>0
            qPDm_b <-  PhD.m.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                 m=m[[i]],q=q,nt = n,cal = cal) %>% as.numeric()
            covm_b = Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]],n)
            qPD_unc_b <- unique(invChatPD_inc(x = ai_B[isn0&tgroup_B=="Tip"],
                                              ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                              q = q,Cs = goalSC,n = n,cal = cal))$qPD
            # btime <- Sys.time()
            # print(paste0("Est boot sample",B,": ",btime-atime))
            return(c(qPDm_b,covm_b,qPD_unc_b))
          }) %>% apply(., 1, sd)
        }else{
          ses <- sapply(1:nboot, function(B){
            # atime <- Sys.time()
            ai_B <- Boots$boot_data[,B]
            isn0 <- ai_B>0
            qPDm_b <-  PhD.m.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                 m=m[[i]],q=q,nt = n,cal = cal) %>% as.numeric()
            covm_b = Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]],n)
            # btime <- Sys.time()
            # print(paste0("Est boot sample",B,": ",btime-atime))
            return(c(qPDm_b,covm_b))
          }) %>% apply(., 1, sd)
        }
      }else{
        if(unconditional_var){
          ses <- rep(NA,length(c(qPDm,covm,qPD_unc$qPD)))
        }else{
          ses <- rep(NA,length(c(qPDm,covm)))
        }
      }
      ses_pd <- ses[1:length(qPDm)]
      ses_cov <- ses[(length(qPDm)+1):(length(qPDm)+length(covm))]
      m_ <- rep(m[[i]],each = length(q)*length(reft))
      method <- ifelse(m[[i]]>n,'Extrapolation',ifelse(m[[i]]<n,'Rarefaction','Observed'))
      method <- rep(method,each = length(q)*length(reft))
      orderq <- rep(q,length(reft)*length(m[[i]]))
      SC_ <- rep(covm,each = length(q)*length(reft))
      SC.LCL_ <- rep(covm-qtile*ses_cov,each = length(q)*length(reft))
      SC.UCL_ <- rep(covm+qtile*ses_cov,each = length(q)*length(reft))
      reft_ = rep(rep(reft,each = length(q)),length(m[[i]]))
      out_m <- tibble(Assemblage = nms[i], nt=m_,Method=method,Order.q=orderq,
                      qPD=qPDm,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
                      SC=SC_,SC.LCL=SC.LCL_,SC.UCL=SC.UCL_,
                      Reference.time = reft_,
                      Type=ifelse(cal=="PD", "PD", "meanPD")) %>%
        arrange(Reference.time,Order.q,nt)
      out_m$qPD.LCL[out_m$qPD.LCL<0] <- 0;out_m$SC.LCL[out_m$SC.LCL<0] <- 0
      out_m$SC.UCL[out_m$SC.UCL>1] <- 1
      if(unconditional_var){
        ses_pd_unc <- ses[-(1:(length(qPDm)+length(covm)))]
        out_C <- qPD_unc %>% mutate(qPD.LCL = qPD-qtile*ses_pd_unc,qPD.UCL = qPD+qtile*ses_pd_unc,
                                    Type=ifelse(cal=="PD", "PD", "meanPD"),
                                    Assemblage = nms[i])
        id_C <- match(c('Assemblage','goalSC','SC','nt', 'Method', 'Order.q', 'qPD', 'qPD.LCL','qPD.UCL','Reference.time',
                        'Type'), names(out_C), nomatch = 0)
        out_C <- out_C[, id_C] %>% arrange(Reference.time,Order.q,nt)
        out_C$qPD.LCL[out_C$qPD.LCL<0] <- 0
      }else{
        out_C <- NULL
      }
      return(list(size_based = out_m, coverage_based = out_C))
    })
  }
  if(unconditional_var){
    ans <- list(size_based = do.call(rbind,lapply(Estoutput, function(x) x$size_based)),
                coverage_based = do.call(rbind,lapply(Estoutput, function(x) x$coverage_based)))
  }else{
    ans <- list(size_based = do.call(rbind,lapply(Estoutput, function(x) x$size_based)))
  }
  return(ans)
}
#=====old version=====
# PhD.m.est = function(aL, m, q,datatype, nt){
#   t_bar <- sum(aL[,1] / nt * aL[,2])
#
#   aL_matrix = as.matrix(aL[,c(1,2)])
#   RPD_m = RPD(aL_matrix, nt, nt-1, q)
#   obs <- RPD(aL_matrix, nt, nt, q)
#   # obs <- PD.qprofile(aL = aL, q = Q, cal="PD" ,datatype = datatype , nforboot = nforboot, splunits = splunits)
#   #asymptotic value
#   asy <- PhD.q.est(aL = aL,q = q, datatype = datatype, nt = nt)
#   asy <- sapply(1:length(q), function(j){
#     max(asy[j],obs[j])
#   })
#   #beta
#   beta <- rep(0,length(q))
#
#   beta0plus <- which(asy != obs)
#   beta[beta0plus] <- (obs[beta0plus]-RPD_m[beta0plus])/(asy[beta0plus]-RPD_m[beta0plus])
#   # if(asy == obs) beta = 0
#   # if(asy != obs) beta =(obs-RPD_m)/(asy-RPD_m)
#   #Extrapolation
#   EPD = function(m,q){
#     m = m-nt
#     out <- sapply(1:length(q), function(i){
#       if( q[i]!=2 ) {
#         obs[i]+(asy[i]-obs[i])*(1-(1-beta[i])^m)
#       }else if( q[i] == 2 ){
#         1/sum( (aL_matrix[,2]/(t_bar)^2)*((1/(nt+m))*(aL_matrix[,1]/nt)+((nt+m-1)/(nt+m))*(aL_matrix[,1]*(aL_matrix[,1]-1)/(nt*(nt-1)))) )
#       }
#     })
#     # if( Q == 0 | Q == 1 ) EPD = obs+(asy-obs)*(1-(1-beta)^m)
#     # if( Q == 2 ) EPD = 1/sum( (aL_matrix[,2]/(t_bar)^2)*((1/(n+m))*(aL_matrix[,1]/n)+((n+m-1)/(n+m))*(aL_matrix[,1]*(aL_matrix[,1]-1)/(n*(n-1)))) )
#     return(out)
#   }
#   Sub = function(m){
#     if(m<nt){
#       RPD(aL_matrix,nt,m,q)
#     }else if(m==nt){
#       obs
#     }else{
#       EPD(m,q)
#     }
#   }
#   sapply(m, Sub)
# }
#=====new version=====
PhD.m.est = function(ai,Lis, m, q, nt, cal){
  t_bars <- as.numeric(t(ai) %*% Lis/nt)
  if(sum(m>nt)>0){
    #Extrapolation
    RPD_m <- RPD(ai,Lis,nt,nt-1,q)
    obs <- RPD(ai, Lis, nt,nt, q)
    EPD = function(m,obs,asy){
      m = m-nt
      out <- sapply(1:ncol(Lis), function(i){
        asy_i <- asy[,i];obs_i <- obs[,i];RPD_m_i <- RPD_m[,i]
        Li <- Lis[,i];t_bar <- t_bars[i]
        asy_i <- sapply(1:length(q), function(j){
          max(asy_i[j],obs_i[j])
        })
        beta <- rep(0,length(q))
        beta0plus <- which(asy_i != obs_i)
        beta[beta0plus] <-(obs_i[beta0plus]-RPD_m_i[beta0plus])/(asy_i[beta0plus]-RPD_m_i[beta0plus])
        outq <- sapply(1:length(q), function(i){
          if( q[i]!=2 ) {
            obs_i[i]+(asy_i[i]-obs_i[i])*(1-(1-beta[i])^m)
          }else if( q[i] == 2 ){
            1/sum( (Li/(t_bar)^2)*((1/(nt+m))*(ai/nt)+((nt+m-1)/(nt+m))*(ai*(ai-1)/(nt*(nt-1)))) )
          }
        })
        outq
      })
      return(out)
    }
    # obs <- PD.qprofile(aL = aL, q = Q, cal="PD" ,datatype = datatype , nforboot = nforboot, splunits = splunits)
    #asymptotic value
    asy <- matrix(PhD.q.est(ai = ai,Lis = Lis,q = q, nt = nt,cal = 'PD'),nrow = length(q),ncol = length(t_bars))
  }else if (sum(m==nt)>0){
    obs <- RPD(ai, Lis, nt,nt, q)
  }
  if(cal=='PD'){
    out <- sapply(m, function(mm){
      if(mm<nt){
        ans <- RPD(ai = ai,Lis = Lis,n = nt,m = mm,q = q)
      }else if(mm==nt){
        ans <- obs
      }else{
        ans <- EPD(m = mm,obs = obs,asy = asy)
      }
      return(as.numeric(ans))
    })
  }else if (cal=='meanPD'){
    out <- sapply(m, function(mm){
      if(mm<nt){
        ans <- RPD(ai = ai,Lis = Lis,n = nt,m = mm,q = q)
      }else if(mm==nt){
        ans <- obs
      }else{
        ans <- EPD(m = mm,obs = obs,asy = asy)
      }
      ans <- sapply(1:length(t_bars), function(i){
        ans[,i]/t_bars[i]
      })
      as.numeric(ans)
    })
  }
  out <- matrix(out,ncol = length(m))
  return(out)
}
Coverage = function(data, datatype, m, nt){
  if(datatype == "incidence_raw") datatype = "incidence"
  ifelse("matrix" %in% class(data) || "data.frame" %in% class(data), type <- "raw", type <- "numeric")
  ifelse(type == "raw", x <- rowSums(data), x <- data )
  if(type=="raw" || datatype=='incidence') u<-sum(x)
  x <- x[x>0]
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (nt - 1) / nt * f1 * (f1 - 1) / 2, (nt - 1) / nt * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, nt*f0.hat/(nt*f0.hat+f1), 1)
  Sub <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < nt) {
      xx <- x[(nt-x)>=m]
      out <- 1-sum(xx / nt * exp(lgamma(nt-xx+1)-lgamma(nt-xx-m+1)-lgamma(nt)+lgamma(nt-m)))
    }
    if(m == nt) out <- 1-f1/nt*A
    if(m > nt) out <- 1-f1/nt*A^(m-nt+1)
    out
  }
  Sub2 <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < nt) {
      xx <- x[(nt-x)>=m]
      out <- 1-sum(xx / u * exp(lgamma(nt-xx+1)-lgamma(nt-xx-m+1)-lgamma(nt)+lgamma(nt-m)))
    }
    if(m == nt) out <- 1-f1/u*A
    if(m > nt) out <- 1-f1/u*A^(m-nt+1)
    out
  }
  sapply(m, FUN = function(i){
    ifelse(datatype!='abundance', Sub2(i), Sub(i) )
  })
}
RE_plot = function(data, type){
  #data <- as.data.frame(data)
  datatype <- ifelse(colnames(data$size_based[,2])=='m','abundance','incidence_raw')
  x <- ifelse(datatype=='incidence_raw', 'sampling units', "individuals")
  x_name <- colnames(data$size_based[,2])
  if(type == 1){
    data <- data$size_based
    id <- match(c(x_name,'Method','qPD','qPD.LCL','qPD.UCL','Assemblage','Order.q',
                  'Reference.time'), names(data), nomatch = 0)
    output <- data[, id]
    xlab_ <- paste0("Number of ", x)
  }else if(type == 2){
    data <- data$size_based %>% filter(Order.q==1)
    id <- match(c(x_name,'Method','SC','SC.LCL','SC.UCL','Assemblage','Order.q',
                  'Reference.time'), names(data), nomatch = 0)
    output <- data[, id]
    xlab_ <- paste0("Number of ", x)
    ylab_ <- "Sample Coverage"
  }else if(type == 3){
    data <- data$coverage_based
    id <- match(c('SC','Method','qPD','qPD.LCL','qPD.UCL','Assemblage','Order.q',
                  'Reference.time'), names(data), nomatch = 0)
    output <- data[, id]
    xlab_ <- "Sample Coverage"
  }
  ylab_ <- unique(data$Type)
  title <- c("Sample-size-based sampling curve", "Sample completeness curve","Coverage-based sampling curve")
  title <- title[type]
  
  
  Assemblage <- unique(data$Assemblage)
  colnames(output) <- c("x", "Method", "y", "LCL", "UCL", "Assemblage","Order.q",'Reference.time')
  output$Method <- as.character(output$Method)
  output$Assemblage <- as.character(output$Assemblage)
  output$Reference.time <- round(output$Reference.time,3)
  output$Reference.time <- factor(paste0('Ref.time = ',output$Reference.time),
                                  levels = paste0('Ref.time = ',unique(output$Reference.time)))
  output_obser <- output %>% filter(Method=="Observed")
  output$Method[output$Method=="Observed"] <- "Rarefaction"
  
  # odr_grp <- as_labeller(c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2"))
  # odr_grp <- list(
  #   '0'="q = 0",
  #   '1'="q = 1",
  #   '2'="q = 2"
  # )
  # refts <- unique(output$Reference.time)
  # reft_grp <- sapply(refts, function(i){
  #   paste0('Ref.time = ',refts[i])
  # })
  # names(reft_grp) <- as.character(refts)
  # plot_labeller <- function(variable,value){
  #   if (variable=='Order.q') {
  #     return(odr_grp[value])
  #   } else {
  #     return(reft_grp[value])
  #   }
  # }
  mylab <- labeller(
    Order.q = c(`0` = "q = 0", `1` = "q = 1",`2` = 'q = 2')
  )
  if(length(Assemblage) == 1){
    outp <- ggplot(output, aes(x = x, y = y))+ theme_bw() +
      geom_ribbon(aes(ymin = LCL, ymax = UCL),fill="#F8766D",alpha=0.2)+geom_line(size=1.5, aes(x = x, y = y, linetype=Method),color="#F8766D")+
      geom_point(size=3, data=output_obser,color="#F8766D")+xlab(xlab_)+ylab(ylab_)+
      scale_linetype_manual(values = c("solid", "dashed"), name="Method",breaks=c("Rarefaction", "Extrapolation"), labels=c("Rarefaction", "Extrapolation"))+
      theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
      ggtitle(title)+guides(linetype=guide_legend(keywidth=2.5))
    if(type!=3) outp <- outp + facet_wrap(Reference.time~Order.q,scales = "free_y",labeller=mylab)
    else if (type==3) outp <- outp + facet_wrap(~Reference.time,scales = "free_y")
    #if(length(unique(output$Reference.time))==1) outp <- outp + theme(strip.background = element_blank(), strip.text.x = element_blank())
  }else{
    outp <- ggplot(output, aes(x = x, y = y, color=Assemblage)) + theme_bw() +
      geom_line(size=1.5, aes(x = x, y = y, color=Assemblage, linetype=Method))+
      scale_color_manual(values = color_nogreen(length(unique(output$Assemblage))))+
      geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Assemblage), alpha=0.2, colour=NA)+
      scale_fill_manual(values = color_nogreen(length(unique(output$Assemblage))))+
      geom_point(size=3, data=output_obser)+xlab(xlab_)+ylab(ylab_)+
      scale_linetype_manual(values = c("solid", "dashed"), name="Method",breaks=c("Rarefaction", "Extrapolation"), labels=c("Rarefaction", "Extrapolation"))+
      theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
      ggtitle(title)+guides(linetype=guide_legend(keywidth=2.5))
    if(type!=2)  outp <- outp + facet_wrap(Reference.time~Order.q,scales = "free_y",labeller=mylab)
    else if (type == 2) outp <- outp + facet_wrap(~Reference.time,scales = "free_y")
    #if(length(unique(output$Reference.time))==1) outp <- outp + theme(strip.background = element_blank(), strip.text.x = element_blank())
  }
  return(outp)
}
#===============EstimatePD===============
invChatPD <- function(datalist, datatype,phylotr, q, reft, cal,level, nboot, conf){
  qtile <- qnorm(1-(1-conf)/2)
  # datalist <- lapply(1:ncol(x), function(i) {tmp = x[, i];names(tmp) = rownames(x);tmp})
  # if (is.null(colnames(x))) {
  #   names(datalist) <- paste0("data", 1:ncol(x))
  # } else {names(datalist) <- colnames(x)}
  # x <- datalist
  # refT <- max(ape::node.depth.edgelength(phylotr))
  if(datatype=='abundance'){
    out <- lapply(datalist,function(x_){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = x_,'abundance',refT = reft)
      # aL$treeNabu$branch.length <- aL$BLbyT[,1]
      # aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      x_ <- x_[x_>0]
      n <- sum(x_)
      #n_sp_samp <- sum(aL_table$tgroup=='Tip')
      est <- invChatPD_abu(x = x_,ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                           q = q,Cs = level, n = n,cal = cal)
      if(nboot>1){
        Boots <- Boots.one(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT,n)
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x_)+f0),rep("Inode",nrow(Li_b)-length(x_)-f0))
        #aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          # atime <- Sys.time()
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          # isn0 <- as.vector(aL_table_b[,1]>0)
          # Li_b_tmp <- Li_b[isn0,]
          # aL_table_b <- aL_table_b[isn0,]
          est_b <- invChatPD_abu(x = ai_B[isn0&tgroup_B=="Tip"],ai = ai_B[isn0],
                                 Lis = Li_b[isn0,,drop=F],q = q,Cs = level,
                                 n = n,cal = cal)$qPD
          
          return(est_b)
        }) %>% matrix(.,nrow = length(q)*length(reft)*length(level)) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,nrow(est))
      }
      est <- est %>% mutate(qPD.LCL=qPD-qtile*ses,qPD.UCL=qPD+qtile*ses)
    }) %>% do.call(rbind,.)
  }else if(datatype=='incidence_raw'){
    out <- lapply(datalist,function(x_){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = x_,'incidence_raw',refT = reft)
      # aL$treeNabu$branch.length <- aL$BLbyT[,1]
      # aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      x_ <- x_[rowSums(x_)>0,colSums(x_)>0]
      n <- ncol(x_)
      est <- invChatPD_inc(x = rowSums(x_),ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                           q = q,Cs = level, n = n,cal = cal)
      if(nboot>1){
        Boots <- Boots.one(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT,n)
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x_)+f0),rep("Inode",nrow(Li_b)-nrow(x_)-f0))
        #aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          # atime <- Sys.time()
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          # isn0 <- as.vector(aL_table_b[,1]>0)
          # Li_b_tmp <- Li_b[isn0,]
          # aL_table_b <- aL_table_b[isn0,]
          est_b <- invChatPD_inc(x = ai_B[isn0&tgroup_B=="Tip"],ai = ai_B[isn0],
                                 Lis = Li_b[isn0,,drop=F],q = q,Cs = level,
                                 n = n,cal = cal)$qPD
          return(est_b)
        }) %>% matrix(.,nrow = length(q)*length(reft)*length(level)) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,nrow(est))
      }
      est <- est %>% mutate(qPD.LCL=qPD-qtile*ses,qPD.UCL=qPD+qtile*ses)
    }) %>% do.call(rbind,.)
  }
  Assemblage = rep(names(datalist), each = length(q)*length(reft)*length(level))
  out <- out %>% mutate(Assemblage = Assemblage,
                        Type=ifelse(cal=="PD", "PD", "meanPD"))
  if(datatype=='abundance'){
    out <- out %>% select(Assemblage,m,Method,Order.q,qPD,qPD.LCL,qPD.UCL,
                          SC,goalSC,Reference.time,Type) %>% arrange(Reference.time,goalSC,Order.q)
  }else if(datatype=='incidence_raw'){
    out <- out %>% select(Assemblage,nt,Method,Order.q,qPD,qPD.LCL,qPD.UCL,
                          SC,goalSC,Reference.time,Type) %>% arrange(Reference.time,goalSC,Order.q)
  }
  out$qPD.LCL[out$qPD.LCL<0] <- 0
  rownames(out) <- NULL
  out
}
#=====old version=====
# invChatPD_abu <- function(aL_table, q, Cs, n){
#   x <- unlist(aL_table$branch.abun[aL_table$tgroup=="Tip"])
#   #n <- sum(x)
#   refC = Coverage(x, 'abundance', n, n)
#   f <- function(m, C) abs(Coverage(x, 'abundance', m, n) - C)
#   mm <- sapply(Cs, function(cvrg){
#     if (refC > cvrg) {
#       opt <- optimize(f, C = cvrg, lower = 0, upper = n)
#       mm <- opt$minimum
#       mm <- round(mm)
#     }else if (refC <= cvrg) {
#       f1 <- sum(x == 1)
#       f2 <- sum(x == 2)
#       if (f1 > 0 & f2 > 0) {
#         A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
#       }
#       if (f1 > 1 & f2 == 0) {
#         A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
#       }
#       if (f1 == 1 & f2 == 0) {
#         A <- 1
#       }
#       if (f1 == 0 & f2 == 0) {
#         A <- 1
#       }
#       mm <- (log(n/f1) + log(1 - cvrg))/log(A) - 1
#       mm <- n + mm
#       mm <- round(mm)
#     }
#     mm
#   })
#   mm[mm==0] <- 1
#   SC <- Coverage(x, 'abundance', mm, n)
#   out <- PhD.m.est(aL = aL_table,m = mm,q = q,datatype = 'abundance',nt = n)
#   out <- as.vector(out)
#   method <- ifelse(mm>n,'Extrapolated',ifelse(mm<n,'Interpolated','Observed'))
#   method <- rep(method,each = length(q))
#   m <- rep(mm,each = length(q))
#   order <- rep(q,length(mm))
#   SC <- rep(SC,each = length(q))
#   tibble(m = m,method = method,Order.q = order,
#          qPD = out,SC=SC)
# }
# invChatPD_inc <- function(aL_table, q, Cs, n){
#   x <- unlist(aL_table$branch.abun[aL_table$tgroup=="Tip"])
#   refC = Coverage(x, 'incidence', n, n)
#   f <- function(m, C) abs(Coverage(x, 'incidence', m, n) - C)
#   mm <- sapply(Cs, function(cvrg){
#     if (refC > cvrg) {
#       opt <- optimize(f, C = cvrg, lower = 0, upper = max(x))
#       mm <- opt$minimum
#       mm <- round(mm)
#     }else if (refC <= cvrg) {
#       f1 <- sum(x == 1)
#       f2 <- sum(x == 2)
#       U <- sum(x)
#       if (f1 > 0 & f2 > 0) {
#         A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
#       }
#       if (f1 > 1 & f2 == 0) {
#         A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
#       }
#       if (f1 == 1 & f2 == 0) {
#         A <- 1
#       }
#       if (f1 == 0 & f2 == 0) {
#         A <- 1
#       }
#       mm <- (log(U/f1) + log(1 - cvrg))/log(A) - 1
#       mm <- n + mm
#       mm <- round(mm)
#     }
#     mm
#   })
#   mm[mm==0] <- 1
#   SC <- Coverage(x, 'incidence', mm, n)
#   out <- PhD.m.est(aL = aL_table,m = mm,q = q,datatype = 'incidence_raw',nt = n)
#   out <- as.vector(out)
#   method <- ifelse(mm>n,'Extrapolated',ifelse(mm<n,'Interpolated','Observed'))
#   method <- rep(method,each = length(q))
#   m <- rep(mm,each = length(q))
#   order <- rep(q,length(mm))
#   SC <- rep(SC,each = length(q))
#   tibble(t_ = m,method = method,Order.q = order,
#          qPD = out,SC=SC)
# }
#=====new version=====
invChatPD_abu <- function(x,ai,Lis, q, Cs, n,cal){
  #x <- unlist(aL_table$branch.abun[aL_table$tgroup=="Tip"])
  refC = Coverage(x, 'abundance', n, n)
  f <- function(m, C) abs(Coverage(x, 'abundance', m, n) - C)
  mm <- sapply(Cs, function(cvrg){
    if (refC > cvrg) {
      opt <- optimize(f, C = cvrg, lower = 0, upper = n)
      mm <- opt$minimum
      mm <- round(mm)
    }else if (refC <= cvrg) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }else if(f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }else if(f1 == 1 & f2 == 0) {
        A <- 1
      }else if(f1 == 0 & f2 == 0) {
        A <- 1
      }
      mm <- ifelse(A==1,0,(log(n/f1) + log(1 - cvrg))/log(A) - 1)
      mm <- n + mm
      mm <- round(mm)
    }
    mm
  })
  mm[mm==0] <- 1
  SC <- Coverage(x, 'abundance', mm, n)
  out <- as.numeric(PhD.m.est(ai = ai,Lis = Lis,m = mm,q = q,nt = n,cal = cal))
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,each = length(q)*ncol(Lis))
  m <- rep(mm,each = length(q)*ncol(Lis))
  order <- rep(q,ncol(Lis)*length(mm))
  SC <- rep(SC,each = length(q)*ncol(Lis))
  goalSC <- rep(Cs,each = length(q)*ncol(Lis))
  reft <- as.numeric(substr(colnames(Lis),start = 2,stop = nchar(colnames(Lis))))
  Reference.time = rep(rep(reft,each = length(q)),length(Cs))
  tibble(m = m,Method = method,Order.q = order,
         qPD = out,SC=SC,goalSC = goalSC,Reference.time = Reference.time)
}
invChatPD_inc <- function(x,ai,Lis, q, Cs, n,cal){ # x is a matrix
  #x <- unlist(aL_table$branch.abun[aL_table$tgroup=="Tip"])
  refC = Coverage(x, 'incidence', n, n)
  f <- function(m, C) abs(Coverage(x, 'incidence', m, n) - C)
  mm <- sapply(Cs, function(cvrg){
    if (refC > cvrg) {
      opt <- optimize(f, C = cvrg, lower = 0, upper = n)
      mm <- opt$minimum
      mm <- round(mm)
    }else if (refC <= cvrg) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      U <- sum(x)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }else if(f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }else if(f1 == 1 & f2 == 0) {
        A <- 1
      }else if(f1 == 0 & f2 == 0) {
        A <- 1
      }
      mm <- ifelse(A==1,0,(log(U/f1) + log(1 - cvrg))/log(A) - 1)
      mm <- n + mm
      mm <- round(mm)
    }
    mm
  })
  mm[mm==0] <- 1
  SC <- Coverage(x, 'incidence', mm, n)
  out <-  as.numeric(PhD.m.est(ai = ai,Lis = Lis,m = mm,q = q,nt = n,cal = cal))
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,each = length(q)*ncol(Lis))
  m <- rep(mm,each = length(q)*ncol(Lis))
  order <- rep(q,ncol(Lis)*length(mm))
  SC <- rep(SC,each = length(q)*ncol(Lis))
  goalSC <- rep(Cs,each = length(q)*ncol(Lis))
  reft <- as.numeric(substr(colnames(Lis),start = 2,stop = nchar(colnames(Lis))))
  Reference.time = rep(rep(reft,each = length(q)),length(Cs))
  tibble(nt = m,Method = method,Order.q = order,
         qPD = out,SC=SC,goalSC = goalSC, Reference.time = Reference.time)
}


