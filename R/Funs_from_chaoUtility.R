check_datatype <- function (datatype)
{
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if ((is.na(pmatch(datatype, TYPE))) | (pmatch(datatype,
                                                TYPE) == -1))
    stop("invalid datatype")
  datatype <- match.arg(datatype, TYPE)
  if (datatype == "incidence_freq")
    datatype <- "incidence"
  return(datatype)
}
phyL_Abu_T_<-function(treeNdata,t_1,rootExtend=T,treeH=0){
  rootlength<-0
  if (t_1>0){
    if(rootExtend==T & t_1>treeH) {
      rootlength<-t_1-treeH
      phyAL<-treeNdata %>% mutate(branch.length.new=case_when(tgroup=="Root" ~ rootlength,
                                                              TRUE~ branch.length)
      )
    }
    else{
      phyAL<-treeNdata %>% mutate(branch.length.new=case_when(tgroup=="Root" ~ rootlength,
                                                              branch.height>=t_1 ~ pmax(0,branch.length-branch.height+t_1),
                                                              TRUE~ branch.length))
    }
  }
  else stop("reference T should >0 ")
  
  branch.length.byT<-phyAL %>% pull(branch.length.new)
  names(branch.length.byT)<-phyAL %>% pull(label)
  return(branch.length.byT)
}
#' @import dplyr
#' @import tidytree
#' @import ape
phyBranchAL_Abu <- function(phylo,data, datatype="abundance",refT=0,rootExtend=T,remove0=T){
  #if(class(phylo) != "phylo")
  if (!inherits(phylo, "phylo"))
    stop("invlid class: only support phylo object")
  
  datatype <- check_datatype(datatype)
  if(datatype=="incidence_freq" | datatype=="incidence")
    stop('only support datatype="incidence_raw"')
  
  labels<-names(data)
  my_match <- match(labels, phylo$tip.label)
  if(sum(is.na(my_match)) > 0) stop("Argument labels and tree Tip not matach")
  
  
  if(datatype=="abundance"){
    
    ###drop Abu=0 tips###
    if(remove0==T){
      dtip = phylo$tip.label[-match(names(data[data>0]), phylo$tip.label)]
      subtree = ape::drop.tip(phylo, dtip)
      subdata = data[data>0]
      
    }
    else{
      subtree<-phylo
      subdata<-data
    }
    
    phylo.root<-length(subtree$tip.label)+1
    
    edgelength<-ape::node.depth.edgelength(subtree)
    #treeH<-max(ape::node.depth.edgelength(subtree))
    treeH<-max(edgelength)
    
    ##change to tibble format: easy to read
    phylo.t <- tidytree::as_tibble(subtree)
    
    
    
    # phylo.t.1<-phylo.t %>% mutate(tgroup=case_when(node<phylo.root ~"Tip",
    #                                                node==phylo.root ~"Root",
    #                                                TRUE ~"Inode"),
    #                               newlabel=case_when(node-phylo.root==0 & label =="" ~ "Root",
    #                                                  label=="" ~paste("I",node-phylo.root,sep=""),
    #                                                  TRUE ~ label)) %>% select(-label) %>% rename(label=newlabel)
    
    phylo.t.1<-phylo.t %>% mutate(branch.length=replace(branch.length, is.na(branch.length),0),
                                  tgroup=case_when(node<phylo.root ~"Tip",
                                                   node==phylo.root ~"Root",
                                                   TRUE ~"Inode"),
                                  newlabel=case_when(node-phylo.root==0 & (label ==""|is.na(label)) ~ "Root",
                                                     label==""|is.na(label) ~paste("I",node-phylo.root,sep=""),
                                                     TRUE ~ label ),
                                  edgelengthv=edgelength,
                                  node.age=case_when(edgelengthv==treeH~0,
                                                     TRUE ~ treeH-edgelengthv),
                                  branch.height=case_when(tgroup=="Tip" ~ branch.length,
                                                          tgroup=="Root"~treeH,
                                                          TRUE ~branch.length+node.age )
    ) %>% select(-label) %>% rename(label=newlabel)
    
    tmp<-tibble(label=names(subdata),x=subdata)
    treeNdata<-full_join(phylo.t.1, tmp, by="label")
    inodelist<-treeNdata %>% filter(tgroup !="Tip") %>% pull(node)
    names(inodelist)<-treeNdata %>% filter(tgroup !="Tip") %>% pull(label)
    inode_x<-sapply(inodelist,function(x){offspring(treeNdata,x,tiponly=T) %>% select(x) %>% sum()})
    
    tmp_all<-bind_rows(tibble(label=names(subdata),branch.abun=subdata),tibble(label=names(inode_x),branch.abun=inode_x))
    treeNdata<-full_join(treeNdata, tmp_all, by="label") %>% select(-x,-edgelengthv,-node.age)
    
    phyL<-sapply(refT,function(y) phyL_Abu_T_(treeNdata,y,rootExtend,treeH))
    colnames(phyL)<-paste("T",refT,sep="")
    treeNdata<-treeNdata %>% select(-branch.height)
    
    
    z <- list("treeNabu"=treeNdata,"treeH"=treeH,"BLbyT"=phyL)
    class(z) <- "Chaophyabu"
    return(z)
  }
  
}
#' @import dplyr
#' @import tidytree
#' @import ape
phyBranchAL_Inc<-function(phylo,data, datatype="incidence_raw",refT=0,rootExtend=T,remove0=T){
  if (!inherits(phylo, "phylo"))
    stop("invlid class: only support phylo object")
  
  datatype <- check_datatype(datatype)
  if(datatype=="incidence_freq" | datatype=="incidence")
    stop('only support datatype="incidence_raw"')
  
  labels<-rownames(data)
  my_match <- match(labels, phylo$tip.label)
  if(sum(is.na(my_match)) > 0) stop("Argument labels and tree Tip not matach")
  
  ###drop Abu=0 tips ###
  if(remove0==T){
    if(datatype=="abundance"){
      sp = unique(names(data)[data>0])
      subdata = data[sp]
      
    }
    else if(datatype=="incidence_raw"){
      sp = unique(rownames(data)[rowSums(data)>0])
      subdata = data[sp, ]
    }
    dtip <- phylo$tip.label[-match(sp,phylo$tip.label)]
    subtree = ape::drop.tip(phylo, dtip)
  }
  else{
    subtree<-phylo
    subdata<-data
  }
  
  labels<-rownames(subdata)
  
  
  if(datatype=="incidence_raw"){
    phylo.root<-length(subtree$tip.label)+1
    
    edgelength<-ape::node.depth.edgelength(subtree)
    #treeH<-max(ape::node.depth.edgelength(subtree))
    treeH<-max(edgelength)
    
    ##change to tibble format: easy to read
    phylo.t <- tidytree::as_tibble(subtree)
    phylo.t.1<-phylo.t %>% mutate(branch.length=replace(branch.length, is.na(branch.length),0),
                                  tgroup=case_when(node<phylo.root ~"Tip",
                                                   node==phylo.root ~"Root",
                                                   TRUE ~"Inode"),
                                  newlabel=case_when(node-phylo.root==0 & (label ==""|is.na(label)) ~ "Root",
                                                     label==""|is.na(label) ~paste("I",node-phylo.root,sep=""),
                                                     TRUE ~ label ),
                                  edgelengthv=edgelength,
                                  node.age=case_when(edgelengthv==treeH~0,
                                                     TRUE ~ treeH-edgelengthv),
                                  branch.height=case_when(tgroup=="Tip" ~ branch.length,
                                                          tgroup=="Root"~treeH,
                                                          TRUE ~branch.length+node.age )
    ) %>% select(-label) %>% rename(label=newlabel)
    
    
    y <- iNEXT::as.incfreq(subdata)
    t <- y[1]
    y <- y[-1]
    names(y) <- labels
    tmp.tip<-data.frame(label=names(y),x=y,stringsAsFactors=F)
    treeNdata<-full_join(phylo.t.1, tmp.tip, by="label")
    
    inode_each<-apply(subdata,2,function(i){
      tmp<-data.frame(label=labels,x=i)
      tmp$label<-as.character(tmp$label)
      tmp.treeNdata<-full_join(phylo.t.1, tmp, by="label")
      inodelist<-tmp.treeNdata %>% filter(tgroup !="Tip") %>% pull(node)
      names(inodelist)<-tmp.treeNdata %>% filter(tgroup !="Tip") %>% pull(label)
      ivalue_each<-sapply(inodelist,function(x){offspring(tmp.treeNdata,x,tiponly=T) %>% select(x) %>% max()})
    })
    inode_x <- rowSums(inode_each)
    tmp.inode<-data.frame(label=names(inode_x),branch.abun=inode_x)
    
    tmp.tip<-tmp.tip %>% rename(branch.abun=x)
    tmp_all<-rbind(tmp.tip,tmp.inode)
    treeNdata<-full_join(treeNdata, tmp_all, by="label") %>% select(-x,-edgelengthv,-node.age)
    
    phyL<-sapply(refT,function(y) phyL_Abu_T_(treeNdata,y,rootExtend,treeH))
    colnames(phyL)<-paste("T",refT,sep="")
    treeNdata<-treeNdata %>% select(-branch.height)
    
    
    z <- list("treeNabu"=treeNdata,"treeH"=treeH,"BLbyT"=phyL)
    class(z) <- "Chaophyabu"
    return(z)
  }
  
}
#' @import dplyr
#' @import tidytree
phylo2phytree<-function(phylo){
  if(class(phylo) != "phylo")
    stop("invlid class: only support phylo object")
  
  phylo.root<-length(phylo$tip.label)+1
  
  ##change to tibble format: easy to read
  phylo.t <- tidytree::as_tibble(phylo)
  
  ###leaves
  phylo.t.leaves<-phylo.t %>% filter(node<phylo.root)
  data.leaves<-phylo.t.leaves$branch.length
  names(data.leaves) <- phylo.t.leaves$label
  
  ###nodes
  phylo.t.nodes<-phylo.t %>% filter(node>=phylo.root)
  phylo.t.nodes<-phylo.t.nodes %>% mutate(Inode=node-phylo.root)
  phylo.t.nodes<-phylo.t.nodes %>% mutate(newlable=ifelse(Inode==0,"Root",paste("I",Inode,sep="")))
  phylo.t.nodes<-phylo.t.nodes %>% mutate(label.new=ifelse(is.na(label)|label=="",newlable,label))
  phylo.t.nodes<-phylo.t.nodes %>% mutate(length.new=ifelse(is.na(branch.length),0,branch.length))
  data.nodes<-phylo.t.nodes$length.new
  names(data.nodes) <- phylo.t.nodes$label.new
  
  ###combine leave node to complete data
  phylo.t.nodes<-phylo.t.nodes %>% select(parent,node,length.new,label.new) %>% rename(branch.length=length.new,label=label.new)
  phylo.t.all<-rbind(phylo.t.leaves,phylo.t.nodes)
  phylo.t.all<-phylo.t.all %>% mutate(tgroup=ifelse(node<phylo.root,"Tip",ifelse(node==phylo.root,"Root","Inode")))
  
  
  ##add treeH
  treeH<-max(ape::node.depth.edgelength(phylo))
  
  
  ##add node.age
  edgelength<-ape::node.depth.edgelength(phylo)
  # node.age<-treeH-ape::node.depth.edgelength(phylo)
  node.age<-sapply(1:length(edgelength),function(x){ifelse(isTRUE(all.equal(edgelength[x],treeH)),0,treeH-edgelength[x])})
  phylo.t.all<-phylo.t.all %>% mutate(node.age=node.age)
  
  z <- list("tips"=data.leaves, "nodes"=data.nodes, "phytree"=phylo.t.all,"treeH"=treeH)
  class(z) <- "chaophytree"
  return(z)
  
}
phy_L_Abu_T_<-function(treeNdata,t_1,rootExtend=T,treeH=0){
  rootlength<-0
  if (t_1>0){
    if(rootExtend==T & t_1>treeH) rootlength<-t_1-treeH
  }
  phyAL<-treeNdata %>% mutate(branch.height=ifelse(tgroup=="Tip",branch.length,node.age))
  phyAL<-phyAL %>% mutate(cumsum.length=ifelse(tgroup=="Tip",branch.length,node.age+branch.length))
  phyAL$refT<-t_1
  phyAL<-phyAL %>% mutate(tmp=cumsum.length-refT)
  phyAL<-phyAL %>% mutate(branch.length.new=ifelse(tgroup=="Root",rootlength,
                                                   ifelse(tmp<0,branch.length,
                                                          ifelse(tgroup=="Tip",branch.length-tmp,
                                                                 ifelse(node.age>refT,0,refT-node.age)))))
  branch.length.byT<-phyAL %>% pull(branch.length.new)
  names(branch.length.byT)<-phyAL %>% pull(label)
  tout <- list("branchL"=branch.length.byT)
  return(tout)
}
#' @import dplyr
#' @import tidytree
#' @import ape
phy_BranchAL_IncBootP<-function(phylo,pdata,refT=0,rootExtend=T,remove0=T){
  #if(class(phylo) != "phylo")
  if (!inherits(phylo, "phylo"))
    stop("invlid class: only support phylo object")
  
  labels<-names(pdata)
  my_match <- match(labels, phylo$tip.label)
  if(sum(is.na(my_match)) > 0) stop("Argument labels and tree Tip not matach")
  
  
  
  ###drop Abu=0 tips###
  if(remove0==T){
    dtip = phylo$tip.label[-match(names(pdata[pdata>0]), phylo$tip.label)]
    subtree = ape::drop.tip(phylo, dtip)
    subdata = pdata[pdata>0]
    
  }
  else{
    subtree<-phylo
    subdata<-pdata
  }
  
  
  ###calculate inode abundance
  chaotr<-phylo2phytree(subtree)
  
  #inodelist<-names(chaotr$nodes)
  tmp<-data.frame(label=names(subdata),x=subdata,stringsAsFactors=F)
  treeNdata<-full_join(chaotr$phytree, tmp, by="label")
  inodelist<-treeNdata %>% filter(tgroup !="Tip") %>% pull(node)
  names(inodelist)<-treeNdata %>% filter(tgroup !="Tip") %>% pull(label)
  inode_x<-sapply(inodelist,function(x){
    tmp<-offspring(treeNdata,x,tiponly=T) %>% mutate(y=1-x) %>% select (y) %>% prod()
    pi<-1-tmp
    return(pi)
  })
  
  
  tmp1<-data.frame(label=names(inode_x),branch.incBP=inode_x)
  tmp2<-tmp %>% rename(branch.incBP=x)
  tmp_all<-rbind(tmp2,tmp1)
  treeNdata<-full_join(treeNdata, tmp_all, by="label") %>% select(-x)
  
  ###calculate inode length
  # treeH and rootlength
  treeH<-chaotr$treeH
  
  phyL<-sapply(refT,function(y) phy_L_Abu_T_(treeNdata,y,rootExtend,treeH))
  names(phyL)<-paste("T",refT,sep="")
  
  z <- list("treeNincBP"=treeNdata,"treeH"=treeH,"BLbyT"=phyL)
  class(z) <- "ChaophyincBP"
  return(z)
  
  
}
