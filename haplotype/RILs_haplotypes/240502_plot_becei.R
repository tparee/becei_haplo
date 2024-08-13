source("/Users/tomparee/Documents/Documents - MacBook Pro de tom/basics_TP.R")
load("~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240502_beceiFounderPhasedHaplotypes_chrII.Rdata")
load("~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240505_rilsAfounderHaplotypeBlocks_chrII.Rdata")

foundercolnames = c("FA.g1", "FA.g2", "FM.g1", "FM.g2")

info=phasedfounderhaplotypes[,1:4]
colnames(info) = c("ID", "chrom", "pos", "cM")

Rilshaplotypes = RIlsA_FounderHaplotypeBlocks_chrII




FuseHapSeparatedByOneSNP = function(HAPLO, posinfo = c("whichsnp", "pos", "cM"), foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2")){
  
  HAPLO = HAPLO[order(HAPLO$whichsnp1,HAPLO$whichsnp2),]
  
  fused = do.call(rbind, lapply(1:nrow(HAPLO), function(i){
    #print(i)
    whichfounder = which(HAPLO[i,foundercolnames]>0)
    founder = matrix(unlist(round(HAPLO[,whichfounder+2])), nrow = nrow(HAPLO))
    founder[founder>0.0] = 1
    
    #Look if there is another haplotype block with corresponding to the same founder
    nextcommon =  sort(unlist( apply(founder, 2, function(x){which(x>0)}) ))
    nextcommon = min(nextcommon[nextcommon >i])
    if(nextcommon==Inf){return(NULL)}
    j = nextcommon 
    
    #Look if there is an overlap between haplotype i & j
    overlapwithnext = HAPLO[j,1] <= HAPLO[i,2]
    
    #gap size between i and j
    gapsize = HAPLO[j,"whichsnp1"] - HAPLO[i,"whichsnp2"] - 1
    
    if(!overlapwithnext & gapsize == 1){
      
      commonfounder = HAPLO[i,foundercolnames] > 0 & HAPLO[j,foundercolnames] > 0
      commonfounder = as.numeric(commonfounder)/sum(commonfounder)
      new = HAPLO[i,]
      new[,foundercolnames] = commonfounder
      new$whichsnp2 = HAPLO[j,]$whichsnp2
      if("pos" %in% posinfo) new$pos2 = HAPLO[j,]$pos2
      if("cM" %in% posinfo) new$cM2 = HAPLO[j,]$cM2
      return(new)
      
    }else{return(NULL)}
  }))
  
  if(is.null(fused)){fused = HAPLO[which(F),]}
  ### Fuse the fused haplotype that overlap
  i = 1
  while(i < nrow(fused)){
    j=i+1
    #pos1=HAPLO[i,1];pos2=HAPLO[i,2]
    overlapwithnext = fused[j,"whichsnp1"] <= fused[i,"whichsnp2"]
    commonfounder = fused[i,foundercolnames]>0 & fused[i,foundercolnames]>0
    
    if(sum(commonfounder)>0 & overlapwithnext){
      new = fused[i,]
      new[,foundercolnames] = commonfounder
      new$whichsnp2 = fused[j,]$whichsnp2
      if("pos" %in% posinfo) new$pos2 = fused[j,]$pos2
      if("cM" %in% posinfo) new$cM2 = fused[j,]$cM2
      
      fused[i,]=new
      fused=fused[-j,]
    }else{i = i+1}
  }
  
  
  
  
  HAPLO = as.data.frame(rbind(HAPLO, fused))
  
  HAPLO = HAPLO[order(HAPLO$whichsnp1,HAPLO$whichsnp2),]
  
  h = 1
  while(h < nrow(HAPLO)){
    
    snp1=HAPLO[h,"whichsnp1"];snp2=HAPLO[h,"whichsnp2"]
    
    isEncased = HAPLO[,1] <= snp1 & HAPLO[,2] >= snp2
    if(sum(isEncased)>1){HAPLO=HAPLO[-h,]}else{h=h+1}
  }
  
  HAPLO
}



putBreakAtMidGeneticDistance=function(HAPLO,info){
  HAPLO = HAPLO[order(HAPLO$whichsnp1,HAPLO$whichsnp2),]
  
  for(i in 1:(nrow(HAPLO)-1)){
    j = i + 1
    
    if(HAPLO[i,"whichsnp2"]>HAPLO[j,"whichsnp1"]){
      cMi = HAPLO[i,"cM2"]
      cMj = HAPLO[j,"cM1"]
      
      newcM = (cMi+cMj)/2
      newpos = approx(info$cM,info$pos, xout = newcM)$y
      newsnp = which.min(abs(info$cM-newcM))
      
      HAPLO[i,"cM2"] = newcM
      HAPLO[j,"cM1"] = newcM
      
      HAPLO[i,"pos2"] = newpos
      HAPLO[j,"pos1"] = newpos
      
      HAPLO[i,"whichsnp2"] = newsnp-1
      HAPLO[j,"whichsnp1"] = newsnp
    }
  }
  
  #HAPLO = HAPLO[,-which(colnames(HAPLO) %in% c("whichsnp1", "whichsnp2"))]
  HAPLO
}





GetBreaks=function(HAPLO,info, mid=T){
  HAPLO = HAPLO[order(HAPLO$whichsnp1,HAPLO$whichsnp2),]
  
  snpsrange = apply(cbind(HAPLO[2:nrow(HAPLO),"whichsnp1"],HAPLO[1:(nrow(HAPLO)-1),"whichsnp2"]), 2, sort)
  cmrange = apply(cbind(HAPLO[2:nrow(HAPLO),"cM1"],HAPLO[1:(nrow(HAPLO)-1),"cM2"]), 2, sort)
  posrange = apply(cbind(HAPLO[2:nrow(HAPLO),"pos1"],HAPLO[1:(nrow(HAPLO)-1),"pos2"]), 2, sort)
  
  snpsrange = matrix(snpsrange, ncol=2)
  posrange = matrix(posrange, ncol=2)
  cmrange = matrix(cmrange, ncol=2)
  
  breaksrange = setNames(as.data.frame(cbind(snpsrange, posrange, cmrange)), c("whichsnp1"," whichsnp2", "pos1","pos2", "cM1", "cM2"))
  breaksrange$rilname = HAPLO$rilname[1]
  
  if(mid){
    newcM = (breaksrange$cM1+breaksrange$cM2)/2
    }else{
    weights = runif(nrow(breaksrange))
    newcM = breaksrange$cM1*weights + breaksrange$cM2*(1-weights)
    }
  
  newpos = approx(info$cM,info$pos, xout = newcM)$y
  newsnp = sapply(newcM, function(x){which.min(abs(info$cM-x))})
  
  breaksrange$interm.cM = newcM
  breaksrange$interm.pos = newpos
  breaksrange$interm.snp = newsnp
  breaksrange
}


Rilshaplotypes = do.call(rbind,lapply(split(Rilshaplotypes, Rilshaplotypes$rilname), function(HAPLO){
  
  print(HAPLO$rilname[1])
  #HAPLO=split(Rilshaplotypes, Rilshaplotypes$rilname)[[1]]
  #HAPLO = subset(Rilshaplotypes,rilname=="B_QG3161")
  HAPLO = FuseHapSeparatedByOneSNP(HAPLO=HAPLO, posinfo = c("whichsnp", "pos", "cM"),
                                   foundercolnames=foundercolnames)
  
  HAPLO
}))

breaks=do.call(rbind,lapply(split(Rilshaplotypes, Rilshaplotypes$rilname), function(HAPLO){
  
  if(nrow(HAPLO)>1){
    mid = GetBreaks(HAPLO,info, mid=T)[,c("rilname", "interm.cM", "interm.pos", "interm.snp")]
    mid$type = "mid"
    randompos = do.call(rbind, lapply(1:20, function(i){
      x=GetBreaks(HAPLO,info, mid=F)[,c("rilname", "interm.cM", "interm.pos", "interm.snp")]
      x$type = paste0("random_", i)
      x
    }))
    
    return(rbind(mid, randompos))
    }
  
}))

breaks = do.call(rbind, lapply(split(breaks, breaks$type), function(x){
  #x=split(breaks, breaks$type)[[1]]
  x=x[order(x$interm.snp),]
  x$cumbreaks = cumsum(rep(1,nrow(x)))
  x$relcumbreaks = x$cumbreaks/max(x$cumbreaks)
  x
}))

info=do.call(rbind, lapply(split(info, info$chrom), function(x){
  x$relcM = x$cM/max(x$cM)
  x
}))


nbreaks = do.call(rbind,lapply(split(Rilshaplotypes, Rilshaplotypes$rilname), function(HAPLO){
  data.frame(HAPLO$rilname[1], nbreaks = nrow(HAPLO))
}))


p1=ggplot()+
  geom_line(data=info, aes(pos/1e6,relcM), alpha=0.4, size=1, color="blue")+
  geom_line(data=subset(breaks, type!="mid"), aes(interm.pos/1e6, relcumbreaks, group=type), alpha=0.4, size=0.5, color="grey")+
  geom_line(data=subset(breaks, type=="mid"), aes(interm.pos/1e6, relcumbreaks), size=1.25)+
  theme_Publication3(base_size = 12)+
  xlab("Physical position (Mb)")+
  ylab("Relative genetic distance")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,20,1))


p2=ggplot()+
  geom_histogram(data=nbreaks, aes(nbreaks), binwidth=1, color="black", fill='grey')+
  theme_Publication2()+coord_cartesian(expand=0)+
  xlab("Number of breaks")

gridExtra::grid.arrange(p1, p2, nrow=1, widths = c(1.5,1))




pcumbreaks=ggplot()+
  #geom_line(data=info, aes(pos/1e6,relcM), alpha=0.4, size=1, color="blue")+
  geom_line(data=subset(breaks, type!="mid"), aes(interm.pos/1e6, cumbreaks, group=type), alpha=0.4, size=0.5, color="grey")+
  geom_line(data=subset(breaks, type=="mid"), aes(interm.pos/1e6, cumbreaks), size=2.5/.pt)+
  theme_Publication3(base_size = 10)+
  xlab("Physical position (Mb)")+
  ylab("Cum. n of breaks")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,20,1))


Rilshaplotypes=do.call(rbind,lapply(split(Rilshaplotypes, Rilshaplotypes$rilname), function(HAPLO){
  
  print(HAPLO$rilname[1])
  #HAPLO=split(Rilshaplotypes, Rilshaplotypes$rilname)[[1]]
  #HAPLO = subset(Rilshaplotypes,rilname=="A_QG3251")
  
  if(nrow(HAPLO)>1){HAPLO = putBreakAtMidGeneticDistance(HAPLO=HAPLO,info=info)}
  
  #HAPLO = HAPLO[,-which(colnames(HAPLO) %in% c("whichsnp1", "whichsnp2"))]
  
  HAPLO
  
}))


divideOverlapHapBlocksByFounderID = function(haplotype, foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2"), info){
  
  haplotype=haplotype[order(haplotype$whichsnp1, haplotype$whichsnp2),]
  
  snpmin = min(haplotype$whichsnp1)
  snpmax = max(haplotype$whichsnp2)
  
  rilname=haplotype$rilname[1]
  
  snphapmatrix=do.call(rbind, lapply(foundercolnames, function(thisfounder){
    
    whichsnphap = 1:snpmax
    snphap = rep(0, length(whichsnphap))
    thisfounderhaplotypes=haplotype[haplotype[,thisfounder]>0,]
    if(nrow(thisfounderhaplotypes)>0){
      for(i in 1:nrow(thisfounderhaplotypes)){
        wsnp1 = thisfounderhaplotypes[i,]$whichsnp1
        wsnp2 = thisfounderhaplotypes[i,]$whichsnp2
        snphap[wsnp1:wsnp2]=1
      }
    }
    
    snphap
  }))
  
  
  snphapid=apply(snphapmatrix, 2, function(x){paste(foundercolnames[x>0], collapse=";")})
  
  breaks = which(diff(as.numeric(factor(snphapid)))!=0)
  breaks=breaks[breaks>snpmin]
  breaks = c(snpmin, breaks, snpmax)
  
  
  haplotype=do.call(rbind,lapply(1:(length(breaks)-1), function(i){
    win=(breaks[i]+1):breaks[i+1]
    id = unique(snphapid[win])
    if(length(id)!=1) stop("multiple id for a single window")
    whichsnp1=win[1]
    whichsnp2=win[length(win)]
    
    nfounder=length(strsplit(id,";")[[1]])
    
    data.frame(whichsnp1= whichsnp1, whichsnp2= whichsnp2, founder=id, nfounder = nfounder)
  }))
  
  
  haplotype$rilname=rilname
  
  haplotype$pos1=info$pos[haplotype$whichsnp1]
  haplotype$pos2=info$pos[haplotype$whichsnp2]
  haplotype$cM1=info$cM[haplotype$whichsnp1]
  haplotype$cM2=info$cM[haplotype$whichsnp2]

  
  haplotype
  
}

#founderscolors = c("#e8b40e","#990000","#cfe2f3","#31496a","#93c71d", "#257403")


founderscolors = c("#BD274B","#422C1B", "#2A6318", "#DEA93F","#F8F2AA", "#3A54A0")
#founderscolors = c("#AA497F","#644329", "#428D4D", "#D47C3E", "#EFE956", "#2E3677")

haploPlotFromat = function(Rilshaplotypes, widthRILs = 0.85, foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2"),
                           founderscolors = c("#AA497F","#644329", "#428D4D", "#D47C3E"), info){
  
  haplop = do.call(rbind, lapply(split(Rilshaplotypes, Rilshaplotypes$rilname), function(haplotype){
    #haplotype = split(Rilshaplotypes, Rilshaplotypes$rilname)[[5]]
    haplotype=divideOverlapHapBlocksByFounderID(haplotype, foundercolnames=foundercolnames,info=info)
    haplotype=do.call(rbind,lapply(1:nrow(haplotype), function(i){
      #print(i)
      x=haplotype[i,]
      if(x$nfounder>1){
        founders=strsplit(x$founder,";")[[1]]
        rely = seq(0,1,1/length(founders))
        rely1 = rely[1:(length(rely)-1)]
        rely2 = rely[2:length(rely)]
        x=x[rep(1,length(rely1),),]
        x$rely1=rely1
        x$rely2=rely2
        x$founder=founders
      }else{
        x$rely1=0;x$rely2=1
      }
      x
    }))
    
    haplotype
  }))
  
  
  haplop$y1 = as.numeric(factor(haplop$rilname))-1 + (haplop$rely1)*widthRILs
  haplop$y2 = as.numeric(factor(haplop$rilname))-1 + (haplop$rely2)*widthRILs
  
  haplop$haplocolor=founderscolors[match(haplop$founder, foundercolnames)]
  
  haplop$haplocolor[is.na(haplop$haplocolor)]="black"
  
  haplop
  
}



#c("#BD274B","#422C1B", "#2A6318", "#DEA93F")

haplop=haploPlotFromat(Rilshaplotypes,
                      widthRILs = 0.85,
                      foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2"),
                      founderscolors = c("#BD274B","#422C1B", "#2A6318", "#DEA93F"),
                       info)

#"#E1E04A", "#2E3677"
#haplop$haplocolor2 = haplop$haplocolor
#haplop$haplocolor = haplop$haplocolor2


#haplop$haplocolor[haplop$haplocolor ==  "#E1E04A"] = "#F8F2AA"
#haplop$haplocolor[haplop$haplocolor ==  "#2E3677"] = "#3A54A0"
#haplop$haplocolor[haplop$haplocolor ==  "#2A6318"] =  "#8FC88C" #"#9BCC8C"



library(ggplot2)
pblocks = ggplot(haplop)+theme_Publication3(base_size = 10)+
  geom_rect(aes(xmin=pos1/1e6,xmax=pos2/1e6, ymin=y1, ymax=y2,fill=haplocolor), color = NA)+
  scale_fill_manual(values = unique(haplop$haplocolor), breaks = unique(haplop$haplocolor))+
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,20,1))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.y =element_blank(),
        axis.line.y = element_line(color = "white"),
        axis.ticks.y = element_line(color = "white"),
        axis.title = element_blank(),
        #axis.text.x = element_text(size = 15),
        panel.background = element_rect(fill = "black", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(10,10,10,10), "pt"))+
  ylab("")+xlab("")+
  xlab("Physical size (Mb)")




#ggplot(info, aes(pos,cM))+geom_point()
#At each snp we caluclate the number of possible founder
#nfounderthere=apply(snphapmatrix, 2, sum)
#We format in a 4 row matrix which match snphapmatrix
#nfounderthere = matrix(rep(nfounderthere,each=nrow(snphapmatrix)), nrow=nrow(snphapmatrix))

#snphapmatrix=snphapmatrix/nfounderthere

#snphapmatrix[,1535:1545]

# while(i<nrow(haplox)){
#   
# 
#   #Haplotype boundaries
#   b1  = haplotype[i,1]
#   b2 = haplotype[i,2]
#   
#   #osearch haplotypes thats overlap with boundaries 1 (b1)
#   # i.e., b1 is conatined between their own boundaries 
#   overlap1 = sort(which(haplotype[,1] < b1 & haplotype[,2] > b1))[1]
#   if(is.na(overlap1)){overlap1=NULL} # Just hapen when i = 1
#   #same for b2
#   overlap2 = max(which(haplotype[,1] < b2 & haplotype[,2] > b2))
#   
#   #sameplace are when two possible haplotype are occupies the same gap wetween two haplotype
#   # example:
#   # 11111111111111111          4444444444444444
#   #                222222222222222222
#   #          33333333333333333333
#   #=> Here both 2 and 3 may be the haplotype between 1 and 4
#   # Note that where the haplotype overlap it means that the founders have the same sequence there
#   # Also, Note that 2 or 3 does not contain the other. When it is the case, it if resolved within haplosearch function
#   # by chosing the (biggest) haplotype which contain the other
#   
#   sameplace = which(haplotype[,1] > haplotype[overlap1,1] & haplotype[,2] < haplotype[overlap2,2])
#   sameplace = sameplace[sameplace>i]
#   
#   if(length(sameplace)>0){print("same place here, function not robust, see what you can do")}
#   
#   #if there is an haplotype that overlap with b2, create a new haplotype
#   if(length(overlap2)>1){print("2 overlaping haplo, function not robust, see what you can do")}
#   
#   if(length(overlap2)==1){
#     
#     
#   }
#   
# }
# 
# 
# colnames()

getFounderFreq = function(Rilshaplotypes,foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2"), df=T){
  snpmin = min(Rilshaplotypes$whichsnp1)
  snpmax = max(Rilshaplotypes$whichsnp2)
  whichsnphap = snpmin:snpmax
  
  foundercount = matrix(0,nrow=length(foundercolnames),ncol=length(whichsnphap))
  totfoundercount = matrix(0,nrow=length(foundercolnames),ncol=length(whichsnphap))
  #haplotype = split(Rilshaplotypes, Rilshaplotypes$rilname)[[1]]
  for(i in 1:nrow(Rilshaplotypes)){
    if(i %% 100 == 0){print(i)}
    x=Rilshaplotypes[i,]
    wsnp1 = x$whichsnp1
    wsnp2 = x$whichsnp2
    foundervalue = matrix(rep(unlist(c(x[,foundercolnames])),length(wsnp1:wsnp2)), ncol=length(wsnp1:wsnp2))
    totfoundervalue = matrix(1, ncol=length(wsnp1:wsnp2),nrow=length(foundercolnames))
    
    foundercount[,wsnp1:wsnp2]=foundercount[,wsnp1:wsnp2]+foundervalue
    totfoundercount[,wsnp1:wsnp2]=totfoundercount[,wsnp1:wsnp2]+totfoundervalue
  }
  
  founderfreq=foundercount/totfoundercount
  rownames(founderfreq)=foundercolnames
  colnames(founderfreq)=whichsnphap
  
  if(df==T){
    founderfreq=reshape2::melt(founderfreq)
    colnames(founderfreq)=c("founder","whichsnp","freq")
  }
 
  founderfreq
  
}


founderfreq=getFounderFreq(Rilshaplotypes,foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2"), df=T)
  
pfreq = ggplot(founderfreq, aes(info$pos[founderfreq$whichsnp]/1e6, freq, color=founder))+
  geom_point(size=1.2/.pt)+ylim(0,0.6)+
  theme_Publication3(base_size=10)+theme(legend.position = "none")+
  scale_color_manual(breaks=c("FA.g1", "FA.g2","FB.g1", "FB.g2", "FM.g1", "FM.g2"),
                     values = c("#BD274B","#422C1B","#E1E04A", "#3A54A0", "#2A6318", "#DEA93F"))+
  xlab("Physical Position (Mb)")+ylab("Frequency")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,20,1))




p=grid.arrange(pblocks+theme(axis.title.x = element_blank(), plot.margin = unit(c(10,5,10,20), "pt")),
             pcumbreaks+theme(axis.title.x = element_blank()),
             pfreq, ncol=1, heights = c(10,2.2,2.2))



ggsave(filename="/Users/tomparee/Desktop/crossA_beceiHaplo.png", plot=p, height=7.4, width=4.6, dpi=600)



nbreaks=unlist(lapply(split(Rilshaplotypes, Rilshaplotypes$rilname), function(x){nrow(x)-1}))
hist(nbreaks)
breaks = do.call(rbind, lapply(split(Rilshaplotypes, Rilshaplotypes$rilname), function(x){
  #print(x$rilname[1])
  #x=split(Rilshaplotypes, Rilshaplotypes$rilname)[[1]]
  #x=subset(Rilshaplotypes, rilname=="A_QG3251")
  x=x[order(x$whichsnp1, x$whichsnp2),]
  
  if(nrow(x)>1){
    #Breaks between founder haplotype blocks
    # => between end and start of two adjacent haplotype blocks
    b = cbind(x$whichsnp2[-length(x$whichsnp2)],x$whichsnp1[-1])
    # => each row correspond to a single break
    # => for a given break b[i,]; if b[i,1] > b[i,2], it means that the two haplotypes overlap
    # (i.e, there is a region where the two founder haplotype are the same)
    #Let's just sort  b[i,1] & b[i,2] so b[i,2 is always higher, and indicate if there is an overlap or no
    b = as.data.frame(t(apply(b,1,function(y){
      newy = sort(y)
      c(newy, ifelse(y[1]>=y[2], "overlap", "gap"))
    })))
    
    b[,1] = as.numeric(b[,1]); b[,2] = as.numeric(b[,2])
    names(b) = c("whichsnp1", "whichsnp2", "type")
    b$rangesize = b$whichsnp2-b$whichsnp1
    b$rilname = x$rilname[1]
    b
  }else{
    NULL
  }
  
}))

#breaks$y=as.numeric(factor(breaks$rilname))


#ggplot(breaks)+theme_minimal()+
#  geom_segment(aes(x=whichsnp1, xend=whichsnp2, y=y,yend=y, color=rangesize))+
#  scale_color_gradient(low="black", high="lightgrey")


breaks=breaks[order(breaks$whichsnp1),]
breaks$cumbreaks = um(rep(1,nrow(breaks)))

breaks$pos = info$pos[breaks$whichsnp1]
ggplot(breaks)+theme_minimal()+
  geom_point(aes(x=pos/1e6, y=cumbreaks))





load("~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240501_beceiFounderPhasedHaplotypes_chrI.Rdata")


foundercolnames=colnames(phasedfounderhaplotypes)[5:10]
phasedfounderhaplotypes.df = reshape::melt(phasedfounderhaplotypes[,c("POS",foundercolnames)], id="POS")

foundercomb=expand.grid(foundercolnames,foundercolnames)
foundercomb = t(apply(foundercomb, 1, sort))
foundercomb=unique(foundercomb)
foundercomb=foundercomb[foundercomb[,1]!=foundercomb[,2],]


founderGenoDivergencePlots=lapply(1:nrow(foundercomb), function(i){
  pair = foundercomb[i,]

  xx=subset(phasedfounderhaplotypes.df, variable %in% pair)
  xx$value[xx$variable==pair[2]]=abs(xx$value[xx$variable==pair[2]]-1)
  xx$y = as.numeric(factor(xx$variable))
  colnames(xx)=c("pos","founder","value", "y")

  p=ggplot()+
    theme_minimal()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(color="black"),
          axis.text.y=element_text(face="bold", color="black"),
          legend.position = "none")+
    geom_point(data=xx,aes(pos/1e6,y,color=founder))+
    geom_line(data=xx,aes(pos/1e6,y, group=paste0(pos,value)), linewidth=0.001)+
    scale_color_manual(breaks=c("FA.g1", "FA.g2","FB.g1", "FB.g2", "FM.g1", "FM.g2"),
                       values = c("#90C98C","darkgreen","#D38C40","#AB4734","#4DBBD5FF","#3C5488FF"))+
    scale_y_continuous(breaks=c(1,2), labels = pair)+
    ylab("")+xlab("")
  
  p
  
})


gridExtra::grid.arrange(grobs=founderGenoDivergencePlots,ncol=3)

ggplot(test2)+
  geom_point(aes(pos,y))+
  geom_line(aes(pos,y, group=paste0(pos,value)), linewidth=0.001)







rils <- read_csv("~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/041724_becei_chr_I_Cross_AB_Founder_RIL_GT_table.csv")
load("~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240501_beceiFounderPhasedHaplotypes_chrI.Rdata")
mm = match(phasedfounderhaplotypes$POS, rils$POS)
rils=rils[mm,]
snps = phasedfounderhaplotypes[,1:4]
colnames(snps)=tolower(colnames(snps))
founders = phasedfounderhaplotypes[,5:10]
foudernames = colnames(founders)
founders = matrix(as.numeric(unlist(c(founders ))), ncol=ncol(founders))
colnames(founders)=foudernames
foundersA = founders[,grepl("A|M", colnames(founders))]
foundersB = founders[,grepl("B|M", colnames(founders))]

rils = rils[,grepl("_QG", colnames(rils))]
rilnames =  colnames(rils)
rils[rils=="1/1"] = "1"
rils[rils=="0/0"] = "0"
rils[rils=="0/1"]=NA
rils = matrix(as.numeric(unlist(c(rils))), ncol=ncol(rils))
colnames(rils) = rilnames

rilsA = rils[, grepl("A_Q",colnames(rils))]
rilsB = rils[, grepl("B_Q",colnames(rils))]



pos1 = 1.5e6
pos2 = 3e6

pos1= 11e6
pos2=12e6

win=which(snps$pos>pos1 & snps$pos <pos2)

founder_r2 = cor(t(phasedfounderhaplotypes[win,c("FA.g1", "FA.g2", "FM.g1", "FM.g2")]),use="pairwise.complete.obs")^2
colnames(founder_r2)=rownames(founder_r2)=win
founder_r2[lower.tri(founder_r2, diag=T)]=NA
founder_r2 = reshape2::melt(founder_r2) 
founder_r2 = founder_r2[!is.na(founder_r2$value),]
founder_r2$type ="founders"

rils_r2 = cor(t(rilsA[win,]),use="pairwise.complete.obs")^2
colnames(rils_r2)=rownames(rils_r2)=win
rils_r2[lower.tri(rils_r2, diag=T)]=NA
rils_r2 = reshape2::melt(rils_r2) 
rils_r2 = rils_r2[!is.na(rils_r2$value),]
rils_r2$type ="rils"

r2 = rbind(founder_r2, rils_r2)

r2$pos1 = snps$pos[r2$Var1]
r2$pos2 = snps$pos[r2$Var2]
r2$pdist=r2$pos2-r2$pos1

r2=do.call(rbind,lapply(split(r2, r2$type), function(x){
  bins = seq(1,max(x$pdist),2e4)
  x = do.call(rbind,lapply(1:(length(bins)-1), function(i){
    mindist = bins[i]
    maxdist = bins[i+1]
    #print(mindist)
    data.frame(dist1=mindist,dist2=maxdist,r2= mean(x$value[x$pdist>=mindist & x$pdist<maxdist]), type=x$type[1])
  }))
}))


ggplot(r2, aes((dist1+dist2)/2e3, r2, color=type))+
  geom_line()+theme_Publication2()+xlab("Physical distance (kb)")+
  ggtitle(paste0(pos1/1e6, " - ", pos2/1e6, " Mb"))+ylim(0,0.85)+
  theme(legend.position = "none")


r2ratio = cbind(r2[r2$type=="rils",1:2], r2_ratio = r2$r2[r2$type=="rils"]/r2$r2[r2$type=="founders"])

ggplot(r2ratio, aes((dist1+dist2)/2e3, r2_ratio))+
  geom_line(color='purple')+theme_Publication2()+xlab("Physical distance (kb)")+
  ggtitle(paste0(pos1/1e6, " - ", pos2/1e6, " Mb"))+ylim(0,1)




