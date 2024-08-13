source("/Users/tomparee/Documents/Documents - MacBook Pro de tom/basics_TP.R")

load("~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240515_beceiFounderPhasedHaplotypes_chrV.Rdata")


foundercolnames=colnames(phasedfounderhaplotypes)[5:10]
phasedfounderhaplotypes.df = reshape::melt(phasedfounderhaplotypes[,c("POS",foundercolnames)], id="POS")

foundercomb=expand.grid(foundercolnames,foundercolnames)
foundercomb = t(apply(foundercomb, 1, sort))
foundercomb=unique(foundercomb)
foundercomb=foundercomb[foundercomb[,1]!=foundercomb[,2],]


splithere = seq(0,max(phasedfounderhaplotypes$POS),1e5)
snpdensity=do.call(rbind, lapply(1:(length(splithere)-1), function(i){
  pos1=(splithere[i]+1)
  pos2=splithere[i+1]
  data.frame(pos1=pos1,pos2=pos2,nsnp=sum(phasedfounderhaplotypes$POS>=pos1 & phasedfounderhaplotypes$POS<pos2))
}))

ggplot(snpdensity)+
  geom_rect(aes(xmin=pos1/1e6,xmax=pos2/1e6,ymin=0,ymax=nsnp))+
  theme_Publication2()+xlab("Physical distance (Mb)")+
  ylab("nSNP")

#ggplot(phasedfounderhaplotypes)+
#  geom_point(aes(x=POS/1e6,y=genotypedistance))+
#  theme_Publication2()+xlab("Physical distance (Mb)")+
#  ylab("Genotype distance with poolseq")


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




