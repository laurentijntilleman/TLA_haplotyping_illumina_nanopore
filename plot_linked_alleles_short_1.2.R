
plot_linked_alleles <- function(networkfile,locus,mutpos=NULL){
    
  network <- read.table(networkfile, header=TRUE, stringsAsFactors = F)
  
  cleaned <- subset(network, network$interaction=="pp")



  haplotype1_links <- subset(network, network$hap==1)
  haplotype2_links <- subset(network, network$hap==2)
  rootsplit <- strsplit(networkfile, split="/")[[1]]
  rootname <- strsplit(rootsplit[length(rootsplit)], split="_network.txt")[[1]][1]
  
  locussplit <- strsplit(locus, split=":")[[1]][2]
  locussplit_coor <- strsplit(locussplit, split="-")[[1]]
  


  nuc.pos <- unique(c(haplotype1_links$pos, haplotype1_links$pos.1, haplotype2_links$pos, haplotype2_links$pos.1))
  #####select rounded 5% and 95% quantile as size
  size <- floor(diff(quantile(nuc.pos, c(0.05,0.95)))/1e3) 
  
  #Take provided range as x-axis to show
  left_edge <- as.numeric(locussplit_coor[1])
  right_edge <- as.numeric(locussplit_coor [2])
  
  N_positions_hap1 <- length(unique(c(haplotype1_links$pos, haplotype1_links$pos.1)))
  N_positions_hap2 <- length(unique(c(haplotype2_links$pos, haplotype2_links$pos.1)))


  distance <- right_edge-left_edge
  dlength <- (distance/100)*5
  arrow_X <- ((left_edge+right_edge)/2)


  plot(range(c(left_edge,right_edge)), c(-3,3), type='n', 
      axes=F, xlab=sprintf("%s bp",format(distance,big.mark=",",scientific=FALSE)), ylab="", main=rootname)
  axis(1,c(left_edge, right_edge))


  arrows(arrow_X-dlength,-3.5,arrow_X+dlength,-3.5,xpd = TRUE, code=3, length=0.15, lwd=1.2) # Arrow for the distance below x axis..
  text(left_edge,1,N_positions_hap1)
  text(left_edge,-1,N_positions_hap2)

  abline(h=0)
  if(!is.null(mutpos)){
    abline(v=as.numeric(mutpos))
  }
  for(i in 1:nrow(haplotype1_links)){
    dat <- haplotype1_links[i,]
    xspline(c(dat$pos, (dat$pos+dat$pos.1)/2, dat$pos.1), c(0, 3,0), c(0,1,0), border=rgb(1,0,0,0.3), lwd=1)
  }	
  
  for(i in 1:nrow(haplotype2_links)){
    dat <- haplotype2_links[i,]
    xspline(c(dat$pos, (dat$pos+dat$pos.1)/2, dat$pos.1),c(0, -3,0), c(0,1,0), border=rgb(0,0,1,0.3), lwd=1)
    
  }
}

