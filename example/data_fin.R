#Analyse localisations Daten
# for paper Petersen et al., in preparation
# Volkmar Liebscher
# 2019-07-06
#
## load library
library("spatstat")
##create pp3 data 
#read in parsed data
dat0<-read.csv2("cell_coordinates.txt",sep="\t",dec=".")
#split different cells
l1<-split(dat0[,1:4],dat0$cell)
#View(l1[[1]])

#this one is ultimate way - every cell separately, split according to types
a1<-lapply(l1, function(v) split(ppx(v,coord.type=c(rep("s",3),"m"),simplify=TRUE),drop=TRUE))

#View(l6[[1]])
#how to plot single cell+type
# i<-1
# for (i in 1:length(a1))
# {pd<-a1[[i]][[1]]$data
# pq<-ppx(pd,coord.type=rep("s",3),domain=list(xrange=range(pd$x),
#                                              yrange=range(pd$y),
#                                              zrange=range(pd$z)),
#         simplify=TRUE)
# pr<-pp3(pd$x,pd$y,pd$z,domain=list(xrange=range(pd$x),yrange=range(pd$y),zrange=range(pd$z)))
# plot(pr,main=paste(i,names(a1[[i]])[1]))
# }
# #estimated  K-function vs. different models, 1 cell+type 
# pdf("a.pdf")
# ke<-K3est(pq,nrval=2^10)
# plot(iso~r,data = ke,type="l")
# lines(I(pi*r^2)~r,data=ke,col=2)
# lines(I(1.333*pi*r^3)~r,data=ke,col=3)
# lines(I(r)~r,data=ke,col=4)
# legend("topleft",col=1:4,lwd=5,legend = c("Daten","2D","3D","1D"))
# dev.off()

# #report gold  K-functions over cells
# pdf("Kacrosscells.pdf",one=TRUE)
# for (i in 1:length(a1)){pd<-a1[[i]]$gold$data
# pq<-ppx(pd,coord.type=rep("s",3),domain=list(xrange=range(pd$x),
#                                              yrange=range(pd$y),
#                                              zrange=range(pd$z)),
#         simplify=TRUE)
# ke<-K3est(pq,nrval=1024)
# plot(iso~r,data = ke,type="l",
#      main=paste("K Funktion Gold für ",names(a1)[[i]]))
# lines(I(pi*r^2)~r,data=ke,col=2)
# lines(I(1.333*pi*r^3)~r,data=ke,col=3)
# legend("topleft",col=1:3,lwd=5,legend = c("Daten","2D","3D"))
# }
# dev.off()

# different models for comparison
#box
runifboxpoint3<-function (n, width=1, nsim = 1, drop = TRUE) 
{
  domain <- as.box3(rep(c(0,1),3)*width)
  result <- vector(mode = "list", length = nsim)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  for (i in 1:nsim) {
    x <- runif(n)
    y <- runif(n)
    z <- runif(n)
    result[[i]] <- pp3(x*width, y*width, z*width, domain)
  }
  if (drop && nsim == 1) 
    return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}
#sphere
runifspherepoint3<-function (n, radius=1, jitter=0,nsim = 1, drop = TRUE) 
{
  domain <- as.box3(rep(c(-1,1),3)*radius)
  result <- vector(mode = "list", length = nsim)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  for (i in 1:nsim) {
    x <- rnorm(n)
    y <- rnorm(n)
    z <- rnorm(n)
    d<-sqrt(x^2+y^2+z^2)
    xx<-x/d*radius+rnorm(n,m=0,sd=jitter/sqrt(3))
    yy<-y/d*radius+rnorm(n,m=0,sd=jitter/sqrt(3))
    zz<-z/d*radius+rnorm(n,m=0,sd=jitter/sqrt(3))
    result[[i]] <- pp3(xx, yy,zz, domain)
  }
  if (drop && nsim == 1) 
    return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}
#halfsphere
runifhalfspherepoint3<-function (n, radius=1, jitter=0,nsim = 1, drop = TRUE) 
{
  domain <- as.box3(c(-1,1,-1,1,0,1)*radius)
  result <- vector(mode = "list", length = nsim)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  for (i in 1:nsim) {
    x <- rnorm(n)
    y <- rnorm(n)
    z <- rnorm(n)
    d<-sqrt(x^2+y^2+z^2)
    xx<-x/d*radius+rnorm(n,m=0,sd=jitter/sqrt(3))
    yy<-y/d*radius+rnorm(n,m=0,sd=jitter/sqrt(3))
    zz<-abs(z/d)*radius*rbinom(n,prob=2/3,size=1)+rnorm(n,m=0,sd=jitter/sqrt(3))
    result[[i]] <- pp3(xx, yy,zz, domain)
  }
  if (drop && nsim == 1) 
    return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

#(double) rings
runifringspoint3<-function (n, radius=1, hdistance=0, jitter=0,nsim = 1, drop = TRUE) 
{
#unit cube, anyway -> comparing K functions
    domain <- as.box3(rep(c(-radius,radius),3))
  result <- vector(mode = "list", length = nsim)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  for (i in 1:nsim) {
    x <- rnorm(n)
    y <- rnorm(n)
    z <- rbinom(n=n,prob=0.5,size=1)
    d<-sqrt(x^2+y^2)
    xx<-x/d*radius+rnorm(n,m=0,sd=jitter/sqrt(3))
    yy<-y/d*radius+rnorm(n,m=0,sd=jitter/sqrt(3))
    zz<- hdistance*(2*z-1)+rnorm(n,m=0,sd=jitter/sqrt(3))
    result[[i]] <- pp3(xx, yy,zz, domain)
  }
  if (drop && nsim == 1) 
    return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}


#plot simulation results for double checking
#scale ratios similar to observed ones 
pdf("phantomsin2.pdf")
#X11()
#layout
par(mfrow=c(2,2),mai=rep(0,4),oma=c(0,0,1,0))
#seed for simulations
set.seed(123)
plot(runifboxpoint3(n=1e3),main=NULL)#"complete randomness")
#plot(runifspherepoint3(n=1e3),main="spherical")
plot(runifspherepoint3(n=1e3,jitter = 1/16),main=NULL)#"spherical,jittered")
#plot(runifhalfspherepoint3(n=1e3),main="half spherical")
#plot(runifringspoint3(n=1e3,hd=.125))#,main="two rings")
plot(runifringspoint3(n=1e3,hd=0,jit=1/16),main=NULL)#"two rings,jittered")
#plot(runifringspoint3(n=1e3,hd=0),main="one ring")
plot(runifringspoint3(n=1e3,hd=0,jit=4/16),main=NULL)#"one ring, jittered")
dev.off()


#pdf("gold_geometry_across_cells.pdf",one=TRUE)
pdf("gold_geometry_across_cells_new%02d.pdf",one=FALSE)
#Parameters for envelope - resolution of distance and No. of simulations
nr<-2^8
#for 5%pointwise confidence
ns<-1599
nrank<-40
#seed for simulations
set.seed(123)

#Names -now anonymous
nam<-c(paste("mutant cell No.",1:8),paste("WT cell No.",1:6))
#jitter Distance(radius of cylinder)
jr<-3.25e-2
for (i in 1:length(a1)){
  #construct point clouds -gold
  pd<-a1[[i]]$gold$data
  #No. points
  ng<-dim(pd)[1]
  #point configuration
  p.g<-ppx(pd,coord.type=rep("s",3),domain=list(xrange=range(pd$x),
                                                yrange=range(pd$y),
                                                zrange=range(pd$z)),
           simplify=TRUE)
  #Box scale
  r.g<-0.5*max(unlist(lapply(list(xrange=range(pd$x),
                                  yrange=range(pd$y),
                                  zrange=range(pd$z)),FUN=diff)))
  r.g2<-0.5*median(unlist(lapply(list(xrange=range(pd$x),
                                      yrange=range(pd$y),
                                      zrange=range(pd$z)),FUN=diff)))
  #envelops - CSR
  env.g<-envelope(p.g,K3est,nrval=nr,nsim=ns,nrank=nrank,VARIANCE = FALSE,correction="trans")
  #spere as comparison
  env.g1<-envelope(p.g,K3est,nrval=nr,nsim=ns,nrank=nrank,
                   simulate=expression(runifspherepoint3(radius=r.g,n=ng,jit=jr)),VARIANCE = FALSE,correction="trans")
  #one ring as comparison
  env.g2<-envelope(p.g,K3est,nrval=nr,nsim=ns,nrank=nrank,
                   simulate=expression(runifringspoint3(radius=r.g,hd=0,n=ng,jit=jr)),VARIANCE = FALSE,correction="trans")
  #more jittering
  env.g3<-envelope(p.g,K3est,nrval=nr,nsim=ns,nrank=nrank,
                   simulate=expression(runifringspoint3(radius=r.g,hd=0,n=ng,jit=2*jr)),VARIANCE = FALSE,correction="trans")
  #even more jittering
  env.g4<-envelope(p.g,K3est,nrval=nr,nsim=ns,nrank=nrank,
                   simulate=expression(runifringspoint3(radius=r.g,hd=0,n=ng,jit=3*jr)),VARIANCE = FALSE,correction="trans")
  
  #K functions for gold particles with envelops - important
  plot(env.g,main=nam[[i]],legend=FALSE,shade=NULL,add=FALSE,xlim=c(0,0.5*r.g),
       xlab=expression(italic(r)~group("[",{mu*m},"]")),
       ylab=expression(K[3](r)))
  plot(env.g1,add=TRUE,col="blue",shadecol=rgb(0,0,1,alpha=0.3))
  plot(env.g4,add=TRUE,col="red",shadecol=rgb(1,0.2,1,alpha=0.3))
  plot(env.g3,add=TRUE,col="red",shadecol=rgb(1,0.2,0.2,alpha=0.3))
  plot(env.g2,add=TRUE,col="orange",shadecol=rgb(1,1,0,alpha=0.3))
  plot(env.g,add=TRUE,col=c("black",rep("darkgreen",3)),shadecol=rgb(0.2,1,0.2,alpha=0.3))
}
dev.off()

