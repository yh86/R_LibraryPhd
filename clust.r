###
### library of modules related to clustering algorithms
###
###

## clust::meanplotSOM
#   plot the mean (centroid) and stand deviation of Self-Organizing maps
#
# PARAMETERS
#   ds:               data frame includes class indicator
#   ci:               swith on computing confidence interval of mean using ysmean.cl.boot
#   ix.class:         index of where class information can be found
#   grid.r:           number of rows of the grid (plotting) system
#   grid.c:           number of columns of the grid (plotting) system
#
# DEPENDENCIES
#   lib::getCentrd
#
# USAGE
#   setwd("C:\\Users\\yongsheng\\Documents\\Yongsheng's Documents\\UM\\research\\hero\\phd\\analysis\\2009.03.02.Draft.FLU.TSA\\2009.05.14.SOM.ggplot")
#   test = read.csv("2009.05.13.Boosting.Test.bsplcls.SCALE.csv",header=T,row.names=1)
#   test.cls=read.csv("2009.05.13.Boosting.Test.bsplcls.SCALE.classSOM.csv",header=T,row.names=1)
#   # make sure they are same order
#   sum(rownames(test)!=rownames(test.cls))
#   test$classSOM = test.cls$classSOM
#   meanplotSOM(ds=test, n.class=8, ix.class=31, grid.r=2,grid.c=4,scale_y_min=-3,scale_y_max=3)
#
## NOTE
##  this does not work
##  the mean and standard deviation is too small relative to the whole data range
##  we will generate a new data set with summarized values only and then plot on that dataset
#   test.plot = melt(test,id="classSOM")
#   test.plot$variable=gsub('H','',test.plot$variable)
#   c = qplot(value, x=variable, data=test.plot, xlab="",ylab="",geom="blank",group="classSOM") + facet_grid(classSOM~.)
#   c = c + stat_sum_dfmX(  fun="mean_cl_boot"
#                          , geom="crossbar",colour=rep("darkgreen",15),width=0.1
#                          , data=test.plot[grep("A",test.plot$variable),]
#                          )
#   c = c + stat_sum_dfmX(  fun="mean_cl_boot"
#                          , geom="crossbar",colour=rep("red",15),width=0.1
#                          , data=test.plot[grep("S",test.plot$variable),]
#                          )
#
meanplotSOM <- function (ds=NULL, ci=FALSE, which.class=0, n.class=0, ix.class=0, grid.r=2, grid.c=4, grid.ht=4, grid.wt=3, scale_y_min=-3, scale_y_max=3) {
  # clust::meanplotSOM

  if(is.null(ds)) stop("Where is the input, you stupid?")
  n.class = length(unique(ds[,ix.class]))
  if((grid.r*grid.c)!=n.class & grid.r!=1 & grid.c!=1) stop("number of class desn't fit grid")
  if(ix.class==0) stop("index of column that contains SOM class identifier")

  # set the b&w theme
  theme_set(theme_bw())
  if(!ci) { # we plot normal mean +/- stdev
    # compute the mean and standard deviation
    test.plot.mean = getCentrd(ds=ds,ixclass=ix.class,fun="mean")
    test.plot.sd = getCentrd(ds=ds,ixclass=ix.class,fun="sd")
    # melt the two data frames
    test.plot.mean = melt(test.plot.mean,id="class")
    test.plot.sd = melt(test.plot.sd,id="class")
    # make sure the order before we merge the two
    test.plot.mean[,'class']==test.plot.sd[,'class'] & test.plot.mean[,'variable']==test.plot.sd[,'variable']
    test.plot = cbind(test.plot.mean, test.plot.sd[,'value'])
    colnames(test.plot)=c("class","x","mean","sd")
  }
  
  if(ci) {  # we plot mean +/- 95% confidence interval
    # compute the mean and 95% CI
    test.plot.mean = getCentrd(ds=ds,ixclass=ix.class,fun="ysmean.cl.boot",what="mean")
    test.plot.meanlo = getCentrd(ds=ds,ixclass=ix.class,fun="ysmean.cl.boot",what="lower")
    test.plot.meanhi = getCentrd(ds=ds,ixclass=ix.class,fun="ysmean.cl.boot",what="upper")
    # melt the data frames
    test.plot.mean = melt(test.plot.mean,id="class")
    test.plot.meanlo = melt(test.plot.meanlo,id="class")
    test.plot.meanhi = melt(test.plot.meanhi,id="class")
    # make sure the order before we merge the two
    tmp = test.plot.mean[,'class']==test.plot.meanlo[,'class'] & test.plot.mean[,'variable']==test.plot.meanlo[,'variable']
    tmp = tmp & test.plot.mean[,'class']==test.plot.meanhi[,'class'] & test.plot.mean[,'variable']==test.plot.meanhi[,'variable']
    if (!tmp) stop("data frames containing summary statistics are not melt or merged correctly")
    # merge the three data frames
    test.plot = cbind(test.plot.mean, test.plot.meanlo[,'value'], test.plot.meanhi[,'value'])
    colnames(test.plot)=c("class","x","mean","lower","upper")
  }
  
  test.plot$Pheno = as.factor(left(test.plot$x,1))

  grid.newpage()
  grid.rect(gp=gpar(fill="white"))
  pushViewport(viewport(layout=grid.layout(grid.r, grid.c)))
  all.class = unique(test.plot$class)
  print("Now i will print the following classes:"); print(all.class)

  # set the error bar correctly according to whether a regular mean+/-stdev or bootstrap mean+/-95%ci
  if(!ci) {limits <- aes(ymax = mean + sd, ymin=mean - sd)}
  if(ci)  {limits <- aes(ymax = upper, ymin=lower)}
  cnt.class = 0
  for (i in 1:grid.r)
    for (j in 1:grid.c) {
      # multiple plots in one grid
      if((grid.r>1 | grid.c>1) & which.class==0) cnt.class = cnt.class+1
      # one plot for one class at a time
      else cnt.class = cnt.class+which.class
      if(cnt.class>n.class) stop("Hrr...you try to plot more than you have...stupid")
      tmp = test.plot[test.plot$class==all.class[cnt.class],]
      tmp.viewplot = viewport(  width=unit(grid.ht,"inches"),height=unit(grid.wt,"inches")
                             ,layout.pos.row=(grid.r-i+1),layout.pos.col=j)
      tmp.gplot = ggplot(tmp, aes(colour=Pheno,y=mean,x=x)) + facet_grid(class~.) + xlab(cnt.class)
      tmp.gplot = tmp.gplot + geom_point() + geom_line() + geom_errorbar(limits,width=1)
      tmp.gplot = tmp.gplot + geom_hline(aes(yintercept=0),colour=I(alpha("black",10/10)))
      tmp.gplot = tmp.gplot + scale_y_continuous(limits=c(scale_y_min,scale_y_max))
      # 2009-09-21
      # added for better visualization
      tmp.gplot = tmp.gplot + opts(   #axis.line=theme_segment(size=1)
                                     panel.border=theme_rect(size=0.7)
                                    , axis.text.y=theme_text(face="bold", size=8, hjust=1)
                                    , axis.text.x=theme_text(face="bold", size=4, hjust=1, vjust=1, angle=45)
                                    , axis.ticks=theme_segment(size=0.6)
                                    , axis.title.x=theme_text(face="bold",size=6,vjust=0)
                                    #, axis.title.y=theme_text(face="bold",size=11,angle=90)
                                    , axis.title.y=theme_blank()
                                    , legend.key=theme_rect(size=0.7)
                                    , legend.text=theme_text(face="bold",size=8)
                                    , plot.margin=unit(rep(0.5,4),"lines")
                                  )
       tmp.gplot = tmp.gplot + scale_colour_manual(value=c("darkgreen","red"))
      
      # plot the color by true RGB color
      # tmp.gplot = tmp.gplot + scale_colour_manual(values=c("blue","red"))
      print(tmp.gplot, vp=tmp.viewplot)
    }
  popViewport()

  return ()

}


##############################################################################################################
# clust::distClass
#   compute Euclidean distance from each observations to class centroid
#
# PARAMETERS
#   ds:                 data set (eg genes by time) with a column indicating class memebership of observations
#   cls.name:           name of the column containing class membership
#   sum.fun:            function to be used to summarize values of each class
#   do.plot:            whether or not to plot the distance
#
# RETURN VALUE(S)
#   ret:                a data frame contains for each observation (gene)
#                       > the distance to its class centeroid
#                       > its original class memebership
#
# DEPENDENCIES
#   1)  lib::getCentrd
#   2)  library(ggplot2)
#
##############################################################################################################
distClass <- function(ds=NULL, cls.name="classSOM", sum.fun="mean", do.plot=FALSE, ... ) {

  ret = NULL

  test = ds
  test.centrd = getCentrd(ds=test,ixclass=which(colnames(test)==cls.name),fun=sum.fun,isRow=FALSE)
  # compute Eucledian distance
  ixclass = which(colnames(test)==cls.name)
  x = apply(    test
              , 1
              , function(x) {
                    centrd = test.centrd[test.centrd$class==x[ixclass],-which(colnames(test.centrd)=='class')]
                    return(distEucl(x[-ixclass],centrd))
                }
            )
  test.distSOM = cbind(distSOM = x
                      ,classSOM = test[match(attr(x,"names"),rownames(test)),cls.name])

  test.distSOM = data.frame(test.distSOM)
  test.distSOM[,cls.name] = as.factor(test.distSOM[,cls.name])
  if (do.plot) {
    theme_set(theme_bw())
    aggplot = ggplot(test.distSOM, aes(distSOM, colour=test.distSOM[,cls.name]))
    aggplot = aggplot + stat_bin( aes(size = ..density..,), binwidth=0.5, geom="point", position="jitter")
    print(aggplot)
  }

  ret = test.distSOM

  return(ret)

}


cat(">>>>>>> Clustering library has been successfully loaded!\n")