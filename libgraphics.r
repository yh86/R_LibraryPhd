#***********************************************************************************************
# libgraphics::plotPtAvg
#    plot data points and the average of the data points as horizontal line
# 
# PAMETERS
#   y:    a vector of values to be plotted
#   x:    a vector of group labels
#
# VALUE
#   None
#
# DEPENDENCIES
#   1) assume the figure is currently active in the console
#
# DBG
#   None
#
# NOTES
#
# USAGE
#    y = c(rnorm(5, mean=1), rnorm(5, mean=5), rnorm(5,mean=20), rnorm(5,mean=10))
#    x = c( rep("a",5), rep("b",5), rep("c",5),rep("d",5))
#    dp.avg(y,x)
#
# KNOWN ISSUES
#
# TIME STAMP
#    2010/07/19
#***********************************************************************************************
plotPtAvg <- function (y, x, lwd=3, x.step=0.5,cex=1.4,jitters=3) {

  #  cex = 1.4
  #  lwd=3
  #  x.step = 0.5
  
  x.dither = rnorm(length(x), mean=0, sd=0.05)
  
  x.tab = table(x)
  xlim = seq(along=x.tab, by=x.step)
  x.new = NULL
  for (i in 1:length(xlim)) {
    x.new = c(x.new, rep(xlim[i],x.tab[i])) 
  }
  
  xlim = c(min(xlim)-x.step/2, xlim, max(xlim)+x.step/2)
  
  
  whatColor = "red"
  par(mar=c(3,3,2,2))
  #plot ( y, x=as.factor(x), type="p")
  plot(jitter(y,factor=jitters),x=x.new+x.dither,type="p", pch=17, col=whatColor, cex=cex
       , xlim=c(min(xlim),max(xlim)), axes=FALSE, main="", xlab="", ylab="")
  
  avg = aggregate(y, by=list(x.new), FUN=mean)
  
  dither.max = 5/8 * (max(x.dither) + abs(min(x.dither)))
  
  for (i in 1:nrow(avg)) {
    x1 = avg[i,1] - abs(dither.max*1.1)
    y1 = avg[i,2]
    y2 = avg[i,2]
    x2 = avg[i,1] + abs(dither.max*1.1)
    segments (x1,y1,x2,y2, lwd=lwd, col=whatColor)
  }
  
  box(lwd=lwd, col="#000000")  
  xlim = xlim[2:(length(xlim)-1)]
  axis(1, at=xlim, labels=names(x.tab), lwd=lwd, lwd.ticks=lwd, font=2,cex.axis=cex)
  axis(2, lwd=lwd, lwd.ticks=lwd,font=2, cex.axis=cex)
  
}



#***********************************************************************************************
# libgraphics::fineFig
#   fine tuning some components of a figure to make it prettier
#
# AUTHOR
#   Yongsheng Huang   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PAMETERS
#   NULL
#
# VALUE
#   None
#
# DEPENDENCIES
#   1) assume the figure is currently active in the console
#
# DBG
#   None
#
# NOTES
#
# USAGE
#  plot(ecdf(rnorm(10000)), col="red", col.01line="grey70", pch="", lwd=3)
#  fineFig()
#
# KNOWN ISSUES
#
# TIME STAMP
#   2009-09-21
#***********************************************************************************************
fineFig <- function () {
 box(lwd=3)
 axis(side=2,lwd=2,font=2)
 axis(side=1,lwd=2,font=2)
}



#***********************************************************************************************
# libgraphics::ys_theme
#   set my default ggplot theme
#   (larger font size, bold border, bold font, bold axis, dotted grid major line)
#
# AUTHOR
#   Yongsheng Huang   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PAMETERS
#   axis.text.x.angle   -   rotation angle of horizonal axis tick label
#
# VALUE
#   None
#
# DEPENDENCIES
#   1) require(ggplot2)
#
# DBG
#   None
#
# NOTES
#
# USAGE
#   ys_theme(axis.text.x.angle=30)
#   qplot(carat, price, data=diamonds)
#
# KNOWN ISSUES
#
# TIME STAMP
#   2009-09-21
#***********************************************************************************************

ys_theme <- function (axis.text.x.angle=0) {
  require(ggplot2)
  theme_set(theme_bw())
  theme_update(
                   panel.border=theme_rect(size=1.5)
                  , strip.background = theme_rect()
                  , panel.grid.major=theme_line(colour="black",size=0.5,linetype="dotted")
                  , axis.text.y=theme_text(face="bold", size=12, hjust=1)
                  , axis.text.x=theme_text(face="bold", size=12, vjust=1, angle=axis.text.x.angle)
                  , axis.ticks=theme_segment(colour="black", size=1, linetype=1.5)
                  , axis.title.x=theme_text(face="bold",size=14)
                  , axis.title.y=theme_text(face="bold",size=14,angle=90)
               )
}



#***********************************************************************************************
# libgraphics::ys.tiff
#   wrapper function of {grDevices}::tiff
#
# AUTHOR
#   Yongsheng Huang   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PAMETERS
#   fname       -   name of outupt (.tif) file
#   fsize       -   size of the output (.tif) file
#
# VALUE
#   none
#
# DEPENDENCIES
#
# DBG
#   setting parmeters for running the function
#
# NOTES
#
# USAGE
#   fname = "test.tif"
#   opt.size = c(1280, 1280)
#   ys.tiff(fname=fname, fsize=opt.size)
#   plot(rnorm(100))
#   dev.off()
#
# KNOWN ISSUES
#
#
# TIME STAMP
#   2009-09-03
#***********************************************************************************************
ys.tiff <- function (fname, fsize) {
  tiff(   filename=fname
        , width=fsize[1], height=fsize[2], units="px"
        , compression = "none"
        , bg = "white"
        , res = 300
        , restoreConsole = TRUE)
}




# use Vennable instead
#library(Vennerable)
#setwd("C:\\Users\\yongsheng\\Documents\\Yongsheng's Documents\\UM\\research\\hess\\ChipSeq\\analysis\\2009.12.30.Draft\\2009.12.30.Venn\\")
#hoxa = Venn (SetNames=c("",""), Weight=c(`01`=103,`11`=696,`10`=2660))
#meis = Venn (SetNames=c("",""), Weight=c(`01`=1851,`11`=1093,`10`=1895))
#pdf("hoxa_meis_venn.pdf")
#plot(meis, doWeights=TRUE, type="circles", show = list(SetLabels = FALSE, Faces = FALSE, DarkMatter = FALSE))
#plot(hoxa, doWeights=TRUE, type="circles", show = list(SetLabels = FALSE, Faces = FALSE, DarkMatter = FALSE))
#dev.off()


# horrible: don't use 
# no proportion
ysVenn <- function(ds=NULL) {
  #***********************************************************************
  # ysVenn
  #   plot Venn diagram given a data frame
  #
  # PARAMETERS
  #   ds  -   data frame with two columns "names" and "group"
  #
  #  USAGE
  #    test = read.csv("venn.csv",header=T,stringsAsFactors=F)
  #    ds = test
  #    ysVenn(ds)
  #***********************************************************************
  if(colnames(ds)!=c("names","group"))
     stop("Give me two columns titled names and group, you stupid!")
  library(gplots)
  ds$group=as.factor(ds$group)
  par(font=2,lty=2,lwd=2,cex=2, mai=rep(0,4))
  venn(ds$names,ds$group,main=NULL)
  return(0)
}

###########################################################################
# graphics::stat_sum_singleX
#   summarize the value at each single X with given function
#
# PARAMETERS
#   fun         -       function to be used to summarize data
#   colour      -       colour for plot
#   geom        -       geometry object to be used in plot
#   size        -       size of geometry object
#
# TIME STAMP
#   2009-05-13
#
# USAGE
#   c = qplot(cyl, mpg, data=mtcars)
#   (c = c + stat_sum_singleX(mean,geom="point",colour="darkblue",size=2))
###########################################################################
stat_sum_singleX <- function(fun, geom="point", colour="red", size=3, ...) {
  stat_summary(fun.y=fun, colour=colour, geom=geom, size = size, ...)
}


###########################################################################
# graphics::stat_sum_dfmX
#   summarize the value on a data.frame with given function
#   (Hmisc package provides a set of useful summary functions)
#
# PARAMETERS
#   fun         -       function to be used to summarize data
#   colour      -       colour for plot
#   geom        -       geometry object to be used in plot
#   size        -       size of geometry object
#
# DEPENDENCIES
#   none
#
# TIME STAMP
#   2009-05-13
#
# NOTE
#   see ?stat_summary for more information
#
# USAGE
#    c = qplot(cyl, mpg, data=mtcars)
#    c = c + stat_sum_dfmX(fun="mean_cl_boot",geom="errorbar",colour="darkgreen",width=0.1)
#    (c = c + stat_sum_dfmX(fun="mean_cl_normal",geom="errorbar",colour="red",width=0.1))
###########################################################################
# Alternatively, you can supply a function that operates on a data.frame.
# A set of useful summary functions is provided from the Hmisc package:

# this will interfere with colour plotting if we want different color at different X tick
#stat_sum_dfmX <- function(fun, geom="crossbar", colour="red", width=0.7, ...) {
#  stat_summary(fun.data=fun, colour=colour, geom=geom, width=width, ...)
#}

stat_sum_dfmX <- function(fun, geom="crossbar", width=0.7, ...) {
  stat_summary(fun.data=fun, geom=geom, width=width, ...)
}


################################################################################
#
# graphics::stat_sum_smooth
#   linear smoothing with 95% CI
#
# TIME STAMP
#     2009-05-13
#
# PARAMETERS
#   fun:        function to be used for summarization
#   geom:       geometry object to be plotted
#   colour:     colour of geometry object
#   degree:     degree of smoother
#   width:      width of geometry object
#   alpha:      opaque of geometry
#   size:       size of geometry object
#
# RETURN VALUES
#   none
#
# USAGE
#   c <- qplot(cyl, mpg, data=mtcars)
#   c + stat_sum_smooth("median_hilow")
################################################################################

stat_sum_smooth <- function(fun, geom="smooth", colour="red", degree=4, width=0.2, alpha=5/10, size=0.4, ...){
  stat_summary(   fun.data=fun
                , method="lm", formula=y~smooth.spline(x,degree)
                , colour=colour
                , geom=geom
                , width=width
                , alpha=alpha
                , size=size
                , ...)
}

################################################################################
#
# libgraphics::circle.cor
#   circles representation of correlation matrix
#
# TIME STAMP
#     2009-06-19
#
# PARAMETERS
#   ds:                 a data frame to be plotted
#   main:               title (main)
#   pch:                point character
#   color:              color to be used
#
# RETURN VALUES
#   none
#
# DEPENDENCIES
#   library(kohonen)
#
# USAGE
# ds = iris[1:4]
# mycolor = c("red", "green3", "blue")[unclass(iris$Species)]
# mypairs (ds, main="Pairs plot", pch=21, color=mycolor)
################################################################################

mypairs <- function(ds=NULL, main="", bg=NULL, pch=NULL, color=NULL) {

  pairs(ds
        , main = main
        , pch = pch
        , bg = bg
        , diag.panel = panel.hist <- function(x, ...)
                                     {
                                        usr <- par("usr"); on.exit(par(usr))
                                        par(usr = c(usr[1:2], 0, 1.5) )
                                        h <- hist(x, plot = FALSE)
                                        breaks <- h$breaks; nB <- length(breaks)
                                        y <- h$counts; y <- y/max(y)
                                        rect(breaks[-nB], 0, breaks[-1], y, col="gray60", ...)
                                      }
       )

}




################################################################################
#
# libgraphics::circle.cor
#   circles representation of correlation matrix
#
# TIME STAMP
#     2009-06-19
#
# PARAMETERS
#   cor:                a correlation matrix
#   tri:                TRUE (default) if plot upper triangular matrix
#   axes:               FALSE (default) whether axes should be plot
#   highlight.cutoff:   a value that determines at what magnitude a circle should be highlighted with red border
#
# RETURN VALUES
#   none
#
# DEPENDENCIES
#   library(kohonen)
#
# USAGE
#   fit = lm(mpg ~ ., mtcars)
#   cor = summary(fit, correlation = TRUE)$correlation
#   circle.cor(cor,tri=TRUE,highlight.cutoff=0.8)
################################################################################
## an example data(mtcars)

circle.cor = function(cor, tri=TRUE, axes = FALSE, xlab = "",
     ylab = "", asp = 1, title = "Correlation matrix circles",
     highlight.cutoff=0.9, font.cex.lab = 0.8, scale=FALSE, ...)
{
  par(lwd=2)

  if (tri) {  # only plot the upper triangular correlation matrix
     # author: Yongsheng Huang
     
     cor[lower.tri(cor)] = 0
     
     varnames = colnames(cor)
     n = nrow(cor)
     par(mar = c(0, 0, 2, 0), bg = "white")
     plot(c(0, n + 0.8), c(0, n + 0.8), axes = axes, xlab = "",
             ylab = "", asp = 1, type = "n")
     ##add grid
     # row
     segments( c(0.5 + (n-1):0, 0.5), 0.5 + 0:(n)
              , rep(n + 0.5, n), 0.5 + 0:(n)
              , col = "gray")

#     segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
#             0.5 + 0:n, col = "gray")
     segments(  0.5 + n:0, rep(n + 0.5, n + 1)
              , 0.5 + n:0, c(0.5, 0:(n-1) + 0.5)
              , col = "gray")

#     segments(0.5,n+0.5, n+0.5,0.5,col="gray")

#     fg = cor
     fg = matrix("black",nrow=nrow(cor),ncol=ncol(cor))
     fg[cor >= highlight.cutoff & cor!=1 ] = "red"
     fg[cor <= -highlight.cutoff & cor!=-1] = "red"

     ##define circles' background color.
     ##black for positive correlation coefficient and white for negative
     bg = cor
     bg[cor > 0] = "black"
     bg[cor <= 0] = "white"   ##plot n*n circles using vector language, suggested by Yihui Xie
     symbols(rep(1:n, each = n), rep(n:1, n), add = TRUE, inches = F
            , circles = as.vector(ifelse(sqrt(abs(cor))/2==0,NA,sqrt(abs(cor))/2))
            , bg = as.vector(bg), fg = as.vector(fg))

     # row names
     text(x=rep(n+1, n), y=n:1, varnames, col = "red", cex=font.cex.lab, font=2)
     # col names
     text(1:n, y=rep(n + 1), varnames, col = "red", cex=font.cex.lab, font=2)
     title(title)

  }
  else { # plot the complete correlation matrix
  
     # author: Taiyun Wei
     varnames = colnames(cor)
     n = nrow(cor)
     par(mar = c(0, 0, 2, 0), bg = "white")
     plot(c(0, n + 0.8), c(0, n + 0.8), axes = axes, xlab = "",
             ylab = "", asp = 1, type = "n")
     ##add grid
     segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
             0.5 + 0:n, col = "gray")
     segments(0.5 + 0:n, rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5,
                     n), col = "gray")

#     segments(0.5,n+0.5, n+0.5,0.5,col="gray")

#     fg = cor
     fg = matrix("black",nrow=nrow(cor),ncol=ncol(cor))
     fg[cor >= highlight.cutoff & cor!=1 ] = "red"
     fg[cor <= -highlight.cutoff & cor!=-1] = "red"

     ##define circles' background color.
     ##black for positive correlation coefficient and white for negative
     bg = cor
     bg[cor > 0] = "black"
     bg[cor <= 0] = "white"   ##plot n*n circles using vector language, suggested by Yihui Xie
     symbols(rep(1:n, each = n), rep(n:1, n), add = TRUE, inches = F
            , circles = as.vector(ifelse(sqrt(abs(cor))/2==0,NA,sqrt(abs(cor))/2))
            , bg = as.vector(bg), fg = as.vector(fg))

     # row names
     text(x=rep(0, n), y=n:1, varnames, col = "red", cex=0.5, font=2)
     # col names
     text(1:n, y=rep(n + 1), varnames, col = "red", cex=0.5, font=2)
     title(title)

#     if(scale==TRUE) {
#       circles = as.vector(ifelse(sqrt(abs(cor))/2==0,NA,cor))
#       cir.legend = seq(from=min(circles,na.rm=T)^2,to=max(circles,na.rm=T)^2,length.out=10)
#       cir.n = length(cir.legend)     
#       symbols(rep(1:cir.n), rep(1, cir.n), add = TRUE, inches = F
#               , circles = cir.legend
#               , bg = as.vector(bg), fg = as.vector(fg), lwd=2)
#     }
  }             
  
}


################################################################################
#
# graphics::my.plot.kohcodes
#   a modified version of plot.kohcodes{kohonen} for generating segments plot
#   when # of variables is > 15
#
# TIME STAMP
#     2009-05-14
#
# PARAMETERS
#   x:                  an kohonen SOM object
#   codeRendering:      "segments","stars",lines"
#   force.segment:      flag indicating whether segment plot be forced if n.var > 15
#
# RETURN VALUES
#   none
#
# DEPENDENCIES
#   library(kohonen)
#
# USAGE
#   my.plot.kohcodes(   test.som, main=NULL, bgcol=NULL, whatmap=NULL
#                     , codeRendering="segments", keepMargins=FALSE
#                     , force.segment=TRUE)
################################################################################

my.plot.kohcodes <- function (x, main, bgcol, whatmap, codeRendering, keepMargins, force.segment=FALSE, ...)
{
    if (!keepMargins) {
        opar <- par(c("mar", "ask"))
        on.exit(par(opar))
    }
    whatmap <- check.whatmap(x, whatmap)
    nmaps <- length(whatmap)
    if (is.list(x$codes)) {
        if (prod(par("mfrow")) < nmaps)
            par(ask = TRUE)
        for (i in 1:nmaps) {
            huhn <- x
            huhn$codes <- huhn$codes[[whatmap[i]]]
            if (length(main) == length(x$codes)) {
                main.title <- main[whatmap[i]]
            }
            else {
                if (length(main) == nmaps) {
                  main.title <- main[i]
                }
                else {
                  if (length(main) == 1) {
                    main.title <- main
                  }
                  else {
                    if (is.null(main)) {
                      if (!is.null(names(x$codes))) {
                        main.title <- names(x$codes)[whatmap[i]]
                      }
                      else {
                        main.title <- "Codes plot"
                      }
                    }
                  }
                }
            }
            if (length(codeRendering) == length(x$codes)) {
                cR <- codeRendering[whatmap[i]]
            }
            else {
                if (length(codeRendering) == nmaps) {
                  cR <- codeRendering[i]
                }
                else {
                  cR <- codeRendering
                }
            }
            plot.kohcodes(huhn, main = main.title, bgcol = bgcol,
                whatmap = NULL, codeRendering = cR, keepMargins = TRUE,
                ...)
        }
    }
    else {
        codes <- x$codes
        nvars <- ncol(codes)
        if (is.null(codeRendering)) {
            if (nvars < 15) {
                codeRendering <- "segments"
                maxlegendcols <- 3
            }
            else {
              codeRendering <- "lines"
            }
        }
        
        # yongsheng
        if(force.segment) {
           codeRendering <- "segments"
           maxlegendcols <- 6
        }
        
        margins <- rep(0.6, 4)
        if (!is.null(main))
            margins[3] <- margins[3] + 2
        par(mar = margins)
        # yongsheng changed
        if (codeRendering == "segments" & (nvars < 15 | force.segment) & !is.null(colnames(codes))) {
        #if (codeRendering == "segments" & nvars < 15 & !is.null(colnames(codes))) {
            plot(x$grid, ylim = c(max(x$grid$pts[, 2]) + min(x$grid$pts[,
                2]), -2))
            current.plot <- par("mfg")
            plot.width <- diff(par("usr")[1:2])
            cex <- 1
            leg.result <- legend(x = mean(x$grid$pts[, 1]), xjust = 0.5,
                y = 0, yjust = 1, legend = colnames(codes), cex = cex,
                plot = FALSE, ncol = min(maxlegendcols, nvars),
                fill = rainbow(nvars))
            while (leg.result$rect$w > plot.width) {
                cex <- cex * 0.9
                leg.result <- legend(x = mean(x$grid$pts[, 1]),
                  xjust = 0.5, y = 0, yjust = 1, legend = colnames(codes),
                  cex = cex, plot = FALSE, ncol = min(maxlegendcols,
                    nvars), fill = rainbow(nvars))
            }
            leg.result <- legend(x = mean(x$grid$pts[, 1]), xjust = 0.5,
                y = 0, yjust = 1, cex = cex, legend = colnames(codes),
                plot = FALSE, ncol = min(maxlegendcols, nvars),
                fill = rainbow(nvars), ...)
            par(mfg = current.plot)
            plot(x$grid, ylim = c(max(x$grid$pts[, 2]) + min(x$grid$pts[,
                2]), -leg.result$rect$h))
            legend(x = mean(x$grid$pts[, 1]), xjust = 0.5, y = 0,
                yjust = 1, cex = cex, plot = TRUE, legend = colnames(codes),
                ncol = min(maxlegendcols, nvars), fill = rainbow(nvars),
                ...)
        }
        else {
            plot(x$grid, ...)
        }
        title.y <- max(x$grid$pts[, 2]) + 1.2
        if (title.y > par("usr")[4] - 0.2) {
            title(main)
        }
        else {
            text(mean(range(x$grid$pts[, 1])), title.y, main,
                adj = 0.5, cex = par("cex.main"), font = par("font.main"))
        }
        if (is.null(bgcol))
            bgcol <- "transparent"
        symbols(x$grid$pts[, 1], x$grid$pts[, 2], circles = rep(0.5,
            nrow(x$grid$pts)), inches = FALSE, add = TRUE, bg = bgcol)
        if (codeRendering == "lines") {
            yrange <- range(codes)
            codes <- codes - mean(yrange)
        }
        else {
            codemins <- apply(codes, 2, min)
            codes <- sweep(codes, 2, codemins)
        }
        switch(codeRendering, segments = {
            stars(codes, location = x$grid$pts, labels = NULL,
                len = 0.4, add = TRUE, col.segments = rainbow(nvars),
                draw.segments = TRUE)
        }, lines = {
            for (i in 1:nrow(x$grid$pts)) {
                if (yrange[1] < 0 & yrange[2] > 0) {
                  lines(seq(x$grid$pts[i, 1] - 0.4, x$grid$pts[i,
                    1] + 0.4, length = 2), rep(x$grid$pts[i,
                    2], 2), col = "gray")
                }
                lines(seq(x$grid$pts[i, 1] - 0.4, x$grid$pts[i,
                  1] + 0.4, length = ncol(codes)), x$grid$pts[i,
                  2] + codes[i, ] * 0.8/diff(yrange), col = "red")
            }
        }, stars = stars(codes, location = x$grid$pts, labels = NULL,
            len = 0.4, add = TRUE))
    }
    invisible()
}


cat(">>>>>>> Graphics library has been successfully loaded!\n")