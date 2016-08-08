########################################################################
# phdgraphics.r
#   graphics utitlity functions and modules specific for phd project
########################################################################



#***********************************************************************************************
# phdgraphics::ggPair
#   plot pairs of data
#
# AUTHOR
#   Yongsheng Huang
#   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PAMETERS
#   ds      -   a data frame
#
# VALUE
#   None
#
# DEPENDENCIES
#   1) library(ggplot2)
#
# DBG
#   setting parmeters for running the function
#
# NOTES
#
#
# USAGE
#   exemplary usage of the function
#
# KNOWN ISSUES
#
#
# TIME STAMP
#   2009-08-13
#***********************************************************************************************
#list.what = loadLW("TNF.H00vsH05.csv")
#test = getDesign(list.what=list.what, list.data=dgep)
#see(test)
#genes = unique(rownames(test))
#n.genes = length(genes)
#
#pdf ("test.ggPair.pdf")
#for (i in 1:n.genes) {
#  a.gene.name = getFeatures(genes[i], annot='Symbol', objannot=dgep$features)
#  a.gene = data.frame(test[genes[i],])
#  a.gene = data.frame(t(a.gene))
#  a.gene$obs = getPheno(rownames(a.gene), annot="Subject", objannot=dgep$obs)
#  a.gene$pheno = getPheno(rownames(a.gene), annot="PhenoNew", objannot=dgep$obs)
#  a.gene$timepoints = getPheno(rownames(a.gene), annot="TimePoints", objannot=dgep$obs)
#  colnames(a.gene)[1] = a.gene.name
#
#  a.gene = reshape ( a.gene
#                        ,v.name=a.gene.name
#                        ,idvar=c("obs")
#                        ,timevar=c("timepoints")
#                        ,direction="wide"
#                      )
#
#  tmp = ggPair(ds=a.gene, obs="obs", group="pheno",vars=c(3,4), main.label=a.gene.name)
#  rm(tmp)
#}
#dev.off()
#
ggPair <- function(ds=NULL, obs="obs", group="pheno", vars=c(3,4), cex.label=4, main.label="") {

  if(is.null(ds)) stop("\n give me a data set to plot, you stupid ... \n")

  suppressMessages(require(ggplot2))
  theme_set(theme_bw())

  ds$Diff = abs(ds[,vars[1]] - ds[,vars[2]])
  ds$Group = as.factor(ds[,group])
  ds$Obs = ds[,obs]
  
  nrange = c(min(ds[,vars]), max(ds[,vars]))

  a = qplot(  x=ds[,vars[1]],y=ds[,vars[2]],data=ds
            ,colour=Group
            ,xlab=colnames(ds)[vars[1]], ylab=colnames(ds)[vars[2]]
            ,xlim=nrange, ylim=nrange
            ,size=Diff, shape=I(19)
            ,geom=c("point")
            ,main=main.label
           )
  a = a + geom_text(aes(label=Obs, size=Diff*cex.label),position="dodge")
  a = a + geom_abline(intercept=0,colour="gray60",size=1,linetype="longdash")

  return(print(a))
  
}






#***********************************************************************************************
# phdgraphics::gplotAS
#   plot temporal signature of a feature (gene/protein) by phenotype
#
# AUTHOR
#   Yongsheng Huang   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PARAMETERS
#   features:       names of features (genes/proteins)
#   obs:            names of observations
#   list.data:      list of data/feature_annotation/obs_annotation to be searched
#   by.what:        by what variable the graphics panel should be organized
#
# VALUE
#   ret     -   return value of this function
#
# DEPENDENCIES
#   1) there are three elements in list.data: data, features, obs
#   2) time variable will always be plotted on the x-axis
#
# DBG
#  dgep = loadData("gep")
#  dprtn = loadData("protein")
#  list.what = loadLW("test.csv")
#    #> list.what = loadLW("test.csv")
#    #List of 2
#    # $ features: chr [1:2] "SIGIRR" "IL18"
#    # $ obs     : chr [1:251] "Z01H005" "Z01H00" "Z01H012" "Z01H021" ...
#  by.what = list("PhenoNew","TimePoints")
#  data.tmp = gplotAS (list.what=list.what, list.data=dgep , by.what=by.what, pdffile="test.pdf")
#  gplotAS (list.what=list.what, list.data=dprtn , by.what=by.what, pdffile="2009.05.26.Prtn.gplotAS.pdf")
#
# NOTE
#  currently only plot by "PhenoNew" and "TimePoints"
#
# USAGE
#   dgep = loadData("gep")
#   fname.list.what = "flu.list.what.asxEDGE.csv"
#   fname.pdf = "2009.09.09.flu.key.genes.Asx.unique.plot.pdf"
#   by.what = list("PhenoNew","TimePoints")
#   list.what = loadLW(fname.list.what)
#   key.genes = gplotAS (list.what=list.what, list.data=dgep , by.what=by.what, pdffile=fname.pdf)
#
#   to print into a jpeg file
#jpeg()
#key.genes = gplotAS (list.what=list.what, list.data=dgep , by.what=by.what, pdffile=fname.pdf)
#dev.off()
#
# KNOWN ISSUES
#
# TIME STAMP
#   2009-09-03
#***********************************************************************************************

gplotAS <- function (list.what=NULL, list.data=NULL, by.what=NULL, try.symbol=TRUE, panel.size.x=6.5,panel.size.y=5, pdffile=NULL) {

  library(ggplot2)

  theme_set(theme_bw())

  # fix this for now (2009-05-25)
  # by.what = list("PhenoNew","TimePoints")


  n.features = length(list.what$features)
  n.obs = length(list.what$obs)

  tmp = getDesign(list.what,list.data=list.data, try.symbol=try.symbol)
  # we would like to have gene symbol printed, instead of probeset names
  if(length(grep("_at$",rownames(tmp)))>nrow(tmp)/2) rownames(tmp) = getFeatures(idfeatures=rownames(tmp),annot='Symbol',objannot=dgep$features)

  # very strict - make sure all features / samples can be mapped
  if(is.null(tmp)) stop("no data was obtained")

  tmp = dmApnd(t(tmp), list.data$obs, what=by.what)
  tmp = sort.df(tmp,by=~PhenoNew+TimePoints)
  # tmp$PT = paste(as.character(tmp$PhenoNew), right(paste('000',gsub("H","",gsub("H0","",as.character(tmp$TimePoints))),sep=''),3),sep="")
  tmp$PT = labelfix(prefix=as.character(tmp$PhenoNew),avec=as.character(tmp$TimePoints))

  if(!is.null(pdffile)) {pdf(pdffile,width=panel.size.x+0.1,height=panel.size.y+0.1)}
  for (i in 1:n.features) {
    print(paste("n=",i," ","now plotting ",list.what$features[i],sep=""))
#    y.lab = list.what$features[i]
    tmp.name.keep = colnames(tmp)[i]
    y.lab = tmp.name.keep
    colnames(tmp)[i] = "value"
    tmp.gplot = qplot(value,x=PT,colour=PhenoNew,data=tmp, geom="blank", ylab=y.lab, xlab="Phenotype + Time")
#    tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_brewer(type="seq", palette=3)
    #tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_manual(value=c("darkgreen","red"))
    tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_manual(value=c("#0099FF","red"))
    # mean_cl_boot: non-parameter bootstrap (Efron) for obtaining 95% CI of population mean without assumption of normality
#    tmp.gplot = tmp.gplot + stat_sum_dfmX("mean_cl_boot",mult=1,geom="crossbar",colour=I("red"),width=0.5)
    tmp.gplot = tmp.gplot + stat_sum_dfmX("mean_cl_boot",mult=1,geom="crossbar",width=0.5)

    # 2009-09-21
    tmp.gplot = tmp.gplot + opts(panel.grid.major=theme_line(size=0.4, linetype="dotted"))
    tmp.gplot = tmp.gplot + opts(   #axis.line=theme_segment(size=1)
                                   panel.border=theme_rect(size=0.7)
                                  , axis.text.y=theme_text(face="bold", size=10, hjust=1)
                                  , axis.text.y=theme_blank()
                                  , axis.text.x=theme_text(face="bold", size=9, hjust=1, vjust=1, angle=45)
                                  , axis.ticks=theme_segment(size=0.6)
                                  , axis.title.x=theme_text(face="bold",size=12)
                                  , axis.title.y=theme_text(face="bold",size=12,angle=90)
                                  , axis.title.y=theme_blank()
                                )
                                
    grid.newpage(recording=TRUE)
    tmp.viewport = viewport(width=unit(panel.size.x,"inches"),height=unit(panel.size.y,"inches"))
    print(tmp.gplot, vp=tmp.viewport)
    # set the column name back to its original name
    colnames(tmp)[i] = tmp.name.keep
  }
  if(!is.null(pdffile)) {dev.off()}

  return (tmp)
  
}


gplotAS_SparseX <- function (list.what=NULL, list.data=NULL, by.what=NULL, try.symbol=TRUE, panel.size.x=6.5,panel.size.y=5, pdffile=NULL) {
# alternative version of gplotAS for STM paper ( sparse X axis with bigger font)
# 2010/08/14

  library(ggplot2)

  theme_set(theme_bw())

  # fix this for now (2009-05-25)
  # by.what = list("PhenoNew","TimePoints")


  n.features = length(list.what$features)
  n.obs = length(list.what$obs)

  tmp = getDesign(list.what,list.data=list.data, try.symbol=try.symbol)
  # we would like to have gene symbol printed, instead of probeset names
  if(length(grep("_at$",rownames(tmp)))>nrow(tmp)/2) rownames(tmp) = getFeatures(idfeatures=rownames(tmp),annot='Symbol',objannot=dgep$features)

  # very strict - make sure all features / samples can be mapped
  if(is.null(tmp)) stop("no data was obtained")

  tmp = dmApnd(t(tmp), list.data$obs, what=by.what)
  tmp = sort.df(tmp,by=~PhenoNew+TimePoints)
  # tmp$PT = paste(as.character(tmp$PhenoNew), right(paste('000',gsub("H","",gsub("H0","",as.character(tmp$TimePoints))),sep=''),3),sep="")
  tmp$PT = labelfix(prefix=as.character(tmp$PhenoNew),avec=as.character(tmp$TimePoints))

  if(!is.null(pdffile)) {pdf(pdffile,width=panel.size.x+0.1,height=panel.size.y+0.1)}
  for (i in 1:n.features) {
    print(paste("n=",i," ","now plotting ",list.what$features[i],sep=""))
#    y.lab = list.what$features[i]
    tmp.name.keep = colnames(tmp)[i]
    y.lab = tmp.name.keep
    colnames(tmp)[i] = "value"
    tmp.gplot = qplot(value,x=PT,colour=PhenoNew,data=tmp, geom="blank", ylab=y.lab, xlab="Phenotype + Time")
#    tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_brewer(type="seq", palette=3)
    #tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_manual(value=c("darkgreen","red"))
    tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_manual(value=c("#0099FF","red"))
    # mean_cl_boot: non-parameter bootstrap (Efron) for obtaining 95% CI of population mean without assumption of normality
#    tmp.gplot = tmp.gplot + stat_sum_dfmX("mean_cl_boot",mult=1,geom="crossbar",colour=I("red"),width=0.5)
    tmp.gplot = tmp.gplot + stat_sum_dfmX("mean_cl_boot",mult=1,geom="crossbar",width=0.5)

    # 2009-09-21
    tmp.gplot = tmp.gplot + opts(panel.grid.major=theme_line(size=0.5, linetype="dotted"))
    tmp.gplot = tmp.gplot + opts(   #axis.line=theme_segment(size=1)
                                   panel.border=theme_rect(size=1.2)
                                  , axis.text.y=theme_text(face="bold", size=16, hjust=0.8)
                                  #, axis.text.y=theme_blank()
                                  , axis.text.x=theme_text(face="bold", size=16, vjust=0.6)
                                  , axis.ticks=theme_segment(size=1)
                                  , axis.title.x=theme_text(face="bold",size=12)
                                  , axis.title.y=theme_text(face="bold",size=12,angle=90,hjust=0.5)
                                  #, axis.title.y=theme_blank()
                                )
    tmp.gplot = tmp.gplot  + scale_x_discrete(breaks=c('-1', 'A012','A036','A060','A084','A108', 'S012','S036','S060','S084','S108'), labels=c(-12,   12,36,60,84,108,   12,36,60,84,108))
    grid.newpage(recording=TRUE)
    tmp.viewport = viewport(width=unit(panel.size.x,"inches"),height=unit(panel.size.y,"inches"))
    print(tmp.gplot, vp=tmp.viewport)
    # set the column name back to its original name
    colnames(tmp)[i] = tmp.name.keep
  }
  if(!is.null(pdffile)) {dev.off()}

  return (tmp)
  
}


#
# 2010-08-23
# plot each data points with different symbols for different subjects
#
# alternative version of gplotAS for STM paper ( sparse X axis with bigger font)
# 
#
gplotAS2_SparseX <- function (list.what=NULL, list.data=NULL, by.what=NULL, try.symbol=TRUE, panel.size.x=6.5,panel.size.y=5, pdffile=NULL) {

  library(ggplot2)

  theme_set(theme_bw())

  # fix this for now (2009-05-25)
  # by.what = list("PhenoNew","TimePoints")

  n.features = length(list.what$features)
  n.obs = length(list.what$obs)

  tmp = getDesign(list.what,list.data=list.data, try.symbol=try.symbol)
  # we would like to have gene symbol printed, instead of probeset names
  if(length(grep("_at$",rownames(tmp)))>nrow(tmp)/2) rownames(tmp) = getFeatures(idfeatures=rownames(tmp),annot='Symbol',objannot=dgep$features)

  # very strict - make sure all features / samples can be mapped
  if(is.null(tmp)) stop("no data was obtained")

  tmp = dmApnd(t(tmp), list.data$obs, what=by.what)
  tmp = sort.df(tmp,by=~PhenoNew+TimePoints)
  # tmp$PT = paste(as.character(tmp$PhenoNew), right(paste('000',gsub("H","",gsub("H0","",as.character(tmp$TimePoints))),sep=''),3),sep="")
  tmp$PT = labelfix(prefix=as.character(tmp$PhenoNew),avec=as.character(tmp$TimePoints))
  tmp$Subject = as.factor(as.integer(gsub("Z","",left(rownames(tmp),3))))
  n.subject = length(unique(tmp$Subject))


  if(!is.null(pdffile)) {pdf(pdffile,width=panel.size.x+0.1,height=panel.size.y+0.1)}
  for (i in 1:n.features) {
    cat(paste("n=",i," ","now plotting ",list.what$features[i],sep=""),"\n")
#    y.lab = list.what$features[i]
    tmp.name.keep = colnames(tmp)[i]
    y.lab = tmp.name.keep
    colnames(tmp)[i] = "value"
    tmp.gplot = qplot(value,x=PT,colour=PhenoNew,data=tmp, geom="blank", ylab=y.lab, xlab="Phenotype + Time")
#    tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_brewer(type="seq", palette=3)
    tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_manual(value=c("#0099ff","#ff0000"))
    # mean_cl_boot: non-parameter bootstrap (Efron) for obtaining 95% CI of population mean without assumption of normality
#    tmp.gplot = tmp.gplot + stat_sum_dfmX("mean_cl_boot",mult=1,geom="crossbar",colour=I("red"),width=0.5)
    #tmp.gplot = tmp.gplot + stat_sum_dfmX("mean_cl_boot",mult=1,geom="crossbar",width=0.5)
    tmp.gplot = tmp.gplot + aes(shape=tmp$Subject) + scale_shape_manual(values=seq(from=0,to=n.subject-1))

    # 2009-09-21
    tmp.gplot = tmp.gplot + opts(panel.grid.major=theme_line(size=0.4, linetype="dotted"))
    tmp.gplot = tmp.gplot + opts(   #axis.line=theme_segment(size=1)
                                   panel.border=theme_rect(size=1.2)
                                  , axis.text.y=theme_text(face="bold", size=16, hjust=0.8)
                                  #, axis.text.y=theme_blank()
                                  , axis.text.x=theme_text(face="bold", size=16, vjust=0.6)
                                  , axis.ticks=theme_segment(size=1)
                                  , axis.title.x=theme_text(face="bold",size=12)
                                  , axis.title.y=theme_text(face="bold",size=12,angle=90,hjust=0.5)
                                  #, axis.title.y=theme_blank()
                                )
    tmp.gplot = tmp.gplot  + scale_x_discrete(breaks=c('-1', 'A012','A036','A060','A084','A108', 'S012','S036','S060','S084','S108'), labels=c(-12,   12,36,60,84,108,   12,36,60,84,108))
                    
    grid.newpage(recording=TRUE)
    tmp.viewport = viewport(width=unit(panel.size.x,"inches"),height=unit(panel.size.y,"inches"))
    print(tmp.gplot, vp=tmp.viewport)
    # set the column name back to its original name
    colnames(tmp)[i] = tmp.name.keep
  }
  if(!is.null(pdffile)) {dev.off()}

  return (tmp)

}


#
# 2010-03-28
# plot each data points with different symbols for different subjects
# 
#
gplotAS2 <- function (list.what=NULL, list.data=NULL, by.what=NULL, try.symbol=TRUE, panel.size.x=6.5,panel.size.y=5, pdffile=NULL) {

  library(ggplot2)

  theme_set(theme_bw())

  # fix this for now (2009-05-25)
  # by.what = list("PhenoNew","TimePoints")

  n.features = length(list.what$features)
  n.obs = length(list.what$obs)

  tmp = getDesign(list.what,list.data=list.data, try.symbol=try.symbol)
  # we would like to have gene symbol printed, instead of probeset names
  if(length(grep("_at$",rownames(tmp)))>nrow(tmp)/2) rownames(tmp) = getFeatures(idfeatures=rownames(tmp),annot='Symbol',objannot=dgep$features)

  # very strict - make sure all features / samples can be mapped
  if(is.null(tmp)) stop("no data was obtained")

  tmp = dmApnd(t(tmp), list.data$obs, what=by.what)
  tmp = sort.df(tmp,by=~PhenoNew+TimePoints)
  # tmp$PT = paste(as.character(tmp$PhenoNew), right(paste('000',gsub("H","",gsub("H0","",as.character(tmp$TimePoints))),sep=''),3),sep="")
  tmp$PT = labelfix(prefix=as.character(tmp$PhenoNew),avec=as.character(tmp$TimePoints))
  tmp$Subject = as.factor(as.integer(gsub("Z","",left(rownames(tmp),3))))
  n.subject = length(unique(tmp$Subject))


  if(!is.null(pdffile)) {pdf(pdffile,width=panel.size.x+0.1,height=panel.size.y+0.1)}
  for (i in 1:n.features) {
    cat(paste("n=",i," ","now plotting ",list.what$features[i],sep=""),"\n")
#    y.lab = list.what$features[i]
    tmp.name.keep = colnames(tmp)[i]
    y.lab = tmp.name.keep
    colnames(tmp)[i] = "value"
    tmp.gplot = qplot(value,x=PT,colour=PhenoNew,data=tmp, geom="blank", ylab=y.lab, xlab="Phenotype + Time")
#    tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_brewer(type="seq", palette=3)
    tmp.gplot = tmp.gplot + geom_jitter(position=position_jitter(width=0.2))  + scale_colour_manual(value=c("#0099ff","#ff0000"))
    # mean_cl_boot: non-parameter bootstrap (Efron) for obtaining 95% CI of population mean without assumption of normality
#    tmp.gplot = tmp.gplot + stat_sum_dfmX("mean_cl_boot",mult=1,geom="crossbar",colour=I("red"),width=0.5)
    #tmp.gplot = tmp.gplot + stat_sum_dfmX("mean_cl_boot",mult=1,geom="crossbar",width=0.5)
    tmp.gplot = tmp.gplot + aes(shape=tmp$Subject) + scale_shape_manual(values=seq(from=0,to=n.subject-1))

    # 2009-09-21
    tmp.gplot = tmp.gplot + opts(panel.grid.major=theme_line(size=0.4, linetype="dotted"))
    tmp.gplot = tmp.gplot + opts(   #axis.line=theme_segment(size=1)
                                   panel.border=theme_rect(size=0.7)
                                  , axis.text.y=theme_text(face="bold", size=10, hjust=1)
                                  , axis.text.y=theme_blank()
                                  , axis.text.x=theme_text(face="bold", size=9, hjust=1, vjust=1, angle=45)
                                  , axis.ticks=theme_segment(size=0.6)
                                  , axis.title.x=theme_text(face="bold",size=12)
                                  , axis.title.y=theme_text(face="bold",size=12,angle=90)
                                  , axis.title.y=theme_blank()
                                )
                    
    grid.newpage(recording=TRUE)
    tmp.viewport = viewport(width=unit(panel.size.x,"inches"),height=unit(panel.size.y,"inches"))
    print(tmp.gplot, vp=tmp.viewport)
    # set the column name back to its original name
    colnames(tmp)[i] = tmp.name.keep
  }
  if(!is.null(pdffile)) {dev.off()}

  return (tmp)

}




#
# T vs Baseline waterplot
#   2D cartesian plot of Baseline vs T time
#
# usage
#
#> str(test)
#'data.frame':   57 obs. of  4 variables:
# $ phenonew: Factor w/ 2 levels "A","S": 1 1 1 1 1 1 1 1 1 1 ...
# $ PC      : num  8.08 8 8.09 8.08 7.12 7.9 8.1 8.09 8.09 5.86 ...
# $ T       : num  8.04 7.84 8.14 7.9 8.01 7.97 8.09 8.13 7.99 7.87 ...
# $ Diff    : num  0.04 0.16 0.05 0.18 0.89 ...

plotBT2D <- function(x){
  test=read.csv("C:\\Users\\yongsheng\\Documents\\Yongsheng's Documents\\UM\\research\\hero\\phd\\analysis\\2009.04.18.EBC.pH\\2009.05.01.EBC.pH.plot.XY.csv",header=T,row.names=1)
  #test$diff = abs(test$T - test$PC)
  #test$T-PC = abs(test$T - test$PC)
  test$Diff = abs(test$T - test$PC)
  #test$phenonew = test$Pheno
  a = qplot(x=PC,y=T,data=test,colour=phenonew,geom="point",xlim=c(4.5,9),ylim=c(4.5,9),size=Diff,shape=I(19))
  a = a + geom_abline(intercept=0,colour="black",size=1,linetype="longdash")
  # aes does not work in the following code, need to understand the ggplot better
  # a = a + geom_abline(aes(intercept=0,colour=I("black"),linetype="twodash"))
  a = a + scale_colour_brewer(pal="Set1")
  a
}

cat(">>>>>>> PHD project graphics library has been successfully loaded!\n")