#############################################################
# libphd.R
#   utitlity functions and modules specific for PHD project
############################################################

#par(font=2,lty=2,lwd=2,cex=2, mai=rep(0,4))
#venn(test$Genes,test$Challenge,main=NULL)
#

#
# make sure generic libraries is loaded
#
chkLib <- function (slib='library.r') {
  flib = slib
  # load yongsheng's own library
  if(slib=='library.r') {
      # specific location on windows machine
      flibwin = "./libraries/library.r"
      # has it been loaded? - check dummie function
       tryCatch( { if(class(flag_ys_library)=="function"){return();print("library exists in environment!")}
                 }
                 ,error = function(ex) {if (is.element(flib,dir())) source(flib)
                                      else source(flibwin)
                                        # 'ex' is an instance of 'simpleError'
                                        #  print(as.character(ex))
                                        #  print(ex$message)
                                       }
                 ,handler=function(){ if (is.element(flib,dir())) source(flib)
                                      else source(flibwin)
                                    }
       )   # end of tryCatch
   }
   
   # load R libraries
   else {
       tryCatch( { if(class(flag_ys_library)==flib){print("library exists in environment!")}
                 }
                 ,error = function(ex) {library(flib)}
                 ,handler=function(){ library(flib)}
       )   # end of tryCatch
   }
}


#
# exec the function to load all necessary libraries
#
chkLib()

#***********************************************************************************************
# phdlib::uni.lme
#   fitting a univariate linear mixed-effects model and extracts either sqrt(R^2) or p-value
#
# AUTHOR
#   Yongsheng Huang   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PAMETERS
#    ds           -   data frame with data for model fitting
#    what.y       -   name of the variable in "ds" to be used as response or y
#    what.x       -   name of the variable in "ds" to be used as response or x
#    what.cluster -   name of the variable in "ds" that indicates cluster/group structure
#    what.ret     -   "p-value" or "cor"; indicates which type of statistics should be returned
#
# VALUE
#   ret     -   return value of this function
#
# DEPENDENCIES
#   1) library {nlme}
#
# DBG
#    ds = test
#    what.y = "symp"
#    what.x = "gep"
#    what.cluster = "subject"
#    what.ret = "p-value"
#    uni.lme(ds=ds, what.y=what.y, what.x=what.x, what.cluster=what.cluster,what.ret=what.ret)
#
#
# NOTES
#
# 1)adding the interaction term will increase R^2
# but it will no longer be univariate linear model
# the R^2 can not be directly interpretated as coefficients
#
# ds.lme = lme(y~x*factor(time), random=~1|subject/time, data=x) # same results > find out "/time" thing
#
# 2) # x$fitted = predict(test.lme, x)
#
#
# USAGE
#   exemplary usage of the function
#
# KNOWN ISSUES
#
#
# TIME STAMP
#   2009-09-30
#***********************************************************************************************

uni.lme <- function (ds=NULL, what.y=NULL, what.x=NULL, what.cluster=NULL, what.ret=NULL) {

  require(nlme)

  idx.y = grep(what.y, colnames(ds))
  idx.x = grep(what.x, colnames(ds))
  idx.c = grep(what.cluster, colnames(ds))
  colnames(ds)[idx.y] = "y"
  colnames(ds)[idx.x] = "x"
  colnames(ds)[idx.c] = "cluster"

  ds.lme = lme(y~x, random=~1|cluster, data=ds)
  # append the fitted value of fixed effects from what.x
  ds$fitted = ds.lme$fitted[,"fixed"]

  # get the coefficients
  ds.coef = ds.lme$coef$fixed[2]

  res.fitted = sum((ds$y - ds$fitted)^2)
  res.expect = sum((ds$y - mean(ds$y))^2)
  res.fitted = ifelse(res.fitted>res.expect, res.expect, res.fitted)
  ds.cor = sqrt(1-res.fitted/res.expect)
  if(ds.coef<0) ds.cor = -ds.cor

  ds.pval = anova(ds.lme)$`p-value`[2]
  ds.pval = ifelse(ds.pval==0, "<0.0001", ds.pval)

  if(what.ret=="cor") ret = ds.cor
  if(what.ret=="p-value") ret = ds.pval

  return(ret)

}




#===================================================================================================
#  phdlib::featurecls
#    function for creating an sub.list.what, containing features based on class labels and obs
#
#  PARAMETERS
#    list.what         an object with features/obs/class labels/ranks
#    n.feature         number of features to be selected in each class
#    rank.name         name of the element in the list.what contains rank information of features
#    cls.name          name of the element in the list.what for class labels
#
#  DEPENDENCIES
#  1) obs        there has to be an element with name "obs" in the list.data
#  2) features   there has to be an element with name "features" in the list.data
#
#  USAGE
#  list.what = loadLW("flu.list.what.csv")
#  tmp = featurecls(list.what=list.what, n.feature=10, rank.name="edgerank",cls.name="classSOM")
#===================================================================================================

labelfix <- function (avec=NULL, prefix=NULL, leading="H", n0first=1, is.baseline=TRUE) {

  #k = c("BL", "H72", "H96", "H108", "H00")
  #tmp = k
  #leading = "H"
  #is.baseline = TRUE

  tmp = avec
  ix.nonbl = which(tmp!="BL")

  if (is.baseline) tmp[tmp=="BL"] = -1
  tmp = as.numeric(gsub(leading, "",tmp))
  n.pudding = max(apply(matrix(tmp),2,nchar))
  pudding = paste(rep("0",n.pudding),sep="",collapse="")
  tmp[ix.nonbl] = right(
                        paste( pudding, tmp[ix.nonbl], sep="")
                       , n.pudding)
   # make time 0 as 0, instead of "000"
   tmp[tmp==pudding] = paste(rep("0",n0first),sep="",collapse="")
   
   # add prefix (such as phenotype A or S)
   # if no prefix given, just prefix the labels with the leading character
   if(is.null(prefix)) tmp[ix.nonbl] = paste(leading, tmp[ix.nonbl], sep="")
   # if only one prefix (e.g., "A") given, append the same prefix  many times
   else if(length(prefix)==1) prefix = rep(prefix, length(tmp))
   # if a vector of prefix given, check the length compatbile with lables and then prefix it upfront
   else if(length(prefix)==length(tmp)) tmp[ix.nonbl] = paste(prefix[ix.nonbl], tmp[ix.nonbl], sep="")
   else stop("check your prefix")

   return (tmp)
   
}

#list.what = loadLW("2009.08.04.key.genes.csv")
#key.genes = gplotAS (list.what=list.what, list.data=dgep , by.what=by.what, pdffile="2009.08.04.key.genes.plot.pdf")
#

#===================================================================================================
#  phdlib::featurecls
#    function for creating an sub.list.what, containing features based on class labels and obs
#
#  PARAMETERS
#    list.what         an object with features/obs/class labels/ranks
#    n.feature         number of features to be selected in each class
#    rank.name         name of the element in the list.what contains rank information of features
#    cls.name          name of the element in the list.what for class labels
#
#  DEPENDENCIES
#  1) obs        there has to be an element with name "obs" in the list.data
#  2) features   there has to be an element with name "features" in the list.data
#
#  USAGE
#  list.what = loadLW("flu.list.what.csv")
#  tmp = featurecls(list.what=list.what, n.feature=10, rank.name="edgerank",cls.name="classSOM")
#===================================================================================================
featurecls <- function(list.what=NULL, n.feature=100, rank.name="edgerank",cls.name="classSOM")
{ # phdlib::featurecls

    #rank.name = "edgerank"
    #cls.name = "classSOM"
    #n.feature=100
    table(list.what[[cls.name]])
    somcls = unique(list.what[[cls.name]])
    somcls = sort(somcls)

    all.features = NULL
    for (i in 1:length(somcls)) {
      tmp = which(list.what[[cls.name]]==somcls[i])
      a.som = data.frame(   features=list.what$features[tmp]
                          , rank=list.what[[rank.name]][tmp]
                          , class=list.what[[cls.name]][tmp]
                          , stringsAsFactors=FALSE)
      a.som = sort.df(a.som, by=~+rank)
      # be very careful here with local variable chaning global variable
      if (nrow(a.som)<n.feature) n.featurecls=nrow(a.som)
      else n.featurecls = n.feature
      tmp = a.som$features[1:n.featurecls]
      msg = paste(" of genes selected from class ", somcls[i], " = ", length(tmp), sep="")
      print(msg)
      all.features = c(all.features,tmp)
    }

    return ( list(features=all.features, obs=list.what$obs) )
}


loadLW <- function(ifile=NULL,obs=NULL,features=NULL) {

  #
  # phdlib::loadLW
  #   utility function to load 'list.what' (to be used with phdlib::getDesign
  #
  # DEPENDENCIES
  #   1) input file should look like the following
  #        obs	        features
  #        Z01H005	    SIGIRR
  #        Z01H00	      IL18
  #        Z01H012
  #
  # RETURN VALUE
  #   list.what   -   list; two elements: features (genes/proteins); obs (observations/samples)
  #
  # DBG
  # list.what = loadLW(ifile="test.csv")
  #
  
  if(is.null(ifile) & !is.null(obs) & !is.null(features)) {
    list.what=list()
    list.what[["obs"]] = obs
    list.what[["features"]] = features
    return(list.what)
  }
  
  test = read.csv(ifile,header=T,stringsAsFactors=F)
  items = colnames(test)

  if (! ("features" %in% items)) stop("we need a column named 'features'")
  if (! ("obs" %in% items)) stop("we need a column named 'obs'")

  list.what = list()
  for (i in 1:length(items)) {
    list.what[[i]] = test[test[,items[i]]!="",items[i]]
    names(list.what)[i]=items[i]
  }

  str(list.what)

  return(list.what)

}


getPair <- function(what=NULL,by.what=NULL,src.annot=NULL,tar.annot=NULL,time.gap=1,fwd=FALSE,bkwd=FALSE) {
  #
  # phdlib::getPair
  #   given a list of samples, find their matching ones in another set, allow time gap
  #
  # PARAMETERS
  #
  #
  # NOTE
  #   a variable Time has to be present in by.what list
  #
  # DBG
  #    obs = dprtn$obs[dprtn$obs$Assay=='A' & dprtn$obs$Challenge=='Z' & dprtn$obs$TimePoints %in% c("BL","H00","H012","H021","H093","H108"),]
  #    what = rownames(obs)
  #    by.what = c("Challenge", "TimePoints", "Subject", "PhenoNew", "Time")
  #    src.annot = dprtn$obs
  #    tar.annot = dgep$obs
  #    time.gap = 0
  #    fwd = F
  #    bkwd = F
  #    x = getPair(what,by.what,src.annot,tar.annot,time.gap=1, fwd=TRUE, bkwd=T)
  #
  if(is.null(what) | is.null(by.what) | is.null(src.annot) | is.null(tar.annot))
      stop("don't you know how to supply the right arguments, stupid ? ")

  cat("Caution >> This function is still primitive. Error may occur. Do check the return values!!! \n")

  ret = NULL
  step = 0.5
  
  src = src.annot[what,by.what]
  src$src = what
  tar = tar.annot[,by.what]
  tar$tar = rownames(tar.annot)
#  if(time.gap==0)
#    merge(src, tar, by = by.what
#        ,all.x=T
#        ,sort = TRUE, suffixes = c(".src",".tar")
#       )

  if (length(what)!=dim(src)[1]) stop("check your input list ... stupid!")

  for (i in 1:length(what)){
    step.cnt = 0
    tmp = src[what[i],]
    # first look for perfect match
    tmp.match = merge(tmp, tar, by=by.what, all.x=T, sort=TRUE,suffixes=c(".src",".tar"))
    # if no perfect match, step through a few time points for closest match
    while(any(is.na(tmp.match))==T & (step.cnt*step<=time.gap)) {

      step.cnt = step.cnt+1
      if (fwd & bkwd) {
        # if search both directions, we do forward first and see if a match can be found
        tmp$Time=tmp$Time + step*step.cnt
        tmp.match = merge(tmp, tar, by=by.what, all.x=T, sort=TRUE,suffixes=c(".src",".tar"))
        if (any(is.na(tmp.match))!=T) next
        # if no match for forard match, we set it to backward match and leave it to outside condition
        else
          tmp$Time=tmp$Time + step*(step.cnt+1)
       }
      else if(fwd)   tmp$Time=tmp$Time + step*step.cnt
      else if(bkwd)  tmp$Time=tmp$Time - step*step.cnt
      tmp.match = merge(tmp, tar, by=by.what, all.x=T, sort=TRUE,suffixes=c(".src",".tar"))
    }

    ret = rbind(ret, tmp.match)
  }

  return(ret)
}




#
# phd::dmApnd
#   append (column-wise) annotations onto a data file according to rownames of the data frame
#
# PARAMETERS
#   data:         data frame to be operated upon (whose rownames to be found in annot)
#   annot:        an annotation object
#   idx.match:    idx of column which to be matched
#   what:         a list of variable names to be appended and (searched in annotation object)
#
# RETURN VALUE
#   data
#
# TIME STAMP
#   2009-05-24
#
# DBG
#   what = list("Assay","Challenge","TimePoints","Subject","PhenoNew","TimePointsT")
#   test = dmApnd(test,dgep$obs,what=what)
dmApnd <- function (data=NULL,objannot=NULL,idx.match=NULL,what=NULL){

  if (is.null(data) | is.null(objannot) | is.null(what)) stop("give me the correct inputs...stupid !!!")
  
  ret = NULL
  data = data.frame(data)
  
  for (i in 1:length(what)) {
    if(is.null(idx.match)) data = cbind(data,objannot[rownames(data),what[[i]]])
    else  data = cbind(data,objannot[data[,idx.match],what[[i]]])
    colnames(data)[dim(data)[2]] = what[[i]]
  }
  ret = data
  return(ret)
}


#
# phd::bsplcls2
#   CRITICAL - Changed fitting over aggregated median of data to direct spline fitting on data
#
#   conditional cubic smoothing spline fitting (conditional on "Time" and "PhenoNew")
#   note: time is the time point when sample is drawn; PhenoNew={"Asx","Sx"}
#
# TIME STAMP
#  2009-05-31
#
# PAMETERS
#  avec:          a vector of values (eg., expression values of one gene from multiple samples)
#  annot:         an object contains annotation for the samples, it has to have at least the following
#                   > Time:       categorical variable on which bspline to be fitted
#                   > PhenoNew:   categorical variable on which separated bspline to be fitted
#  degfree:       degree of spline fitting
#
# DEPENDENCIES
#
# RETURN VALUE
#   ret            -  a dataframe of design matrix (row: features; column: observations)
#
# DBG
#
# USAGE
#  test = NULL; test.a=NULL; test.bsplfit2=NULL;
#  test = read.table("plasma.AC.val.txt", sep="\t",header=T,row.names=1)
#  test.a = read.table("plasma.AC.cov.txt", sep="\t",header=T,row.names=1)
#  test.a = t(test.a)
#  test = t(scale(t(test),center=T,scale=T))
#  test.bsplcls=apply(test, 1, bsplcls, annot=dgep$obs,degfree=4)
#  write.csv(test.bsplcls, "test.bsplcls.csv")
#
# NOTE
# 1)some error checking mechanism is needed
#   in one case, we had "Z01H00", and when we substring it by right("Z01H00",5), we end up with
#   "1H00", and then we gsub("H","","1H00"), we have a time knot at 100, instead of the correct 0
#   we eventually fix it by as.numeric(gsub("H","",substr(rownames(test.a),4,10)))
#   also, in the annotaiton file, TimePoints are like "H00", not the sample name.
#   Nonetheless, we need a summary printout of what was fitted to make sure no stupid things being done
#
# 2) 2009-05-31: THIS IS NOT TRUE ANY MORE IN BSPLCLS2 FITTING, which is closer to implementation of Storey's EDGE fitting method
#    data points are first aggregated using MEDIAN at each time point; then cubic spline fitting of median values
#    this is an implmentation of what Storey did in EDGE
#
# 3)  2009-05-31: IMPLEMENTED
#     we could fit cubic spline directly to the data points
#
#     in one example FLU>6868_at>H00-H108, it doesn't seem to be too much differences.
#     need to investigate later again to make sure
#

bsplcls2 <- function (avec=NULL, annot=NULL, degfree=0) {

 if(degfree==0) stop("fitting a spline with df=0, are you crazy?")
 if(is.null(annot)) stop("where do you want me to find pheno and time info, you stupid?")

 avec.bsplfit = NULL
 test = avec

 gname = attr(avec,"names")

 # Time = as.numeric(annot[match(gname, rownames(annot),nomatch=NA),'Time'])
 # critical
 # this makes sure we don't get tripped by time like 69.5 which is coded as H069.
 # critical
 Time = as.numeric(gsub("H","",annot[match(gname, rownames(annot),nomatch=NA),'TimePoints']))
 PhenoNew = annot[match(gname, rownames(annot),nomatch=NA),'PhenoNew']

 # fitting a smoothing spline (Hastie) for individual PhenoNew (namely, 'A' and 'S')
 tmp = sort(unique(PhenoNew))
 for (i in 1:length(tmp)) {
  ix = which(PhenoNew==tmp[i])
  apheno.spl=NULL; apheno.spl = smooth.spline(x=Time[ix], y=avec[ix], df=degfree)

  # equivalent way of getting fitted value, other than apheno.spl$y
  sppred=NULL; sppred = predict(object=apheno.spl,x=sort(unique(Time[ix])))
  # make sure this is the same as the fitted
  if(any(sppred$y!=apheno.spl$y)) stop("critical error: smoothing spline fitting does not match")

  tmpfit = sppred$y
  names(tmpfit)= paste(tmp[i],"H",apheno.spl$x,sep="")
  avec.bsplfit = c(avec.bsplfit,tmpfit)
 }

 return (avec.bsplfit)

}

#
# phd::bsplcls
#   conditional cubic smoothing spline fitting (conditional on "Time" and "PhenoNew")
#   note: time is the time point when sample is drawn; PhenoNew={"Asx","Sx"}
#
# TIME STAMP
#  2009-05-13
#
# PAMETERS
#  avec:          a vector of values (eg., expression values of one gene from multiple samples)
#  annot:         an object contains annotation for the samples, it has to have at least the following
#                   > Time:       categorical variable on which bspline to be fitted
#                   > PhenoNew:   categorical variable on which separated bspline to be fitted
#  degfree:       degree of spline fitting
#
# DEPENDENCIES
#
# RETURN VALUE
#   ret            -  a dataframe of design matrix (row: features; column: observations)
#
# DBG
#
# USAGE
#  test = NULL; test.a=NULL; test.bsplfit2=NULL;
#  test = read.table("plasma.AC.val.txt", sep="\t",header=T,row.names=1)
#  test.a = read.table("plasma.AC.cov.txt", sep="\t",header=T,row.names=1)
#  test.a = t(test.a)
#  test = t(scale(t(test),center=T,scale=T))
#  test.bsplcls=apply(test, 1, bsplcls, annot=dgep$obs,degfree=4)
#  write.csv(test.bsplcls, "test.bsplcls.csv")
#
# NOTE
# 1)some error checking mechanism is needed
#   in one case, we had "Z01H00", and when we substring it by right("Z01H00",5), we end up with
#   "1H00", and then we gsub("H","","1H00"), we have a time knot at 100, instead of the correct 0
#   we eventually fix it by as.numeric(gsub("H","",substr(rownames(test.a),4,10)))
#   also, in the annotaiton file, TimePoints are like "H00", not the sample name.
#   Nonetheless, we need a summary printout of what was fitted to make sure no stupid things being done
#
# 2) data points are first aggregated using MEDIAN at each time point; then cubic spline fitting of median values
#    this is an implmentation of what Storey did in EDGE
#
# 3) we could fit cubic spline directly to the data points
#     in one example FLU>6868_at>H00-H108, it doesn't seem to be too much differences.
#     need to investigate later again to make sure
#
bsplcls <- function (avec=NULL, annot=NULL, degfree=0, what.fun="median") {
 avec.bsplfit = NULL
 test = avec
 gname = attr(avec,"names")

 # Time = as.numeric(annot[match(gname, rownames(annot),nomatch=NA),'Time'])
 # critical
 # this makes sure we don't get tripped by time like 69.5 which is coded as H069.
 # critical
 Time = as.numeric(gsub("H","",annot[match(gname, rownames(annot),nomatch=NA),'TimePoints']))
 PhenoNew = annot[match(gname, rownames(annot),nomatch=NA),'PhenoNew']

# # Splits the data into subsets,
# # computes summary statistics for each,
# # and returns the result in a convenient form.
 test.avg = aggregate(x=test, by=list(Time=Time, PhenoNew=PhenoNew), FUN=what.fun)

# colnames(test.avg)[which(colnames(test.avg)=="x")]=gname
 tmp = unique(test.avg$PhenoNew)
 for (i in 1:length(tmp)) {
  apheno.avg = test.avg[test.avg$PhenoNew==tmp[i],]
  # fitting a B-spline
  apheno.spl = smooth.spline(x=apheno.avg$Time, y=apheno.avg$x, df=degfree)
#  # plot original data points
#  plot(apheno.avg$Time, apheno.avg$x)
#  lines(apheno.spl,col="blue")
#  # plot the fitted line
#  lines(pp <- predict(apheno.spl,apheno.avg$Time), col = "blue")
#  # plot fitted value
#  points(apheno.avg$Time, pp$y, pch=3, col="dark red")
  sppred = predict(object=apheno.spl,x=apheno.avg$Time)

  tmpfit = sppred$y
  names(tmpfit)= paste(tmp[i],"H",apheno.avg$Time,sep="")
  avec.bsplfit = c(avec.bsplfit,tmpfit)
 }

 return (avec.bsplfit)
 
}




getDesign <- function(list.data=NULL, list.what=NULL, try.symbol=TRUE) {
  #
  # phdlib::getDesign
  #   construct the design matrix for downstream analysis
  #
  # TIME STAMP
  #  2009-05-07
  #
  # PAMETERS
  #  list.data       -  list of three elements 
  #                     > a) data matrix; b) feature annotation; c) observation annotation
  #  list.what       -  list of two elements
  #                     > a) features; b) observations
  #
  # DEPENDENCIES
  #
  # RETURN VALUE
  #   ret            -  a dataframe of design matrix (row: features; column: observations)
  #
  # DBG
  #
  # USAGE
  #
  #   test.design = getDesign(list.data=dprtn, list.what=test)
  #   (test.design = na.omit(test.design))
  #
  if(is.null(list.data) | is.null(list.what)) stop("\n Oops, arguments are not correctly specified")

  n.features = length(list.what$features)
  n.obs = length(list.what$obs)
  
  if(length(unique(list.what$features))!=n.features | length(unique(list.what$obs))!=n.obs )
    stop("check your obs/features for duplicates...stupid")

  ix.features = match(list.what$features,rownames(list.data$data))
  # if features are not completely matched, let's try see if supplied are "Symbol"
  if(try.symbol & any(is.na(ix.features))) {
    ix.features = match(list.what$features,list.data$data$Symbol)
  }
  ix.obs = match(list.what$obs,colnames(list.data$data))
  
  if(length(ix.features)!=n.features | length(ix.obs)!=n.obs) stop("\n Oops, not a complete match, you moron")
  
  if(any(is.na(ix.features),is.na(ix.obs))) {
   cat("not all features/observations can be matched in provided data \n")
   cat("features not found  >>  ")
   cat(list.what$features[which(is.na(ix.features))])
   cat("\n")
   cat("observations not found  >>  ")
   cat(list.what$obs[which(is.na(ix.obs))])
   cat("\n")
   stop("check your data...moron!")
  }

  ret = data.frame(list.data$data)[ix.features,ix.obs]
  
  return (ret)
}


#
# func::venn elem
#
# USAGE
# input
#
#Challenge	Symbol
#H	CRP
#H	IL10
#H	IL5
#H	IL1R1
#R	B2M
#R	CRP
#R	IL10
#R	TNFRSF1B
#R	CSF3
#R	IL5
#Z	B2M
#Z	CRP
#Z	APCS
#Z	TNFRSF1B
#Z	HP
#Z	IL18
#Z	IL1R1
#Z	VCAM1
#Z	FGA
#Z	VWF
#
venn.elem <- function (x) {
  test = read.csv("C:\\Users\\yongsheng\\Desktop\\prtncm.csv",header=T,stringsAsFactor=FALSE)
  groups = unique(test$Challenge)
  n.groups = length(groups)
  all.groups = list(); for (i in 1:n.groups) all.groups[[i]] = test[test$Challenge==groups[i],]$Symbol;
  tmp = list()
  cnt=1

  # all intersect
  common = all.groups[[cnt]]
  for (i in 2:n.groups) common = intersect(common,all.groups[[i]]);
  tmp[[cnt]] = common; names(tmp)[cnt]='common'

  # pair-wise intersect
  for (i in 1:(n.groups-1))
    for (j in (i+1):n.groups) {
      cnt=cnt+1;
      tmp[[cnt]]=setdiff(intersect(all.groups[[i]],all.groups[[j]]), common)
      names(tmp)[cnt]=paste(groups[i],"vs",groups[j],sep="")
  }

  # one group
  for (i in 1:n.groups) {
    cnt=cnt+1;

    others=NULL;
    for (j in 1:n.groups)
      if (j!=i) {
        others=union(others, all.groups[[j]])
      }

    tmp[[cnt]]=setdiff(all.groups[[i]],others);
    names(tmp)[cnt]=paste(groups[i],sep="")
  }

  return (tmp)
}


saplstats <- function(list.data=NULL) {

  #
  # phdlib::saplstats
  #  compute basic stats for each challenge study
  #
  # TIME STAMP
  #  2009-10-06
  #
  # PAMETERS
  #  list.data   -  list; three elements a) data matrix; b) feature annotation; c) observation annotaiton
  #
  # DEPENDENCIES
  #
  # RETURN VALUE
  #   ret        -  list; each element (per assay) is a list of elements for individual challenges
  #
  # DBG
  #
  # USAGE
  #
  # need much more work to make this more elegant   (2009-05-03)
  #

  if(is.null(list.data)) stop("Incorrect arguments specified")
  
  ret = list()
  
  list.data$obs = list.data$obs[list.data$obs$Use != "N",]
  
  assay = sort(unique(list.data$obs$Assay))
  print(assay)
  chlng = sort(unique(list.data$obs$Challenge))
  print(chlng)
  for (i in 1:length(assay)) {
  
    one.assay = list()

    for (j in 1:length(chlng)) {
        tmp = list.data$obs[list.data$obs$Assay==assay[i] & list.data$obs$Challenge==chlng[j],]

        pre.bl = sort(rownames(tmp)[grep("BL",rownames(tmp))])
        pre.pc = sort(rownames(tmp)[grep("H00$",rownames(tmp))])
        pre = list.data$data[,c(pre.bl,pre.pc)]
        pre.N = ncol(pre)
        pre.mean = apply(pre,1,mean)
        pre.stdev = sqrt(apply(pre,1,var))

        if(length(pre.bl)>0 & length(pre.pc)>0){
          ix.bl=match(left(pre.bl,4),left(pre.pc,4))
          ix.pc=match(left(pre.pc,4),left(pre.bl,4))
          pairpre = list.data$data[,c(pre.bl[ix.bl],pre.pc[ix.pc])]
          pairpre.N = ncol(pairpre)/2
          pairpre.meandiff = apply(pairpre, 1, function(x){mean(x[(1+pairpre.N):(2*pairpre.N)] - x[1:pairpre.N])})
          pairpre.stdevdiff = apply(pairpre, 1, function(x){sqrt(var(x[(1+pairpre.N):(2*pairpre.N)] - x[1:pairpre.N]))})
        }
        one.assay[[j]]=list(pre.N=pre.N
                            ,pre.mean=pre.mean
                            ,pre.stdev=pre.stdev
                            ,pairpre.N=pairpre.N
                            ,pairpre.meandiff=pairpre.meandiff
                            ,pairpre.stdevdiff=pairpre.stdevdiff)
        names(one.assay)[j]=paste(assay[i],chlng[j],sep="")
      }
      ret[[i]] = one.assay; names(ret)[i]=paste(assay[i],sep="")
  }
  return(ret)
}



saplpair <- function(xtab=NULL, what=NULL) {
  #
  # func::saplpair
  #  obtain samples paired at diff time points
  #
  # TIME STAMP
  #  2009-05-03
  #
  # PAMETERS
  #  list.data   -  list; three elements a) data matrix; b) feature annotation; c) observation annotaiton
  #  what        -  a vector of two time points
  #
  # DEPENDENCIES
  #
  # RETURN VALUE
  #   ret     -   a list with two elements for two time points' samples
  #
  # DBG
  #
  #
  # USAGE
  #
  #
  if(is.null(xtab) | is.null(what)) stop("arguments are not correctly specified")
  
  ret = list()
  
  assay = unique(xtab$assay)
  
  tmp = apply(xtab[,what], 1, function(x) {return(!is.na(all(x==1)))})
  tmp = xtab[tmp,what]
  for(i in 1:nrow(tmp)) for(j in 1:ncol(tmp)) tmp[i,j]=paste(assay,rownames(tmp)[i],colnames(tmp)[j],sep="")
  
  for(i in 1:ncol(tmp)) {ret[[i]] = tmp[,i]; names(ret)[i]=colnames(tmp)[i]}
  
  return(ret)
}

#
# phd::sortData
#  sort a data frame (or matrix) by one (or multiple) criteria (eg, "challenge" and "pheno")
#
# TIME STAMP
#  2009-05-16
#
# PAMETERS
#  data:      a data frame (or matrix)
#  objannot:  an object contains annotation where sorting criteria to be searched
#  what:      a list of critria (eg., list("challenge","phenonew","subject")
#  by:        1 = "row", 2 = "column"; indicating how operation should be carried out
#
# DEPENDENCIES
#  none
#
# RETURN VALUE
#   ret:      a re-ordered data frame in the same dimension/format as input
#
# DBG
#   lib::sort.df
#
# USAGE
#   test = read.csv("ebc.data.csv",header=T,row.names=1)
#   what = list("Challenge","PhenoNew","TimePoints","Subject")
#   test.sort = sortData(test, what=what, objannot=dprtn$obs, by=2)
#
sortData <- function (data=NULL, objannot, what=NULL, by=1, ...) {

  if(is.null(data) | is.null(objannot) | is.null(what)) stop("get me the right arguments input, stupid!")

  ret = NULL;
  
  # for convenience, make sure we always operate on rows
  # flip the sorted data back as before when we are done
  if (by==2) data = data.frame(t(data))

  n.criteria = length(what)
  for (i in 1:n.criteria) {
    data = cbind(data,objannot[match(rownames(data),rownames(objannot)),what[[i]]])
    colnames(data)[dim(data)[2]] = what[[i]]
  }

  criteria = NULL
  for (i in 1:n.criteria) {
    criteria = paste(criteria, " +", what[[i]], sep=" ")
  }
  criteria = as.formula(paste(" ~ ", criteria, sep=""))

  data = sort.df(data, by = criteria)

  if (by==2) data = data.frame(t(data))
  
  ret = data
  
  return (ret)
}


#
# phd::prepData
#  pre-process data for downstream analysis (such as zero out a protein/gene etc)
#
# TIME STAMP
#  2009-05-16
#
# PAMETERS
#  data:      a data frame (or matrix)
#  what:      what operation should be done
#             > "zero" - wipe out data entries with 0 if more than half of values are 0s
#             >
#             >
#  by:        1 = "row", 2 = "column"; indicating how operation should be carried out
#
# DEPENDENCIES
#  none
#
# RETURN VALUE
#   ret:      a data frame in the same dimension/format as input
#
# DBG
#
#
# USAGE
#   test.prep = prepData(test, what="zero",by=1); see(test.prep)
#
# NOTE
#   as of 2009/05/16, only "zero" option is implemented for preprocessing protein amino assay data
prepData <- function (data=NULL, what="zero", by=1, print=F...) {

  ret = NULL;
  msg = NULL;
  nperc = 0.5
  affected = NULL;
  
  if (what=="zero") {
     ix.zero = apply(data, by, function(x){sum(x==0)>length(x)*nperc})
     if(by==1) {data[ix.zero,]=0; affected = rownames(data)[ix.zero==T];}
     else if (by==2) {data[,ix.zero] = 0; affected = colnames(data)[ix.zero==T];}
     msg = paste("total of ", sum(ix.zero), " has ", nperc*100, "% of cells equal to 0", sep="")
     ret = data
  }
  print(msg)
  if(print) print(affected)
  return (ret)
}



#
# phd::saplxtab
#  counting samples by assay (mRNA/proteomics); time; subject
#
# TIME STAMP
#  2009-05-02
#
# PAMETERS
#  list.data       -  a list with three elements a) data matrix; b) feature annotation; c) observation annotaiton
#
# DEPENDENCIES
#
# RETURN VALUE
#   ret     -   a list with each elements for an assay, whereas each element being a list contain counts for individual challenges
#
# DBG
#
#
# USAGE
#   dprtn.xtab = saplxtab(dprtn); sink("dprtn.xtab.txt"); print(dprtn.xtab); sink();
#   dprtn.xtab.name = saplxtab(dprtn, names.fill=TRUE); sink("dprtn.xtab.name.txt"); print(dprtn.xtab.name); sink();
#   dgep.xtab = saplxtab(list.data=dgep,names.fill=FALSE); sink("dgep.xtab.txt"); print(dgep.xtab); sink();
#
saplxtab <- function(list.data=NULL, names.fill=FALSE, replicates.rm=TRUE) {

  ret = list()
  idata = list.data$data
  ifeature = list.data$features
  iobs = list.data$obs

  if(replicates.rm) { # remove all replicates from counting
      idata = idata[,- grep("R$",colnames(idata))]
      iobs = iobs[-grep("R$", iobs$IDShort),]
  }

  #
  # later
  #
  # some quality control is needed here
  #
  # 1. make sure the column names are there
  # 2. exclude those samples with "N" under "Use" column
  #
  #
  assay = unique(iobs$Assay)

  for (i in 1:length(assay)) {

    lassay = list()

    # get all annotation for one assay
    iobs.assay = iobs[iobs$Assay==assay[i],]
    chalng = sort(unique(iobs.assay$Challenge))

    for(j in 1:length(chalng)) {
        timepoints = sort(unique(iobs.assay[iobs.assay$Challenge==chalng[j],]$TimePoints))
        subj = sort(unique(iobs.assay[iobs.assay$Challenge==chalng[j],]$Subject))
        tmp=matrix(NA,nrow=length(subj),ncol=length(timepoints))
        rownames(tmp)=subj
        colnames(tmp)=timepoints
        # obtain infomation about pheno, subject etc
        lassay[[j]]=data.frame( assay=rep(assay[i],length(subj))
                               ,phenonew=getPheno(idsapls=subj,objannot=list.data$obs,by.subj=TRUE)
                               ,tmp);
        #
        # check in list.data$data for data points and populate the sample matrix
        #
        # define an on-the-fly function to do so  (may separate this into a stand-alone function at future)
        f <- function(x) {
                      for(i in 1:nrow(x))
                       for(j in 3:ncol(x)) {
                        prefix.assay = ifelse(x[i,]$assay=="mRNA","",as.character(x[i,]$assay))
                        x[i,j]=ifelse(length(grep(paste(prefix.assay,rownames(x)[i], colnames(x)[j], sep=""),colnames(idata)))==0
                                      ,NA
                                      ,1)
                       }
                      return(x) }
        # apply the function
        lassay[[j]] = f(lassay[[j]])
        #lassay[[j]] = f(idata)

        if(names.fill){
          for (ir in 1:nrow(lassay[[j]]))
            for (ic in 1:ncol(lassay[[j]]))
              if(!is.na(lassay[[j]][ir,ic]) & lassay[[j]][ir,ic]==1)
                lassay[[j]][ir,ic] = paste(ifelse(lassay[[j]][ir,'assay']=="mRNA", "", lassay[[j]][ir,'assay'])
                                            ,rownames(lassay[[j]])[ir]
                                            ,colnames(lassay[[j]])[ic]
                                            ,sep="")
        }

        names(lassay)[j]=paste(assay[i],chalng[j],sep='')

        rm(tmp)
    }

    ret[[i]] = lassay; names(ret)[i]=assay[i]

  }

  return(ret)
}


#   saplxtab2(Subject~TimePoints, data=x)

saplxtab2 <- function(formula=~.,data,exclude = c(NA, NaN), drop.unused.levels = FALSE) {

  if (missing(formula) && missing(data)) 
      stop("must supply either 'formula' or 'data'")
  if (!missing(formula)) {
      formula <- as.formula(formula)
      if (!inherits(formula, "formula")) 
          stop("'formula' missing or incorrect")
  }

  if (any(attr(terms(formula, data = data), "order") > 1)) 
      stop("interactions are not allowed")
      
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
      m$data <- as.data.frame(data)
  #m$... <- m$exclude <- m$drop.unused.levels <- m$sparse <- NULL

  m[[1L]] <- as.name("model.frame")
  mf <- eval(m, parent.frame())

  if(length(formula)==3L) {
    i <- attr(attr(mf, "terms"), "response")
    by <- mf[-i]
    y <- mf[[i]]
  } 
  else {
    return()
  }         

  by <- unlist(lapply(by, function(u) {
          if (is.factor(u))  u <- levels(u) 
          else u
        }))
 
  y <- unlist(lapply(y, function(u) {
          if (is.factor(u)) u <- as.character(u)
          else u
        }))

  tmp = cbind(by, y)
  
  obs = unique(y)
  times = unique(by)
  ret = matrix(0, nrow=length(times), ncol=length(obs))
  colnames(ret) = obs; rownames(ret)=times;

  apply(tmp, 1, function(u) {ret[u[1], u[2]] <<- 1})
  rm(tmp)
  return(ret)
}



#***********************************************************************************************
# phdlib::loadData
#   loading phd challenge data into environment
#
# PAMETERS
#   what        -   what type of dataset to load
#
# VALUE
#   ldata       -   a list with three elements
#                   a) data matrix; b) feature annotation; c) observation annotaiton
#
# DEPENDENCIES
#   1) assuming specific file location (windows machine)
#   2) assuming all input file has HEADER and 1st column contains row.names
#   3) 2010-02-07 Changes: for ianot, only the usable samples (Use=="Y") will be assigned row.names
#   4) {library}::countComments
#
# DBG
#   setting parmeters for running the function
#
# USAGE
#   exemplary usage of the function
#
# TIME STAMP
#   2010-10-05
#***********************************************************************************************

loadData <- function( what=NULL, lnprint=0,ddir=NULL
                     ,nameAnot=NULL,nameGene=NULL,nameClinics=NULL,nameGep=NULL) {

  if(is.null(what)) stop("What type of data needs to be load: protein or gep?")

  ldata = list()
  ldata$data = NULL           # data matrix (assuming to be observations by features)
  ldata$features = NULL       # annotaiton of features (genes/proteins)
  ldata$obs = NULL            # annotation of observations
  
  msg = "data loading failed"

  #
  # loading protein data
  #
  if(what=="protein") {
    # set up the data file names
    ddata = "C:\\YS_Documents\\UM\\research\\hero\\phd\\analysis\\data\\proteomics\\"
    if(length(grep("unix",.Platform))>0) ddata=NULL
    fdata = "darpa.HRF.Proteomics.PWUE._CURRENT_.data.csv"
    fprtn = "darpa.HRF.Proteomics.PWUE._CURRENT_.proteins.csv"
    fsapl = "darpa.HRF.Proteomics.PWUE._CURRENT_.samples.csv"
    # read in files
    pmdata = read.csv(cstr(c(ddata,fdata)),header=TRUE,row.names=1,stringsAsFactors=FALSE)
    pmprtn = read.csv(cstr(c(ddata,fprtn)),header=TRUE,row.names=1,stringsAsFactors=FALSE)
    pmsapl = read.csv(cstr(c(ddata,fsapl)),header=TRUE,row.names=1,stringsAsFactors=FALSE)

    # look at the data
    if(lnprint>0){
      pmdata[1:lnprint,1:lnprint]
      pmprtn[1:lnprint,]
      pmsapl[1:lnprint,1:lnprint]
    }
    
    ldata$data = pmdata
    ldata$features = pmprtn
    ldata$obs = pmsapl

    msg = "protein data has been sucessfully loaded ..."
  }

  #
  # loading gene expression data
  #
  if(what=="gep") {
#    ddata = "C:\\Users\\yongsheng\\Documents\\Yongsheng's Documents\\UM\\research\\hero\\phd\\Analysis\\data\\joint\\"
    ddata = "C:\\YS_Documents\\UM\\research\\hero\\phd\\Analysis\\data\\gep\\"
    if(!is.null(ddir)) ddata = ddir
#    if(length(grep("unix",.Platform))>0) ddata=NULL       # if it is on unix server, default to current folder
#    fanot = "darpa.HRF.annotation.2009.06.21.csv"
#    fgep = "darpa.HRF.joint_hrv_rsv_flu.rma.(n=730).ccdf.adj.OUTPUT.csv"
#    fgene = "darpa.HRF.HG-U133A_2.Entrez.Annotation.csv"
#    fanot = "darpa.HRFW.annotation.2009.10.22.csv"
#    fgene = "darpa.HRF.HG-U133A_2.Entrez.Annotation.csv"
#    fgep =  "darpa.HRFW.joint_hrv_rsv_flu_h1n1.rma.(n=1108).ccdf.adj.OUTPUT.csv"
#    fclinics = "darpa.HRFW.clinics.2009.09.28.csv"
    fanot = "darpa.h3n2.annotation.csv"
    fgene = "darpa.HRF.HG-U133A_2.Entrez.Annotation.csv"
    #fgep =  "darpa.HRFW.joint_hrv_rsv_flu_h1n1.rma.(n=1108).ccdf.adj.OUTPUT.csv"
    fgep = "darpa.gep.rma.processed.data.ComBat.Adjust_2010-Oct-03.csv"
    fclinics = "darpa.HRFW.clinics.2009.09.28.csv"

    if(!is.null(nameAnot)) fanot = nameAnot
    if(!is.null(nameGene)) fgene = nameGene
    if(!is.null(nameGep))  fgep = nameGep
    if(!is.null(nameClinics)) fclinics = nameClinics

#    ianot = read.csv(cstr(c(ddata,fanot)),header=T,row.names=1,stringsAsFactors=FALSE)
#    igene = read.csv(cstr(c(ddata,fgene)),header=T,row.names=1,stringsAsFactors=FALSE)
#    igep = read.csv(cstr(c(ddata,fgep)),header=T,row.names=1,stringsAsFactors=FALSE)
#    iclinics = read.csv(cstr(c(ddata,fclinics)),header=T,row.names=1,stringsAsFactors=FALSE)
    ianot = read.csv(file.path(ddata,fanot),header=T,stringsAsFactors=FALSE)        # 2010-02-07 chg for duplicate rownames
    if (length(grep("Use",names(ianot)))>0) {
      tmpIx = which(ianot[,"Use"]!="Y"); if(length(tmpIx)!=0) {ianot = ianot[-tmpIx, ]}
    }
    rownames(ianot) = ianot[,1]
    # ianot = ianot[-which(ianot[,"Use"]!="Y"), ]; rownames(ianot) = ianot[,1]      # 2010-02-07 chg for duplicate rownames
    
    igene = read.csv(file.path(ddata,fgene),header=T,row.names=1,stringsAsFactors=FALSE)
    n_skip=0; n_skip = countComments(fname=file.path(ddata,fgep),patn_comments="#",nrows=6)
    igep = read.csv(file.path(ddata,fgep),header=T,row.names=1,stringsAsFactors=FALSE,skip=n_skip)
    iclinics = read.csv(file.path(ddata,fclinics),header=T,row.names=1,stringsAsFactors=FALSE)

    #
    # add a column of Symbol for easy view of genes
    #
    if (is.na(match("symbol",tolower(colnames(igep))))) igep$Symbol = igene$Symbol[match(rownames(igep),igene$Probesets)]
    
    if(lnprint>0){
      print(igep[1:lnprint,1:lnprint])
      print(ianot[1:lnprint,])
      print(igene[1:lnprint,])
      print(iclinics[1:lnprint,])
    }

    ldata$data = igep
    ldata$features = igene
    ldata$obs = ianot
    ldata$clinics = iclinics

    msg = "gene expression data has been sucessfully loaded ..."
  }
  
  cat(msg,"\n")
  
  return(ldata)
  
}






#
# func::cstrSample
#   step-wise constrained sampling
#
# TIME STAMP
#     2009-03-21
#
# PARAMETERS
#   indata      - data frame of input data contains all info needed for sample selection
#   proportion  - proportion of samples should be selected (# of samples * proportion) will be rounded
#   ifcstr      - boolean value that signals whether a constrained sampling should be done
#   ixcstr      - index of column upon which constrained sampled should be conditioned
#
# DEPENDENCIES
#  annotation file should be either in current folder or (windows machine) at specific location
#
# DBG
#
# USAGE
#
# TO DO (2009-03-23)
# ADD class balance (Balancing via resampling)
#

cstrSample <- function (indata=NULL,proportion=NULL,ifcstr=FALSE,ixcstr=NULL) {
#  if((indata==NULL)||(proportion==NULL)) stop("Arguments are missing")

  cstrspl=indata[,ixcstr[1]]
  for (i in 2:length(ixcstr)) cstrspl = paste(cstrspl,indata[,ixcstr[i]],sep="")

  indata$cstrspl = cstrspl

  cnt.cstrspl = table(cstrspl)
  cnt.cstrspl = round(cnt.cstrspl * proportion,0)
  n.cstrspl = length(cnt.cstrspl)
  

  ixspls = NULL
  for (i in 1:n.cstrspl) {
    whatspl = names(cnt.cstrspl)[i]
    howmany = (cnt.cstrspl)[i]
    one.sampling = sample(which(indata$cstrspl==whatspl),howmany)
    ixspls = c(ixspls, one.sampling)
  }

  return(indata[ixspls,])
}
#y = cstrSample(x,0.3,T,c(2,4))


#
# func::getPheno
#   to get pheno labels given a vector of sample labels
#
# DEPENDENCIES
#  annotation file should be either in current folder or (windows machine) at specific location
#
#getPheno <- function (idsapls=NULL,annot='',objannot=NULL,prnt=FALSE, by.subj=FALSE) {
getPheno <- function (idsapls=NULL,annot='',objannot=NULL,by.subj=FALSE) {

  if(is.null(objannot)) stop("Please provide annotation object")

  ret = NULL

  # start>> 2009-05-02 added support for getPheno with Subject
  if(by.subj==TRUE) return(unlist(lapply(idsapls,function(x){unique(objannot[objannot$Subject==x,]$PhenoNew)})))
  # end>> 2009-05-02 added support for getPheno with Subject

  # load the library if it's not already been loaded
  # note there is NO quote next to flag_ys_library, R will look for this object in global environment
  # if(class(flag_ys_library)=="function"){} else {source("C:\\Users\\yongsheng\\Documents\\Yongsheng\'s Documents\\Develop\\R\\library.r")}

  # note the following code does not work since "flag_ys_library"
  # will be considered as a local variable taking on a string form.
  # the check will always end up with source the library
  #  if(is.element("flag_ys_library",objects())){} else {source("C:\\Users\\yongsheng\\Documents\\Yongsheng\'s Documents\\Develop\\R\\library.r")}

#  if(annot=='') colPheno='PhenoNew'
#  else colPheno=annot
#
#  colPheno = 'PhenoNew'
#  ddata = "C:\\Users\\yongsheng\\Documents\\Yongsheng's Documents\\UM\\research\\hero\\phd\\Analysis\\data\\joint\\"
#  fanot = "darpa.HRF.annotation.2009.03.21.csv"
  
#  if(is.element(fanot,dir())) ianot=read.csv(fanot,header=T,row.names=1,stringsAsFactors=FALSE)
#  else ianot = read.csv(cstr(c(ddata,fanot)),header=T,row.names=1,stringsAsFactors=FALSE)

#  if(prnt) {print(str(ianot)); stop()}

  # assign the annotation object if supplied
#  if(!is.null(objannot)) ianot=objannot

#  ret = ianot[match(idsapls,rownames(ianot)),colPheno]

  ret = objannot[match(idsapls,rownames(objannot)),annot]
  
  return(ret)
  
}

#
# func::getFeatures
#   to get meta info about features (genes/proteins) given a vector of features
#
# DEPENDENCIES
#  annotation file should be either in current folder or (windows machine) at specific location
#
# USAGE
#  getFeatures(gsub("`","",attr(which(abs(coef(x.test[[2]][[1]]))>0)[-1],"names")),annot='Symbol',objannot=dgep$features)
getFeatures <- function (idfeatures=NULL,annot='',objannot=NULL) {

  if(is.null(objannot)) stop("Please provide annotation object")

  ret = NULL

  if(annot=='') ret = objannot[match(idfeatures,rownames(objannot)),]
  else ret = objannot[match(idfeatures,rownames(objannot)),annot]

  return(ret)

}


#
# func:getGene
#   obtain annotation for a list of features (gene or protein)
#
# PARAMETER
#   idgenes:        a list of input genes/proteins
#   feature.annot:  an object of annotation
#   what:           a list of annotations to be obtained
#   by.symbol:      whether to search by symbol (gene)
#
# TIME STAMP
#   2009-05-24
#
# DEPENDENCIES
#  annotation file should be either in current folder or (windows machine) at specific location
#
getGene <- function (idgenes=NULL,feature.annot=NULL,what=NULL,by.symbol=FALSE) {

  ret = NULL

  if(!is.null(feature.annot)) {

    ix = NULL

    if(by.symbol) ix = match(idgenes,feature.annot$Symbol)
    else ix = match(idgenes, rownames(feature.annot))

    if(is.null(ix)|length(ix)!=length(idgenes)) stop("some genes did not match")

    if (is.null(what)) ret = feature.annot[ix,]
    else ret = feature.annot[ix,what]

    return (data.frame(ret))

  }

  # load the library if it's not already been loaded
  # note there is NO quote next to flag_ys_library, R will look for this object in global environment
  if(class(flag_ys_library)=="function"){} else {source("C:\\Users\\yongsheng\\Documents\\Yongsheng\'s Documents\\Develop\\R\\library.r")}
  ddata = "C:\\Users\\yongsheng\\Documents\\Yongsheng's Documents\\UM\\research\\hero\\phd\\Analysis\\data\\joint\\"
  fgene = "darpa.HRF.HG-U133A_2.Entrez.Annotation.csv"
  if(is.element(fgene,dir())) igene=read.csv(fgene,header=T,row.names=1,stringsAsFactors=FALSE)
  else igene = read.csv(cstr(c(ddata,fgene)),header=T,row.names=1,stringsAsFactors=FALSE)

  ret = NULL
  ret = igene[match(idgenes,rownames(igene)),]

  print("This annotation is provided for backward compatibility and might be outdated")
  print("To get most recent annotation, specify an annotation object")

  return(ret)
}




#***********************************************************************************************
# phdlib::splsel
#    function for sample selection
#    - current implmentation only includes absolute time points
#    - relative clinic time point should be added
#
# PAMETERS
#    idata     - data frame of data
#    ianotsapl - data frame of sample annotation
#    ianotvar  - data frame of feature (variable/biomarkers) annotation
#
# VALUE
#
# DEPENDENCIES
#   1) assuming specific file location (windows machine)
#   2) assuming all input file has HEADER and 1st column contains row.names
#
#  DBG:
#   dgep = loadData("gep")
#   list.data=dgep; ta = c(0,0); vir = c('Z','Z'); npheno=2; ord='pheno'  # flu:asx-sx: H00
#   do.log=FALSE; raffy=FALSE; rrep=TRUE;
#
# USAGE
#   exemplary usage of the function
#
# TIME STAMP
#     2009-08-12
#***********************************************************************************************

#
# proteomic data
#   note proteomic data is assumed to be  (observation x protein-biomarkers)
#   a transposed version of gene expression data
#   ta = c(-1,0,12,21,93,108); vir = rep('AZ',6); npheno=2; ord='pheno';
#    x = splsel(ta,vir,idata=t(pmdata),ianotspl=pmsapl,ianotvar=pmprtn)

#ta = c(72,117,108); vir = c('C','R','Z'); npheno=2; ord='pheno'
#ta = c(-1,0); vir = c('Z','Z'); npheno=2; ord='pheno'  # flu:asx: baseline and H00
#
#ta = c(-1,0); vir = c('Z','Z'); npheno=-1; ord='pheno'  # flu:asx: baseline and H00
#ta = c(-1,21); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H012 for flu
#ta = c(0,12); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H012 for flu
#ta = c(0,21); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H012 for flu
#ta = c(0,45); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H045 for flu
#ta = c(0,60); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H060 for flu
#ta = c(0,69); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H069 for flu
#ta = c(0,84); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H069 for flu
#
#ta = c(0,0); vir = c('Z','Z'); npheno=2; ord='pheno'  # flu:asx-sx: H00

#ta = c(5); vir = c('Z','Z'); npheno=2; ord='pheno'  # flu:asx-sx: H005
#
#ta = c(-1,0); vir = c('R','R'); npheno=2; ord='pheno'  # flu:asx: baseline and H00
#
#ta = c(-1,0); vir = c('C','C'); npheno=2; ord='pheno'  # flu:asx: baseline and H00
#
#
#ta = c(0,-1,0,-1); vir = c('R','R','Z','Z'); npheno=2; ord='pheno'  # HRV-and-FLU: H00
#splsel(ta,vir,npheno,ord,raffy=TRUE,rrep=TRUE)

# ta = c(0,5,0,5); vir = c('Z','Z'); npheno=2; ord='pheno'  # FLU: H00 and H05
# ta = c(0,12,0,12); vir = c('Z','Z'); npheno=2; ord='pheno'  # FLU: H00 and H12

# test = splsel(ta,vir,list.data=dgep, npheno,ord,raffy=TRUE,rrep=TRUE)

splsel <- function (ta=NULL, vir=NULL, list.data=NULL, npheno=2, ord='pheno', do.log=FALSE, raffy=FALSE, rrep=TRUE) {

  ret = NULL
  
  igep = list.data$data           # make a local copy of gep
  ianot = list.data$obs           # make a local copy of annotation
  igene = list.data$features      # make a local copy of gene annotation

  # specific to gene expression
  # remove problematic samples (C14 and Unknwn)
  if (TRUE) {
    ixrmv = which(ianot$Use=='N')
    if (length(ixrmv>0)) ianot=ianot[-ixrmv,]
  }

  # specific to gene expression
  # remove technical replicates
  if (rrep==TRUE) {
    ixrmv = grep("R$",rownames(ianot))                  # n=14
    if (length(ixrmv>0)) ianot=ianot[-ixrmv,]
  }

  # remove affy control probes
  if (raffy==TRUE) {
    ixaffy = grep('AFFX',rownames(igep))
    if(length(ixaffy>0)) igep=igep[-ixaffy,]

    ixaffy = grep('AFFX', rownames(igene))
    if(length(ixaffy>0)) igene = igene[-ixaffy,]
  }

  # explicitely assign pheno to avoid confusion
  if(npheno==-1) pheno='A'
  if(npheno==1) pheno='S'
  if(npheno==2) pheno=c('A','S')

  # number of viral challenges
  nvir = length(vir)

  ix = NULL
  sapl = NULL
  if (ord=='pheno') {
    for(i in 1:npheno) {
      for (j in 1:nvir) {
        ix =  which(ianot$Time>(ta[j]-1) & ianot$Time<(ta[j]+1) & ianot$PhenoNew ==pheno[i] & paste(ifelse(ianot$Assay=="mRNA","",ianot$Assay),ianot$Challenge,sep="") %in% vir[j])
        sapl = c(sapl, sort(rownames(ianot[ix,])))
      }
    }
  }

  if (ord=='challenge') {
    for (j in 1:nvir) {
      for(i in 1:npheno) {
        ix =  which(ianot$Time>(ta[j]-1) & ianot$Time<(ta[j]+1) & ianot$PhenoNew ==pheno[i] & paste(ifelse(ianot$Assay=="mRNA","",ianot$Assay),ianot$Challenge,sep="") %in% vir[j])
        sapl = c(sapl, sort(rownames(ianot[ix,])))
      }
    }
  }

#  x = ifelse(do.log==TRUE, log(igep[,sapl],10), igep[,sapl])
  x = igep[,sapl]

  # add gene annotation info
  x = cbind(igene[rownames(x),],x)
  x = rbind(PhenoNewInt=as.factor(ianot[colnames(x),'PhenoNew']), x)
  x = rbind(PhenoNew=ianot[colnames(x),'PhenoNew'], x)

  sname = paste("darpa.hrf",".",cstr(vir),".",cstr(pheno),".csv",sep="")
  write.csv(x, sname)

  cat("Sample gep has been written into file:","\n",sname,"\n",sep="")
  cat("Dimension:"," ", "nrows=",nrow(x)," by ", "ncols=",ncol(x),"\n",sep="")
  cat("Total number of samples = ", length(sapl),sep=" ","\n")

  ret = x
  
  return(ret)
  
}


cat(">>>>>>> PHD project library has been successfully loaded!\n")










#
# the following code was copied from phd.R
# it might be duplicated code that can be deleted later on
# but need to verify that is true
#
# 2009-03-12
#
##
##  func::splsel
##    function for sample selection
##    - current implmentation only includes absolute time points
##    - relative clinic time point should be added
##
##  TIME STAMP
##     2009-03-05
##
##  PARAMETERS
##
##  DEPENDENCIS
##
##  DBG:
##  ta = c(72,117,108)
##  vir = c('C','R','Z')
##  npheno=2
##  nvir = 3
##  ord='pheno'
##
#
#ta = c(72,117,108); vir = c('C','R','Z'); npheno=2; ord='pheno'
#ta = c(-1,0); vir = c('Z','Z'); npheno=2; ord='pheno'  # flu:asx: baseline and H00
#
#ta = c(-1,0); vir = c('Z','Z'); npheno=-1; ord='pheno'  # flu:asx: baseline and H00
#ta = c(-1,21); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H012 for flu
#ta = c(0,12); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H012 for flu
#ta = c(0,21); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H012 for flu
#ta = c(0,45); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H045 for flu
#ta = c(0,60); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H060 for flu
#ta = c(0,69); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H069 for flu
#ta = c(0,84); vir = c('Z','Z'); npheno=-1; ord='pheno'  # get asx H00 and H069 for flu
#
#ta = c(0,0); vir = c('Z','Z'); npheno=2; ord='pheno'  # flu:asx-sx: H00
#ta = c(5); vir = c('Z','Z'); npheno=2; ord='pheno'  # flu:asx-sx: H005
#
#ta = c(-1,0); vir = c('R','R'); npheno=2; ord='pheno'  # flu:asx: baseline and H00
#
#ta = c(-1,0); vir = c('C','C'); npheno=2; ord='pheno'  # flu:asx: baseline and H00
#
#
#ta = c(0,-1,0,-1); vir = c('R','R','Z','Z'); npheno=2; ord='pheno'  # HRV-and-FLU: H00
#splsel(ta,vir,npheno,ord,raffy=TRUE,rrep=TRUE)
#
#
#
#splsel <- function (ta=NULL, vir=NULL, npheno=2, ord='pheno', raffy=FALSE, rrep=TRUE) {
#
#  igep = igep       # make a local copy of gep
#  ianot = ianot     # make a local copy of annotation
#  igene = igene     # make a local copy of gene annotation
#
#  # remove problematic samples (C14 and Unknwn)
#  if (TRUE) {
#    ixrmv = which(ianot$Use=='N')
#    ianot=ianot[-ixrmv,]
#  }
#
#  # remove technical replicates
#  if (rrep==TRUE) {
#    ixrmv = grep("R$",rownames(ianot))                  # n=14
#    ianot=ianot[-ixrmv,]
#  }
#
#  # remove affy control probes
#  if (raffy==TRUE) {
#    ixaffy = grep('AFFX',rownames(igep))
#    igep=igep[-ixaffy,]
#
#    ixaffy = grep('AFFX', rownames(igene))
#    igene = igene[-ixaffy,]
#  }
#
#
#  # explicitely assign pheno to avoid confusion
#  if(npheno==-1) pheno='A'
#  if(npheno==1) pheno='S'
#  if(npheno==2) pheno=c('A','S')
#
#  # number of viral challenges
#  nvir = length(vir)
#
#  ix = NULL
#  sapl = NULL
#  if (ord=='pheno') {
#    for(i in 1:npheno) {
#      for (j in 1:nvir) {
#        ix =  which(ianot$Time>(ta[j]-1) & ianot$Time<(ta[j]+1) & ianot$PhenoNew ==pheno[i] & ianot$Challenge %in% vir[j])
#        sapl = c(sapl, sort(rownames(ianot[ix,])))
#      }
#    }
#  }
#
#  if (ord=='challenge') {
#    for (j in 1:nvir) {
#      for(i in 1:npheno) {
#        ix =  which(ianot$Time>(ta[j]-1) & ianot$Time<(ta[j]+1) & ianot$PhenoNew ==pheno[i] & ianot$Challenge %in% vir[j])
#        sapl = c(sapl, sort(rownames(ianot[ix,])))
#      }
#    }
#  }
#
#  x = igep[,sapl]
#
#  # add gene annotation info
#  x = cbind(igene[rownames(x),],x)
#  x = rbind(PhenoNewInt=as.factor(ianot[colnames(x),'PhenoNew']), x)
#  x = rbind(PhenoNew=ianot[colnames(x),'PhenoNew'], x)
#
#  sname = paste("darpa.hrf",".",cstr(vir),".",cstr(pheno),".csv",sep="")
#  write.csv(x, sname)
#
#  cat("Sample gep has been written into file:","\n",sname,"\n",sep="")
#  cat("Dimension:"," ", "nrows=",nrow(x)," by ", "ncols=",ncol(x),"\n",sep="")
#  cat("Total number of samples = ", length(sapl),sep=" ","\n")
#  print(rownames(ianot[sapl,]))
#
#}
#
#
#
#
#
#
# old version of code for selecting proteomic data
# dis
#
#selsapl <- function (ds=NULL,what=''){
#  if(is.null(ds)) stop("need name of dataset object")
#  if(what=='') stop("three letters for samples to be searched - HRV/RSV/FLU")
#  spatn = ifelse(what=='HRV','^AC',ifelse(what=='RSV','^AR','^AZ'))
#  ixflu = grep(spatn,rownames(ds))
#  if(length(ixflu)>0) return(ds[ixflu,])
#}
#
#
#pmdataflu = selsapl(pmdata,'FLU')
#pmdataflu = scale(pmdataflu)        #z-transformation of data (has been verified)
