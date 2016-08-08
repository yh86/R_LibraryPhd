#*******************************************************************
# Library.R
#   utitlity functions and modules
#
# AUTHOR
#   Yongsheng Huang   {huangys}@umich.edu
#   Ph.D Candidate
#   Bioinformatics Graduate Program
#   Department of Statistics
#   University of Michigan
#   2017 Palmer Commons, 100 Washtenaw Ave, Ann Arbor, MI 48109-2218
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# DATE (Last Updated)
#   2009-12-07
#*******************************************************************


#
# dummie function for verification of library has been loaded or not
#
flag_ys_library <- function () {return(1)}
#


flag_ys_unix=FALSE; if(Sys.getenv("OSTYPE")=="linux") flag_ys_unix=TRUE

getDfRowCol <- function(index=NULL,ds=NULL,doprnt=FALSE) {
  #**********************************************************************************
  # library::getDfRowCol
  #   given an index of a matrix/dataframe, find corresponding row and column indices
  #
  # PARAMETERS
  #   index  - integer; index of the cell
  #   ds     - df/matrix; data matrix or data frame
  #   doprnt - boolean; print out results
  # 
  # VALUE
  #   a list of two elements: (1) irow --- index of row; (2) icol --- index of column
  #
  # DEPENDENCIES
  #  
  # DBG
  #
  # NOTES
  #
  # USAGE
  #  num=100; a=matrix(rnorm(25),nrow=5); a[3,2]=num;  
  #  getDfRowCol (index=which(a==num),ds=a,doprnt=TRUE)
  #
  # KNOWN ISSUES
  #
  # TIME STAMP
  #   2010-10-18
  #**********************************************************************************  
  
  ir = index %% nrow(ds)
  ic = (index - ir) / nrow(ds) + 1
  if(doprnt==T) cat("row index: ",ir,"; ", "column index: ",ic,";\t"
                    ,"value: ", ds[ir,ic],sep="","\n\n")

  return(list(irow=ir, icol=ic))
}


slice <- function (data=NULL, howmany_cols=200,align=c(0.5,0.5),do_print=FALSE) 
{
  #***********************************************************************************************
  # library::slice
  #   Get a slice from a matrix or data frame
  #
  # PAMETERS
  #   data         -   a matrix or data frame 
  #   howmany_cols -   int; how many columns of the matrix to be obtained
  #   align        -   a vector of two int; specifying how the slice to be obtained
  #                    (0.5,0.5) equal number of columns from center of the matrix
  #                    (0.66,0.34) remove 2/3 left of center and 1/3 right of center
  #   do_print     -   boolean; do nothing but print out a summary of no. of columns sliced off
  #
  # VALUE
  #   ret          -   a matrix slice of the original data matrix 
  #
  # DEPENDENCIES
  #
  # DBG
  #
  # NOTES
  #
  # USAGE
  #   ncol_tot=20; ncol_get=10; xtest=matrix(rnorm(10*ncol_tot),ncol=ncol_tot); 
  #   for(i in 1:ncol_tot) xtest[,i] = i;
  #   slice(xtest,ncol_get,align=c(0.8,0.2),do_print=FALSE)
  #
  # KNOWN ISSUES
  #
  # TIME STAMP
  #   2010-10-14
  #***********************************************************************************************
  
  if(class(data)!="matrix"&class(data)!="data.frame") stop("data should be a matrix or data frame")

  n_cols = ncol(data)
  if(do_print) {
    cat("The total number of elements in this dataset are: \n\n"); 
    print(n_cols)
    cat("\n")  
    return(NULL);
  }
  
  ndiff_cols = max(n_cols) - howmany_cols     # number of cols to be removed
  colsDiscard = c( floor(ndiff_cols*align[1]), ndiff_cols-floor(ndiff_cols*align[1]) )
  cat("Total of ", colsDiscard[1], " columns were removed from left tail (including 0 valued). \n")
  cat("Total of ", colsDiscard[2], " columns were removed from right tail (including 0 valued). \n") 
    
  ret = data

  if(ndiff_cols!=0) ret = ret[,c((colsDiscard[1]+1) : (ncol(ret)-colsDiscard[2]))]
  n.col = nchar(ncol(ret))
  colnames(ret) = paste("p",right(paste(paste(rep("0",n.col),collapse="")
                      ,seq(1:ncol(ret)),sep=""), n.col),sep="")
   
  return(ret)

}




countComments <- function (fname = "test.csv", patn_comments = "#",nrows=100) {
  #
  # library::countComments
  #   get a head count of total number of comments lines in an input file 
  #
  #
  # fname         -   name of the input file
  # patn_comments -   pattern of leading character(s) in a comment line
  # nrows         -   number of lines to check; 
  #                   automatically double 'nrows' if all lines are found to be comments, 
  #    
  # value         -   integer; no. of lines that are comments
  #
  tmp = read.delim(fname,stringsAsFactors=FALSE,comment.char="",header=FALSE,nrows=nrows)
  nrow_comments = length(unlist(apply(tmp, 1, grep, pattern=patn_comments)))
  if(nrow_comments == nrows) nrow_comments = countComments(fname, patn_comments, nrows*2)
  else 
    cat("\nTotal of ", nrow_comments, " lines of comments are found in ",fname,".\n\n" )
  return(nrow_comments)
  
}



#
# function to get a vector by slicing a row out of a data frame
#
# DBG
#L3 <- LETTERS[1:3]
#(d <- data.frame(cbind(x=1, y=1:10), fac=sample(L3, 10, replace=TRUE)))
#ixR = sample(1:nrow(d), 1)
#class(oneRowDf(d, ixR, retVector=TRUE))     # a vector
#class(d[ixR,])                              # data frame
#class(oneRowDf(d, ixR, retVector=FALSE))    # data frame
oneRowDf <- function (dat, ixR, retVector=TRUE) {

  ret = NULL

  if (retVector) {
  dat = as.matrix(dat)
  ret = dat[ixR,]
  }

  else
    ret = dat[ixR, ]

  return(ret)

}



#
# apply t.test on one vector instead of x, y vectors, esp useful for paired t-test
#
tTestVec <- function(vec,type="paired",whatRet="p.value") {

  n = length(vec)

  if(type=="paired" & n%%2!=0) stop("request paired test by vector is asymmetric")

  n = n/2

  tmp = NULL
  options(show.error.messages = FALSE)
  if(type=="paired") tmp = try(t.test(vec[1:n], vec[(n+1):(2*n)], paired=TRUE))
  else cat("other t-test \n")
  options(show.error.messages = TRUE)

  if (class(tmp)=="try-error") return(NA)

  return(switch(whatRet
                ,p.value = tmp$p.value
                ,statistic = tmp$statistic
                ,estimate = tmp$estimate
        )
  )

}



#
# get the first occurence of "what" in a string "string"
#
first<-function(what=NULL,string=NULL) {
  if(is.null(what) | is.null(string)) stop("check your argument passing \n")
  string = rawToChar(charToRaw(string),multiple=TRUE)
  ret = grep(what,string)
  ret = ifelse(length(ret)>0, ret[1], 0)
  return(ret)
}
#
# get the last occurence of "what" in a string "string"
#
last<-function(what=NULL,string=NULL) {
  if(is.null(what) | is.null(string)) stop("check your argument passing \n")
  string = rawToChar(charToRaw(string),multiple=TRUE)
  ret = grep(what,string)
  ret = ifelse(length(ret)>0, ret[length(ret)], 0)
  return(ret)
}


#
# loading external libraries that will be used
#
if(!flag_ys_unix) suppressMessages(library(ggplot2))
# overloading theme_bw function (only necessary for ggplot version 0.8.2)
#theme_bw <- function (base_size = 12)
#{
#    structure(list(axis.line = theme_blank(), axis.text.x = theme_text(size = base_size *
#        0.8, lineheight = 0.9, vjust = 1), axis.text.y = theme_text(size = base_size *
#        0.8, lineheight = 0.9, hjust = 1), axis.ticks = theme_segment(colour = "black",
#        size = 0.2), axis.title.x = theme_text(size = base_size,
#        vjust = 1)
#        , axis.title.y = theme_text(size = base_size,angle = 90, vjust = 0.5)
#        , axis.ticks.length = unit(0.3,"lines")
#        , axis.ticks.margin = unit(0.5, "lines")
#        , legend.background = theme_rect(colour = NA)
#        , legend.key = theme_rect(colour = "grey80")
#        , legend.key.size = unit(1.2, "lines"), legend.text = theme_text(size = base_size *
#            0.8), legend.title = theme_text(size = base_size *
#            0.8, face = "bold", hjust = 0), legend.position = "right",
#        panel.background = theme_rect(fill = "white", colour = NA),
#        panel.border = theme_rect(fill = NA, colour = "grey50"),
#        panel.grid.major = theme_line(colour = "grey90", size = 0.2),
#        panel.grid.minor = theme_line(colour = "grey98", size = 0.5),
#        panel.margin = unit(0.25, "lines"), strip.background = theme_rect(fill = "grey80",
#            colour = "grey50"), strip.label = function(variable,
#            value) value, strip.text.x = theme_text(size = base_size *
#            0.8), strip.text.y = theme_text(size = base_size *
#            0.8, angle = -90), plot.background = theme_rect(colour = NA),
#        plot.title = theme_text(size = base_size * 1.2), plot.margin = unit(c(1,
#            1, 0.5, 0.5), "lines")), class = "options")
#}



##-------- Showing all the extra & some char graphics symbols ---------
# Author
# R System
# Usage
#pchShow()
#pchShow(c("o","O","0"), cex = 2.5)
pchShow <-
  function(extras = c("*",".", "o","O","0","+","-","|","%","#"),
           cex = 3, ## good for both .Device=="postscript" and "x11"
           col = "red3", bg = "gold", coltext = "brown", cextext = 1.2,
           main = paste("plot symbols :  points (...  pch = *, cex =",
                        cex,")"))
  {
    nex <- length(extras)
    np  <- 26 + nex
    ipch <- 0:(np-1)
    k <- floor(sqrt(np))
    dd <- c(-1,1)/2
    rx <- dd + range(ix <- ipch %/% k)
    ry <- dd + range(iy <- 3 + (k-1)- ipch %% k)
    pch <- as.list(ipch) # list with integers & strings
    if(nex > 0) pch[26+ 1:nex] <- as.list(extras)
    plot(rx, ry, type="n", axes = FALSE, xlab = "", ylab = "",
         main = main)
    abline(v = ix, h = iy, col = "lightgray", lty = "dotted")
    for(i in 1:np) {
      pc <- pch[[i]]
      ## 'col' symbols with a 'bg'-colored interior (where available) :
      points(ix[i], iy[i], pch = pc, col = col, bg = bg, cex = cex)
      if(cextext > 0)
          text(ix[i] - 0.3, iy[i], pc, col = coltext, cex = cextext)
    }
  }




#
# library::llist
#   function for getting the length of each element in a list
#
# PARAMETER
#   obj   -   an object of class list
#
llist <- function (obj=NULL) {
  stopifnot(class(obj)=="list")
  tmp = unlist(lapply(obj,length))
  return(tmp)
}

# USAGE
#   aArrayList=list()
#   for (i in 1:6) aArrayList[[i]] = rnorm(i+10)
#   list2array.test = list2array(aArrayList)
#   list2array.test = list2array(aArrayList,retArray=T)
#   apply(list2array.test, 1, mean, na.rm=TRUE)
#
# DEPENDENCY
#   library::llist()
list2array <- function(obj=NULL, align="center", retArray=FALSE, filling=NA) {
  stopifnot(class(obj)=="list")
  ret = NULL
  maxlen = max(llist(obj))
  minlen = min(llist(obj))
  switch( align
         ,center = {
                    ret=  lapply( obj
                                 ,function(x){
                                              thislen = length(x)
                                              difflen = maxlen-thislen
                                              if(difflen%%2 == 0) {leftlen = difflen/2; rightlen = difflen/2;}
                                              else {leftlen=ceiling(difflen/2); rightlen=floor(difflen/2)}
                                              return(c(rep(filling,leftlen), x, rep(filling,rightlen)))
                                            }
                                )
                   }
         ,left = {
                    ret=  lapply( obj
                                 ,function(x){
                                              thislen = length(x)
                                              difflen = maxlen-thislen
                                              return(c(x, rep(filling,difflen)))
                                            }
                                )
                  }
         ,right = {
                    ret=  lapply( obj
                                 ,function(x){
                                              thislen = length(x)
                                              difflen = maxlen-thislen
                                              return(c(rep(filling,difflen),x))
                                            }
                                )
                  }
         );

  # use sapply() to get an array if so desire
  if (retArray) while(class(ret)=="list") ret = sapply(ret, function(x){return(x)})

  return (ret)

}





str2word <- function (onestring=NULL) {

  stopifnot(!is.null(onestring))
  
  ret = paste(toupper(left(onestring,1)),tolower(substring(onestring,first=2,last=10000)),sep="")

  return(ret)
}

#
# library::gettd
#   get the name of a Target directory (in contrarory to getwd()
#
# PARAMETER
#   look      -   to which direction that we should look in the file system hierarchy
#                 "up" for parent directory; "down" for child sub-directory
#   level     -   how many levels should we recursive step up or down
#   name.dir  -   vector of name(s) of directories to be searched  for the target directory
#
gettd <-function (look="up", level=1, name.dir=NULL) {

  stopifnot(!is.null(name.dir))
  stopifnot(level==length(name.dir))

  prefix=ifelse(look=="up", "..",".")
  path = NULL

  for (i in 1:level) {
    path = file.path(prefix,name.dir[i])
  }

  return (path)

}





#
# list objects in current environment with more details than ls()
# author: Tim Beissbarth, Anja von Heydebreck, Florian Markowetz, Wolfgang Huber
#
# DBG
#   property="class"
#   member="%!in%"
#   envir<-parent.frame();
#   substitute( {
#                                     fcn <-get(property, mode="function")
#                                     fcn(member)
#                                    }
#                                    , list = list(
#                                                   member=as.name(member),
#                                                   property=property,
#                                                   fcn=as.name(property)
#                                                  )
#                                 )
##      {
##          class <- get("class", mode = "function")
##          class(`%!in%`)
##      }

ll <- function() {
  df <-NULL;
  envir<-parent.frame();
  for (member in ls(envir=envir)) {
    oneRow <-list(member=member)
    
    for (property in c("class", "mode", "dim", "length")) {

      value <-eval( substitute( {
                                  fcn <-get(property, mode="function")
                                  fcn(member)
                                 }
                                 , list = list(
                                                member=as.name(member),
                                                property=property,
                                                fcn=as.name(property)
                                               )
                              ),
                    envir=envir
                  );
      if (is.null(value))
        value <-"NULL"
      else if (is.vector(value) && length(value) > 1)
        value <-sprintf("c(%s)", paste(value, collapse=","))
      else if (is.list(value))
        value <-unlist(value);

      if (length(value) > 0)
        value <-value[1];

      value <-as.character(value);

      oneRow[[property]] <- value;

    }
    
    df <-rbind(df, as.data.frame(oneRow));

  }

  df
}

#
# reverse selection of one vector that is not in another
#
"%!in%" <- function(x,y) x[!x %in% y] #--  x without y


#
# remove rows or columns of NAs from a dataframe or matrix
#
df.na <- function(df=NULL, by.row=TRUE) {

  ret = NULL
  
  by.what = ifelse(by.row, 1, 2)
  
  idx.na = which(apply(is.na(df), by.what, any)==T)
  
  if(by.what==1) ret = df[-idx.na, ]
  else ret = df[, - idx.na]
  
  return(ret)
  
}

#
# loverloaded Hmisc::sdl function
# won't be plotting 2*stdev in ggplot2
#
smean.sdl  <- function (x, mult = 1, na.rm = TRUE) {
  if (na.rm)
      x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0)
      return(c(Mean = NA, Lower = NA, Upper = NA))
  xbar <- sum(x)/n
  sd <- sqrt(sum((x - xbar)^2)/(n - 1))
  c(Mean = xbar, Lower = xbar - mult * sd, Upper = xbar + mult * sd)
}
    
#
# take top n rows from a data frame
#
top <- function (ds=NULL, by=NULL, n=1)
{ # library::top
    if(by[[1]] !="~") stop("Argument 'by' must be a one-sided formula, stupid")
    ds = sort.df(ds, by=by)
    return(ds[1:n,])
}


#
# drop unused levels in a data frame
#
# currently only works with columns
# need to add support of rows
#
k <- function (a=NULL) {
x=NULL;
for (i in 1:ncol(a)) {
  x=c(x,class(a[,i])=="factor")
}

for (j in 1:ncol(a)) { if (x[j]==TRUE) {tmp=a[,j]; tmp=tmp[drop=T]; a[,j]=tmp;}}
return(a)
}


#  lib:: randCorr
#     compute p-value of correlation statistics given two vectors
#
#  PARAMETERS
#    x, y:           two vectors of numeric values
#    n.perm:         number of permutations to be run (default=1000)
#    what.method:    what type of correlation test to be ran c("pearson","kendall","spearman")
#
#  RETURN VALUES
#    ret:            p.value based on permutation
#
#  TIME STAMP
#    2009-05-25
#
#  DBG
#    x = seq(1:10)
#    y = seq(2:11)
#    randCorr(x,y)

randCorr <- function(x, y, n.perm=1000, what.method="pearson") {

  if(length(x)!=length(y)) stop("two vectors are not equal length")

  n = length(x)
  # need to take care unix version
#  corr.obs = cor(x,y,use="everything",method=what.method)
  corr.obs = cor(x,y,method=what.method)
  cnt = 0
  for (i in 1:n.perm) {
    ix.rand = sample(1:n)
    iy.rand = sample(1:n)
    corr.rand = cor(x[ix.rand],y[iy.rand],method=what.method)
    cnt = cnt + ifelse(abs(corr.obs)>abs(corr.rand), 0, 1)
  }

  ret = cnt/n.perm

  return(ret)
}


#  lib:: mfind
#   loacate a value in a data frame or matrix
#
#  PARAMETERS
#    data:           a data frame or matrix
#    index:          index of a value to be searched in the provided data frame
#
#  RETURN VALUES
#    a list of two elements > row and column name of where value is located
#
#  TIME STAMP
#    2009-05-25
#
#  DBG
#    n = 100
#    n.c = 20
#    n.r = n/n.c
#    x = matrix(rnorm(n),ncol=n.c)
#    rownames(x)=seq(1:n.r)
#    colnames(x)=seq(1:n.c)
#    mfind(x,which.max(x))

mfind <- function (data=NULL, index=NULL) {

  retls = NULL

  if(length(index)>1) print("later")

  ix = index %% nrow(data)
  iy = (index - ix) / nrow(data) + 1

  retls = list(c(rownames(data)[ix],colnames(data)[iy]))
  
  return(retls)

}


#
# sort a data frame
#
sort.df <- function(x, by)
{ # library::sortdf

    # Author: Kevin Wright
    # with some ideas from Andy Liaw
    # http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html

    # x: A data.frame
    # by: A one-sided formula using + for ascending and - for escending
    #     Sorting is left to right in the formula

    # Useage is:
    # library(nlme);
    # data(Oats)
    # sort.df(Oats, by= ~nitro-Variety)

    if(by[[1]] != "~")
        stop("Argument 'by' must be a one-sided formula.")

    # Make the formula into character and remove spaces
    formc <- as.character(by[2])
    formc <- gsub(" ", "", formc)
    # If the first character is not + or -, add +
    if(!is.element(substring(formc, 1, 1), c("+", "-")))
        formc <- paste("+", formc, sep = "")

    # Extract the variables from the formula
    vars <- unlist(strsplit(formc, "[\\+\\-]"))
    vars <- vars[vars != ""] # Remove any extra "" terms

    # Build a list of arguments to pass to "order" function
    calllist <- list()
    pos <- 1 # Position of + or -
    for(i in 1:length(vars)){
        varsign <- substring(formc, pos, pos)
        pos <- pos + 1 + nchar(vars[i])
        if(is.factor(x[, vars[i]])){
            if(varsign == "-") {
                calllist[[i]] <- -rank(x[, vars[i]])
            } else {
                calllist[[i]] <- rank(x[, vars[i]])
            }
        } else {
            if(varsign == "-") {
                calllist[[i]] <- -x[, vars[i]]
            } else {
                calllist[[i]] <- x[,vars[i]]
            }
        }
    }
    return(x[do.call("order", calllist), ])
}

#
# lib::see
#
#   look at a small piece of a big data matrix
#
see <- function (ds=NULL, rev=FALSE, n.rows=6, n.cols=6) {
    dimr = dim(ds)[1]
    dimc = dim(ds)[2]

    if(dimr<n.rows) n.rows=dimr
    if(dimc<n.cols) n.cols=dimc
    if(!rev) print(ds[1:n.rows,1:n.cols])
    else print(ds[(dimr-n.rows):dimr,(dimc-n.cols):dimc])
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## MODULE::B-Spline Fitting
##  to perform a gene-wise B-spline fitting (gep~time) and reconstruct the
##  a profile for each gene with fitted value
##
## functions
##  1) func::bspltrfm
##  2) func::bsplfit
##
## TIME STAMP
##   2009-03-09
##
## to be complete
##   xyplot of genes
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#test = read.csv("C:\\temp\\test.csv",header=T,row.names=1)
#ixgep = c(8:258)
#ixedge = c(1:7)
#test.gep = test[,ixgep]
#gbspl = NULL
#for (i in 1:nrow(test.gep)){
#  gbspl = cbind(gbspl,bspltrfm(test.gep[i,],splannt=ianot, degfree=4))
##  colnames(gbspl)[i]=rownames(test)[i]
#}
##colnames(gbspl) = igene[match(rownames(test.gep),rownames(igene)),'Symbol']
## do NOT use symbol - not reliable; some symbols might be NULL
#colnames(gbspl) = rownames(igene)[match(rownames(test.gep),rownames(igene))]
#gbspl=t(gbspl)
## append original annotation and EDGE statistics back
#ix = match(rownames(gbspl),rownames(test))
#gbspl = cbind(test[ix,ixedge],gbspl)
#
#
##
## func::bspltrfm
##   prepare data and calls func::bsplfit to fit a bspline
##
## PARAMETERS
##   agene       -     a vector of numeric values
##   splannt     -     annotation files for the name of each element in agene
##
## RETURN VALUE
##
## DEPENDENCIES
##   func::bsplfit
##
## TIME STAMP
##   2009-05-09
##
#bspltrfm <- function(agene=NULL, splannt=ianot, degfree=4){
## 2009-03-08
## PROBLEM OF APPLY CALL NEED TO BE FIXED
##  if(any(is.null(splannt),is.null(agene))) stop("Input sample annotation file and GEP for a gene")
#  ix = match(colnames(agene),rownames(splannt))
#  # THE FOLLOWING LINE KEEPS THIS FUNCTION TO BE CALLED BY APPLY
#  tobspl = cbind(t(agene),splannt[ix,c('Time','PhenoNew')])
#  tobspl$PhenoNew = as.factor(tobspl$PhenoNew)
#  ret=bsplfit(agene=tobspl,ixgep=1,degfree=degfree)
#  return(ret)
#}
#
#
#
##
## func::bsplfit
##   function for a) fitting b-spline b) obtain fitted values
##
## PARAMETER
##   a dataframe containing one gene's gep values with at least three fields
##   1) Time     the x values for spline fitting
##   2) GEP      (column name of genename) gep intensities for at least one gene
##   3) PhenoNew: class labels for phenotype
##
## RETURN VALUE
##   agene.bspltrfm - a vector of b-spline fitted values
##
## TIME STAMP
##   2009-03-09
##
#bsplfit <- function (agene=NULL, ixgep=1, degfree=4) {
# agene.bsplfit = NULL
# test = agene
# gname = colnames(test)[ixgep]
## print(gname)
#
# # Splits the data into subsets,
# # computes summary statistics for each,
# # and returns the result in a convenient form.
# #???????????????????????????????
# # don't know why the name of the summarized values
# # are not automatically assign, an 'x' is given at
# # run time instead.
# # need to figure it out later
# #???????????????????????????????
# test.avg = aggregate(x=test[,ixgep], by=list(Time=test$Time, PhenoNew=test$PhenoNew), FUN="mean")
# colnames(test.avg)[which(colnames(test.avg)=="x")]=gname
# tmp = unique(test.avg$PhenoNew)
# for (i in 1:length(tmp)) {
#  apheno.avg = test.avg[test.avg$PhenoNew==tmp[i],]
#  # fitting a B-spline
#  apheno.spl = smooth.spline(x=apheno.avg$Time, y=apheno.avg[,gname], df=degfree)
###  plot(apheno.avg$Time, apheno.avg[,gname])
###  lines(apheno.spl,col="blue")
###  lines(pp <- predict(apheno.spl,apheno.avg$Time), col = "blue")
###  points(apheno.avg$Time, pp$y, pch=3, col="dark red")
#  sppred = predict(object=apheno.spl,x=apheno.avg$Time)
#  tmpfit = sppred$y
#  names(tmpfit)= paste(tmp[i],"H",apheno.avg$Time,sep="")
#  agene.bsplfit = c(agene.bsplfit,tmpfit)
##  agene.bsplfit = c(agene.bsplfit,sppred$y)
#  # name the rows as "AH26" - Asx Hour 26
##  names(agene.bsplfit)=paste("H",tmp,apheno.avg$Time,sep="")
# }
#
# return (agene.bsplfit)
#}
#
#

############################################################################
# htslib::isOverlap
#   check if two vectors are overlapping with each other
#
# PARAMETERS
#   avec, bvec  - two vectors to be checked
#
# TIME STAMP
#   2009-04-06
############################################################################

isOverlap <- function(avec=NULL,bvec=NULL) {
  return(ifelse(length(intersect(avec,bvec))==0,FALSE,TRUE))
}

#
# permCor
#   correlation test with permutation for p-value
#
# PARAMETERS
#   avec:   first vector of numbers
#   bvec:   second vector of numbers
#   n.perm: # of permutations the vector should be shuffled
#
# RETURN VALUES
#   a list of two values
#     > perm.cor:   correlation coefficient with
#     > p-value:    defined as p-value = sum{(# cor.rand) >= (# cor.obs)} / n.perm
#
permCor <- function (avec=NULL,bvec=NULL,n.perm=100) {

  n.false=0
  len = length(avec)

  rol.obs = cor(avec,bvec)
  rol.rand = NULL
  for (i in 1:n.perm) {
    rol.rand = cor(sample(avec,len),sample(bvec,len))
    if (rol.rand >= rol.obs) n.false = n.false+1
  }
  pval = n.false/n.perm
  return(list(c(rol.obs,pval)))
  
}




##
##
## lib::getCentrd
##   compute class mean (centroid) based on class labels
##
## PARAMETERS
##   ds:    data frame for which function to be applied on
##   ixclass:   index of column that cotains class labels
##   what:  "mean"/"lower"/"upper", only used if "ysmean.cl.boot" method is specified
##
## RETURN VALUE
##   ds.centroid: data frame of centroids of each class
##
##

getCentrd <- function (ds=NULL,ixclass=1,fun="mean",what=NULL, isRow=FALSE) {

  if (fun=="ysmean.cl.boot") {
    # specify to compute mean with confidence interval using bootstrap method
    ds.centroid = aggregate(x=ds[,-ixclass],by=list(class=ds[,ixclass]), FUN=fun, ret.what=what)
  }
  else {
    ds.centroid = aggregate(x=ds[,-ixclass],by=list(class=ds[,ixclass]), FUN=fun)
  }
  
  return(ds.centroid)
}


##
##
## lib::ysmean.cl.boot
##   a modification of Hmisc::smean.cl.boot - change return value from a vector to a scalar
##   so that this function can be used with aggregate
##
## PARAMETERS
##   x:    an input vector
##   conf.int:    desired confidence interval to be computed
##   B:           number of bootstraps
##   ret.what:    what values to be returned >> "mean", "lower", "upper"
##
## DEPENDENCIES
##   {Hmisc}
##
## RETURN VALUE
##   res:         the exact return value is dependent on parameter ret.what
##
## USAGE
#    x = rnorm(100)
#    ysmean.cl.boot(x,ret.what="mean")
#    ysmean.cl.boot(x,ret.what="lower")
#    ysmean.cl.boot(x,ret.what="upper")

ysmean.cl.boot <- function (x, conf.int = 0.95, B = 1000, ret.what=NULL, na.rm = TRUE, reps = FALSE) {
  # lib::ysmean.cl.boot
  
    if (is.null(ret.what)) stop("please specify what you are looking for: mean/lower/upper")
  
    library(Hmisc)
    if (na.rm)
        x <- x[!is.na(x)]
    n <- length(x)
    xbar <- mean(x)
    if (n < 2)
        return(c(Mean = xbar, Lower = NA, Upper = NA))
    z <- unlist(lapply(1:B, function(i, x, N) sum(x[.Internal(sample(N,
        N, TRUE, NULL))]), x = x, N = n))/n
    quant <- quantile(z, c((1 - conf.int)/2, (1 + conf.int)/2))
    names(quant) <- NULL
    # ys comment
    # res <- c(Mean = xbar, Lower = quant[1], Upper = quant[2])
    if(ret.what=="mean") res = xbar
    else if (ret.what=="lower") res=quant[1]
    else if (ret.what=="upper") res=quant[2]

    if (reps)
        attr(res, "reps") <- z
    res
}


##
## distEucl
##   compute Eucledian distance (square it to get squared Eucledian distance)
##
## PARAMETERS
##   avec:    first vector
##   bvec:    second vector
##
## RETURN VALUE
##   ret:  a scalar form of Eculedian distance between two vectors
##
##

distEucl <- function (avec, bvec) {

   diffvec = avec - bvec
   ret = (sqrt(sum(diffvec * diffvec)))

   return (ret)
   
}

##
## tscale
##   transposed version of scale function
##
## PARAMETERS
##   ds:      a data frame upon which scale should be performed ROW-WISE
##   ixcol:   colname index that subsetting the data
##
## RETURN VALUE
##   ret:  a matrix of row-wise scaled values with same dimension as input
##
## DEBUG
##

tscale <- function (ds=NULL, ixcol=NULL) {

  ret = NULL

  y = ds[,ixcol]
  y = t(y)
  y.scale = scale(y)
  ret = t(y.scale)
  return (ret)
}



# get left-most part of a string
left <- function (x=NULL,n=0){
  return(substring(x,1,n))
}

# get right-most part of a string
right <- function (x=NULL,n=0){
  slen = nchar(x)
  return(substring(x,slen-n+1,slen))
}


#
# cstr
#   function for concatenating multiple strings with no white space in between
#
# PARAMETERS
#   avec: a vector of strings to be concatenated
#
# RETURN VALUE
#   ret:  a string that is concatenation of all strings in the input vector
#
# DEBUG
#   x = "kiss"
#   y = "ass"
#   z = "yes"
#   a = cstr(c(x,y,z))
#   print(cstr(c(x,y,z)))
#   [1] "kissassyes"
#
cstr <- function (avec = NULL) {
  ret = NULL
  for (i in 1:length(avec))
    ret = paste(ret,avec[i],sep="")

  return(ret)
}

scat <- function (...) {
   tmp = list(...)
   tmp = sapply(tmp, StrTrim)
   ret = paste(tmp, collapse="")
   return(ret)
}


#
# findClosest
#   function to find the element in a vector that
#   is the closest to a reference element indicated 
#   by its index in that same vector
#
# PARAMETERS
#    avec:    a vector of numbers to be searched
#    aref:    the reference number
# 
# RETURN VALUE
#   the element that was found to be closest to reference
#
findClosest <- function (avec = NULL, arefpos=NULL) {
    #arefpos=1
    #avec = c(36370211, 90,91,163,172)
    #aref = round(199/2)      
    pos = 0                   
    adist = 9999999
    avec = as.numeric(avec)
    #aref = round(avec[arefpos]/2,0)     # center of the sequence
    aref = avec[arefpos]
    avec = avec[-arefpos]

    for (i in 1:length(avec)) {
      if (abs(avec[i]-aref)<adist) {
        pos = i
        adist=abs(avec[i]-aref)
      }
    }
    # return 0 if the search vector is 0 on all elements
    return (ifelse(pos==0,0,avec[pos]))
}



# RMSE function
rmse <- function(x,y) sqrt(mean((x-y)^2))
 

#
# mtxC2N
#   function to convert a matrix with characters into numeric
#
# PARAMETERS
#    mtx:    a matrix with elements as characters
#
# RETURN VALUE
#   a matrix with numeric elements
#
mtxC2N <- function (mtx=NULL) {
  return (mtx = matrix(as.numeric(mtx),nrow=nrow(mtx),ncol=ncol(mtx),byrow=TRUE,dimnames=dimnames(mtx)))
}


##
## fillSqrMatrix
##   function to fill in square matrix 
##
## PARAMETERS
##    sqrmtx:    a square matrix (symmetric)
##    upper:    1 - upper triangular; 0 - lower triangular
##
##
##  note (2008-12-17)
##  two more functionality need to be implemented
##  1. add function to flip a triangular matrix, ie, upper->lower and lower->upper
##  2. add function to fill in lower triangular matrix (currently works only for upper triangular)
fillSqrMatrix <- function (sqrmtx = NULL, upper=1, val=-1) {

#  ifelse(upper %in% c(0,1)), next(), print('need to know if upper or lower triangular')

  x = sqrmtx

  if (val!=(-1)) {
    idx = NULL
    if(upper==1) {idx=upper.tri(x)} 
    else {idx=lower.tri(x)}
    x[idx] = val
    return (x)
  }

  for (ixcol in 1:ncol(sqrmtx))
    for (ixrow in 1:nrow(sqrmtx))
      #   sqrmtx[j][i] = sqrmtx[i][j]
      #   if(ixcol<ixrow) x[ixrow,ixcol] = sqrmtx[ixcol,ixrow]
      if (upper==1 & ixcol<ixrow) x[ixrow,ixcol] = sqrmtx[ixcol,ixrow]   
      else if (upper==0 & ixcol>ixrow) x[ixrow,ixcol] = sqrmtx[ixcol,ixrow]   

      
  return (x)
}

# comment this on 12/22/2008
fillSqrMatrix.old <- function (sqrmtx = NULL, upper=1, val=-1) {

  x = sqrmtx

  if (val!=(-1)) {
    for (ixcol in 1:ncol(sqrmtx))
      for (ixrow in 1:nrow(sqrmtx))
        #   sqrmtx[j][i] = sqrmtx[i][j]
        #   if(ixcol<ixrow) x[ixrow,ixcol] = sqrmtx[ixcol,ixrow]
        if (upper==1 & ixcol<ixrow) x[ixrow,ixcol] = val
    return (x)
  }

  for (ixcol in 1:ncol(sqrmtx))
    for (ixrow in 1:nrow(sqrmtx))
      #   sqrmtx[j][i] = sqrmtx[i][j]
      #   if(ixcol<ixrow) x[ixrow,ixcol] = sqrmtx[ixcol,ixrow]
      if (upper==1 & ixcol<ixrow) x[ixrow,ixcol] = sqrmtx[ixcol,ixrow]   
  return (x)
}

##
## mysetwd
##   function to set working directory
##
mysetwd <- function (strpath=NULL,strname=NULL){
  path = NULL
  aname = NULL
  newdirname = NULL
  
  msg = NULL
  
  if(is.null(strname)) {
    print("Please specify name of the directory to be created");
    msg = paste("Current directory is: ", getwd(), sep="")
    print(msg)
    return (FALSE)
  }
  
  path = ifelse(is.null(strpath),"C:\\temp\\r.working\\", strpath)

  newdirname = paste(path,Sys.Date(),"\\.",strname,sep="")

  dir.create(newdirname,recursive=TRUE)
  
  setwd(newdirname)
 
  msg = paste("Working directory is set to: ", getwd(), sep="")
  print(msg)
  return (TRUE)
}


mynorm <- function (avector = NULL) {
  ret = sqrt(sum(avector^2))  
  return (ret)  
}


# function 
# unnorcorrcoef = sum(aPeakReads * aConservationScore) / (sqrt(sum(aPeakReads^2))*sqrt(sum(aConservationScore^2)))      

# a wrapper to time the code
mytimer <- function (afunction=NULL) {
  if (is.null(afunction)) print("USAGE: mytimer(a function)")
  print(system.time(afunction))
}

# function for sort a matrix by specific column or by its row name
# use both sort and sort.list just to demonstrate how to use each
# individual one

# the major difference is that sort only works for matrix

msort <- function (amatrix, acolumn){ 
  # sort by the row name
  if (acolumn==0)    
    amatrix[sort.list(rownames(amatrix)),] 
  else 
    # amatrix[sort(amatrix[,acolumn],index.return=T)$ix,] 
    amatrix[sort.list(amatrix[,acolumn]),]
}


# stack based on regexpression pattern
# x = grep("pval",colnames(test))
# not completed yet and need to add two more things
# a. stack on rows
# b. stack based on regexp 
ArrayStack <- function (aFrame, vIdx, diM=2) {
  aRet = NULL
  k = 0
  if (diM==2) { 
    for (i in x) {
      k = k + 1
      aRet = rbind(aRet,cbind(colnames(aFrame)[i],rep(k,length(aFrame[,i])),aFrame[,i]))
    }
  }
  return (aRet)
}  



#
# search a substring within another string
#
StrSearch <- function(aPattern, aString, isCaseSensitive=T) {
  idx = 1
  for (idx in 1:nchar(aString)) {
    if(is.na(grep(aPattern, substring(aString,idx,nchar(aString)))[1])==T) break;
  }
  return (StrSearch = idx -1)
}

StrSearch2 <- function(aPattern, aString, isCaseSensitive=T, what="last") {
  idx = 1
  len = nchar(aString)
  pos = 0
  if(what=="last") {
    for (idx in 1:nchar(aString)) {
      if(is.na(grep(aPattern, substring(aString,idx,nchar(aString)))[1])==T) return(idx-1);
    }
  }

  if(what=="all") {
    tmp = StrSearch2(aPattern,substr(aString,idx,len),what="last")
    pos = tmp
    while(tmp!=0) {
        idx = tmp - 1
        tmp = StrSearch2(aPattern,substr(aString,1,idx),what="last")
        if(tmp!=0) pos = c(pos, tmp)
    }
    pos = rev(pos)
  }

  return(pos)

}

#
# strgsub
#   replace a substring inside another string
strgsub <- function(sStart, sStop, replacement="", x){

  iStart = StrSearch2(sStart, x, what="last")
  iStop =  StrSearch2(sStop, x, what="all")
  ret = gsub(substr(x, iStart, iStop[which(iStop>iStart)[1]]), replacement,x)
  if(is.na(ret)==TRUE) ret = x

  return(ret)

}

# search a string pattern in an array of strings
# return the index of the element in the array where patter is found
# ret: -1 if not found
#      >0 otherwise
StrSearchArray <- function(aPattern, aStringArray, isCaseSensitive=T) {
  idx = -1
  for (i in 1:length(aStringArray)) {
    if (length(grep(aPattern, aStringArray[i]))!=0) {
      idx = i
      break
    }
  }
  return (StrSearch = idx)
}

# remove leading and trailing white spaces
StrTrim <- function (aString) {
 return (gsub("^[[:space:]]+|[[:space:]]+$","",aString))
}


# convert a list into a vector
# by default, stack() returns a dataframe, we convert it into a matrix for easy manipulation
List2Matrix <- function (aList,keep.listname=FALSE,newMethod=FALSE) {
  
  ret = NULL

  nelements = length(aList)
  
  # assign name for each element in the list if there is no names
  # otherwise, the stack function will complain
  if(is.null(names(aList))) names(aList) = seq(1:nelements)

  if (newMethod) {
    ret = matrix(unlist(aList),nrow=nelements,byrow=TRUE)
    return (ret)
  }

  if (keep.listname) {return (as.matrix(stack(aList)))}
  else {as.matrix(stack(aList)[,-2,drop=FALSE])}
}


#
# func::drop.levels
#   Drop unused factor levels from all factors in a data.frame
#   Author: Kevin Wright.  Idea by Brian Ripley.
#
drop.levels <- function(dat){
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}



###########################################################################
# library::scattersmooth
#
#   Author: P. H. C. Eilers and J. J. Goeman
#
#   A function for better visualisation of scatterplots involving thousands
#   or millions of dots.
#
#   Scatterplots of microarray data generally contain a very large number
#   of dots, making it difficult to get a good impression of their
#   distribution in dense areas. We present a fast and simple algorithm
#   for two-dimensional histogram smoothing, to visually enhance scatterplots.
#
############################################################################
scattersmooth <- function(x,y, nbin = 200, lambda = 10, ndot = 500, ...)
{

    fillhist <- function(xb, yb, nb)
    {
        H <- matrix(rep(0, prod(nb)), nb[1], nb[2])
        for (i in 1:length(xb))
        {
            H[xb[i],yb[i]] <- H[xb[i],yb[i]] + 1
        }
        H
    }

    # check correct input
    if (length(x) != length(y))
        stop("lengths of x and y do not match")
    if ( (length(x) < 2) | (length(y) < 2) )
        stop("x and y should be vectors")
    if ( !is.numeric(x) | !is.numeric(y) )
        stop("x and y should contain numeric values")
    if ( all(!is.numeric(lambda)) | all(lambda < 0) | (length(lambda) > 2) )
        stop("lambda should be numeric and positive")
    if ( length(lambda) == 1 )
        lambda <- c(lambda, lambda)
    if (all(!is.numeric(nbin)) | all(nbin < 1) | (length(nbin) > 2) )
        stop("nbin should be a strictly positive integer")
    if ( length(nbin) == 1 )
        nbin <- c(nbin, nbin)
    if (!is.numeric(ndot) | (ndot < 0) | (length(ndot) > 1) )
        stop("ndot should be a strictly positive integer")
    ndot <- floor(ndot)
    m <- length(x)

    # Put the x-values into bins
    xmin <- min(x)
    xmax <- max(x)
    dx <- (xmax - xmin) / (nbin[1] - 1)
    xbin <- floor(1 + (x - xmin) / dx)
    xscale <- xmin + (1:nbin[1] - 0.5) * dx

    # Put the y-values into bins
    ymin <- min(y)
    ymax <- max(y)
    dy <- (ymax - ymin) / (nbin[2]-1)
    ybin <- floor(1 + (y - ymin) / dy)
    yscale <- ymin + (1:nbin[2] - 0.5) * dy

    # Create the unsmoothed histogram
    H <- fillhist(xbin, ybin, nbin)

    # Calculate the smoothing matrix
    D1x <- diff(diag(nbin[1]))
    D1y <- diff(diag(nbin[2]))
    D2x <- diff(D1x)
    D2y <- diff(D1y)
    Qx <- diag(nbin[1]) + lambda[1]^2 * t(D2x) %*% D2x + 2 * lambda[1] * t(D1x) %*% D1x
    Qy <- diag(nbin[2]) + lambda[2]^2 * t(D2y) %*% D2y + 2 * lambda[2] * t(D1y) %*% D1y

    # Smooth
    H <- t(solve(Qy, t(solve(Qx, H))))

    # Plot coloured image
    image(x = xscale, y = yscale, z = -H, xlab = "", ylab = "", col = heat.colors(100), ...)
    #image(x = xscale, y = yscale, z = -H, xlab = "", ylab = "", col = gray(100), ...)

    # Plot selection of dots
    if (ndot > 0) {
      ndot <- min(m, ndot)
      sel <- sort.list(rnorm(m))[1:ndot]
      points(x[sel], y[sel], cex = 0.1)
    }
}


cat(">>>>>>> Yongsheng's library has been successfully loaded!\n")