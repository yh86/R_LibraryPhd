###############################################################################
# phdpath                                                                     #
#   library for temporal pathway analysis for PHD project                     #
#                                                                             #
# AUTHOR                                                                      #
#   Yongsheng Huang   {huangys}@eecs.umich.edu                                #
#   Ph.D Candidate                                                            #
#   Bioinformatics Graduate Program and Department of Statistics              #
#   University of Michigan                                                    #
#   2017 Palmer Commons, 100 Washtenaw Ave, Ann Arbor, MI 48109-2218          #
#   http://www-personal.umich.edu/~huangys                                    #
#   http://www.hyperfocal.org                                                 #
#                                                                             #
# DATE (Last Updated)                                                         #
#   2010-08-10                                                                #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#dbgKid="04620"; print(getKeggName(dbgKid))
getKeggName <- function(kid=NULL,n=5) {

  require(KEGG.db)

  if(class(kid)!="list") kid = as.list(kid)
  kid = lapply(kid, function(x){right(paste(paste(rep(0,n),collapse=""),x,sep=""),n)})
  # apparantly, the mget function takes care of input of a list
  # therefor, if kid is a list, it's fine without explicit casting
  ret <- mget(as.character(kid), KEGGPATHID2NAME, ifnotfound=list(NA))
  if(class(ret)=="list") ret = unlist(ret)
  
  return(ret)
  
}

#print(unifyKeggId(kid))
unifyKeggId <- function(kid=NULL,n=5){
  if(class(kid)!="list") kid = as.list(kid)
  ret = unlist(lapply(kid, function(x){right(paste(paste(rep(0,n),collapse=""),x,sep=""),n)}))
  return(ret)
}


#***********************************************************************************************
# library::function
#   given a vector pathway indices, plot pair-wise (e.g., baseline vs T) plot for each component gene
#
# AUTHOR
#   Yongsheng Huang
#   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys; http://www.hyperfocal.org
#
# PAMETERS
#   list.what       -   list of samples to use
#   gmt.idxr        -   indices of pathways
#   gmt.broad       -   pathway databases (e.g., Broad C2)
#   prefix          -   prefix in the output files
#   save.data       -   c(TRUE, FALSE); whether to save the expression data
#
# VALUE
#   None
#
# DEPENDENCIES
#   1) phdgraphics::ggPair
#
# DBG
#  dgep = loadData("gep")
#  fname = "Flu.H00vsH05.csv"
#  list.what = loadLW(fname)
#  gmt.idxr = c(606, 1306, 1587, 1641, 444, 914, 1413, 1423, 1672, 1725, 1730, 1735, 1740, 1769)
#  gmt.broad = loadBroadGMT()
#
# NOTES
#   1) In the future, we can extend this to more than pairwise plot (with a function name)
#
# USAGE
#  dgep = loadData("gep")
#  fname = "Flu.H00vsH05.csv"
#  list.what = loadLW(fname)
#  gmt.idxr = c(606, 1306, 1587, 1641, 444, 914, 1413, 1423, 1672, 1725, 1730, 1735, 1740, 1769)
#  gmt.broad = loadBroadGMT()
#  plotSigPathwys(list.what=list.what, list.data=dgep, gmt.broad=gmt.broad, gmt.idxr=gmt.idxr, prefix=fname, save.data=TRUE)
#
# KNOWN ISSUES
#
#
# TIME STAMP
#   2009-08-11
#***********************************************************************************************
plotSigPathwys <- function (list.what=NULL, list.data=NULL, gmt.broad=NULL, gmt.idxr=NULL, prefix=fname, save.data=TRUE) {

  pathways = getBroadGenes(gmt.idxr=gmt.idxr, gmt.broad=gmt.broad, what="entrez")
  pathways = lapply(pathways, function(x){return(paste(x, "_at", sep=""))})

  for (i in 1:length(pathways)) {

    tmp = getFeatures(pathways[[i]], annot="Symbol",objannot=list.data$features)
    idx.na = which(is.na(tmp)==T)
    if(length(idx.na)!=0) {
      list.what$features=pathways[[i]][-idx.na]
      list.what$symbol=tmp[-idx.na]
    }
    else {
      list.what$features=pathways[[i]]
      list.what$symbol=tmp
    }

    # remove the duplicates - may need more elegant solution than this
    list.what$features = unique(list.what$features)
    list.what$symbol = unique(list.what$symbol)
    
    test = getDesign(list.what=list.what, list.data=list.data)
    genes = unique(rownames(test))
    n.genes = length(genes)

    prefix = gsub(".csv", "", prefix)
    output.name = paste(prefix, ".", names(pathways)[i], ".ggPair",sep="")
    # save a copy of the data
    if(save.data) write.csv(test, paste(output.name,".csv",sep=""))
    
    pdf (paste(output.name,".pdf",sep=""))
    for (j in 1:n.genes) {
      a.gene.name = getFeatures(genes[j], annot='Symbol', objannot=list.data$features)
      a.gene = data.frame(test[genes[j],])
      a.gene = data.frame(t(a.gene))
      a.gene$obs = getPheno(rownames(a.gene), annot="Subject", objannot=list.data$obs)
      a.gene$pheno = getPheno(rownames(a.gene), annot="PhenoNew", objannot=list.data$obs)
      a.gene$timepoints = getPheno(rownames(a.gene), annot="TimePoints", objannot=list.data$obs)
      colnames(a.gene)[1] = a.gene.name

      a.gene = reshape ( a.gene
                            ,v.name=a.gene.name
                            ,idvar=c("obs")
                            ,timevar=c("timepoints")
                            ,direction="wide"
                          )

      tmp = ggPair(ds=a.gene, obs="obs", group="pheno",vars=c(3,4), main.label=a.gene.name)
      rm(tmp)
    }
    dev.off()
  }

}



#***********************************************************************************************
# library::function
#   what does this library do
#
# AUTHOR
#   Yongsheng Huang
#   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PAMETERS
#   gmt.idxr    -   row index of pathways in the GMT (Broad GSEA) file, e.g., output from SAM GSA analysis
#   what        -   c("entrez","symbol"); what type of gene information to be returned
#
# VALUE
#   ret         -   return value of this function
#
# DEPENDENCIES
#   1) condition A
#   2) condition B
#
# DBG
#   setting parmeters for running the function
#
# NOTES
#
# USAGE
#    # load Broad molecular database
#    gmt.broad = loadBroadGMT()
#    # define which pathways to look at (with a vector of numbers indicating the indices)
#    gmt.idxr = c(1413,1423)
#    getBroadGenes(gmt.idxr=gmt.idxr, gmt.broad=gmt.broad, what="symbol")
#    getBroadGenes(gmt.idxr=gmt.idxr, gmt.broad=gmt.broad, what="entrez")
#
# KNOWN ISSUES
#   1) some genes may not be mapped due to the nature of gene nomenclature
#      for instance, C10ORF7 can't be found, but C10orf7 is valid.
#      nothing we can do now and hopefully the annotation becomes better
#
# TIME STAMP
#   2009-08-08
#
#***********************************************************************************************
getBroadGenes <- function(gmt.idxr=NULL, gmt.broad=NULL, what="entrez") {

  if(is.null(gmt.idxr)) stop("\n >>> need row index of gmt file, you idiot")
  if(is.null(gmt.broad)) stop("\n >>> give me the broad database (.gmt), you stupid")

  suppressMessages(require(GSEABase,quietly=TRUE))

  # get the pathways based on row index of pathways
  pathwy = gmt.broad[gmt.idxr]
  pathwy.names = unlist(lapply(pathwy, setName))
  pathwy.geneIds = lapply(pathwy, geneIds)

  if(what=="entrez") {
    suppressMessages(require(org.Hs.eg.db))
    #   mget("C10orf7", org.Hs.egALIAS2EG, ifnotfound=NA)
    #   mget("KNTC2", org.Hs.egALIAS2EG)
    #   mget(msba, revmap(org.Hs.egSYMBOL))
    pathwy.geneIds = lapply(  pathwy.geneIds
                            , function(x){ x = gsub("ORF", "orf", x);
                                          return( unlist(mget(x, envir=org.Hs.egALIAS2EG, ifnotfound=NA)) )
                                         }
                            )
  }

  # assign pathway names
  names(pathwy.geneIds) = pathwy.names

  return (pathwy.geneIds)

}

#***********************************************************************************************
# phdpath::loadBroadGMT
#   load .gmt molecular database (Broad GSEA)
#
# AUTHOR
#   Yongsheng Huang
#   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PAMETERS
#   gmt.broad     -     molecular signature data base (Broad GSEA), e.g., c1, c2, etc.
#                       if null, will default to "c2.all.v2.5.symbols.gmt"
#
# VALUE
#   ret     -   an object contains a collection of Broad molecular database
#
# DEPENDENCIES
#   1) {GSEABase}
#
# DBG
#   setting parmeters for running the function
#
# NOTES
#
# USAGE
#   gmt.broad = loadBroadGMT()
#
# KNOWN ISSUES
#
#
# TIME STAMP
#   2009-08-13
#
#***********************************************************************************************
loadBroadGMT <- function(gmt.broad=NULL) {

  if (is.null(gmt.broad)) {
      cat("No Broad database location specified, load the c2 database automatically \n")
      gmt.broad = "C:\\Users\\yongsheng\\Documents\\Yongsheng's Documents\\Develop\\R\\database\\c2.all.v2.5.symbols.gmt"
  }

  suppressMessages(require(GSEABase, quietly=TRUE))
  
  # read in the GMT database
  ret = getGmt(gmt.broad)

  show(ret)

  return(ret)
}


#
# globaltest expects usual statistical data-format (row by column = subjects by covariates)
#
# we change this into the common genomic dataset format
#

#***********************************************************************************************
# phdpath::pathgt
#    pathway group testing function for performing pathway analysis
#    at multiple time points, independently
#
# AUTHOR
#   Yongsheng Huang
#   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PARAMETERS
#   para_1            kdajf
#
# DEPENDENCIES
#   1) if Broad database to be used but not specified, assuming its location at
#     "C:\\Users\\yongsheng\\Documents\\Yongsheng's Documents\\Develop\\R\\database\\msigdb_v2.5.xml"
#   2) assuming data is in typical statistical data format (row by column = obs by genes)
#      NOTE that globaltest assumes common genomic data format (row by col = genes by observations)
#      Use gt.options(transpose=TRUE/FALSE) to change default, as shown in the function
#   3) library("org.Hs.eg.db")
#   4) library(globaltest)
#   5) assuming gene names are automatically changed by R (prefixing "100_at" with an "X")
#   6) functions {gtKEGG}
#
# VALUE
#   ret       a list of three elements
#             - a list of pathway analysis results;
#             - gene symbol that can be used in covariates() plot
#             - and (for validation purpose) a list of samples selected at each time point
#
# USAGE
#   dgep = loadData("gep")
#   list.what = loadLW("flu.list.what.csv")
#   test = data.frame(t(getDesign(list.data=dgep, list.what=list.what)))
#   test.KEGG = pathgt (ds=test, dgep=dgep)
#   test.Broad = pathgt (ds=test, dgep=dgep, db.what="Broad", category.broad=c("c2"))
#
# DBG
#     n.permutation = 10000
#     db.what = "Broad"
#     xml.broad = NULL
#
# NOTES
#   1) when running Broad database analysis, it's always a good idea to run only a category at a time
#      otherwise, the number of pathways identified can be huge and less useful to problem at hand
#       c1 (n=386)   positional gene sets
#       c2	(n=1892)  curated gene sets
#       c3	(n=837)   motif gene sets
#       c4 (n=883)   computational gene sets
#       c5	(n=1454)  GO gene sets
#
# KNOWN ISSUES
#
#
# TIME
#   2010-11-02
#===================================================================================================
pathgt <- function (ds=NULL, dgep=NULL, n.permutation=10000, entrez=NULL, db.what="KEGG", xml.broad=NULL, category.broad=c("c2"),baseline=FALSE) {

  require("org.Hs.eg.db") || stop("please install org.Hs.eg.db first")
  require(globaltest) || stop("please install globaltest first")

  gt.options(transpose=T) # critical for matrix as data input
  
  if(is.null(entrez)) gene.mapping = as.list(gsub("_at", "", substring(colnames(ds),2,100)))
  else gene.mapping = entrez                # 2010-03-08 >> if probeset ids are not "entrez_at", we need a list of entrez genes
  names(gene.mapping) = colnames(ds)

  genes.symbol = gsub("^X", "", colnames(ds))
  genes.symbol = getFeatures(genes.symbol,annot='Symbol',objannot=dgep$features)

  test.pheno = dgep$obs[rownames(ds),]
  v.timepoints = sort(unique(test.pheno$TimePoints))
  n.timepoints = length(v.timepoints)

  ret = list()
  sample.keep = list()

  if(baseline){ # contrast each timepct versus baselines

    for (i in 2:n.timepoints){
  
      tmp =  test.pheno[test.pheno$TimePoints == v.timepoints[1] | test.pheno$TimePoints == v.timepoints[i], ]

      a.pheno = as.factor(test.pheno[rownames(tmp),]$TimePoints)
      a.data = as.matrix(ds[rownames(tmp),])
      sample.keep[[i]] = cbind(samples=rownames(tmp),phenotype=a.pheno)
      rm(tmp)
  
      tmp = NULL
  
  
      if(db.what=="Broad") {
  
        require(GSEABase) || stop("get package{GSEABase} installed first!")
  
        if (is.null(xml.broad)) xml.broad = "C:\\YS_Documents\\Develop\\R\\database\\msigdb_v2.5.xml"
        broad = getBroadSets(xml.broad)
  
        a.Broad = gtBroad (   a.pheno
                            , t(a.data)
                            , collection=broad
                            , category = category.broad
                            , probe2entrez=gene.mapping
                            , annotation="org.Hs.eg.db"
                            , permutations=n.permutation
                            )
        tmp = a.Broad
      }
  
      if (db.what=="KEGG") {
        a.KEGG = gtKEGG(  a.pheno
                        , t(a.data)
  #                      , alias = gene.symbol
                        , probe2entrez=gene.mapping, annotation="org.Hs.eg.db", permutations = n.permutation)
        tmp = a.KEGG
      }
  
      if (db.what=="GO") {
        cat("GO not implemeted yet")
        tmp = NULl
      }
  
      ret[i-1] = tmp
      names(ret)[i-1] = paste(db.what, v.timepoints[i], sep="")
    }  
  }




  else{
    for (i in 1:n.timepoints){
  
      tmp =  test.pheno[test.pheno$TimePoints == v.timepoints[i], ]
      a.pheno = as.factor(test.pheno[rownames(tmp),]$PhenoNew)
      a.data = as.matrix(ds[rownames(tmp),])
      sample.keep[[i]] = cbind(samples=rownames(tmp),phenotype=a.pheno)
      rm(tmp)
  
      tmp = NULL
  
  
      if(db.what=="Broad") {
  
        require(GSEABase) || stop("get package{GSEABase} installed first!")
  
        if (is.null(xml.broad)) xml.broad = "C:\\YS_Documents\\Develop\\R\\database\\msigdb_v2.5.xml"
        broad = getBroadSets(xml.broad)
  
        a.Broad = gtBroad (   a.pheno
                            , t(a.data)
                            , collection=broad
                            , category = category.broad
                            , probe2entrez=gene.mapping
                            , annotation="org.Hs.eg.db"
                            , permutations=n.permutation
                            )
        tmp = a.Broad
      }
  
      if (db.what=="KEGG") {
        a.KEGG = gtKEGG(  a.pheno
                        , t(a.data)
  #                      , alias = gene.symbol
                        , probe2entrez=gene.mapping, annotation="org.Hs.eg.db", permutations = n.permutation)
        tmp = a.KEGG
      }
  
      if (db.what=="GO") {
        cat("GO not implemeted yet")
        tmp = NULl
      }
  
      ret[i] = tmp
      names(ret)[i] = paste(db.what, v.timepoints[i], sep="")
    }
  }
  # this only gets pvalues not the multiple testing corrected
  #pvals = lapply(ret, p.value)
  #lapply(pvals, function(x) {print(which(x<0.05))})

  ret = list(path.ret=ret, symbol=genes.symbol, samples=sample.keep)
  
  return (ret)

}

# 
# replaced this version with the one above on Nov 02, 2010
# 
#pathgt <- function (ds=NULL, dgep=NULL, n.permutation=10000, entrez=NULL, db.what="KEGG", xml.broad=NULL, category.broad=c("c2")) {
#
#  require("org.Hs.eg.db") || stop("please install org.Hs.eg.db first")
#  require(globaltest) || stop("please install globaltest first")
#
#  gt.options(transpose=T) # critical for matrix as data input
#  
#  if(is.null(entrez)) gene.mapping = as.list(gsub("_at", "", substring(colnames(ds),2,100)))
#  else gene.mapping = entrez                # 2010-03-08 >> if probeset ids are not "entrez_at", we need a list of entrez genes
#  names(gene.mapping) = colnames(ds)
#
#  genes.symbol = gsub("^X", "", colnames(ds))
#  genes.symbol = getFeatures(genes.symbol,annot='Symbol',objannot=dgep$features)
#
#  test.pheno = dgep$obs[rownames(ds),]
#  v.timepoints = sort(unique(test.pheno$TimePoints))
#  n.timepoints = length(v.timepoints)
#
#  ret = list()
#  sample.keep = list()
#
#  for (i in 1:n.timepoints){
#
#    tmp =  test.pheno[test.pheno$TimePoints == v.timepoints[i], ]
#    a.pheno = as.factor(test.pheno[rownames(tmp),]$PhenoNew)
#    a.data = as.matrix(ds[rownames(tmp),])
#    sample.keep[[i]] = cbind(samples=rownames(tmp),phenotype=a.pheno)
#    rm(tmp)
#
#    tmp = NULL
#
#
#    if(db.what=="Broad") {
#
#      require(GSEABase) || stop("get package{GSEABase} installed first!")
#
#      if (is.null(xml.broad)) xml.broad = "C:\\Users\\yongsheng\\Documents\\Yongsheng's Documents\\Develop\\R\\database\\msigdb_v2.5.xml"
#      broad = getBroadSets(xml.broad)
#
#      a.Broad = gtBroad (   a.pheno
#                          , t(a.data)
#                          , collection=broad
#                          , category = category.broad
#                          , probe2entrez=gene.mapping
#                          , annotation="org.Hs.eg.db"
#                          , permutations=n.permutation
#                          )
#      tmp = a.Broad
#    }
#
#    if (db.what=="KEGG") {
#      a.KEGG = gtKEGG(  a.pheno
#                      , t(a.data)
##                      , alias = gene.symbol
#                      , probe2entrez=gene.mapping, annotation="org.Hs.eg.db", permutations = n.permutation)
#      tmp = a.KEGG
#    }
#
#    if (db.what=="GO") {
#      cat("GO not implemeted yet")
#      tmp = NULl
#    }
#
#    ret[i] = tmp
#    names(ret)[i] = paste(db.what, v.timepoints[i], sep="")
#  }
#
#  # this only gets pvalues not the multiple testing corrected
#  #pvals = lapply(ret, p.value)
#  #lapply(pvals, function(x) {print(which(x<0.05))})
#
#  ret = list(path.ret=ret, symbol=genes.symbol, samples=sample.keep)
#  
#  return (ret)
#
#}

#***********************************************************************************************
# phdpath::pathgtPathwy
#   1) combine global test statistics (z-score, p-values, or statistic) of pathways accross
#      multiple comparisons (e.g., multiple time points)
#   2) adjust nominal pvalues with specified (default "BH") multiple testing correction
#   3) if desired, plot the pvalues across multiple time points
#
# AUTHOR
#   Yongsheng Huang
#   {huangys}@umich.edu
#   PhD Candidate, Bioinformatics, U. of Michigan
#   http://www-personal.umich.edu/~huangys
#   http://www.hyperfocal.org
#
# PARAMETERS
#   obj.pathgt  -   a globaltest object containing results of pathway analysis from multiple time points
#   what        -   c("p-value", "z-score", "statistic")
#   multtest    -   type of multiple testing correction; c("holm","fdr", "BH", "BY", "Bonferoni", "none")
#   na.rm       -   whether to remove the NA pathways (usually no input genes are present in such pathways)
#   fpdf        -   set this parameter if a pdf file needs to be plotted
#   fpdf.single -   (default TRUE), whether to plot one pathway at a single pdf page
#
# VALUE
#   ret         -   a vector holding counts of signifciant pathways in each pathway analysis
#
# DEPENDENCIES
#
# DBG
#
# NOTES
#   1) -log(p-values, 10)  are used to transform all the p-values on the plotting
#
# USAGE
#   dgep = loadData("gep")
#   list.what = loadLW("flu.list.what.csv")
#   test = data.frame(t(getDesign(list.data=dgep, list.what=list.what)))
#   test.KEGG = pathgt (ds=test, dgep=dgep)
#   test.Broad = pathgt (ds=test, dgep=dgep, db.what="Broad", category.broad=c("c2"))
#   test.KEGG.pathway = pathgtPathwy(test.KEGG, what="p-value", multtest="BH", na.rm=F);     # get p-values
#       test.KEGG.pathway = test.KEGG.pathway[sort(rownames(test.KEGG.pathway)),]
#   test.KEGG.pathwayz = pathgtPathwy(test.KEGG, what="z-score", na.rm=T);      # get z-scores
#       test.KEGG.pathwayz = test.KEGG.pathwayz[sort(rownames(test.KEGG.pathwayz)), ]
#   test.KEGG.pathway = pathgtPathwy(test.KEGG, what="p-value", multtest="BH", na.rm=T, fpdf.single=TRUE, fpdf="test.KEGG.Pathway.pvaladj.pdf");
#   test.KEGG.pathway = pathgtPathwy(test.KEGG, what="p-value", multtest="BH", na.rm=T, fpdf.single=F, fpdf="test.KEGG.Pathway.pvaladj.pdf");
#
#   test.Broad.pathway = pathgtPathwy(test.Broad, what="p-value", multtest="BH", na.rm=T, fpdf.single=F, fpdf="test.Broad.pvaladj.pdf");
#   test.Broad.pathway = pathgtPathwy(test.Broad, what="p-value", multtest="BH", na.rm=T, fpdf.single=TRUE, fpdf="test.Broad.pvaladj.single.pdf");
#
# KNOWN ISSUES
#   known bugs needed to be fixed
#
# TIME STAMP
#   2009-10-03
#
# DBG
#   setting parmeters for running the function
#
# NOTES
#
# 1) 2009-10-02
# plot p-values of a list of sepcified pathways (parameter >> plot.subset)
#
# this portion of code is not finished because it turned out to be
# less interesting when you plot multiple pathways at the same time you
# only see lines cross over with each other (mesh-like)
# therefore we do not continue pursuing this plotting scheme
# pathgtPathwy(test.Broad, what="p-value", multtest="BH", na.rm=T, fpdf.single=FALSE, fpdf="test.Broad.pvaladj.single.xset.pdf", plot.subset=x.test);
#***********************************************************************************************
pathgtPathwy <- function (obj.pathgt=NULL, what="p-value", alpha=0.05, multtest="holm", na.rm=FALSE, fpdf=NULL, plot.subset=NULL, fpdf.single=TRUE, fpdf.width=11, fpdf.height=6) {

  pathgt = obj.pathgt[[1]]

  # set the lowest bound of p-value that we will plot
  p.lower = 0.00001

  # much easier if we only need z-score of the pathways
  if (what=="z-score") {
    ret = sapply(pathgt, z.score)
    ret = df.na (df=ret, by.row=T)               # remove NA rows from return values
    return(ret)
  }

 pathgt = lapply(pathgt, result)

  # only needed if it is Broad database, we add an "alias" column so that it is consistent to KEGG database output
  if (!(length(grep("alias",colnames(pathgt[[1]])))>0)) {
    pathgt = lapply(pathgt, function(x){x$alias = rownames(x); return(x);})
  }

  # now we sort each dataframe by the names of the pathways
  pathgt = lapply(pathgt, sort.df, by=~alias)

  ret = pathgt[[1]][,c("alias","alias")]      # this duplicate of "alias" is done so on purpose
                                              # so that we can get a dataframe easily

  if(what=="p-value")  ret = cbind(ret, sapply(pathgt, function(x){return(p.adjust(x[,what],method=multtest))}))
  else ret = cbind(ret, sapply(pathgt, function(x){return(x[,what])}))

  ret = ret[,-2]                              # now it's good time to get rid of the duplicated "alias"

  # remove the NA from data set
   if(na.rm) ret = ret[-which(apply(ret, 1, function(x){sum(is.na(x))})>0),]

  #
  # if so desired, we plot the p-value of each individual pathway in the database
  #
  if(!is.null(fpdf) & fpdf.single) {
    pdf(fpdf, width=fpdf.width, height=fpdf.height)
    for (i in 1:dim(ret)[1]){
      tmp = ret[i,]
      tmp.name = tmp$alias
      tmp.pval = t(tmp[,-1])
      tmp.pval[tmp.pval==0] = p.lower
      tmp.pval = -log(tmp.pval,10)
      plot(tmp.pval, ylim=c(0,ceiling(max(tmp.pval))),type="l", lwd=3, xlab="Time", ylab="-log(p-value)",font.lab=2, font.axis=2, main=tmp.name, xaxt="n", axes=FALSE)

#      x.labels = gsub("Broad", "", gsub("KEGG","",rownames(tmp.pval)))
      x.labels = gsub("^KEGG","",rownames(tmp.pval))
      x.labels = gsub("^Broad","",x.labels)
      x.labels = gsub("^GO","",x.labels)
      axis( side=1, at = 1:dim(tmp.pval)[1], labels = x.labels
            , cex.axis=1, lwd=2, lwd.ticks=2, font=2)
#      axis( side=2, at = seq(0,(ceiling(max(tmp.pval)+1)),by=1), labels = seq(0,(ceiling(max(tmp.pval)+1)),by=1)
      axis( side=2, at = seq(0,ceiling(max(tmp.pval)),by=1)#, labels = seq(0,ceiling(max(tmp.pval)),by=1)
            , cex.axis=1, lwd=2, lwd.ticks=2, font=2, las=1)
      points(tmp.pval, col="red", pch=19)
      abline(h=-log(alpha,10),col="gray60", lty=2, lwd=2)
    }
    dev.off()
  }

#  #
#  # 2009-10-02
#  # plot p-values of a list of sepcified pathways
#  #
#  # note: this portion of code is not finished because it turned out to be
#  # less interesting when you plot multiple pathways at the same time you
#  # only see lines cross over with each other (mesh-like)
#  # therefore we do not continue pursuing this plotting scheme
#
#  if(!is.null(fpdf) & !is.null(plot.subset)) {
#    pdf(fpdf, width=11, height=6)
#
#    tmp.set = ret[plot.subset,]
#    tmp.set.pval = (tmp.set[,-1])
#    tmp.set.pval[tmp.set.pval==0] = p.lower
#    tmp.set.pval = -log(tmp.set.pval)
#    # keep the max and set the ylim for plot
#    ylims = c(0, max(tmp.set.pval))
#    # update p-values with -log(p-values)
#    tmp.set[,-1] = tmp.set.pval
#    for (i in 1:dim(tmp.set)[1]){
#      tmp = tmp.set[i,]
#      tmp.name = tmp$alias
#      tmp.pval = t(tmp[,-1])
##      tmp.pval[tmp.pval==0] = p.lower
##      tmp.pval = -log(tmp.pval)
#      if (i==1){
#        plot( tmp.pval, type="l", lwd=3, xlab="Time"
##              , yaxt="n", ylim=ylims, ylab="-log(p-value)"
#              , ylim=ylims, ylab="-log(p-value)"
#              ,font.lab=2, font.axis=2
#              #, main=tmp.name
#              , xaxt="n", axes=FALSE)
#      }
#      else {
#        par(new=T)
#        plot( tmp.pval, type="l", lwd=3, xlab="Time"
##              , yaxt="n", ylim=ylims, ylab="-log(p-value)"
#              , ylim=ylims, ylab="-log(p-value)"
#              ,font.lab=2, font.axis=2
#              #, main=tmp.name
#              , xaxt="n", axes=FALSE)
#
#      }
##      x.labels = gsub("Broad", "", gsub("KEGG","",rownames(tmp.pval)))
#      x.labels = gsub("^KEGG","",rownames(tmp.pval))
#      x.labels = gsub("^Broad","",x.labels)
#      x.labels = gsub("^GO","",x.labels)
#      axis( side=1, at = 1:dim(tmp.pval)[1], labels = x.labels
#            , cex.axis=1, lwd=2, lwd.ticks=2, font=2)
#      axis( side=2, at = seq(0,max(tmp.pval),by=1), labels = seq(0,max(tmp.pval),by=1)
#            , cex.axis=1, lwd=2, lwd.ticks=2, font=2, las=1)
#      points(tmp.pval, col="red", pch=19)
#      abline(h=-log(0.05),col="gray60", lty=2, lwd=2)
#    }
#    dev.off()
#  }


  #
  # if so desired, we plot the p-value for all pathways in the database
  #
  if(!is.null(fpdf) & !fpdf.single) {
    pdf(fpdf, height=7, width=9.8)
    for (i in 1:dim(ret)[1]){
      tmp = ret[i,]
      tmp.name = tmp$alias
      tmp.pval = t(tmp[,-1])
      tmp.pval[tmp.pval==0] = p.lower
      tmp.pval = -log(tmp.pval,10)
      if(i!=1) par(new=T)
      point.colour = ifelse(tmp.pval>(-log(alpha, 10)),"red","gray70")
#      plot( jitter(tmp.pval)~jitter(1:length(tmp.pval))
      tmp.pval = rbind(rbind(0,tmp.pval), 0)
      point.colour = c("white",point.colour, "white")
      a = seq(1:dim(tmp.pval)[1]) + rnorm(dim(tmp.pval)[1])/10
      plot( tmp.pval~a
            , type="p", pch=17, cex=0.7, col=point.colour
            , xlab="Time", ylab="-log(p-value)"
            , ylim=c(0,round(-log(p.lower)))
            , yaxt="n", xaxt="n")

    }

    # get rid of the first and last fake elements
#    tmp.pval = tmp.pval[-c(1, dim(tmp.pval)[1]),]
    # x.labels = gsub("Broad", "", gsub("KEGG","",rownames(tmp.pval)))
      x.labels = gsub("^KEGG","",rownames(tmp.pval))
      x.labels = gsub("^Broad","",x.labels)
      x.labels = gsub("^GO","",x.labels)
#      axis(side=1, at = 1:length(tmp.pval), labels = x.labels, cex.axis=0.8, font=2)
      axis(side=1, at = seq(1:dim(tmp.pval)[1]), labels = x.labels, cex.axis=0.8, font=2)
      axis(side=2, at = 0:(round(-log(p.lower))), labels=0:(round(-log(p.lower))), cex.axis=0.8, font=2)
      abline(h=-log(alpha,10),col="gray60", lty=2, lwd=2)

    dev.off()
  }

  return (ret)

}




#===================================================================================================
# phdpath::pathgtCnt
#    function for counting the number of significant pathways in each analysis
#    among analyses of multiple time points, with multiple testing correction
#
# PARAMETERS
#   obj.pathgt      a globaltest object containing results of pathway analysis from multiple time points
#   p.cut           cutoff value for calling a p-value significant after multiple testing correction
#   what.mult       c("holm", "BH", "BY")
#
# DEPENDENCIES
#
# VALUE
#   ret             a vector holding counts of signifciant pathways in each pathway analysis
#
# USAGE
#   dgep = loadData("gep")
#   list.what = loadLW("flu.list.what.csv")
#   test = data.frame(t(getDesign(list.data=dgep, list.what=list.what)))
#   test.KEGG = pathgt (ds=test, dgep=dgep)
#   test.Broad = pathgt (ds=test, dgep=dgep, db.what="Broad", category.broad=c("c2"))
#   test.KEGG.cnt = pathgtCnt(obj.pathgt=test.KEGG, p.cut=0.05, what.mult="BH"); print(test.KEGG.cnt)
#   test.Broad.cnt = pathgtCnt(obj.pathgt=test.Broad, p.cut=0.01, what.mult="none"); print(test.Broad.cnt)
#
# DBG
#   debugging code here
#
# NOTES
#   1) "BH" will give exactly the same results as "fdr" - need to check `multtest` function later
#
# KNOWN ISSUES
#
#
# TIME
#   2009-08-10
#===================================================================================================
pathgtCnt <- function (obj.pathgt=NULL, p.cut=0.05,what.mult="holm"){

  pathgt = obj.pathgt[[1]]
  pathgt = lapply(pathgt, result)
  pval.adj = lapply(  pathgt
                    , function(x) { ret=p.adjust(x[,"p-value"], method=what.mult); return(ret) }
                   )
  ret = unlist(lapply(pval.adj, function(x){ret=sum(x<=p.cut, na.rm=TRUE);return(ret)}))
  
  return (ret)
  
}

#===================================================================================================
# phdpath::pathgtCut
#    selecting significant pathways and their adjusted-p-values from multiple analyses
#    and assemble them into one data frame (to facilitate visualization)
#
# PARAMETERS
#   obj.pathgt      a globaltest object containing results of pathway analysis from multiple time points
#   p.cut           cutoff value for calling a p-value significant after multiple testing correction
#   what.mult       c("holm", "BH", "BY")
#
# DEPENDENCIES
#   1) phdpath::pathgtCnt
#
# VALUE
#   ret             a dataframe holding all sig. pathways and their adjusted pvalues
#
# USAGE
#   dgep = loadData("gep")
#   list.what = loadLW("flu.list.what.csv")
#   test = data.frame(t(getDesign(list.data=dgep, list.what=list.what)))
#   test.KEGG = pathgt (ds=test, dgep=dgep)
#   test.Broad = pathgt (ds=test, dgep=dgep, db.what="Broad", category.broad=c("c2"))
#   test.pval =  pathgtCut(test.KEGG,p.cut=0.05,what.mult="BH"); print(test.pval);
#
# DBG
#   debugging code
#
# NOTES
#   1) "BH" will give exactly the same results as "fdr" - need to check `multtest` function later
#
# KNOWN ISSUES
#
#
# TIME
#   2009-08-10
#===================================================================================================
pathgtCut <- function (obj.pathgt=NULL, p.cut=0) {

  ret = NULL
  
  pathgt = obj.pathgt[[1]]
  pathgt.cntsig =  pathgtCnt(obj.pathgt,p.cut=0.05,what.mult=what.mult)

  for (i in 1:length(pathgt)) {
    n.sig = pathgt.cntsig[i]

    if(n.sig==0) {
      ret = rbind(ret, data.frame(adj.pval=NA, category=names(pathgt)[i]))
      next
    }
    
    tmp = pathgt[[i]]
    tmp.pval = p.adjust(p.value(tmp),method=what.mult)[1:n.sig]
    names(tmp.pval) = names(tmp)[1:n.sig]
    
    ret = rbind(ret, data.frame(adj.pval=tmp.pval, category=names(pathgt)[i]))
    
   }

  return (ret)
}



#===================================================================================================
# phdpath::pathgtPlotCov
#    plot significant pathways and their adjusted-p-values from multiple analyses
#    and assemble them into one data frame (to facilitate visualization)
#
# PARAMETERS
#   obj.pathgt      a globaltest object containing results of pathway analysis from multiple time points
#   p.cut           cutoff value for calling a p-value significant after multiple testing correction
#   what.mult       c("holm", "BH", "BY")
#   what.plot       what type of covariates plot should be produced
#                   what.plot = c("p-value", "statistic", "z-score", "weighted")
#   do.cluster      defalut FALSE; whether clustering of covariates should be conducted
#                   > will cause issue if a pathway only have one gene
#
# DEPENDENCIES
#   1) phdpath::pathgtCnt
#   2) phdpath::ys.covariates.pval / phdpath::ys.covariates (the first one returns p-values of each gene)
#
# VALUE
#   None; One pdf file is generated for each time point
#
# USAGE
#   dgep = loadData("gep")
#   list.what = loadLW("flu.list.what.csv")
#   test = data.frame(t(getDesign(list.data=dgep, list.what=list.what)))
#   test.KEGG = pathgt (ds=test, dgep=dgep)
#   pathgtPlotCov(obj.pathgt=test.KEGG, p.cut=0.05, what.plot="statistic", do.cluster=FALSE)
#
#   pathgtPlotCov(obj.pathgt=test.Broad, p.cut=0.05, what.mult="BH", what.plot="p-value", do.cluster=FALSE)
#
# DBG
#
#
# NOTES
#   1) "BH" will give exactly the same results as "fdr" - need to check `multtest` function later
#   2) if `do.cluster=TRUE`, it will give an error when a pathway only contains one gene
#   3) for each analysis (at a time point), a PDF file is generated at current working directory
#   4) if a pathway pvalue=0, then the plot of "p-value" will not work
#
# KNOWN ISSUES
#
#
# TIME
#   2009-08-10
#===================================================================================================
pathgtPlotCov <- function(obj.pathgt=NULL, p.cut=0.05, y.las=2, what.mult="holm", what.plot="z", do.cluster=FALSE) {
   ret = NULL
   pathgt = obj.pathgt[[1]]
   pathgt.cntsig =  pathgtCnt(obj.pathgt, p.cut=p.cut, what.mult=what.mult)
#   pathgt.cntsig =  pathgtCnt(pathgt, p.cut=p.cut, what.mult=what.mult)

   pathgt.names = names(pathgt)
   genes.symbol = obj.pathgt[["symbol"]]
   
   for (i in 1:length(pathgt)) {
    n.sig = pathgt.cntsig[i]
    if(n.sig==0) next

    tmp = pathgt[[i]]
    fname = paste(pathgt.names[i],".pdf",sep="")
    if(do.cluster) do.cluster=(result(tmp)[,"#Cov"]>1)[1:n.sig]         # doesn't work now
    #covariates(tmp[1:n.sig], what=what.plot, cluster=do.cluster, pdf=fname, alias=genes.symbol)
    # ys.covariates(tmp[1:n.sig], what=what.plot, cluster=do.cluster, pdf=fname, alias=genes.symbol)
    ret = ys.covariates.pval(tmp[1:n.sig], what=what.plot, y.las=y.las, cluster=do.cluster, pdf=fname, alias=genes.symbol, colors=c("red","green"))
   }
#   return(ret)
}


#===================================================================================================
# phdpath::ys.covariates.pval
#   wrap function for covariates with fix of p-value plot and output p-value
#
# VALUE
#   p-value of individual genes on this pathway
#
# PARAMETER
#   what.ret    -    either genes or the whole pathway; what values to return
#   y.las       -    default 2;  numeric in {0,1,2,3}; the style of axis labels.
#                    0:
#                    always parallel to the axis [default],
#                    1:
#                    always horizontal,
#                    2:
#                    always perpendicular to the axis,
#                    3:
#                    always vertical.
#
# TIME
#   2009-10-04
#===================================================================================================
ys.covariates.pval <- function (object, what = c("p-value", "statistic", "z-score",
    "weighted"), cluster = "average", alpha = 0.05, sort = TRUE,
    legend = TRUE, colors, alias, cex.labels = 0.6, pdf, trace,
    what.ret = "genes",
    y.las =2 )
{
    if ((length(object) > 1) && missing(pdf))
        stop("length(object) > 1. Please reduce to a single test result or specify an output file.")
    if (missing(trace))
        trace <- gt.options()$trace
    if (missing(alias))
        alias <- NULL
    if (!missing(pdf)) {
        if (tolower(substr(pdf, nchar(pdf) - 3, nchar(pdf))) !=
            ".pdf")
            pdf <- paste(pdf, ".pdf", sep = "")
        pdf(pdf)
    }
    what <- substr(match.arg(tolower(what), c("p-value", "statistic",
        "z-score", "weighted")), 1, 1)
    if (is.character(legend)) {
        object@legend$cov <- legend
        legend <- TRUE
    }
    if (!is.null(alias)) {
        if (is.environment(alias))
            alias <- as.list(alias)
        if (is.list(alias))
            alias <- unlist(alias)
        if (length(alias) == object@functions$df()[3] && is.null(names(alias)))
            names(alias) <- object@functions$cov.names()
    }
    for (jj in 1:length(object)) {
        obj <- object[jj]
        if (is.null(obj@weights))
            weights <- rep(1, size(obj))
        else weights <- obj@weights[[1]]
        if (is.null(obj@subsets)) {
            subset <- seq_len(size(obj))
            obj@subsets[[1]] = subset
        }
        else subset <- obj@subsets[[1]]
        ttl <- names(obj)
        if (!is.null(alias(obj)))
            ttl <- paste(ttl, "-", alias(obj))
        test <- function(set) {
            obj@functions$test(subset[set], weights[set])
        }
        leaves <- t(sapply(1:size(obj), function(i) {
            test(i)
        }))
        rownames(leaves) <- obj@functions$cov.names(subset)
        if (is.null(rownames(leaves)))
            rownames(leaves) <- subset
        if (what == "w") {
            leaves[, c("S", "ES", "sdS")] <- leaves[, c("S",
                "ES", "sdS")] * matrix(weights(obj), size(obj),
                3)
        }
        pps <- -log10(leaves[, "p"])
        #
        # yongsheng modified here
        #
        if (length(pps)==1) pps[is.infinite(pps)]=-log(0.0001)
        else if(all(is.infinite(pps))) pps[is.infinite(pps)]=-log(0.0001)
        else pps[is.infinite(pps)]=max(pps[!is.infinite(pps)],na.rm=T)
        #
        #
        #
        bars <- switch(what, p = pps, z = (leaves[, "S"] - leaves[,
            "ES"])/leaves[, "sdS"], w = , s = leaves[, "S"])
        names(bars) <- rownames(leaves)
        if (!is.null(alias)) {
            if (length(alias) == length(bars) && (is.null(names(alias)) ||
                is.null(names(bars))))
                names(bars) <- alias
            else {
                names(bars) <- alias[names(bars)]
            }
        }
        if (sort)
            order.bars <- -pps
        else order.bars <- 1:length(bars)
        margins <- par("mai")
        dendrogram <- ((!is.logical(cluster)) || (cluster)) &&
            (substr(cluster, 1, 1) != "n")
        if (is.logical(cluster) && dendrogram)
            cluster <- "average"
        if (dendrogram) {
            cors <- obj@functions$cor(subset)
            if (obj@directional)
                dd <- as.dist(1 - cors)
            else dd <- as.dist(1 - abs(cors))
            hc <- as.dendrogram(hclust(dd, method = cluster))
            hc <- reorder(hc, wts = order.bars, agglo.FUN = min)
            sorter <- unlist(hc)
            obj@result = rbind(obj@result, leaves)
            obj@subsets = c(list(obj@functions$cov.names(obj@subsets[[1]])),
                unlist(obj@functions$cov.names(subset)))
            obj@extra <- NULL
            obj@structure <- NULL
            obj = inheritance(obj, sets = hc, trace = trace,
                stop = 1)
            sigcol <- function(branch, sig, top) {
                setlist = obj@subsets
                labels = unlist(dendrapply(branch, function(n) attributes(n)$label))
                newlabels = unlist(branch)
                labels[newlabels] = labels
                branchset <- unlist(branch)
                selected <- names(setlist)[which(lapply(setlist,
                  setequal, labels[branchset]) == T)]
                if (sig) {
                  sig <- obj@extra$inheritance[which(names(obj) ==
                    selected)] <= alpha
                }
                uit <- branch
                attr(uit, "edgePar") <- list(col = ifelse(sig,
                  1, gray(0.8)), lwd = ifelse(sig, 2, 1))
                if (sig && top)
                  attr(uit, "nodePar") <- list(pch = 20)
                if (!is.leaf(branch)) {
                  for (i in 1:length(branch)) {
                    uit[[i]] <- sigcol(branch[[i]], sig, FALSE)
                  }
                }
                return(uit)
            }
            hc <- sigcol(hc, sig = TRUE, top = TRUE)
            par(mai = c(0, max(1, margins[2]), margins[3:4]))
            layout(as.matrix(1:2), heights = c(1, 2))
            ylab <- ifelse(obj@directional, "correlation", "absolute correlation")
            plot(hc, leaflab = "none", yaxt = "n", ylab = ylab,
                mgp = c(4, 1, 0))
            axis(2, at = seq(0, 2, by = 0.2), labels = 1 - seq(0,
                2, by = 0.2), las = 2)
        }
        else {
            sorter <- sort.list(order.bars)
        }
        leaves <- leaves[sorter, , drop = FALSE]
        bars <- bars[sorter]
        labwidth <- max(strwidth(names(bars), "inches", cex.labels)) +
            0.2
        par(mai = c(max(margins[1], labwidth * 1.3), max(1, margins[2]),
            if (dendrogram) 0 else margins[3], margins[4]))
        positive <- obj@functions$positive(subset)[sorter]
        if (missing(colors))
            if (max(positive) <= 2)
                colors <- 3:2
            else colors <- rainbow(max(positive), start = 0,
                end = 1/2)
        if (all(positive %in% 0:1))
            cols <- ifelse(positive, colors[1], colors[2])
        else cols <- colors[positive]
        ylab <- switch(what, p = "p-value", z = "z-score", s = "test statistic",
            w = "weighted test statistic")
        ylims <- switch(what, z = , p = range(bars), w = , s = range(c(bars,
            leaves[, "ES"])))
        if (ylims[1] > 0)
            ylims[1] <- 0
        if (legend) {
            nbars <- length(bars)
            room <- (ylims[2] - max(bars[trunc(nbars * 0.6):nbars]))/diff(ylims)
            ylims[2] <- ylims[2] + diff(ylims) * max(0, 0.1 *
                length(colors) - room)
        }
        #
        # yongsheng modified here
        #
        mids <- drop(barplot(bars, xaxt="n", yaxt = "n", las = 2, ylab = ylab,
            ylim = ylims, mgp = c(4, 1, 0), col = cols, cex.names = cex.labels))
        axis(side=1,at=mids,labels=names(bars),lwd=2,lwd.ticks=2,font=2,las=2,cex=2)
        #
        #
        #
        if (dendrogram) {
            mb <- max(bars)
            for (i in 1:length(bars)) lines(c(mids[i], mids[i]),
                c(max(0, bars[i]) + 0.01 * mb, mb * 1.2), col = gray(0.8),
                lty = 3)
        }
        if (what == "p") {
            abline(h=(-log(alpha,10)),col="grey60",lty=2,lwd=2)
            maxlogp <- max(bars, na.rm = TRUE)
            labs <- seq(0, maxlogp, by = max(1, maxlogp%/%5))
            if (length(labs) == 1)
                labs <- log10(c(1, 2, 10/3, 5))
            else if (length(labs) <= 2)
                labs <- outer(log10(c(1, 2, 5)), labs, "+")
            else if (length(labs) <= 4)
                labs <- outer(log10(c(1, 10/3)), labs, "+")
            #
            # yongsheng modified here
            #
            #axis(2, at = labs, labels = 10^-labs, las = 2, font=2, lwd=2,lwd.ticks=2)
            axis(2, at = labs, labels = 10^-labs, las = y.las, font=2, lwd=2,lwd.ticks=2)
        }
        else axis(2, las = 2)
        if (what %in% c("s", "w")) {
            sapply(seq_along(mids), function(i) {
                lines(c(mids[i] - 0.5, mids[i] + 0.5), rep(leaves[i,
                  "ES"], 2), lwd = 3)
                sapply(seq_len(max(0, (bars[i] - leaves[i, "ES"])/leaves[i,
                  "sdS"])), function(k) lines(c(mids[i] - 0.5,
                  mids[i] + 0.5), rep(leaves[i, "ES"] + k * leaves[i,
                  "sdS"], 2)))
            })
        }
        abline(0, 0)
        if (legend)
            legend("topright", obj@legend$cov, fill = colors,
                bg = "white")
        layout(1)
        par(mai = margins)
        if (!missing(pdf)) {
            title(ttl)
        }
    }

    if (!missing(pdf))
        dev.off()
    if (length(object) == 1) {
        out <- obj
    }
    else out <- NULL

    if(what.ret=="pathway") ret = 10^(-pps)
    else if(what.ret=="genes") ret = leaves

#    return(10^(-pps))
    return (ret)
}



##===================================================================================================
## phdpath::ys.covariates
##   wrap function for covariates with fix of p-value plot
##
## TIME
##   2009-08-10
##===================================================================================================
#ys.covariates <- function (object, what = c("p-value", "statistic", "z-score",
#    "weighted"), cluster = "average", alpha = 0.05, sort = TRUE,
#    legend = TRUE, colors, alias, cex.labels = 0.6, pdf, trace)
#{
#    if ((length(object) > 1) && missing(pdf))
#        stop("length(object) > 1. Please reduce to a single test result or specify an output file.")
#    if (missing(trace))
#        trace <- gt.options()$trace
#    if (missing(alias))
#        alias <- NULL
#    if (!missing(pdf)) {
#        if (tolower(substr(pdf, nchar(pdf) - 3, nchar(pdf))) !=
#            ".pdf")
#            pdf <- paste(pdf, ".pdf", sep = "")
#        pdf(pdf)
#    }
#    what <- substr(match.arg(tolower(what), c("p-value", "statistic",
#        "z-score", "weighted")), 1, 1)
#    if (is.character(legend)) {
#        object@legend$cov <- legend
#        legend <- TRUE
#    }
#    if (!is.null(alias)) {
#        if (is.environment(alias))
#            alias <- as.list(alias)
#        if (is.list(alias))
#            alias <- unlist(alias)
#        if (length(alias) == object@functions$df()[3] && is.null(names(alias)))
#            names(alias) <- object@functions$cov.names()
#    }
#    for (jj in 1:length(object)) {
#        obj <- object[jj]
#        if (is.null(obj@weights))
#            weights <- rep(1, size(obj))
#        else weights <- obj@weights[[1]]
#        if (is.null(obj@subsets)) {
#            subset <- seq_len(size(obj))
#            obj@subsets[[1]] = subset
#        }
#        else subset <- obj@subsets[[1]]
#        ttl <- names(obj)
#        if (!is.null(alias(obj)))
#            ttl <- paste(ttl, "-", alias(obj))
#        test <- function(set) {
#            obj@functions$test(subset[set], weights[set])
#        }
#        leaves <- t(sapply(1:size(obj), function(i) {
#            test(i)
#        }))
#        rownames(leaves) <- obj@functions$cov.names(subset)
#        if (is.null(rownames(leaves)))
#            rownames(leaves) <- subset
#        if (what == "w") {
#            leaves[, c("S", "ES", "sdS")] <- leaves[, c("S",
#                "ES", "sdS")] * matrix(weights(obj), size(obj),
#                3)
#        }
#        pps <- -log10(leaves[, "p"])
#        #
#        # yongsheng modified here
#        #
#        if (length(pps)==1) pps[is.infinite(pps)]=-log(0.0001)
#        else pps[is.infinite(pps)]=max(pps[!is.infinite(pps)],na.rm=T)
#        #
#        #
#        #
#        bars <- switch(what, p = pps, z = (leaves[, "S"] - leaves[,
#            "ES"])/leaves[, "sdS"], w = , s = leaves[, "S"])
#        names(bars) <- rownames(leaves)
#        if (!is.null(alias)) {
#            if (length(alias) == length(bars) && (is.null(names(alias)) ||
#                is.null(names(bars))))
#                names(bars) <- alias
#            else {
#                names(bars) <- alias[names(bars)]
#            }
#        }
#        if (sort)
#            order.bars <- -pps
#        else order.bars <- 1:length(bars)
#        margins <- par("mai")
#        dendrogram <- ((!is.logical(cluster)) || (cluster)) &&
#            (substr(cluster, 1, 1) != "n")
#        if (is.logical(cluster) && dendrogram)
#            cluster <- "average"
#        if (dendrogram) {
#            cors <- obj@functions$cor(subset)
#            if (obj@directional)
#                dd <- as.dist(1 - cors)
#            else dd <- as.dist(1 - abs(cors))
#            hc <- as.dendrogram(hclust(dd, method = cluster))
#            hc <- reorder(hc, wts = order.bars, agglo.FUN = min)
#            sorter <- unlist(hc)
#            obj@result = rbind(obj@result, leaves)
#            obj@subsets = c(list(obj@functions$cov.names(obj@subsets[[1]])),
#                unlist(obj@functions$cov.names(subset)))
#            obj@extra <- NULL
#            obj@structure <- NULL
#            obj = inheritance(obj, sets = hc, trace = trace,
#                stop = 1)
#            sigcol <- function(branch, sig, top) {
#                setlist = obj@subsets
#                labels = unlist(dendrapply(branch, function(n) attributes(n)$label))
#                newlabels = unlist(branch)
#                labels[newlabels] = labels
#                branchset <- unlist(branch)
#                selected <- names(setlist)[which(lapply(setlist,
#                  setequal, labels[branchset]) == T)]
#                if (sig) {
#                  sig <- obj@extra$inheritance[which(names(obj) ==
#                    selected)] <= alpha
#                }
#                uit <- branch
#                attr(uit, "edgePar") <- list(col = ifelse(sig,
#                  1, gray(0.8)), lwd = ifelse(sig, 2, 1))
#                if (sig && top)
#                  attr(uit, "nodePar") <- list(pch = 20)
#                if (!is.leaf(branch)) {
#                  for (i in 1:length(branch)) {
#                    uit[[i]] <- sigcol(branch[[i]], sig, FALSE)
#                  }
#                }
#                return(uit)
#            }
#            hc <- sigcol(hc, sig = TRUE, top = TRUE)
#            par(mai = c(0, max(1, margins[2]), margins[3:4]))
#            layout(as.matrix(1:2), heights = c(1, 2))
#            ylab <- ifelse(obj@directional, "correlation", "absolute correlation")
#            plot(hc, leaflab = "none", yaxt = "n", ylab = ylab,
#                mgp = c(4, 1, 0))
#            axis(2, at = seq(0, 2, by = 0.2), labels = 1 - seq(0,
#                2, by = 0.2), las = 2)
#        }
#        else {
#            sorter <- sort.list(order.bars)
#        }
#        leaves <- leaves[sorter, , drop = FALSE]
#        bars <- bars[sorter]
#        labwidth <- max(strwidth(names(bars), "inches", cex.labels)) +
#            0.2
#        par(mai = c(max(margins[1], labwidth * 1.3), max(1, margins[2]),
#            if (dendrogram) 0 else margins[3], margins[4]))
#        positive <- obj@functions$positive(subset)[sorter]
#        if (missing(colors))
#            if (max(positive) <= 2)
#                colors <- 3:2
#            else colors <- rainbow(max(positive), start = 0,
#                end = 1/2)
#        if (all(positive %in% 0:1))
#            cols <- ifelse(positive, colors[1], colors[2])
#        else cols <- colors[positive]
#        ylab <- switch(what, p = "p-value", z = "z-score", s = "test statistic",
#            w = "weighted test statistic")
#        ylims <- switch(what, z = , p = range(bars), w = , s = range(c(bars,
#            leaves[, "ES"])))
#        if (ylims[1] > 0)
#            ylims[1] <- 0
#        if (legend) {
#            nbars <- length(bars)
#            room <- (ylims[2] - max(bars[trunc(nbars * 0.6):nbars]))/diff(ylims)
#            ylims[2] <- ylims[2] + diff(ylims) * max(0, 0.1 *
#                length(colors) - room)
#        }
#        mids <- drop(barplot(bars, yaxt = "n", las = 2, ylab = ylab,
#            ylim = ylims, mgp = c(4, 1, 0), col = cols, cex.names = cex.labels))
#        if (dendrogram) {
#            mb <- max(bars)
#            for (i in 1:length(bars)) lines(c(mids[i], mids[i]),
#                c(max(0, bars[i]) + 0.01 * mb, mb * 1.2), col = gray(0.8),
#                lty = 3)
#        }
#        if (what == "p") {
#            maxlogp <- max(bars, na.rm = TRUE)
#            labs <- seq(0, maxlogp, by = max(1, maxlogp%/%5))
#            if (length(labs) == 1)
#                labs <- log10(c(1, 2, 10/3, 5))
#            else if (length(labs) <= 2)
#                labs <- outer(log10(c(1, 2, 5)), labs, "+")
#            else if (length(labs) <= 4)
#                labs <- outer(log10(c(1, 10/3)), labs, "+")
#            axis(2, at = labs, labels = 10^-labs, las = 2)
#        }
#        else axis(2, las = 2)
#        if (what %in% c("s", "w")) {
#            sapply(seq_along(mids), function(i) {
#                lines(c(mids[i] - 0.5, mids[i] + 0.5), rep(leaves[i,
#                  "ES"], 2), lwd = 3)
#                sapply(seq_len(max(0, (bars[i] - leaves[i, "ES"])/leaves[i,
#                  "sdS"])), function(k) lines(c(mids[i] - 0.5,
#                  mids[i] + 0.5), rep(leaves[i, "ES"] + k * leaves[i,
#                  "sdS"], 2)))
#            })
#        }
#        abline(0, 0)
#        if (legend)
#            legend("topright", obj@legend$cov, fill = colors,
#                bg = "white")
#        layout(1)
#        par(mai = margins)
#        if (!missing(pdf)) {
#            title(ttl)
#        }
#    }
#    if (!missing(pdf))
#        dev.off()
#    if (length(object) == 1) {
#        out <- obj
#    }
#    else out <- NULL
#    return(invisible(out))
#}
#



cat(">>>>>>> Pathway analysis library (PHD project) has been successfully loaded!\n")





#
#Reference: Bioconductor Workshop Apr 2009 (Seattle)
#http://www.bioconductor.org/workshops/2009/SeattleApr09/annotate/annotationBasics.R
#
####################################################
#### chunk number 1: orgDemo
####################################################
###load the package
#library("org.Hs.eg.db")
#
###look what we just loaded
#ls(2)
#
###Just like user accessible functions, all Mappings have a manual page which
###will show you what to expect as well as where the data came from
## ?org.Hs.egCHRLOC
#
###Have a peak:
#as.list(org.Hs.egCHRLOC[1:4])
#
###for the stop locations use:
#as.list(org.Hs.egCHRLOCEND[1:4])
#
###or can use get, mget etc. with the entrez gene ID
#EGs = c("10","100","1000")
#mget(EGs, org.Hs.egCHRLOC, ifnotfound=NA)
#mget(EGs, org.Hs.egCHRLOCEND, ifnotfound=NA)
#
###You can also retrieve ENSEMBL IDs using this package
#mget(EGs, org.Hs.egENSEMBL, ifnotfound=NA)
#
###And GO IDs
#mget(EGs[1], org.Hs.egGO, ifnotfound=NA)
#
###And KEGG pathway IDs etc.
#mget(EGs, org.Hs.egPATH, ifnotfound=NA)
#
#
###Other convenient functions
###Lkeys, RKeys, mappedLkeys(),mappedRkeys()
#Lkeys(org.Hs.egENZYME)[110:112]
#Rkeys(org.Hs.egENZYME)[110:112]
#
###Left keys and right keys can be mapped or un-mapped.
#length(Lkeys(org.Hs.egPATH))
#length(Rkeys(org.Hs.egPATH))
#length(mappedLkeys(org.Hs.egPATH))
#length(mappedRkeys(org.Hs.egPATH))
#
###keys() and mappedkeys() both return the left keys
#length(keys(org.Hs.egPATH))
#length(mappedkeys(org.Hs.egPATH))
#
#
###revmap() can USUALLY be used to reverse the direction of a mapping.
#PATHIDs = unlist(mget(EGs[1], org.Hs.egPATH, ifnotfound=NA))[[1]]
#PATHIDs
#mget(as.character(PATHIDs), revmap(org.Hs.egPATH), ifnotfound=NA)
#
#
###toTable
#toTable(revmap(org.Hs.egPATH))[1:4,]
#
#
###special symbols: packagename(), _dbfile(), _dbinfo(), and _dbconn()
#ls(2)
#
###package info
#org.Hs.eg()
#
###location of the database file
#org.Hs.eg_dbfile()
#
###Data frame with Information
#org.Hs.eg_dbInfo()
#
###Connection object
#org.Hs.eg_dbconn()
#
#
#
####################################################
#### chunk number 2: dbExamples
####################################################
#
###simple examples of using the DBI interface
#library(org.Hs.eg.db)
#dbconn <- org.Hs.eg_dbconn()
#sql <- "SELECT * FROM genes LIMIT 4;"
#result <- dbGetQuery(dbconn, sql)
#result
#
#
###Here is a simple join of the type that generates the mappings
#sql <- "SELECT * FROM genes,kegg WHERE genes._id=kegg._id;"
#result <- dbGetQuery(dbconn, sql)
#result[1:4,]
#
#
#
###simple example of a join across DBs
###get the entrez genes in humans and in mouse which share a kegg pathway ID.
#
### 1st we must get the pathway to the database file
#path = system.file("extdata", "org.Mm.eg.sqlite", package = "org.Mm.eg.db" )
#
###then we have to attach the database
#sql <- paste("ATTACH '",path,"' AS Mm;",sep="")
#dbSendQuery(dbconn, sql)
#
###Then we can make our query
#sql <- "SELECT g.gene_id, k.path_id, mk.path_id, mg.gene_id
#        FROM genes AS g, kegg AS k, Mm.kegg AS mk, Mm.genes AS mg
#        WHERE g._id=k._id AND mg._id=mk._id AND k.path_id=mk.path_id
#        limit 1000;"
#result <- dbGetQuery(dbconn, sql)
#result[1:4,]
#
#
###can even use SQLite specific stuff
#sql = "SELECT * FROM sqlite_master;"
#result = dbGetQuery(dbconn, sql)
#head(result)
#
#
#
#
####################################################
#### chunk number 3: chipExamples
####################################################
###Things work very similarly to an org package
#library(hgu95av2.db)
#ls(2)
#
###Gene symbols and aliases
#probes = keys(hgu95av2SYMBOL)[1:4]
#probes
#
###You can use the probes in the same way you would have used Entrez Genes
###Here is a mapping that retrieves NCBIs official gene symbols
#mget(probes, hgu95av2SYMBOL, ifnotfound=NA)
#
###And here is a mapping that retrieves all known gene symbols
#mget(probes, revmap(hgu95av2ALIAS2PROBE), ifnotfound=NA)
#
#
###Be careful with Aliases as they are NOT unique
###This is easier to demonstrate with an org package.
###But its even more of a problem in a chip package.
#library(org.Hs.eg.db)
#EG2AliasList = as.list(org.Hs.egALIAS2EG)
#EG2AliasList["KAT"]
#
###We can calculate out how many things match each Alias like this:
#lengths = unlist(lapply(EG2AliasList, length))
#lengths["KAT"]
#table(lengths)
#
###And just out of curiosity:
#lengths[lengths==36]
#EG2AliasList["VH"]
#
#
#
####################################################
#### chunk number 4: customChipExample
####################################################
###1st you need the appropriate DBO package
#library(AnnotationDbi)
#available.db0pkgs()
#
###Then you need to get that package
###But Don't actually do this step (if you are copy/pasting along)
## source("http://bioconductor.org/BiocLite.R")
## biocLite("human.db0")
#
###Then you can get a tab-delimited file that has your probes paired with IDs
#hcg110_IDs = system.file("extdata", "hcg110_ID", package="AnnotationDbi")
#head(read.delim(hcg110_IDs,header=FALSE))
#
###For this example lets not actually write anything to the file sys.
#tmpout = tempdir()
#
###Then you can make the package
#makeHUMANCHIP_DB(affy=FALSE,
#                 prefix="hcg110",
#                 fileName=hcg110_IDs,
#                 baseMapType="gb",
#                 outputDir = tmpout,
#                 version="1.0.0",
#                 manufacturer = "Affymetrix",
#                 chipName = "Human Cancer G110 Array",
#                 manufacturerUrl = "http://www.affymetrix.com")
#
#
#
#
####################################################
#### chunk number 5: GOExamples
####################################################
###You may have already noticed that the organism packages have some GO
###information in them already.  This mapping represents the relationship
###between these EG IDs and the GO IDs.  For example: org.Hs.egGO,
###org.Hs.egGO2EG, and org.Hs.egGO2ALLEGS
#ls("package:org.Hs.eg.db")[22:24]
###There are two types of such mappings. org.Hs.egGO will map GO terms to
###entrez gene IDs, while org.Hs.egGO2ALLEGS maps GO terms and relevant child
###terms to specific entrez gene IDs The man pages will help you remember which
###is which.
#
#
###All the other GO information is found in GO.db
#library(GO.db)
#ls("package:GO.db")
#
###The mapping that you usually want is the one that describes all the terms
#keys = keys(GOTERM[1:500])
#x = mget(as.character(keys), GOTERM, ifnotfound=NA)
#x[30]
#
###GO is a directed acyclic graph, so there are many parent and child
###relationships among terms.
###Therefore GO.db also has mappings to tell you about the parent and child
###terms as well as all the ancestor or offspring terms
#mget(as.character(names(x[30])),GOBPCHILDREN, ifnotfound=NA)
#
#
#
#
####################################################
#### chunk number 6: biomaRtExamples
####################################################
###Getting the data from biomaRt:
#
#library("biomaRt")
###Choose a database
#listMarts()[1:5,]
#
###Get the current ensembl database.
#ensembl = useMart("ensembl")
#
###List the datasets therein
#listDatasets(ensembl)[1:10,]
###Then set up so that you use that for this session
###(we will choose the mouse one from NCBI build 37.1):
#ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
#
###List attributes
#attributes = listAttributes(ensembl)
#attributes[1:10,]
#
###And filters
#filters = listFilters(ensembl)
#filters[1:10,]
#
###Some entrez gene IDs
#EGs = c("18392","18414","56513")
#
###1st a Simple example to just get some gene names:
#getBM(attributes = "external_gene_id",
#      filters = "entrezgene",
#      values = EGs,
#      mart=ensembl)
#
#
#
####################################################
#### chunk number 7: biomartDemoContinued
####################################################
###Transcript starts and ends:
#getBM(attributes = c("entrezgene","transcript_start","transcript_end"),
#      filters = "entrezgene",
#      values = EGs,
#      mart=ensembl)
#
#
#
####################################################
#### chunk number 8: SessionInfo
####################################################
#sessionInfo()
#
#
#
