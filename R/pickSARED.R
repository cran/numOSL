#####
pickSARED <- 
function(obj, rsdED.limit=NULL, rcy.interval=NULL, 
         rcp1.limit=NULL, rcp2.limit=NULL, method=NULL, 
         fom.limit=NULL, rcs.limit=NULL, outfile=NULL) {
    UseMethod("pickSARED")
} #
### 2016.07.14.
pickSARED.default <- 
function(obj, rsdED.limit=NULL, rcy.interval=NULL, 
         rcp1.limit=NULL, rcp2.limit=NULL, method=NULL, 
         fom.limit=NULL, rcs.limit=NULL, outfile=NULL) {
    ### Stop if not.
    stopifnot(class(obj)=="calSARED", is.list(obj), length(obj)==10L,
              all(names(obj) %in% c("LMpars","N","failFit.NO","saturate.NO","failED.NO",
              "extrapolate.NO","largeRcy.NO","largeRcp5.NO","largeRcp10.NO","tab")),
              is.null(rsdED.limit) || (is.numeric(rsdED.limit) && length(rsdED.limit)==1L),
              is.null(rcy.interval) || (is.numeric(rcy.interval) && length(rcy.interval)==2L),
              is.null(rcp1.limit) || (is.numeric(rcp1.limit) && length(rcp1.limit)==1L),
              is.null(rcp2.limit) || (is.numeric(rcp2.limit) && length(rcp2.limit)==1L),
              is.null(method) || (method %in% c("Inter", "Extra") && length(method)==1L),
              is.null(fom.limit) || (is.numeric(fom.limit) && length(fom.limit)==1L),
              is.null(rcs.limit) || (is.numeric(rcs.limit) && length(rcs.limit)==1L),
              is.null(outfile) || (is.character(outfile) && length(outfile)==1L))
    ###
    ###
    originGrainNumber <- as.numeric(substr(rownames(obj$tab), start=9L, stop=10000L))
    ###
    if (!is.null(rsdED.limit)) {
        indexValue <- obj$tab[,"rsdED",drop=TRUE]
        selected <- which(indexValue<rsdED.limit)
        tab <- obj$tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
    } else {
        tab <- obj$tab
    } # end if.
    ###
    if (!is.null(rcy.interval)) {
        indexValue <- tab[,"RecyclingRatio",drop=TRUE]
        selected <- which(indexValue>rcy.interval[1L] & indexValue<rcy.interval[2L])
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
    } # end if.
    ###
    if (!is.null(rcp1.limit)) {
        indexValue <- tab[,"Recuperation1",drop=TRUE]
        selected <- which(indexValue<rcp1.limit)
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
    } # end if.
    ###
    if (!is.null(rcp2.limit)) {
        indexValue <- tab[,"Recuperation2",drop=TRUE]
        selected <- which(indexValue<rcp2.limit)
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
    } # end if.
    ###
    if (!is.null(method)) {
        indexValue <- tab[,"Method",drop=TRUE]
        selected <- which(indexValue==method)
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
    } # end if.
    ###
    if (!is.null(fom.limit)) {
        indexValue <- tab[,"FOM",drop=TRUE]
        selected <- which(indexValue<fom.limit)
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
    } # end if.
    ###
    if (!is.null(rcs.limit)) {
        indexValue <- tab[,"RCS",drop=TRUE]
        selected <- which(indexValue<rcs.limit)
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
    } # end if.
    ###
    if (!is.null(outfile)) {
        write.csv(tab, file=paste(outfile,".csv",sep=""))
    } # end if.
    ###
    ###
    rejectGrainNumber <- originGrainNumber[which(!originGrainNumber %in% 
        as.numeric(substr(rownames(tab), start=9L, stop=10000L)))]
    if (length(rejectGrainNumber)==0L) rejectGrainNumber <- NULL
    ###
    ###
    output <- list("sarED"=tab[,c(1L,2L),drop=FALSE],
                   "reject.NO"=rejectGrainNumber)
    return(output)
    ###
} # end function pickSARED.default.
