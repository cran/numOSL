#####
pickBINdata <-
function(obj_loadBIN, Position=NULL, Grain=NULL, Run=NULL, Set=NULL,
         DType=NULL, IRRTime=NULL, NPoints=NULL, LType=NULL, Low=NULL,
         High=NULL, Rate=NULL, Temperature=NULL, Delay=NULL, On=NULL, 
         Off=NULL, LightSource=NULL, AnTemp=NULL, TimeSinceIrr=NULL, 
         view=TRUE, manual.select=FALSE, force.matrix=FALSE) {
    UseMethod("pickBINdata")
} #
### 2017.01.11.
pickBINdata.default <- 
function(obj_loadBIN, Position=NULL, Grain=NULL, Run=NULL, Set=NULL,
         DType=NULL, IRRTime=NULL, NPoints=NULL, LType=NULL, Low=NULL,
         High=NULL, Rate=NULL, Temperature=NULL, Delay=NULL, On=NULL, 
         Off=NULL, LightSource=NULL, AnTemp=NULL, TimeSinceIrr=NULL, 
         view=TRUE, manual.select=FALSE, force.matrix=FALSE) {
    ### Stop if not.
    stopifnot(class(obj_loadBIN)=="loadBIN", 
              names(obj_loadBIN)==c("records","tab"),
              is.null(Position) || is.numeric(Position),
              is.null(Grain) || is.numeric(Grain),
              is.null(Run) || is.numeric(Run),
              is.null(Set) || is.numeric(Set),
              is.null(DType) || is.character(DType),    
              is.null(IRRTime) || is.numeric(IRRTime), 
              is.null(NPoints) || is.numeric(NPoints),
              is.null(LType) || is.character(LType),
              is.null(Low) || is.numeric(Low),
              is.null(High) || is.numeric(High),
              is.null(Rate) || is.numeric(Rate),
              is.null(Temperature) || is.numeric(Temperature),
              is.null(Delay) || is.numeric(Delay),
              is.null(On) || is.numeric(On),
              is.null(Off) || is.numeric(Off),                  
              is.null(LightSource) || is.character(LightSource), 
              is.null(AnTemp) || is.numeric(AnTemp),
              is.null(TimeSinceIrr) || is.numeric(TimeSinceIrr),            
              is.logical(view), length(view)==1L, 
              is.logical(manual.select), length(manual.select)==1L,
              is.logical(force.matrix), length(force.matrix)==1L)

    ###
    search_fun <- function(x,y) any(abs(x-y)<=.Machine$double.eps^0.25)
    ###
    ### Search disc Position (numeric vector).
    tab <- obj_loadBIN$tab
    if (!is.null(Position))  {
        all_value <- tab[,"Position",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, Position)) 
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search Grain number (numeric vector).
    if (!is.null(Grain))  {
        all_value <- tab[,"Grain",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, Grain)) 
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search Run number (numeric vector).
    if (!is.null(Run))  {
        all_value <- tab[,"Run",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, Run)) 
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search Set number (numeric vector).
    if (!is.null(Set))  {
        all_value <- tab[,"Set",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, Set))
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search data Type  (character vector).
    if (!is.null(DType))  {
        all_value <- tab[,"DType",drop=TRUE]
        select_index <- which(all_value %in% DType)
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search irradiation Time(s) (numeric vector).
    if (!is.null(IRRTime))  {
        all_value <- tab[,"IRRTime",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, IRRTime))
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search number of points (numeric vector).
    if (!is.null(NPoints))  {
        all_value <- tab[,"NPoints",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, NPoints)) 
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search luminescence Type (characters vector).
    if (!is.null(LType))  {
        all_value <- tab[,"LType",drop=TRUE]
        select_index <- which(all_value %in% LType)
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search Low values (numeric vector).
    if (!is.null(Low))  {
        all_value <- tab[,"Low",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, Low)) 
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search High values (numeric vector).
    if (!is.null(High))  {
        all_value <- tab[,"High",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, High)) 
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search Rate values (numeric vector).
    if (!is.null(Rate))  {
        all_value <- tab[,"Rate",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, Rate))
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search Temperature values (numeric vector).
    if (!is.null(Temperature))  {
        all_value <- tab[,"Temperature",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, Temperature))
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search Delay values (numeric vector).
    if (!is.null(Delay))  {
        all_value <- tab[,"Delay",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, Delay))
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search On values (numeric vector).
    if (!is.null(On))  {
        all_value <- tab[,"On",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, On))
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search Off values (numeric vector).
    if (!is.null(Off))  {
        all_value <- tab[,"Off",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, Off))
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    
    ### Search light Source (character vector).
    if (!is.null(LightSource))  {
        all_value <- tab[,"LightSource",drop=TRUE]
        select_index <- which(all_value %in% LightSource)
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search AnTemp (numeric vector).
    if (!is.null(AnTemp))  {
        all_value <- tab[,"AnTemp",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, AnTemp)) 
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###

    ### Search TimeSinceIrr (numeric vector).
    if (!is.null(TimeSinceIrr))  {
        all_value <- tab[,"TimeSinceIrr",drop=TRUE]
        select_index <- which(apply(as.matrix(all_value), MARGIN=1L, search_fun, TimeSinceIrr))
        if (length(select_index)>0L) {
            tab <- tab[select_index,,drop=FALSE]
        } else {
            stop("Error: no BIN data satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    ###
    ###
    ID_name <- newTab <- c()
    Position_level <- sort(as.numeric(levels(factor(tab[,"Position",drop=TRUE]))))
    Grain_level <- sort(as.numeric(levels(factor(tab[,"Grain",drop=TRUE]))))
    row_name_tab <- as.numeric(row.names(tab))
    nloop <- 0L
    ###
    for (i in Position_level) {
        for (j in Grain_level) {
            Index <- which(tab[,"Position",drop=TRUE]==i & 
                           tab[,"Grain",drop=TRUE]==j)
            length_Index <- length(Index)
            ###
            if (length_Index>0L) {
                nloop <- nloop + 1L
                newTab <- rbind(newTab, 
                                cbind(Selected=TRUE,tab[Index,,drop=FALSE]), 
                                rep(x=" ", times=21L), 
                                rep(x=" ", times=21L))
                ###            
                ID_name <- c(ID_name, 
                             row_name_tab[Index], 
                             -9L*nloop, 
                             -9L*nloop-1L)
            } # end if.
            ###
        } # end for. 
    } # end for.
    ###
    ###
    row.names(newTab) <- ID_name
    ###
    if (manual.select==TRUE) {
        newTab <- edit(name=newTab, title="Select Records")
    } else if (view==TRUE) {
        View(x=newTab, title="Selected Records")
    } # end if.
    ###
    ###
    newTab <- newTab[newTab[,"Selected",drop=TRUE] %in% c("T","TRUE"),,drop=FALSE]
    if (nrow(newTab)==0L) stop("Error: no BIN data satisfies the given conditions!")
    ###
    row_name_newTab <- as.numeric(row.names(newTab))
    newRecords <- obj_loadBIN$records[row_name_newTab[row_name_newTab>0L]]
    ###
    ###
    newPosition_level <- sort(as.numeric(levels(factor(newTab[,"Position",drop=TRUE]))))
    newGrain_level <- sort(as.numeric(levels(factor(newTab[,"Grain",drop=TRUE]))))
    ###
    ###
    BINdata <- list()
    agID <- c()
    fail_matrix_agID <- c()
    nloop <- 0L
    ###
    for (i in newPosition_level) {
        for (j in newGrain_level) {
            newIndex <- which(as.numeric(newTab[,"Position",drop=TRUE])==i & 
                              as.numeric(newTab[,"Grain",drop=TRUE])==j)
            length_newIndex <- length(newIndex)
            ###============================================================
            if (length_newIndex>0L) {
                ###-------------------------------------------------
                if (force.matrix==TRUE) {
                    NPoints_vec <- sapply(newRecords[newIndex], attr, "NPoints")
                    Low_vec <- sapply(newRecords[newIndex], attr, "Low")
                    High_vec <- sapply(newRecords[newIndex], attr, "High")                 
                    Delay_vec <- sapply(newRecords[newIndex], attr, "Delay")
                    Off_vec <- sapply(newRecords[newIndex], attr, "Off")
                    ###
                    IRRTime_vec <- sapply(newRecords[newIndex], attr, "IRRTime")
                    ###
                    all_equal_NPoints <- all(abs(NPoints_vec-mean(NPoints_vec))<=.Machine$double.eps)
                    all_equal_Low <- all(abs(Low_vec-mean(Low_vec))<=.Machine$double.eps)
                    all_equal_High <- all(abs(High_vec-mean(High_vec))<=.Machine$double.eps)
                    all_equal_Delay <- all(abs(Delay_vec-mean(Delay_vec))<=.Machine$double.eps)
                    all_equal_Off <- all(abs(Off_vec-mean(Off_vec))<=.Machine$double.eps)
                    ###
                    if (all_equal_NPoints==TRUE &&
                        all_equal_Low==TRUE && 
                        all_equal_High==TRUE && 
                        all_equal_Delay==TRUE &&
                        all_equal_Off==TRUE) {
                        ###
                        nloop <- nloop + 1L
                        v_step <- (High_vec[1L]-Low_vec[1L])/NPoints_vec[1L]
                        x_values <- seq(from=v_step, to=High_vec[1L], by=v_step)
                        ij_matrix <- x_values
                        ###
                        for (k in seq(length_newIndex)) {
                            ij_matrix <- cbind(ij_matrix, newRecords[newIndex][[k]])
                        } # end for.                       
                        colnames(ij_matrix) <- c("x", paste("y.", seq(length_newIndex), sep=""))
                        ###
                        BINdata[[nloop]] <- ij_matrix
                        attr(BINdata[[nloop]], "IRRTime") <- IRRTime_vec
                        attr(BINdata[[nloop]], "agID") <- c("NO"=nloop, "Position"=i, "Grain"=j)
                        attr(BINdata[[nloop]], "Delay") <- Delay_vec[1L]
                        attr(BINdata[[nloop]], "Off") <- Off_vec[1L]
                    } else {
                        fail_matrix_agID <- rbind(fail_matrix_agID, c(nloop,i,j))
                    } # end if. 
                    ###
                } else {
                    nloop <- nloop + 1L
                    BINdata[[nloop]] <- newRecords[newIndex]
                    agID <- rbind(agID, c(nloop,i,j))
                } # end if.
                ###-------------------------------------------------
            } # end if.
            ###============================================================
        } # end for. 
    } # end for. 
    ###
    ###
    if (length(BINdata)==0L) BINdata <- NULL
    ###
    NPG <- function(x) paste("[NO=",x[1L],",Position=",x[2L],",Grain=",x[3L],"]",sep="")
    ###
    if(!is.null(fail_matrix_agID)) {
        rownames(fail_matrix_agID) <- NULL
        colnames(fail_matrix_agID) <- c("NO", "Position", "Grain")
        cat("ID failed in matrix transformation (unequal 'NPoints','Low','High','Delay',or 'Off'):\n")
        print(apply(fail_matrix_agID, MARGIN=1L, NPG))
    } # end if.
    ###        
    if (!is.null(agID)) {
        rownames(agID) <- NULL
        colnames(agID) <- c("NO","Position", "Grain")
    } # end if.
    output <- list("BINdata"=BINdata, "agID"=agID)
    class(output) <- "pickBIN"
    invisible(output)
    ###
} # end function pickBINdata.default.
#####
