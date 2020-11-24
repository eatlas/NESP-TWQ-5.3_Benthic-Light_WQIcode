WQ_reportCardColors = c('#00734D','#B0D235','#F0C918','#F47721','#ED1C24')
WQ_reportCardColorPalette = colorRampPalette(WQ_reportCardColors)

WQ_gradeBoundaries.mmp = c(0.5 + (2/3*0.5), 0.5+1/3*0.5, 0.5, 0.25,0)
WQ_gradeBoundaries.ghhp = c(0.85, 0.65, 0.5, 0.25,0)
WQ_gradeBoundaries.uniform = c(0.8, 0.6, 0.4, 0.2,0)

WQ_gradeBoundaries = function(type='MMP') {
    gradeBoundaries=NA
    if (type=='MMP') gradeBoundaries=c(1,WQ_gradeBoundaries.mmp)
    if (type=='GHHP') gradeBoundaries=c(1,WQ_gradeBoundaries.ghhp)
    if (type=='Uniform') gradeBoundaries=c(1,WQ_gradeBoundaries.uniform)
    return(gradeBoundaries)
}
WQ_generateGrades <- function(x,type='MMP') {
    if (!type %in% c('Uniform','MMP','GHHP')) warning('type must be either MMP, GHHP or Uniform')
    if (type=='MMP')
        g=ifelse(is.na(x),'NA',ifelse(x>=0.5 + (2/3*0.5), 'A', ifelse(x>=0.5+1/3*0.5, 'B', ifelse(x>=0.5, 'C',  ifelse(x>=0.25, 'D', 'E')))))
    if (type=='GHHP')
        g=ifelse(is.na(x),'NA',ifelse(x>=0.85, 'A', ifelse(x>=0.65, 'B', ifelse(x>=0.5, 'C',  ifelse(x>=0.25, 'D', 'E')))))
    if (type=='Uniform')
        g=ifelse(is.na(x),'NA',ifelse(x>=0.8, 'A', ifelse(x>=0.6, 'B', ifelse(x>=0.4, 'C',  ifelse(x>=0.2, 'D', 'E')))))
    return(g)
}

WQ_gradeMids = function(type='MMP') {
    gradeBoundaries=NA
    if (type=='MMP') gradeBoundaries=zoo:::rollmean(c(1,WQ_gradeBoundaries.mmp),2)
    if (type=='GHHP') gradeBoundaries=zoo:::rollmean(c(1,WQ_gradeBoundaries.ghhp),2)
    if (type=='Uniform') gradeBoundaries=zoo:::rollmean(c(1,WQ_gradeBoundaries.uniform),2)
    return(gradeBoundaries)
}

##################################################################
## The following function generates a water year vector         ##
## The water year is defined as 1st Oct through to 30 September ##
## input:                                                       ##
##    Dt: a Date vector                                         ##
## output:                                                      ##
##    waterYear: a integer representing the water year          ##
##################################################################
WQ_waterYear <- function(Dt) {
    as.numeric(as.character(format(Dt+(as.Date("1970-12-31")-as.Date("1970-10-01")+1), format="%Y")))
}

############################################################################
## The following function generates a season vector                       ##
## The Wet Season is defined as Nov-Apr                                   ##
## input:                                                                 ##
##    Dt: a Date vector                                                   ##
## output:                                                                ##
##    Season: a character vector representing the Season ("Wet" or "Dry") ##
############################################################################
WQ_season <- function(Dt) {
    factor(ifelse(lubridate:::month(Dt,label=TRUE) %in% c("Nov","Dec","Jan","Feb","Mar","Apr"), "Wet","Dry"))
    }

