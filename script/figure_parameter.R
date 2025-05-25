library(colorspace)

fontType <- "Helvetica"

fontSizeL <- 10 # pt
fontSizeM <- 8
fontSizeS <- 6

lineThick <- 0.75 # pt
lineMedium <- 0.5
lineThin <- 0.25

panelUnit <- 30 # mm
panelMargin <- 1.5

mmToInch <- 0.03937007874
mmToLineUnit <- 1/2.13
mmToLinePlotgarden <- 1/0.75
ptToMM <- 1/2.845


strong_red <- "#CB333A"
strong_blue <- "#4851A0"
weak_red <- lighten(strong_red, amount = 0.4)   # FF7D81
weak_blue <- lighten(strong_blue, amount = 0.4) # 8A91DD
no_grey <- "#A8A8A8"

strong_teel <- "#0892A5"
strong_green <- "#23CE6B" # A485
strong_darkgreen <- "#054A29"
strong_yellow <- "#FFBA49"
strong_orange <- "#F18F01" # dTAG
strong_lightpurple <- "#BD93D8"
strong_purple <- "#9E33CB" # Epi

panelSize <- function(num, unit = panelUnit, margin = panelMargin){
  return(num*unit - 2*margin)
}

genome <- "mm10"

label_kb_mb <- function(x) {
  ifelse(x >= 1000000, paste0(x / 1000000, "Mb"), paste0(x / 1000, "kb"))
}