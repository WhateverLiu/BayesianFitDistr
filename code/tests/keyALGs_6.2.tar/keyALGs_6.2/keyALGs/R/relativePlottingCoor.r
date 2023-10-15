

findPlotCoor = function(relaX = 0.5, relaY = 0.5, realX = NULL, realY = NULL)
{
  if(!is.null(relaX)) return(c(relaX, relaY))
  x = diff(par("usr")[1:2]) * relaX + par("usr")[1]
  y = diff(par("usr")[3:4]) * relaY + par("usr")[3]
  c(x, y)
}




















