plotIRanges<- function(ir){
  bins <- disjointBins(IRanges(start(ir), end(ir) + 1))
  dat <- cbind(as.data.frame(ir), bin = bins)
  ggplot(dat) +
  geom_rect(aes(xmin = start, xmax = end,ymin = bin, ymax = bin + 0.5)) +
  theme_bw()
}
