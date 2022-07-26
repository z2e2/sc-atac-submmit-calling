library(rtracklayer)
library(GenomicRanges)
library(plyr)
library(tidyverse)
library(data.table)

library(biosignals)


root = '/Genomics/pritykinlab/zzhao/sc-atac-submmit-calling/sc-atac-submmit-calling'
motif_regression_out = "results/motif_regression_output/"
motif_regression_out = paste(root, motif_regression_out, sep="/")
# load peaks from archR, convert the format to bed first.
all.atac.peaks <- rtracklayer::import(paste(motif_regression_out, "MACS2_peaks.bed", sep="/"))

names(all.atac.peaks) <- sapply(1:length(all.atac.peaks), function(i) paste("peak", i, sep=""))
                                
                                
                                
atac.signal <- rtracklayer::import(paste(motif_regression_out, "sc_atac_fragments.bw", sep="/"), as = "RleList")
                                
peaks.signal <- atac.signal[all.atac.peaks]
names(peaks.signal) <- names(all.atac.peaks)
gen_kernels <- function(bw) {
    g0 <- generateKernel("gaussian", bandwidth = bw, deriv = 0)
    g1 <- generateKernel("gaussian", bandwidth = bw, deriv = 1)
    g2 <- generateKernel("gaussian", bandwidth = bw, deriv = 2)
    list("g0"=g0, "g1"=g1, "g2"=g2)
}
# summit calling using biosignals package from Leslie Lab: https://bitbucket.org/leslielab/biosignals
findVectorMax <- function(v, g0, g1, g2, min.dist = 150) {
    # pad zeros at the beginning and the end of the vector
    bw = as.integer(length(g0)/2)
    v = c(rep(0, bw), as.numeric(v), rep(0, bw))
    v.conv <- convolve1d(v, g0)[(bw+1):(length(v)-bw)]
    v.conv.d1 <- convolve1d(v, g1)[(bw+1):(length(v)-bw)]
    v.conv.d2 <- convolve1d(v, g2)[(bw+1):(length(v)-bw)]
    v.conv.d1.zero <- zeroCrossings(v.conv.d1)
    points.max <- v.conv.d1.zero[v.conv.d2[v.conv.d1.zero] < 0]
    if (length(points.max) == 0) {
        return(c())
    }
    points.max <- data.table(x = points.max)
    points.max$value <- v.conv[points.max$x]
    points.max <- points.max[order(-value)]
    vector.max <- c()
    for (i in 1:nrow(points.max)) {
        x <- points.max[i, ][, x]
        # instead of compare against the 0 and length of v, 
        if (min(abs(x - c(vector.max, -min.dist, (length(v.conv)+min.dist)))) >= min.dist) {
            vector.max <- c(vector.max, x)
        }
    }
    vector.max
}
kernels <- gen_kernels(150)
peaks.summits <- sapply(peaks.signal, function(peak) findVectorMax(peak, kernels$g0, kernels$g1, kernels$g2))
summit.count <- sapply(peaks.summits, length)
peaks.s0.signal <- peaks.signal[summit.count == 0]
length(peaks.s0.signal)
saveRDS(peaks.summits, paste(motif_regression_out, "SUMMITS-TMP-all-atac-peaks-summits-lists.rds", sep="/"))
                        
peaks.summits.final <- peaks.summits
print("prepare final set of regions around peak summits")
summit.width <- 150
peaks.summits.dt <- data.table(ldply(peaks.summits.final, matrix, .id = "peak"))
colnames(peaks.summits.dt) <- c("peak", "summit")
peaks.summits <- all.atac.peaks[peaks.summits.dt$peak, ]                        

                        
stopifnot(start(peaks.summits) + peaks.summits.dt$summit <= end(peaks.summits))
start(peaks.summits) <- start(peaks.summits) + peaks.summits.dt$summit
end(peaks.summits) <- start(peaks.summits)                                
                                
peaks.summits$peak <- peaks.summits$name
peaks.summits <- sort(peaks.summits)
peaks.summits$name <- paste0("summit", seq_along(peaks.summits))

names(peaks.summits) <- peaks.summits$name
peaks.summits <- resize(peaks.summits, width = summit.width, fix = "center")

write.table(as.data.frame(peaks.summits), paste(motif_regression_out, "SUMMITS-peak_summits.more_info.bed", sep="/"), quote= F, row.names = F, sep ="\t")


export(peaks.summits, paste(motif_regression_out, "SUMMITS-all-atac-summits.bed", sep="/"))

                                
                                
                                
                                
                                