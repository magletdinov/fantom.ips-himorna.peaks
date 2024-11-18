options(warn=-1)

library("ChIPpeakAnno")
library("biomaRt")
library("GenomicFeatures")

peaksAnno<-function(filename, output_filename) {
    txdb <- makeTxDbFromGFF('/export/home/public/agletdinov_shared/annotations/FANTOM_CAT.lv2_permissive.gtf.gz')
    
    anno <- toGRanges(txdb, format='gene')
    gr1 <- toGRanges(filename, format="BED", skip=1)

    overlaps.anno <- annotatePeakInBatch(gr1, AnnotationData=anno, output="overlapping", maxgap=1000L)
    print("count annotation, write to file")
    write.table(overlaps.anno, file=output_filename, quote=FALSE, sep="\t") 
}

peaksAnno("/export/home/agletdinov/work/git_projects/fantom.ips-himorna.peaks/preprocessing/methyl_bed/methyl_concat/methyl_concat.bed", "/export/home/agletdinov/work/git_projects/fantom.ips-himorna.peaks/preprocessing/methyl_bed/methyl_concat/peaks_anno.csv")