setwd('C:\\Users/xbw-1/Documents/GitHub/QLS_Projects/assignment1B/')
library('Rsamtools')
bamFile=scanBam('17.41000000-42000000.HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam')

# names in the bam
names(bamFile[[1]])

lst <- lapply(names(bamFile[[1]]), function(elt) {
  do.call(c, unname(lapply(bamFile, "[[", elt)))
})
names(lst) <- names(bamFile[[1]])

bamFile_df=do.call("DataFrame", lst)