
qc.snpstats = function( snpstats, single.chr = T ) {
	s = read.table(snpstats, header = T)
	par(mfcol = c(1,4+2*single.chr))
	if(single.chr) plot(s$position, s$MAF, xlab = "position", ylab = "MAF", main = "SNP MAF")
	plot(s$position, s$information, xlab = "position", ylab = "Information", main = "SNP Information" )
	hist(s$MAF, n = 100, xlab = "MAF", main = "Histogram of MAF")
	bb = c(seq(0, 0.1, 0.01), 0.15, 0.3, 0.5, 0.75, 1)
	d = hist(s$missing,breaks = bb , plot = F)
	d1 = d$counts
	plot(0:1, 0:1, type ="n", axes = F, xlab = "", ylab = "", main = "SNP Missing rate counts")
	for(i in 1:10) text(0.3, 0+(i-1)*.1, paste("[",bb[i],"-",bb[i+1],"] : ", d1[i], sep=""))
	  if(single.chr) plot(s$position, s$HWE, xlab = "position", ylab = "HWE -log10 p-value", main = "SNP HWE")
	hist(s$HWE, n = 20, xlab = "HWE -log10 p-value", main = "Histogram of HWE -log10 p-value")
}

qc.stats = function(snpstats, samplestats, single.chr = T) {
	s = read.table(snpstats, hea = T)
	s1 = read.table(samplestats, hea = T)

	par(mfcol = c(3,2+2*single.chr))

	if(single.chr) plot(s$position, s$MAF, xlab = "position", ylab = "MAF", main = "SNP MAF")
	plot(s$position, s$information, xlab = "position", ylab = "Information", main = "SNP Information" )
	hist(s$MAF, n = 20, xlab = "MAF", main = "Histogram of MAF")

	if(single.chr) plot(s$position, s$missing, xlab = "position", ylab = "SNP Missing proportion", main = "SNP missing rate", ylim = c(0,1))
	plot(s$position, s$information, xlab = "position", ylab = "SNP Information", main = "SNP Information")
	bb = c(seq(0, 0.1, 0.02), 0.15, 0.3, 0.5, 0.75, 1)
	d = hist(s$missing,breaks = bb , plot = F)
	d1 = d$counts
	plot(0:1, 0:1, type ="n", axes = F, xlab = "", ylab = "", main = "SNP Missing rate counts")

	for(i in 1:10) text(0.3, 0+(i-1)*.1, paste("[",bb[i],"-",bb[i+1],"] : ", d1[i], sep=""))
	  if(single.chr) plot(s$position, s$HWE, xlab = "position", ylab = "HWE -log10 p-value", main = "SNP HWE")
	hist(s$HWE, n = 20, xlab = "HWE -log10 p-value", main = "Histogram of HWE -log10 p-value")

	plot(s1$missing, s1$heterozygosity, xlab = "Missing proportion", ylab = "Heterozygosity", main = "Missing vs Heterozygosity", xlim = c(0,1), ylim = c(0,1))
}


# Example usage:
# qc.stats("stats.txt", "stats1.txt", single.chr = T)

# Or
# qc.stats("stats.txt", "stats1.txt", single.chr = F)

