qc.stats = function(snpstats, samplestats, single.chr = T) {
	s = read.table(snpstats, hea = T)
	s1 = read.table(samplestats, hea = T)

	par(mfcol = c(2,2+2*single.chr))

	if(single.chr) plot(s$position, s$MAF, xlab = "position", ylab = "MAF", main = "SNP MAF")
	hist(s$MAF, n = 20, xlab = "MAF", main = "Histogram of MAF")

	if(single.chr) plot(s$position, s$missing, xlab = "position", ylab = "SNP Missing proportion", main = "SNP missing rate", ylim = c(0,1))
	bb = c(seq(0, 0.1, 0.02), 0.15, 0.3, 0.5, 0.75, 1)
	d = hist(s$missing,breaks = bb , plot = F)
	d1 = d$counts
	plot(0:1, 0:1, type ="n", axes = F, xlab = "", ylab = "", main = "SNP Missing rate counts")

	for(i in 1:10) text(0.3, 0+(i-1)*.1, paste("[",bb[i],"-",bb[i+1],"] : ", d1[i], sep=""))
	  if(single.chr) plot(s$position, -log10(s$HWE), xlab = "position", ylab = "HWE -log10 p-value", main = "SNP HWE")
	hist(-log10(s$HWE), n = 20, xlab = "HWE -log10 p-value", main = "Histogram of HWE -log10 p-value")

	plot(s1$missing, s1$het, xlab = "Missing proportion", ylab = "Heterozygosity", main = "Missing vs Heterozygosity", xlim = c(0,1), ylim = c(0,1))
}


# Example usage:
# qc.stats("stats.txt", "stats1.txt", single.chr = T)

# Or
# qc.stats("stats.txt", "stats1.txt", single.chr = F)

