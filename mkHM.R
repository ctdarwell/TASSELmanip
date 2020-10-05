library(LDheatmap)
library(vcfR)
library(snpStats)
require(grid)

#loop cycles thru VCFS for 12 rice chromosomes
for (i in c(1,2,3,4,5,6,7,8,9,10,11,12)){
  vcf <- read.vcfR( paste('chrom_', i, '.vcf', sep=''), verbose = FALSE )
  vcf <- vcfR2SnpMatrix(vcf)
  colnames(vcf$data) <- vcf$genetic.distance

  MyHeatmap <- LDheatmap(vcf$data, vcf$genetic.distance,
          color=heat.colors(100), flip = T, title = '', geneMapLocation=0.2, geneMapLabelX=0) #NB geneMapLabelX=5 is a fudge to rid the scale bar

  write.csv(MyHeatmap$LDmatrix, file = paste('chrom_', i, '.csv', sep=''), quote = F)
}
