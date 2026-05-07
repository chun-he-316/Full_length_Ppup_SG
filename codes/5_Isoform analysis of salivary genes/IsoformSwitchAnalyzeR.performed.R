library(IsoformSwitchAnalyzeR)
salmonQuant <- importIsoformExpression(parentDir = "salmon_results", addIsofomIdAsColumn = TRUE)
str(salmonQuant)
samples <- c("PpCAa", "PpCAb", "PpCAc", "PpSgGDa", "PpSgGDb", "PpSgGDc","PpSgPDa", "PpSgPDb", "PpSgPDc")
condictions <- c(rep_len("PpCA", 3), rep_len("PpSgGD", 3),rep_len("PpSgPD", 3))
designM <- data.frame(sampleID=samples, condition=condictions)
designM

switchList <- importRdata(isoformCountMatrix=salmonQuant$counts, 
                          isoformRepExpression=salmonQuant$abundance, 
                          designMatrix=designM, 
                          isoformExonAnnoation="Ppup.final.gtf", 
                          isoformNtFasta="transcripts.fasta",
                          fixStringTieAnnotationProblem = TRUE)
sar2 <- preFilter(switchList)
# default parameters
# geneExpressionCutoff = 1 
# isoformExpressionCutoff = 0 
# IFcutoff=0.01
# dIFcutoff = 0.1
# removeSingleIsoformGenes = TRUE
sar3 <- isoformSwitchTestDEXSeq(sar2)

ggplot(data=sar3$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.4 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()

write.csv(
  sar3$isoformFeatures,
  file = "isoformFeatures.csv",
  row.names = FALSE
)
summary(switchList)
