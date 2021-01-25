##################################################################
## Jesse Engreitz
## January 23, 2021
## Script to plot statistics for prime editor pegRNAs:
##   distance from nick to edit
##   RT template lengths
##   pegRNAs per variant
##   pegRNAs per position
##   GC content of RT template
##   GC content of primer binding site

suppressPackageStartupMessages(library("optparse"))

## parse options
option.list <- list(
  make_option("--input", type="character", help="pegRNA design table"),
  make_option("--edits", type="character", help="Variant table"),
  make_option("--outfile", type="character", help="Basename for output files"),
  make_option("--figureLabel", type="character", default="", help="Sample name or other label to include in figure title")
)
opt <- parse_args(OptionParser(option_list=option.list))


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))


pegs <- read.delim(opt$input, check.names=F, stringsAsFactors=F)
variants <- read.delim(opt$edits, check.names=F, stringsAsFactors=F)


mytheme <- theme_classic() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15))


pdf(file=paste0(opt$outfile), width=4, height=4)

## Distance
p <- pegs %>% ggplot(aes(x=editPositionRelativeToNick0Based)) + geom_histogram(binwidth=1, na.rm=T) 
p <- p + ggtitle(opt$figureLabel) + xlim(0,NA) + xlab("Distance from nick to variant (bp)") + mytheme + ylab("# pegRNAs")
print(p)

## RT template lengths
p <- pegs %>% ggplot(aes(x=rtTemplateLength)) + geom_histogram(binwidth=1, na.rm=T) 
p <- p + ggtitle(opt$figureLabel) + xlim(0,NA) + xlab("Length of RT template (nt)") + mytheme + ylab("# pegRNAs")
print(p)

## PBS lengths
p <- pegs %>% ggplot(aes(x=pbsLength)) + geom_histogram(binwidth=1, na.rm=T) 
p <- p + ggtitle(opt$figureLabel) + xlim(0,NA) + xlab("Length of primer binding site (nt)") + mytheme + ylab("# pegRNAs")
print(p)

## GC content of RT templates
p <- pegs %>% ggplot(aes(x=rtGCpct)) + geom_histogram(binwidth=1, na.rm=T) 
p <- p + ggtitle(opt$figureLabel) + xlim(0,100) + xlab("GC of RT template (%)") + mytheme + ylab("# pegRNAs")
print(p)

## GC content of PBS
p <- pegs %>% ggplot(aes(x=pbsGCpct)) + geom_histogram(binwidth=5, na.rm=T) 
p <- p + ggtitle(opt$figureLabel) + xlim(0,100) + xlab("GC of primer binding site (%)") + mytheme + ylab("# pegRNAs")
print(p)

## Spacer usage
p <- pegs %>% group_by(spacer) %>% summarise(n=n()) %>% ggplot(aes(x=n)) + geom_histogram(binwidth=1, na.rm=T) 
p <- p + ggtitle(opt$figureLabel) + xlim(0,NA) + xlab("pegRNAs per spacer") + mytheme + ylab("# spacers")
print(p)

## Variant coverage
pegs <- pegs %>% mutate(variantName=factor(variantName, levels=unique(variants$name)))
p <- table(pegs$variantName) %>% as.data.frame() %>% ggplot(aes(x=Freq)) + geom_histogram(binwidth=1, na.rm=T) 
p <- p + ggtitle(opt$figureLabel) + xlab("pegRNAs per variant") + mytheme + ylab("# variants")
print(p)

## Position coverage
variants <- variants %>% rowwise() %>% mutate(positionName=paste0(chr,":",start,"-",end))
pegs <- pegs %>% rowwise() %>% mutate(positionName=paste0(chr,":",variantStart,"-",variantEnd)) %>% as.data.frame() 
pegs <- pegs %>% mutate(positionName=factor(positionName, levels=unique(variants$positionName)))
p <- table(pegs$positionName) %>% as.data.frame() %>% ggplot(aes(x=Freq)) + geom_histogram(binwidth=1, na.rm=T) 
p <- p + ggtitle(opt$figureLabel) + xlab("pegRNAs per variant position") + mytheme + ylab("# positions")
print(p)

## Gaps per region
regions <- table(pegs$variantName) %>% as.data.frame() %>% rename(variantName=Var1) %>% merge(variants %>% select(name,region), by.x="variantName", by.y="name", all.y=TRUE)
if (length(unique(regions$Freq)) > 1) {
  regions <- regions %>% mutate(Covered=factor(Freq, labels=c("No","Yes")))
  p <- regions %>% ggplot(aes(x=region, fill=Covered)) + geom_bar(position="stack", na.rm=T) 
  p <- p + ggtitle(opt$figureLabel) + xlab("Design regions") + mytheme + ylab("# variants")
  p <- p + theme(axis.text.x = element_text(angle = 45, size=8, vjust=0.90, hjust=1))
  p <- p + scale_fill_manual(values=c(Yes="black",No="red"))
  print(p)
}

## Counts by priority group
if ("priority" %in% colnames(pegs)) {
  regions <- merge(regions, pegs %>% group_by(variantName) %>% summarise(priority=min(priority)), by="variantName", all.x=TRUE)
  regions <- regions %>% mutate(SourcePriority=factor(replace_na(priority,"None")))
  p <- regions %>% ggplot(aes(x=region, fill=SourcePriority)) + geom_bar(position="stack", na.rm=T) 
  p <- p + ggtitle(opt$figureLabel) + xlab("Design regions") + mytheme + ylab("# variants")
  p <- p + theme(axis.text.x = element_text(angle = 45, size=8, vjust=0.90, hjust=1))
  print(p)
}

invisible(dev.off())


save.image(paste0(opt$outfile,".RData"))

