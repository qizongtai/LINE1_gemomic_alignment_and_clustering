#!/usr/bin/Rscript
# use simplified cytoband.txt (cytoband_hg19_2), which doesn't have grey regions in chr plot
# plot genome-wide distribution of features
# single-color chr (all white), can be changed in line 58-60
# ====++============++======+==+=========
# =======+=======+==+=+==============
# ========+===============++===
# =++=====+=========+====
# ======+============

library(getopt);

opt.spec = matrix(ncol=5, byrow=T, data=c(
    'input',      'i', '1', 'character', 'input features in BED format',
    'outpdf',     'o', '1', 'character', 'output PDF plot',
    'cytoband',   'c', '2', 'character', 'cytoband',
    'help',       'h', '0', 'logical',   'print this message'
));
script.name = sub('^.*/', '', strsplit(commandArgs(FALSE)[4],"=")[[1]][2]);
usage = getopt(spec=opt.spec, command=script.name, usage=T);
opt = getopt(spec=opt.spec);

# handle command line options
if (!is.null(opt$help)) {
    cat(usage)
    q(status=0)
}

if (is.null(opt$input)) {
    cat('missing --input\n');
    q(status=1);
}

if (is.null(opt$outpdf)) {
    cat('missing --outpdf\n');
    q(status=1);
}

if (is.null(opt$cytoband)) {
    cat('missing --cytoband\n');
    q(status=1);
}

cytobands = read.table(opt$cytoband, header=F, sep='\t', as.is=T, quote='');
colnames(cytobands) = c('chr', 'start', 'end', 'band');
nband = nrow(cytobands);
cytobands$chr = sub('^chr', '', cytobands$chr);
autosome = grep('^[0-9]*$',cytobands$chr, value =T);
max.autosome = max (as.integer(autosome));
num.x = max.autosome + 1;
num.y = max.autosome + 2;
num.chr = max.autosome + 2;
cytobands$chr[cytobands$chr == 'X'] = num.x;
cytobands$chr[cytobands$chr == 'Y'] = num.y;
cytobands$chr = as.integer(cytobands$chr);

kneg = cytobands$band == 'gneg';
kpos = grepl('^gpos', cytobands$band);
colors = rep('black', nband);
colors[kneg] = 'white';
colors[kpos] = 'white';
width = rep(0.1, nband);
width[kneg] = 0.15;
width[kpos] = 0.15;

features = read.table(opt$input, header=F, sep='\t', as.is=T, quote='')[,1:3];
colnames(features) = c('chr', 'start', 'end');
features$chr = sub('^chr', '', features$chr);
features$chr[features$chr == 'X'] = num.x;
features$chr[features$chr == 'Y'] = num.y;
features$chr = as.integer(features$chr); 
#NAs will be introduced if features exist in regions [neither autosomes(chr1-22) or sex chr(X/Y)]
features = features[complete.cases(features[,1]),];
#remove NAs in the dataframe

pdf(file=opt$outpdf);
par(las=1, tcl=0.3, mar=c(2,4.5,2,2), cex=0.8);
plot(c(0,max(cytobands$end)), c(1,num.chr), ylim=c(num.chr,1), type='n', axes=F, ann=F, bty='n');
axis(side=2, at=1:num.chr, labels=c(1:max.autosome,'X','Y'), las=1, col='0');
rect(xleft=cytobands$start, xright=cytobands$end, ytop=cytobands$chr+width, ybottom=cytobands$chr-width, col=colors, lwd=0.3);
rect(xleft=features$start, xright=features$end, ytop=features$chr+0.13, ybottom=features$chr-0.13, border=NA, col='#FF0000FF',lwd=0.1);
dev.off();
