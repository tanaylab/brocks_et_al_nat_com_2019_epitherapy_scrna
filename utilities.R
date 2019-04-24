#####################
#####################
# OPTIONS
#####################
options(dplyr.width = 200)
options(width=200)
options(stringsAsFactors = FALSE)
options(gmax.data.size = 1e9)
setwd('/home/davidbr/tmp')

#####################
# LIBRARIES
#####################
# getLoadedDLLs()
# dyn.unload()
shhh <- suppressPackageStartupMessages # It's a library, so shhh!

shhh(library("devtools"))
load_all("/home/davidbr/src/metacell/")
scdb_init("/home/davidbr/tmp/", force_reinit = TRUE)
scfigs_init("/home/davidbr/tmp/figs/")
shhh(require(tgconfig))
tgconfig::override_params("~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/metacell/tes/h1299_dacsb.yaml",package="metacell")

shhh(library(tidyverse))
shhh(library(ggsignif))
shhh(library(ape))
shhh(library(scales))
shhh(library(cowplot))
shhh(library(PWMEnrich))
shhh(library(pheatmap))
shhh(library(fields))
shhh(library(magrittr))
shhh(library(dplyr))
shhh(library(GenomicAlignments))
shhh(library(gpatterns))
shhh(library(ggthemes))
shhh(library(glue))
shhh(require(gridExtra))
shhh(library(Matrix))
shhh(library(Rsamtools))
shhh(library(ggseqlogo))
shhh(library(ggrepel))
shhh(library(rtracklayer))
shhh(require(devtools))
shhh(library(reshape2))
shhh(library(ggplot2))
shhh(library(ggrepel))
shhh(library(tglkmeans))
shhh(library(ggrastr))
shhh(library(Biostrings))
shhh(library(broom))
shhh(library("BSgenome.Hsapiens.UCSC.hg38"))
#library(ggbio)
#library("BSgenome.Mmusculus.UCSC.mm10")
shhh(library(RColorBrewer))
shhh(library(tgstat))
shhh(require(tgconfig))
gdb.init('/net/mraid14/export/tgdata/db/tgdb/hg38/trackdb')
rename <- dplyr::rename
count <- dplyr::count
slice <- dplyr::slice
ggsave <- ggplot2::ggsave
ensembl= as_tibble(fread('/net/mraid14/export/data/users/davidbr/general/data/ensembl/ensembl_gene_ids.txt')) 
te_anno_global = fread('~/hg38/rawdata/repeatmasker/repeat_anno.txt')
theme_set(theme_classic())

#####################
# FUNCTIONS
#####################
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

sparseToTidy <- function(sparse_mat, rowid = 'i', colid = 'j', valid = 'n'){
  tidy_tibble <- melt(as.matrix(sparse_mat), varnames = c(rowid, colid), value.name = valid) %>%
    as_tibble
  return(tidy_tibble)
}

tidyToSparse = function(df_tidy)
{
	colnames(df_tidy) = c('a', 'b', 'c')
	data = df_tidy %>%
		mutate_at(c('a', 'b'), funs(factor(.)))
		
	data_sparse = sparseMatrix(as.integer(data$a), as.integer(data$b), x = data$c)
	colnames(data_sparse) = levels(data$b)
	rownames(data_sparse) = levels(data$a)
	
	return(data_sparse)
}

calcCovFromFile <- function(filepaths, regions)
{
  f <- function(filepath, regions)
  {
    coverage <- import.bw(filepath)
    if (max(coverage$score, na.rm = TRUE) != 1) { coverage$score <- percent_rank(coverage$score) }
    res <- calcCov(regions, coverage, summary_fun = 'mean')
    return(res)
  }  
    
  commands <- list()
  for (filepath in filepaths)
  {
    commands[[filepath]] <- glue("f('{filepath}', regions)")
  }
    
  res <- gpatterns::gcluster.run2(command_list = commands, io_saturation = TRUE, jobs_title = 'importbigWigs')
  return(res)
}

feature.selection = function(x, y, binsize = 100, selection_threshold = 2, larger = TRUE)
{	
	#x = x[x != 0]
	#y = na.omit(y)
	x = as.numeric(x)
	y = as.numeric(y)
	ord = order(x)
	d = y[ord]
	d_trend = zoo::rollmedian(d, k = (binsize + 1), fill='extend')
	if (larger == TRUE) 
	{
		selected = d - d_trend > selection_threshold
	} else 
	{
		selected = d - d_trend < selection_threshold
	}
	selected = selected[order(ord)]
	return(selected)
}

tidyToMatrix = function(tidy_counts)
{
	m = tidy_counts %>% 
		select(unique_id, CB, n) %>% 
		spread(CB, n) %>%
		tibble.to.matrix
		
	m[is.na(m)] = 0
	return(m)
}

mafft = function(DNAStringSet, name = 'te', save_alignment = FALSE, threads = 1, type = 'fast')
{
	
	writeXStringSet(DNAStringSet, '/home/davidbr/tmp/mafft_input.fasta', format = 'fasta')
	# RUN MAFFT #
	if(save_alignment == TRUE)
	{
		outfile = paste0('/home/davidbr/davidbr/general/data/multiple_sequence_alignments/', name, '_mafft_aligned.fasta')
	} else {
		outfile = '/home/davidbr/tmp/mafft_output.fasta'
	}
	if (type == 'fast') 
	{
		print('fast alignment')
		cmd = glue('/home/davidbr/davidbr/tools/mafft/bin/mafft --ep 0 --op 1.1 --retree 1 --maxiterate 0 --nomemsave --thread {threads} /home/davidbr/tmp/mafft_input.fasta > {outfile}')
	}
	if (type == 'accurate')
	{
		print('accurate alignment')
		cmd = glue('/home/davidbr/davidbr/tools/mafft/bin/mafft --ep 0 --op 1.1 --retree 2 --maxiterate 1000 --nomemsave --thread {threads} /home/davidbr/tmp/mafft_input.fasta > {outfile}')
	}
	system(cmd)
	if(save_alignment == FALSE)
	{
		alignment = readBStringSet(outfile)
		return(DNAStringSet(toupper(alignment)))
	}
}

summariseMetaAndModuleUMIs = function(umi_mat, modules, metagrouping)
{
	metagrouping_filt = na.omit(metagrouping[colnames(umi_mat)])
	
	m_template = matrix(0, nrow = max(modules), ncol = max(metagrouping_filt, na.rm = TRUE), dimnames = list(as.character(1:max(modules)), as.character(1:max(metagrouping_filt, na.rm = TRUE))))
	for (mc_id in unique(metagrouping_filt))
	{
		for (module in unique(modules))
		{
			row_ids = names(modules[modules == module])
			barcodes = names(metagrouping_filt[metagrouping_filt == mc_id])
			m_template[module, mc_id] = mean(umi_mat[row_ids, barcodes], na.rm = TRUE)
			#print('done')
		}
		#print('done')
	}
	return(m_template)
}

addSCMat = function(mat, id)
{
	scmat = tgScMat(mat, stat_type = 'umi', cell_metadata = data.frame(0, row.names = colnames(mat), type = '10x', 
		batch_set_id = rep(NA, ncol(mat)), amp_batch_id = NA, seq_batch_id = NA, spike_count = 0))
	scdb_del_mat(id)
	scdb_add_mat(id = id, mat = scmat)	
}

sparse.matrix.to.tidy = function(sparse_matrix)
{
	tib = as_tibble(tidy(sparse_matrix)) %>% rename(repname_full = row, id = column, umi_counts = value)
	tib = tib %>% separate(repname_full, c('repname', 'bin_id', 'rep_size', 'repfamily', 'rep_id', 'chrom' ,'start', 'end', 'strand', 'copy_number', 'bin_size'), sep='\\.')
	return(tib)	
}

plotMSA=function(x, cluster = FALSE) 
{	
	if (cluster) {
		d = dist.dna(as.DNAbin(x), model = 'indel')
		hc = hclust(d, 'ward.D')
		ordering = hc$order
	} else {
		ordering = 1:length(x)
	}
	
	x=as.matrix(x)
	x=gsub("-", 0, x); 
	x=gsub("A", 1, x); 
	x=gsub("T", 2, x); 
	x=gsub("G", 3, x); 
	x=gsub("C", 4, x);
	x=apply(x, 2, as.numeric)
	image(y = 0:nrow(x), x = 0:ncol(x), z = t(x[ordering, ]), col=c("white", "#5DA731", "#D63317", "#DFD93C", "#2F4C9B"), yaxt='n', xaxt='n')
}

calcCov = function(regions, covs, norm = FALSE, summary_fun = sum)
{
	if (class(regions) != 'GRanges')
  {
    regions = makeGRangesFromDataFrame(regions)
  }
	regions$score = 0
	seqlevelsStyle(covs) = seqlevelsStyle(regions)
	overlaps = findOverlaps(regions, covs)
	signal = covs[to(overlaps)]$score
	regions[unique(from(overlaps))]$score = as.numeric(tapply(signal, from(overlaps), summary_fun, na.rm = TRUE))
	if(length(unique(from(overlaps))) != length(regions))
	{
		regions[-unique(from(overlaps))]$score = NA
	}
	if (norm)
	{
		region_widths = width(regions)
		regions$score = as.integer(round(regions$score / region_widths * median(region_widths)))
	}
	return(regions$score)
}

calcTELoad = function(tes_tidy, genes_tidy)
{
	te_load_df = full_join(tes_tidy %>% group_by(CB) %>% summarise(te_umis = sum(n)), genes_tidy %>% group_by(CB) %>% summarise(gene_umis = sum(n))) %>% 
		replace_na(list(te_umis = 0, gene_umis = 0)) %>%
		mutate(te_load = te_umis / (te_umis + gene_umis)) %>%
		arrange(CB)
		
	te_load = structure(te_load_df$te_load, names = te_load_df$CB)
	return(te_load)
}

calcCellSize = function(genes_tidy)
{
	df = genes_tidy %>% group_by(CB) %>% summarise(umis = sum(n)) %>%
		arrange(CB)
	res = structure(df$umis, names = df$CB)
	return(res)
}

calcTSSEnr = function(genes_tidy)
{
	cell_stats = genes_tidy %>% 
		group_by(CB, tss) %>% 
		summarise(umis = sum(n)) %>% 
		ungroup %>% 
		spread(tss, umis, sep = '_') %>%
		mutate(tss_perc = tss_TRUE / (tss_FALSE + tss_TRUE) * 100) %>%
		arrange(CB)
		
	res = structure(cell_stats$tss_perc, names = cell_stats$CB)
	return(res)
}

rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

norm.umis = function(umi_counts, factor_by_median = TRUE, norm_size = TRUE, k_scale_umi = 7)
{
	umis = umi_counts
	n_median = 1
	if(factor_by_median) {
		n_median = median(colSums(umis))
		# cat("..use median norm ", n_median)
	}
	if(norm_size) {
		# cat(", norm size ")
		umis = n_median * t(t(umis)/(1+colSums(umis)))
	}
	cell_fps = log2(1+k_scale_umi*umis)
	return(cell_fps)
}

image2 = function(m, cluster_rows = FALSE, cluster_cols = FALSE, nclust = 10)
{
	m = as.matrix(m)
	if (cluster_rows == TRUE) 
	{ 
		d = tgs_dist(m)
		hc = hclust(d, 'ward.D2')
		ct_rows = cutree(hc, nclust)
		#cluster clusters
		ordering = as.character(unlist(tapply(1:nrow(m), as.numeric(ct_rows), function(x) 
		{
			if(length(x) > 1)
			{
				d = tgs_dist(m[x,])
				hc = hclust(d, 'ward.D2')
				ordered_rownames = rownames(m[x,])[hc$order]
			} else { ordered_rownames = rownames(m)[x] }
			
			return(ordered_rownames)		
		})))
	} else { ordering = 1:nrow(m) }
	if (cluster_cols == TRUE) 
	{ 
		d = tgs_dist(t(m))
		hc = hclust(d, 'ward.D2')
		ct_cols = cutree(hc, nclust)
	} else { ct_cols = 1:ncol(m) }
	
	image.plot(t(m[ordering, order(ct_cols)]), col = c(colorRampPalette(c('darkgreen', 'orange', 'darkred'))(256)), yaxt='n', xaxt ='n')
	
	if (cluster_rows == TRUE) { 
		steps = cumsum(table(ct_rows))/max(cumsum(table(ct_rows)))
		axis(2, at = steps, labels=1:nclust) 
		segments(x0=-1, y0=steps, x1=2, y1=steps, lwd = 2)
	}
	if (cluster_cols == TRUE) { 
		steps = cumsum(table(ct_cols))/max(cumsum(table(ct_cols)))
		axis(1, at = steps, labels=1:nclust) 
		segments(x0=steps, y0=0, x1=steps, y1=1)
	}
	axis(3, at = seq(0, 1, 1/(ncol(m)-1)), labels=colnames(m) , las = 2)
	return(ordering)
}

ensembl.to.gene.symbol  = function(ensembl_ids)
{
	ensembl_ids = as_tibble(data.frame('Gene_stable_ID' = ensembl_ids))
	#ensembl = ensembl %>% rename(Gene_stable_ID = Gene_ID) %>% distinct
	#ensembl = ensembl %>% filter(!duplicated(ensembl$Gene_name)) # remove duplicated Genes
	#ensembl = filter(ensembl, Gene_ID %in% rownames(m))
	joined = left_join(ensembl_ids, ensembl, by = c('Gene_stable_ID'))
	return(pull(joined, 2))
}

annotateTE = function(repnames, column = 'repclass')
{
	repnames = as_tibble(repnames) %>% rename(repname = value)
	annos = left_join(repnames, te_anno_global, by = 'repname')
	return(pull(annos, column))
}

meta.aggr = function(sm, metas, aggr_func = sum)
{
	
	sm = sm[, metas$id]
	if (identical(colnames(sm), metas$id) == FALSE) { print('Attention: Barcodes do not match!') }
	meta_agr = t(apply(sm, 1, function(x) { tapply(x, metas$clust, aggr_func) }))
	return(meta_agr)
}

read.clust.fp = function(filepath)
{
	clustfp = as_tibble(fread(filepath, header=F, skip=3))
	colnames(clustfp) = c('Gene.name', 1:(ncol(clustfp) -1))
	return(clustfp)
}

read.xy.file = function(filepath) #metacell cluster file
{
	clusters = as_tibble(fread(filepath, header=F, drop=1))
	colnames(clusters) = c('id', 'x', 'y', 'clust', 'color', 'type', 'batch')
	return(clusters)
}

sep.te.rownames = function(m)
{
	ids = rownames(m)
	tib = as_tibble(data.frame('repname_full' = ids))
	tib = tib %>% separate(repname_full, c('repname', 'bin_id', 'rep_size', 'repfamily', 'rep_id', 'chrom' ,'start', 'end', 'strand', 'copy_number', 'bin_size'), sep='\\.')
	return(tib)
}

sep.full.te.name = function(tib)
{
	counts = tib %>% separate(1, c('repname', 'bin_id', 'repfamily', 'rep_id', 'chrom' ,'start', 'end', 'strand'), sep='\\.')
	return(counts)
}

genome.coverage = function(reads, track_dir = 'path_to_dir')
{
	seqlevelsStyle(reads) = 'UCSC'
	covs_plus = as(coverage(reads[strand(reads) == '+']), 'GRanges')
	covs_minus = as(coverage(reads[strand(reads) == '-']), 'GRanges')
	strand(covs_plus) = '+'
	strand(covs_minus) = '-'
	covs_plus = covs_plus[covs_plus$score > 0]
	covs_minus = covs_minus[covs_minus$score > 0]
	covs = c(covs_plus, covs_minus)
	return(covs)
}

strsplit2 = function(string, sep = '\\.', column = 1)
{
	string = as_tibble(string)
	string = separate(string, value, as.character(1:column), sep = sep)
	if (length(column) > 1)
	{
		res = select(string, column)
	} else {
		res = pull(string, column)
	}
	return(res)
}

bin.TEs = function(TEs, binsize) {
	bins=slidingWindows(TEs, width=binsize, step=binsize)
	names(bins) = 1:length(bins)
	bins=unlist(bins)
	indeces = as.numeric(names(bins))
	bins$name = TEs[indeces]$name
	bins$unique_name = make.unique(TEs[indeces]$unique_name)
	matches=grep('.[0-9]$', bins$unique_name)
	bins[-matches]$unique_name = paste0(bins[-matches]$unique_name, '.0')
	return(bins)
}

chunk.values = function(data_range, steps)
{
	chunks=lapply(split(data_range, ceiling(seq_along(data_range)/steps)), range)
	return(chunks)
}

convertStrand = function(df) {
	if (sum(grepl('1$', df$strand)) == 0)
	{
		df$strand = as.character(df$strand)
		df[df[,'strand']=='+','strand']=1
		df[df[,'strand']=='-','strand']=-1
		df[df[,'strand']=='*','strand']=0
		df$strand = as.numeric(df$strand)
		if(sum(colnames(df) == 'seqnames') > 0) 
		{
			df = df %>% rename(chrom = seqnames)
		}
	} else {
		df[df[,'strand']==1,'strand']='+'
		df[df[,'strand']==-1,'strand']='-'
		df[df[,'strand']==0,'strand']='*'
	}
	return(df)
}

reads.to.ctss = function(reads) {
	reads=granges(reads, use.mcols=TRUE)
	end(reads[strand(reads)=="+"])=start(reads[strand(reads)=="+"])
	start(reads[strand(reads)=="-"])=end(reads[strand(reads)=="-"])
	return(reads)
}

tibble.to.matrix = function(tib) # takes as input a tibble with first column as rownames
{
	m = as.matrix(tib[,-1])
	rownames(m) = pull(tib, 1)
	return(m)
}

read.ctss.from.bam = function(bam_path, which_regions = NULL, uniquely = FALSE) 
{
	print('reading bam file')
	bam_index = gsub('.bam$', '.bai', bam_path)
	if (!is.null(which_regions)) 
	{
		seqlevelsStyle(which_regions) = 'NCBI'
		reads=readGAlignments(bam_path, bam_index,
		use.names=TRUE, param=ScanBamParam(which = which_regions, tag = c("NH", "CB"),
		flag=scanBamFlag(isPaired=TRUE, isProperPair=TRUE, isFirstMateRead=TRUE)))
		if(length(reads) ==0) 
		{ 
		reads = readGAlignments(bam_path, bam_index,
			use.names=TRUE, param=ScanBamParam(which = which_regions, tag = c("NH", "CB"),
			flag=scanBamFlag(isPaired=FALSE)))
		}
	} else {
		reads=readGAlignments(bam_path, 
			use.names=TRUE, param=ScanBamParam(tag = c("NH", "CB"),
			flag=scanBamFlag(isPaired=TRUE, isProperPair=TRUE, isFirstMateRead=TRUE)))
		if(length(reads) ==0) 
		{ 
		reads = readGAlignments(bam_path, 
			use.names=TRUE, param=ScanBamParam(tag = c("NH", "CB"),
			flag=scanBamFlag(isPaired=FALSE)))
		}
	}
	print('finished')
	reads.gr = granges(reads, use.mcols=TRUE)
	if(uniquely == TRUE)
	{
		unique_reads = reads.gr[reads.gr$NH == 1]
		reads.gr = unique_reads	
	}
	tss = reads.to.ctss(reads.gr)
	tss = sort(tss)
	return(tss)
}

get.tss = function(genes)
{
	transcripts = genes[genes$type == 'transcript']
	end(transcripts[strand(transcripts) == '+']) = start(transcripts[strand(transcripts) == '+'])
	#start(transcripts[strand(transcripts) == '+']) = start(transcripts[strand(transcripts) == '+']) - extend
	
	start(transcripts[strand(transcripts) == '-']) = end(transcripts[strand(transcripts) == '-'])
	#end(transcripts[strand(transcripts) == '-']) = end(transcripts[strand(transcripts) == '-']) + extend
	tss = unique(transcripts)
	return(tss)
}

tes.2d = function (TE_counts, clusters, TE_grep, combine) {
	color_scale = colorRampPalette(c('darkgreen', 'orange', 'darkred'))(10)
	specific_TE_matrix = TE_counts[grep(TE_grep, rownames(TE_counts), fixed=TRUE), intersect(clusters$id, colnames(TE_counts))]
	specific_TE_matrix = specific_TE_matrix[rowSums(specific_TE_matrix)>1, ]
	if (combine == TRUE) { 
		specific_TE_matrix = colSums(as.matrix(specific_TE_matrix)) 
		specific_TE_counts_log = log2(specific_TE_matrix+1)
		color_scale_sel = color_scale[as.numeric(cut(c(as.numeric(specific_TE_counts_log), max(as.numeric(specific_TE_counts_log))), breaks = 10))]
		plot(x = clusters$x, y = clusters$y, col = color_scale_sel, pch=20, cex = 0.5, main=TE_grep, yaxt='n', xaxt='n')
	} else {
		specific_TE_counts_log = log2(specific_TE_matrix+1)
		par(mfrow = c(round(sqrt(nrow(specific_TE_matrix))), round(sqrt(nrow(specific_TE_matrix))+1)), mai = c(0,0,0,0))
		for (i in 1:nrow(specific_TE_matrix))
		{
			color_scale_sel = color_scale[as.numeric(cut(c(as.numeric(specific_TE_counts_log[i,]), max(as.numeric(specific_TE_counts_log[i,]))), breaks = 10))]
			plot(x = clusters$x, y = clusters$y, col = color_scale_sel, pch=20, cex = 0.5, yaxt='n', xaxt='n')
		}
	}
}

genes.2d = function (gene_counts, metas, genes) {
	gene_counts = gene_counts[, metas$id]
	color_scale = colorRampPalette(c('darkgreen', 'orange', 'darkred'))(10)
	counts = gene_counts[intersect(genes, rownames(gene_counts)), ]
	if(length(intersect(genes, rownames(gene_counts))) > 1) 
	{
		counts = colSums(as.matrix(counts)) 
	}
	counts = log2(counts + 1)
	color_scale_sel = color_scale[as.numeric(cut(c(counts, max(counts)), breaks = 10))]
	if (length(intersect(genes, rownames(gene_counts))) > 10)
	{
		plot(x = metas$x, y = metas$y, col = color_scale_sel, pch=20, cex = 0.5, yaxt='n', xaxt='n')
	} else {
		plot(x = metas$x, y = metas$y, col = color_scale_sel, pch=20, cex = 0.5, yaxt='n', xaxt='n', main = genes)
	}
}

read.10x.matrix = function(filepath, tidy = TRUE, min_umi_count = 0) {
	matrix_10x = as.matrix(tgutil::fread_mm(filepath))
	genes = fread(paste0(dirname(filepath), '/genes.tsv'), header=F)
	barcodes = fread(paste0(dirname(filepath), '/barcodes.tsv'), header=F)
	rownames(matrix_10x) = pull(genes, 2)
	colnames(matrix_10x) = pull(barcodes, 1)
	matrix_10x = matrix_10x[rowSums(matrix_10x) >= min_umi_count, ]
	matrix_10x = as(matrix_10x, 'dgCMatrix')
	if (tidy == TRUE) { return(as_tibble(tidy(matrix_10x))) } else {	return(matrix_10x) }
}

vectToTib = function(vect)
{
	res = data.frame('id' = names(vect), value = vect) %>%
		as_tibble
	return(res)
}

writeSparseMatrix = function(matrix_object, filedir) {
	matrix_10x = as(matrix_object, 'dgCMatrix')
	genes = rownames(matrix_10x)
	barcodes = colnames(matrix_10x)
	writeMM(matrix_10x, paste0(filedir, 'matrix.mtx'))
	write(barcodes, paste0(filedir, 'barcodes.tsv'))
	write(genes, paste0(filedir, 'features.tsv')) 
	#print('Attention: Genes file needs 2 columns: ENSEMBL and GENE SYMBOL')
}