source('/net/mraid14/export/data/users/davidbr/general/scripts/utilities.R')

#########################
# FUNCTIONS
#########################
plotCorMatrix = function(mat_cor, ct, plot_diag = TRUE)
{
	#diag(mat_cor) = NA
	n = nrow(mat_cor)
	
	#png('test.png', w=160+n*15, h=120+n*15)
	mat_cor_ordered = mat_cor[names(ct), names(ct)]
	#mat_cor_ordered[upper.tri(mat_cor_ordered)] = NA
	
	mat_cor_ordered[mat_cor_ordered > 0.15] = 0.16
	mat_cor_ordered[mat_cor_ordered < -0.15] = - 0.16
	
	if(!plot_diag)
	{
		diag(mat_cor_ordered) = NA
		mat_cor_ordered[upper.tri(mat_cor_ordered)] = NA
	}
	
	cols = c('black', colorRampPalette(c("darkblue", 'blue', 'white', 'orange', "red"))(1000), 'darkred')
	
	image(y = 0:nrow(mat_cor_ordered), x = 0:ncol(mat_cor_ordered), t(mat_cor_ordered), 
		col = cols, 
		breaks = c(-0.16, c(seq(-0.15, 0.15, l = 1001)), 0.16), yaxt = 'n', xaxt = 'n', xlab = '', ylab = '')
	
	steps = c(0, cumsum(table(ct)))
	
	abline(h = steps, v = steps, lwd = 0.5)
}

clusterCorMatrix = function(mat_cor, ct_height = NULL, nclust = NULL)
{
	d =tgs_dist(mat_cor)
	hc = hclust(d, 'ward.D')
	#hc = hclust(as.dist(1 - mat_cor), "ward.D2")

	if (!is.null(ct_height))
	{	
		ct = cutree(hc, h = ct_height)
	} else {
		ct = cutree(hc, nclust)
	}
	#names(ct) = strsplit2(names(ct), '\\|', 1)
	ct_sorted = sort(ct)
	return(ct_sorted)
}

getGeneticDist = function(tes_full_length, mat_cor)
{
	tes_sel_gr = tes_full_length %>% 
		filter(row_id %in% strsplit2(rownames(mat_cor), '\\|', 3)) %>% 
		select(chrom, start, end, strand, row_id) %>% 
		makeGRangesFromDataFrame(keep.extra.columns = TRUE)
	seqlevelsStyle(tes_sel_gr) = 'UCSC'
	seqs = getSeq(Hsapiens, tes_sel_gr)
	names(seqs) = tes_sel_gr$row_id
	msa = mafft(seqs, save = FALSE, threads = 18)
	
	d = as.matrix(dist.dna(as.DNAbin(msa), model = 'indelblock'))
	#d = tgs_cor(as.matrix(d), spearman = TRUE)
	return(d)
}

plotGeneticDistHeatmap = function(dist_m, ct)
{

	zlim_max = max(max(dist_m), 350)
	row_ids = strsplit2(names(ct), '\\|', 3)
	d_ordered = as.matrix(dist_m)[row_ids, row_ids]
	d_ordered[lower.tri(d_ordered)] = NA
	# pheatmap(d_ordered, cluster_cols = F, cluster_rows = F,
		# show_rownames = FALSE, show_colnames = FALSE, col = colorRampPalette(c("darkgreen", 'orange', 'darkred'))(1000),
		# annotation_col = annotation)
	diag(d_ordered) = NA
	#image(t(d_ordered), col = rev(colorRampPalette(c(brewer.pal(11, 'Spectral')))(1000)), 
	image(y = 0:nrow(d_ordered), x = 0:ncol(d_ordered), z = t(d_ordered), col = colorRampPalette(rev(c('white', 'orange', 'darkred', 'black')))(1000),
		zlim = c(0, zlim_max), yaxt = 'n', xaxt = 'n', ylab = '', xlab = '', add = TRUE)
	
	steps = c(0, cumsum(table(ct)))
	
	abline(h = steps, v = steps, lwd = 0.5)	
}

plotClustTEFraction = function(ct)
{
	ct_mod = ct
	names(ct_mod) = strsplit2(names(ct_mod), '\\|', 1)
	#names(ct_mod) = annotateTE(names(ct_mod), 'repclass')

	cor_clusts = tibble(repfamily = names(ct_mod), clust = as.numeric(ct_mod))
	clust_counts = count(cor_clusts, repfamily, clust)
	m = clust_counts %>% 
		spread(clust, n, fill = 0) %>% 
		tibble.to.matrix
	m_perc = t(t(m) / colSums(m))
	m_perc_expanded = m_perc[,rep(1:ncol(m), table(ct_mod))]

	barplot(m_perc_expanded[, colnames(m_perc_expanded)], border = NA, horiz = TRUE, xlab = '', axisnames = FALSE, axes = FALSE)
}

annotateTECustom = function(ids)
{
	TE_loci = ids
	TE_families = strsplit2(TE_loci, '\\|', 1)
	TE_families_raw = TE_families
	
	TE_families[annotateTE(TE_families, 'repclass') == 'DNA'] = 'DNA'
	TE_families[annotateTE(TE_families, 'repclass') == 'SINE'] = 'SINE'
	TE_families[annotateTE(TE_families, 'repclass') == 'LINE'] = 'LINE'
	TE_families[grepl('ERVL', annotateTE(TE_families, 'repfamily'))] = 'ERVL-MaLR/ERVL'
	TE_families[grepl('LTR12', TE_families)] = 'ERV9'
	TE_families[grepl('ERV1|ERVK', annotateTE(TE_families, 'repfamily'))] = 'ERV1/ERVK'
	
	TE_annotation = tibble(repname = TE_families_raw, custom = TE_families, unique_id = TE_loci) 

	return(TE_annotation)
}

plotTEClassAnno = function(TE_annotation, ct)
{
	cols_filt = cols[TE_annotation$custom %>% unique]
	
	num_levels = TE_annotation$custom %>% as.factor %>% as.numeric
	names(num_levels) = TE_annotation$custom
	num_levels_m = as.matrix(num_levels)
	
	image(t(num_levels_m), col = cols_filt[order(names(cols_filt))], yaxt = 'n', xaxt = 'n')
}

clusterMetas = function(mat, mat_all, mc, ct)
{
	mcs = as.factor(mc@mc)
	
	d = tgs_dist(t(log2(mc@mc_fp[names(ct),])))
	hc = hclust(d, 'ward.D')
	
	mcs_reordered = factor(mcs, levels(mcs)[hc$order])
	mcs_reordered_filt = mcs_reordered[colnames(mat)] %>% sort
	
	meta_clustering = as.data.frame.factor(mcs_reordered_filt) %>% 
		as_tibble %>% 
		rownames_to_column %>% 
		rename(CB = rowname, mc = mcs_reordered_filt) %>%
		mutate(mc = as.numeric(as.character(mc))) %>%
		left_join(. , as.data.frame(colSums(mat_all)) %>% 
						as.tibble %>% 
						rownames_to_column %>% 
						rename(CB = rowname, size = `colSums(mat_all)`), by = 'CB')
	
	meta_clustering = meta_clustering %>% group_by(mc) %>% mutate(CB = CB[order(size)], size = size[order(size)]) %>% ungroup
	
	return(meta_clustering)
}

plotVarMean = function(gstat)
{
	total = gstat$tot
	ds_var = gstat$ds_var
	ds_mean = gstat$ds_mean
	
	plot(ds_mean, ds_var/ds_mean, log = 'xy', pch = 20)
}

plotTotalTEs = function(mat_all, meta_clustering)
{
	total_te_umis = colSums(mat_all[, meta_clustering$CB])
	
	cols = tol21rainbow[meta_clustering$mc]
	cols = 'grey'
	
	plot(total_te_umis, type = 'h', col = cols, xaxt = 'n', xlab = '', ylab = '', ylim = c(0,5000))
	abline(h = 5000)
}

plotMetaCell = function(meta_clustering, mat_ds, ct, smooth_bin_size = 5)
{
	meta_clustering_filt = meta_clustering %>% filter(CB %in% colnames(mat_ds))
	mat_ds_filt = mat_ds[names(ct), meta_clustering_filt$CB]
	mat_trans = log2(1 + 7 * mat_ds_filt)
    lus = apply(mat_trans - apply(mat_trans, 1, median), 2, function(x) pmax(x,0))
    lus_smoo = t(apply(lus, 1, function(x) zoo::rollmean(x, smooth_bin_size, fill= 'extend')))
	
	#lus_smoo = as.matrix(mat_trans)

    image(y = 0:nrow(lus_smoo), x = 0:ncol(lus_smoo), z = t(lus_smoo), 
		col = colorRampPalette(c('white', 'darkblue', 'orange', 'darkred'))(256),
		xaxt='n', yaxt='n', xlab = '', ylab = '')
		
	steps_v = c(0, cumsum(table(meta_clustering$mc)[meta_clustering$mc %>% unique]))
	steps_h = c(0, cumsum(table(ct)))
	
	abline(v = steps_v, lwd = 0.5, h = steps_h)
}

plotTEFractionPerClass = function(tes_tidy, meta_clustering, TE_annotation)
{
	tes_tidy_anno_filt = inner_join(tes_tidy, TE_annotation %>% select(custom, unique_id), by = 'unique_id')
	
	cell_stats = tes_tidy_anno_filt %>% 
		group_by(CB, custom) %>% 
			summarise(umis = sum(n)) %>%
		group_by(CB) %>% 
			mutate(te_frac = umis / sum(umis)) %>% 
			ungroup
			
	m = cell_stats %>%
		select(CB, custom, te_frac) %>%
		spread(custom, te_frac, fill = 0) %>%
		tibble.to.matrix %>%
		t
		
	m_filt = m[, meta_clustering$CB]

	barplot(m_filt[c('DNA', 'LINE', 'ERV9', 'ERVL-MaLR/ERVL', 'ERV1/ERVK', 'SINE'), ], border = cols[c('DNA', 'LINE', 'ERV9', 'ERVL-MaLR/ERVL', 'ERV1/ERVK', 'SINE')], axes = FALSE, axisnames = FALSE)
}

plotFigure = function(mat_all, mat_cor, dist_m, meta_clustering, tes_tidy, ct, output_dir = '~/davidbr/proj/epitherapy/manuscript/figures/main/te_meta/')
{
	old_wd = getwd()
	setwd(output_dir)
	
	ct_1 = ct[ct %in% 1:2]
	ct_2 = ct[!ct %in% 1:2]
	
	# PLOT FULL COR MAT
	png('cor_matrix_full.png', w = 5000, h = 5000, res = 1200)
	par(mar = c(rep(0.25, 4)))
	plotCorMatrix(mat_cor, ct, plot_diag = TRUE)
	dev.off()
	
	# PLOT HIGHLIGHT COR MAT + GENETIC DIST
	png('cor_matrix_full_highlight1.png', w = 5000, h = 5000, res = 1200)
	par(mar = c(rep(0.25, 4)), yaxs="i", xaxs = 'i')
	layout(matrix(c(1,2), ncol = 2), heights = c(1,1), widths = c(2, 0.2))
	
	TE_annotation = annotateTECustom(names(ct_1))
	plotCorMatrix(mat_cor, ct_1, plot_diag = FALSE)
	plotGeneticDistHeatmap(dist_m, ct_1)
	plotTEClassAnno(TE_annotation, ct_1)
	dev.off()
	
	# PLOT HIGHLIGHT COR MAT2 + GENETIC DIST2
	png('cor_matrix_full_highlight2.png', w = 5000, h = 5000, res = 1200)
	par(mar = c(rep(0.25, 4)), yaxs="i", xaxs = 'i')
	layout(matrix(c(1,2), ncol = 2), heights = c(1,1), widths = c(2, 0.2))
	
	TE_annotation = annotateTECustom(names(ct_2))
	plotCorMatrix(mat_cor, ct_2, plot_diag = FALSE)
	plotGeneticDistHeatmap(dist_m, ct_2)
	plotTEClassAnno(TE_annotation, ct_2)
	dev.off()
	
	#PLOT META HEATMAP UMI COUNTS
	
	png('umi_matrix_heatmap_1.png', w = 10000, h = 5000, res = 1200)
	par(mar = c(rep(0.25, 4)), yaxs="i", xaxs = 'i', lwd = 0.5)
	plotMetaCell(meta_clustering, mat_ds, ct_1)
	dev.off()
	
	png('umi_matrix_heatmap_2.png', w = 10000, h = 5000, res = 1200)
	par(mar = c(rep(0.25, 4)), yaxs="i", xaxs = 'i', lwd = 0.5)
	plotMetaCell(meta_clustering, mat_ds, ct_2)
	dev.off()
	
	png('total_TE_umi_counts.png', w = 10000, h = 1000, res = 1200)
	par(mar = c(rep(0.25, 4)), yaxs="i", xaxs = 'i', lwd = 0.5)
	plotTotalTEs(mat_all, meta_clustering)
	dev.off()
	
	TE_annotation = annotateTECustom(names(ct))
	
	png('TE_class_fraction.png', w = 10000, h = 1000, res = 1200)
	par(mar = c(rep(0.25, 4)), yaxs="i", xaxs = 'i', lwd = 0.5)
	plotTEFractionPerClass(tes_tidy, meta_clustering, TE_annotation)
	dev.off()
	
	if(cell_line == 'h1299')
	{
	p = data.frame(mod_9 = colMeans(mat_ds[names(ct[ct == 9]),]), 
	          mod_7 = colMeans(mat_ds[names(ct[ct == 7]),])) %>% 
		as_tibble %>%  
		ggplot(aes(mod_7, mod_9)) + 
			geom_point_rast(size = 2, alpha = 0.1)
	ggsave('antagonistic_te_mods_scatterplots.pdf', p, device = 'pdf', width = 4, height = 4)
	}
	
	setwd(old_wd)
}
		   

###############
# LOAD
###############
te_full_length = readRDS('~/davidbr/general/data/repeatmasker/hg38/full_length_te_coords.RDS')

cell_line = 'h1299'
tes_tidy = readRDS(paste0("~/davidbr/proj/epitherapy/data/", cell_line, "/10x/dacsb/counts/te_tss_bin_counts_tidy.RDS"))
mc = readRDS(paste0('~/davidbr/proj/epitherapy/data/', cell_line, '/10x/dacsb/metacell/tes/output/metacell_object.RDS'))
mat_all = readRDS(paste0('~/davidbr/proj/epitherapy/data/', cell_line, '/10x/dacsb/metacell/tes/output/umi_matrix_all.RDS'))
mat_ds = readRDS(paste0('~/davidbr/proj/epitherapy/data/', cell_line, '/10x/dacsb/metacell/tes/output/umi_matrix_marker_ds.RDS'))
mat = readRDS(paste0('~/davidbr/proj/epitherapy/data/', cell_line, '/10x/dacsb/metacell/tes/output/umi_matrix_marker_raw.RDS'))
###############

# COR MATRIX
mat_cor = tgs_cor(t(log2(1 + 7 * as.matrix(mat_ds))), spearman = FALSE)
dist_m = getGeneticDist(te_full_length, mat_cor)
ct = clusterCorMatrix(mat_cor, nclust = 10)#ct_height = 4)
saveRDS(ct, paste0('~/davidbr/proj/epitherapy/workspaces/', cell_line, '_te_modules.RDS'))

ct_1 = ct[ct %in% 1:2]
ct_2 = ct[!ct %in% 1:2]
TE_annotation = annotateTECustom(names(ct))
#clust_stats = clusterEnrichment(ct)
meta_clustering = clusterMetas(mat_ds, mat_all, mc, ct_2)

tol21rainbow = c("#984EA3", "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 	"#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", 'darkred', 'darkorange2', 'orange', 'gold','pink')
names(tol21rainbow) = 1:length(tol21rainbow)

cols = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'ERV9' = "#4DAF4A", 'ERVL-MaLR/ERVL' = "#004529", 'ERV1/ERVK' = "#FFFFE5", 'SINE' = "#984EA3")

plotFigure(mat_all, mat_cor, dist_m, meta_clustering, tes_tidy, ct, output_dir = paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/te_meta/', cell_line, '/'))




