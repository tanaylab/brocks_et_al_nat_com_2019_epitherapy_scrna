#########################
# FUNCTIONS FOR FIG 4
#########################
mcell_mat_rpt_cor_anchors_custom = function(mat_id, gene_anchors, gene_anti = c(), downsample_n = NA, ignore_cells = NULL)
{
        mat = scdb_mat(mat_id)
        if(is.null(mat)) {
                stop("missing mat ", mat_id)
        }
	
        if(!is.null(ignore_cells)) {
                umis = mat@mat[, !colnames(mat@mat) %in% ignore_cells]
        } else {
                umis = mat@mat
        }

        if(is.na(downsample_n)) {
                downsample_n = quantile(colSums(umis), 0.05)
        }

        mat_ds = scm_downsamp(umis, downsample_n)
        mat_ds = mat_ds[rowSums(mat_ds)>10,]
		mat_ds = as.matrix(log2(1 + 7 * mat_ds))
        csize = colSums(mat@mat[,colnames(mat_ds)])
        gcors = data.frame(sz_cor = apply(mat_ds, 1, cor, csize))
        for(g in gene_anchors) {
                gcors[,g] = apply(mat_ds, 1, cor, mat_ds[g,])
        }
        for(g in gene_anti) {
                gcors[,g] = apply(mat_ds, 1, cor, mat_ds[g,])
        }
        N1 = length(gene_anchors) + 1
        N2 = N1 + length(gene_anti)

        if(length(gene_anchors) == 1) {
                gcors$max = gcors[,2]
        } else {
                gcors$max = apply(gcors[,2:N1], 1, max)
        }
        if(length(gene_anti) == 0) {
                gcors$neg_max = 0
        } else if(length(gene_anti) == 1) {
                gcors$neg_max = gcors[,N2]
        } else {
                gcors$neg_max = apply(gcors[,(N1+1):N2], 1, max)
        }
		res = gcors %>% as.tibble %>% rownames_to_column %>% rename(final_id = rowname)
		
		return(res)
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

filterGeneMatrix = function(genes_m, genes_tidy, exclude_mt = TRUE)
{
	if(exclude_mt)
	{
	genes_whitelist = genes_tidy %>% 
		filter(grepl('protein_coding', gene_type)) %>% 
		filter(n_tss == 0 | tss) %>%
		filter(!grepl('RP11', gene_name)) %>%
		filter(seqnames != 'MT') %>%
		#filter(!grepl('_', gene_name, fixed = TRUE)) %>%
		pull(final_id) %>% 
		unique
	} else {
	genes_whitelist = genes_tidy %>% 
		filter(grepl('protein_coding', gene_type)) %>% 
		filter(n_tss == 0 | tss) %>%
		filter(!grepl('RP11', gene_name)) %>%
		#filter(!grepl('_', gene_name, fixed = TRUE)) %>%
		pull(final_id) %>% 
		unique
	}	
	genes_m_filt = genes_m[genes_whitelist, ]
	
	return(genes_m_filt)
}

diff_expr = function(mcs, mat_ds, mcs1, mcs2, reg=5, min_max_umi=50, nms1=NULL, nms2=NULL)
{
        if (is.null(nms1)) {
                nms1 = names(mcs)[mcs %in% mcs1]
        }
        if (is.null(nms2)) {
                nms2 = names(mcs)[mcs %in% mcs2]
        }
        nms1 = intersect(colnames(mat_ds), nms1)
        nms2 = intersect(colnames(mat_ds), nms2)
        message(sprintf("comparing %d vs %d cells", length(nms1), length(nms2)))

        df = data.frame(row.names=rownames(mat_ds), gene=rownames(mat_ds), tot1=rowSums(mat_ds[, nms1]), tot2=rowSums(mat_ds[, nms2]))
        norm_by = min(sum(df$tot1), sum(df$tot2))
        df$tot1 = df$tot1 / sum(df$tot1) * norm_by
        df$tot2 = df$tot2 / sum(df$tot2) * norm_by
        df = df[pmax(df$tot1, df$tot2) >= min_max_umi, ]
        df$enr = log2( (df$tot1 + reg) / (df$tot2 + reg))
        df = df[order(df$enr, decreasing=T), ] %>% as_tibble
		
		return(df)
}

createMeta = function(sparse_matrix, id, mcs)
{
	scmat = tgScMat(sparse_matrix, stat_type = 'umi', cell_metadata = data.frame(0, row.names = colnames(sparse_matrix), type = '10x', 
					batch_set_id = rep(NA, ncol(sparse_matrix)), amp_batch_id = NA, seq_batch_id = NA, spike_count = 0))
	scdb_add_mat(id = id, mat = scmat)				
					
	outliers = setdiff(colnames(sparse_matrix), names(mcs))
	mcell_new_mc(mc_id = id, mc = mcs, outliers = outliers, scmat = scmat)

	# mcell_gset_from_mc_markers(gset_id = 'markers', mc_id = id)
	# mcell_mc_plot_marks(mc_id = id, gset_id = 'markers', mat_id = id, plot_cells = TRUE)
}

plotUMIsOn2D = function(mc2d, umi_vector)
{
	# par(mfrow = c(3,3), mar = rep(0.5,4))
	# for(i in unique(cc_genes))
	# {
	#umi_vector = colSums(genes_m_filt_ds[names(cc_genes[cc_genes == i]), ])
	umi_vector_norm = pmax(umi_vector - median(umi_vector), 0)
	umi_vector_norm = umi_vector / mean(umi_vector)
	
	shared_barcodes = intersect(names(mc2d@sc_x), names(umi_vector))
	
	x = mc2d@sc_x[shared_barcodes]
	y = mc2d@sc_y[shared_barcodes]
	
	umi_vector_norm = umi_vector_norm[shared_barcodes]
	
	#cols = colorRampPalette(c("#FAFAFA", "orange","red","purple","black"))(1000)[as.numeric(cut(umi_vector_norm, breaks=1000))]
	cols = c('lightgrey', colorRampPalette(c("#FAFAFA", "orange","red","purple","black"))(1000))[as.numeric(cut(umi_vector_norm, 
		breaks = c(-Inf, seq(1, max(umi_vector_norm), l = 1000))))]
	
	plot(x, y, pch = 20, cex = 1, xaxt = 'n', xlab = '', yaxt = 'n', ylab = '', col = alpha(cols, alpha = 0.25))
	# plot(x, y, pch = 20, cex = 1, xaxt = 'n', xlab = '', yaxt = 'n', ylab = '', col = alpha(cols, alpha = 0.5), main = i)
	# }
}

calcUMIsPerMod = function(mat_ds, mods, summary_stat = sum)
{
	if(length(setdiff(names(mods), rownames(mat_ds))) > 0)
	{
		overlapping_features = intersect(names(mods), rownames(mat_ds))
		print('attention, mat and modules are not compatible')
		print(paste0('working on: ', length(overlapping_features)))
		mods = mods[overlapping_features]
		mat_ds = mat_ds[overlapping_features,]
	}
	
	
	mat_ds_sorted = mat_ds[names(mods),]
	
	res = apply(mat_ds_sorted, 2, function(cell_counts) 
	{
		res = tapply(cell_counts, as.numeric(mods), summary_stat)
		return(res)
	})
	return(res)
}

plotMarkerGeneHeatmap = function(genes_tidy, te_full_length, nclust = 8)
{
	mcell_gset_from_mc_markers(gset_id = 'markers', mc_id = 'genes_te_meta')
	mc_genic = scdb_mc('genes_te_meta')
	markers = names(scdb_gset('markers')@gene_set)
	
	mat_to_plot = log2(mc_genic@mc_fp[markers, ])
	mat_to_plot[mat_to_plot > 1.5] = 1.6
	mat_to_plot[mat_to_plot < -1.5] = -1.6
	
	te_overlaps = geneTSSOverlapTE(genes_tidy, te_full_length, markers)
	te_overlap_m = matrix(0, nrow = nrow(mat_to_plot), ncol = 1, dimnames = list(rownames(mat_to_plot), 'TE'))
	te_overlap_m[rownames(te_overlaps),] = 1
	
	chromosomes = structure(genes_tidy %>% filter(final_id %in% markers) %>% select(final_id, seqnames) %>% distinct %>% pull(seqnames), 
	                        names = genes_tidy %>% filter(final_id %in% markers) %>% select(final_id, seqnames) %>% distinct %>% pull(final_id))
	
	########################
	# PLOTTING
	########################
	annotation_rows = data.frame(row.names = rownames(mat_to_plot), 
								 'TE' = te_overlap_m[rownames(mat_to_plot),], 
								 'chrom' = chromosomes[rownames(mat_to_plot)])
	
	p = pheatmap(mat_to_plot[rev(rownames(mat_to_plot)),],
			col = c('darkblue', colorRampPalette(get_param("mcp_heatmap_fp_shades", package = 'metacell'))(256), 'black'), 
			breaks = c(-1.6, seq(-1.5, 1.5, l = 257), 1.6),
			clustering_method = 'ward.D',
			cutree_rows = nclust,
			border_color = NA,
			treeheight_row = 0,
			treeheight_col = 0,
			annotation_row = annotation_rows)
	
	print(p)
			
	return(p$tree_col$labels[p$tree_col$order])
}

plotDiffExpr = function(mcs1 = 1, mcs2 = 2:3, mcs, mat_ds, min_max_umi = 25, reg = 5)
{
	diff_genes = diff_expr(mcs, mat_ds, mcs1, mcs2, min_max_umi, reg) %>% 
		mutate(color = ifelse(enr > 1, 'darkred', 'black'), 
		       color = ifelse(enr < -1, 'darkblue', color),
			   rank1 = rank(enr),
			   rank2 = rank(-enr)) %>%
		as_tibble
		
	p_pdf = diff_genes %>% 
		ggplot(aes(x = tot2 + 1, y = tot1 + 1, col = color, label = ifelse(rank1 <= 5 | rank2 <= 5, gene, ''))) + 
		geom_point(size = 0.25) + 
		scale_x_log10(limits = c(1, 5e4), breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + 
		scale_y_log10(limits = c(1, 5e4), breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) +
		scale_color_manual(values = c('black', 'blue', 'red')) + 
		geom_text_repel() +
		theme(legend.position="none") + #, axis.title.x = element_blank(), axis.title.y = element_blank()) + 
		labs(x = paste(sort(mcs2), collapse = '|'), y =  paste(sort(mcs1), collapse = '|'))
		
	return(list(p_pdf, diff_genes))
	#return(diff_genes)
}

selectIFNsUp = function(dmso_tidy, dacsb_tidy, min_max_umi = 25, reg = 5)
{
	ifns = as_tibble(fread('~/davidbr/proj/epitherapy/data/interferome/interferome_all.txt'))

	ifns_stat = ifns %>% 
		filter(Inteferome_Type == 'I') %>%
			group_by(Gene_Name) %>% 
				summarise(nup = sum(Fold_Change >= 2), 
				          ndown = sum(Fold_Change <= -2), 
						  mean_fc = mean(Fold_Change), 
						  nstudies = length(unique(Dataset_ID)), 
						  consistent_up = nup > 0 & ndown == 0) %>% 
				ungroup

	sel_ifns = filter(ifns_stat, consistent_up & mean_fc >= 2 & nstudies > 1) %>%
		pull(Gene_Name) %>% 
		unique %>% 
		sort
		
	dmso = dmso_tidy %>% group_by(gene_name) %>% summarise(tot1 = sum(n))
	dacsb = dacsb_tidy %>% group_by(gene_name) %>% summarise(tot2 = sum(n))
		
	df = full_join(dmso, dacsb, by = 'gene_name') %>% replace_na(list(tot1 = 0, tot2 = 0))
	norm_by = min(sum(df$tot1), sum(df$tot2))
    df$tot1 = df$tot1 / sum(df$tot1) * norm_by
    df$tot2 = df$tot2 / sum(df$tot2) * norm_by
    df = df[pmax(df$tot1, df$tot2) >= min_max_umi, ]
    df$enr = log2( (df$tot1 + reg) / (df$tot2 + reg))
    df = df[order(df$enr, decreasing=T), ] %>% 
		mutate(ifn_gene = gene_name %in% sel_ifns) %>%
		as_tibble
	
	p = df %>% 
		ggplot(aes(tot1, tot2, color = enr < -2 & ifn_gene, alpha = enr < -2 & ifn_gene)) + 
			geom_point_rast(size = 1) + 
			scale_x_log10() + 
			scale_y_log10() +
			scale_alpha_discrete(range = c(0.05, 1)) +
			scale_color_manual(values = c('black', 'red')) +
			theme(legend.position = 'none')
	
	print(p)
		
	ifns_up = df %>% 
		filter(enr < -2 & ifn_gene) %>%
		pull(gene_name)
		
	fwrite(list(ifns_up), paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/ifn_genes_selected.txt'), sep = '\t')
		
	return(ifns_up)
}

plotTEmodMC = function(tes_ds, te_mods, mc_object, mc_ordering, mod_ordering, cell_size, s_score, m_score, smoothed = TRUE, mcs1, mcs2)
{
	tes_tidy_ds = tes_ds %>% as.matrix %>% as.data.frame.table %>% as.tibble %>% rename(unique_id = Var1, CB = Var2, n = Freq)
	te_mods = vectToTib(te_mods) %>% rename(unique_id = id, mod = value)
	mc_df = vectToTib(mc_object@mc) %>% rename(CB = id, mc = value)
	csize_df = vectToTib(cell_size / 1000) %>% rename(CB = id, csize = value)
	#te_load_df = vectToTib(te_load * 100) %>% rename(CB = id, load = value)
	s_score_df = vectToTib(s_score) %>% rename(CB = id, g1s = value)
	m_score_df = vectToTib(m_score) %>% rename(CB = id, g2m = value)
	
	te_mod_df = inner_join(tes_tidy_ds, te_mods) %>% 
		inner_join(mc_df) %>%
		mutate(mc = factor(mc, mc_ordering),
		       mod = factor(mod, mod_ordering)) %>% 
		group_by(mod, CB) %>% 
			summarise(mod_umis = sum(n), mc = unique(mc)) %>% 
		group_by(mod) %>% 
			mutate(fc = log2((mod_umis + 1) / mean(mod_umis + 1))) %>% 
			ungroup
		# group_by(mc, mod) %>% 
			# summarise(mean_fc = mean(fc), 
			          # first_quartile = quantile(fc, 0.25),
					  # third_quartile = quantile(fc, 0.75)) %>% 
		# ungroup
			   
	anno_df = full_join(csize_df, s_score_df) %>%
		full_join(m_score_df) %>%
		full_join(mc_df) %>% 
		replace_na(list(csize = 0, g1s = 0, g2m = 0)) %>%
		filter(!is.na(mc)) %>%
		mutate(mc = factor(mc, mc_ordering)) %>%
		gather(type, 'value', csize, g1s, g2m) %>% 
		# group_by(mc, type) %>%
			# summarise(mean_val = mean(value), sd_val = sd(value)) %>%
			ungroup
			
	sig_levels = te_mod_df %>% 
		group_by(mod) %>% 
			summarise(pval = wilcox.test(fc[mc %in% mc_ordering[mcs1]], fc[mc %in% mc_ordering[mcs2]])$p.value) %>% 
		mutate(qval = p.adjust(pval), sig = ifelse(qval < 0.05, '*', 'ns'), sig = ifelse(qval < 0.005, '**', sig), sig = ifelse(qval < 0.0005, '***', sig)) %>%
		select(mod, sig)
		
	if (smoothed)
	{
		p1 = te_mod_df %>%
			ggplot(aes(y = fc, x = as.numeric(mc))) +
				#geom_line() +
				geom_smooth(method = 'loess', span = 0.4, se = FALSE) +
				#geom_ribbon(aes(ymin = first_quartile, ymax = third_quartile), alpha = 0.1, colour=NA) +
				#coord_cartesian(ylim = c(-1.5,0.5)) +
				facet_wrap(~mod, ncol = 1, strip.position = 'right', drop = FALSE) + 
				scale_x_continuous(breaks = seq(1, length(mc_ordering), 1), labels = mc_ordering) +
				theme(strip.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), panel.spacing = unit(0.5, "lines"), 
					  legend.position = 'none',  strip.text.y = element_blank()) +
				geom_abline(intercept = 0, slope = 0, lty = 2, col = 'grey') +
				geom_text(data = sig_levels, mapping = aes(x = -Inf, y = -Inf, label = sig), x = 5.5, y = 0)
			
		p2 = anno_df %>%
			ggplot(aes(y = value, x = as.numeric(mc))) +
				geom_smooth(method = 'loess', span = 0.4, se = FALSE) +
				#geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), alpha = 0.1, colour=NA) +
				facet_wrap(~type, ncol = 1, scales = 'free_y', strip.position = 'right') + 
					scale_x_continuous(breaks = seq(1, length(mc_ordering), 1), labels = mc_ordering) +
					theme(strip.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), panel.spacing = unit(0.5, "lines"), 
						  legend.position = 'none',  strip.text.y = element_blank())
	} else {
		p1 = te_mod_df %>%
			ggplot(aes(y = fc, x = as.numeric(mc), group = mc)) +
				#geom_boxplot(outlier.shape = NA) +
				geom_line() + 
				stat_summary(fun.y = 'mean') +
				#geom_ribbon(aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc), alpha = 0.1, colour=NA) +
				#ylim(-1.5, 1.5) +
				facet_wrap(~mod, ncol = 1, strip.position = 'right') + 
				scale_x_continuous(breaks = seq(1, length(mc_ordering), 1), labels = mc_ordering) +
				theme(strip.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), panel.spacing = unit(0.05, "lines"), 
					  legend.position = 'none',  strip.text.y = element_blank())
			
		p2 = anno_df %>%
			ggplot(aes(y = value, x = as.numeric(mc), group = mc)) +
				geom_boxplot(outlier.shape = NA) +
				#geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), alpha = 0.1, colour=NA) +
				facet_wrap(~type, ncol = 1, scales = 'free_y', strip.position = 'right') + 
				scale_x_continuous(breaks = seq(1, length(mc_ordering), 1), labels = mc_ordering) +
				theme(strip.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), panel.spacing = unit(0.05, "lines"), 
					  legend.position = 'none',  strip.text.y = element_blank())
	}
	
	p_final = plot_grid(p1, p2, rel_heights = c(10, 3), ncol = 1, align = 'hv', axis = 'l')
						  
	return(p_final)
}

getCCGenes = function(mat_id = 'genes_m_filt', gene_anchors = h1299_dacsb_cc_anchors, gene_anti = h1299_dacsb_cc_anti, cor_threshold = 0.2, downsamp_n = 7500, 
                      ignore_cells = NULL)
{
	cc_genes_df = mcell_mat_rpt_cor_anchors_custom(mat_id, gene_anchors, gene_anti, downsample_n = downsamp_n, ignore_cells = ignore_cells)
	cc_genes = cc_genes_df %>% filter(max > neg_max & max >= cor_threshold) %>% pull(final_id)
	
	names(cc_genes) = cc_genes
	cc_gset = tgGeneSets(cc_genes, "cell cycle")
	scdb_add_gset("cc_genes", cc_gset)

	# set_param(param = 'scm_n_downsamp_gstat', value = downsamp_n, package = 'metacell')
	# mcell_gset_split_by_dsmat(gset_id = 'cc_genes', mat_id = mat_id, K = 12)

	# mcell_plot_gset_cor_mats(gset_id = 'cc_genes', scmat_id = mat_id, downsamp = TRUE)
	# dev.off()
	# mcell_gset_remove_clusts(gset_id = 'cc_genes', filt_clusts = keep_clusts, new_id = "cc_genes_filt", reverse=T)
	
	res = scdb_gset('cc_genes')@gene_set
	return(res)
}

runMeta = function(mat_id = 'genes_m_filt', gset_id = 'cc_genes_final')
{
	n_resamp = 500
	gstat_id = paste0(cell_line, '_', mat_id, '_cc_gstat')
	graph_id = paste0(cell_line, '_', mat_id, '_cc_graph')
	coc_id = paste0(cell_line, '_', mat_id, '_cc_coc')
	mc_id = paste0(cell_line, '_', mat_id, '_cc_mc')
	mc2d_id = paste0(cell_line, '_', mat_id, '_cc_mc2d')
	#mc_id_filt = 'h1299_genes_dacsb_filt'
	
	mcell_add_gene_stat(mat_id = mat_id, gstat_id = gstat_id, force = TRUE)
	#gene_whitelist = scdb_gstat(gstat_id) %>% filter(tot >= 25 & ds_top3 >= 2 & ds_vm_norm >= 0.1 & name %in% names(scdb_gset(gset_id)@gene_set)) %>% pull(name)
	#scdb_add_gset(gset_id, gset_new_restrict_nms(scdb_gset(gset_id), gene_whitelist, inverse = FALSE, desc = 'cc_genes_filt'))

	mcell_add_cgraph_from_mat_bknn(graph_id = graph_id, mat_id = mat_id, gset_id = gset_id, K = 100, dsamp = TRUE)
	mcell_coclust_from_graph_resamp(coc_id = coc_id, graph_id = graph_id, min_mc_size = 100, p_resamp = 0.75, n_resamp = n_resamp) # creates metas min mc_size = 100 
	mcell_mc_from_coclust_balanced(mc_id = mc_id, coc_id = coc_id, mat_id = mat_id, K = 30, min_mc_size = 100, alpha = 2) # min_mc_size = 80
	#mcell_mc_split_filt(new_mc_id = mc_id_filt, mc_id = mc_id, mat_id = mat_id, T_lfc = 5, plot_mats = F) # IF ERROR TRY LOAD DBSCAN
	set_param('mcell_mc2d_K', value = 20, package = 'metacell') # hct 20 
	set_param('mcell_mc2d_max_confu_deg', value = 2.5, package = 'metacell') # hct 2.5 
	set_param('mcell_mc2d_T_edge', value = 0.05, package = 'metacell') # 0.05 hct 
	set_param('mcell_mc2d_proj_blur', value = 0.01, package = 'metacell') # 0.01
	
	mcell_mc2d_force_knn(mc2d_id = mc2d_id, mc_id = mc_id, graph_id = graph_id)
	# par(mfrow = c(2,1), mar = rep(0,4))
	# plotUMIsOn2D(mc2d_id =  mc2d_id, colSums(genes_m_filt_ds[names(cc_mods[cc_mods== 2]),]))
	# plotUMIsOn2D(mc2d_id =  mc2d_id, colSums(genes_m_filt_ds[names(cc_mods[cc_mods == 4]),]))
	
	#plotting
	mcell_gset_from_mc_markers(gset_id = 'markers', mc_id = mc_id)
	mcell_mc_plot_marks(mc_id = mc_id, gset_id = gset_id, mat_id = mat_id, plot_cells = TRUE)
	mcell_mc_plot_marks(mc_id = mc_id, gset_id = gset_id, mat_id = mat_id, plot_cells = FALSE)
	mcell_mc2d_force_knn(mc2d_id = mc2d_id, mc_id = mc_id, graph_id = graph_id)
	mcell_mc2d_plot(mc2d_id = mc2d_id)
	mcell_mc_plot_subheats(mc_id = mc_id, mat_id = mat_id)
}

plotCorGeneTE = function(gene_mat, te_load, top_n = 15)
{
	if (identical(colnames(gene_mat), names(te_load)))
	{
		cors = apply(gene_mat, 1, function(x) cor(x, te_load))
	} else {
		print ('warning: incompatible barcode names')
	}
	
	top_genes = c(head(sort(cors), top_n), tail(sort(cors), top_n))
	
	# plotting
	barplot(top_genes, 
		las =2, ylim = c(-0.6, 0.6), 
		col = rep(c('darkblue', 'darkred'), times = c(top_n, top_n)), 
		ylab = 'Correlation to TE load')
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

plotTEFPs = function(mc_te, col_order = NULL, te_modules = ct_te)
{
	mat_to_plot = log2(mc_te@mc_fp[rev(names(te_modules)),])
	mat_to_plot[mat_to_plot > 1.5] = 1.6
	mat_to_plot[mat_to_plot < -1.5] = -1.6
	
	if(!is.null(col_order))
	{
		mat_to_plot = mat_to_plot[, col_order]	
	}
		
	gaps_row = as.numeric(cumsum(rev(table(te_modules))))
	
	te_anno = annotateTECustom(rownames(mat_to_plot))
	
	annotation_rows = data.frame(row.names = rownames(mat_to_plot), 
								 'TE' = te_anno$custom)
	annotation_cols = list(TE = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'ERV9' = "#4DAF4A", 'ERVL-MaLR/ERVL' = "#004529", 'ERV1/ERVK' = "#FFFFE5", 'SINE' = "#984EA3"))
								 

	cols = c('darkblue', rev(colorRampPalette(c(brewer.pal(11, 'BrBG')))(256)), 'black') 
	
	pheatmap(mat_to_plot, 
			 cluster_cols = FALSE,
			 cluster_rows = FALSE,
			 gaps_row = gaps_row,
			 show_rownames = FALSE,
			 col = cols, 
			 breaks = c(-1.6, seq(-1.5, 1.5, l = 257), 1.6), 
			 clustering_method = 'ward.D',
			 treeheight_row = 0,
			 border = NA,
			 annotation_row = annotation_rows,
			 annotation_colors = annotation_cols)
}

plotMCTEload = function(tes_tidy, genes_tidy, mcs_te = mcs_te, mc_order = gene_clustering)
{
	mcs = data.frame(mcs_te) %>% rownames_to_column() %>% rename(CB = rowname, mc = mcs_te) %>% as_tibble

	te_umis = left_join(mcs, tes_tidy) %>% mutate(custom = annotateTECustom(unique_id)$custom, custom = ifelse(custom == repname, repclass, custom)) %>% group_by(mc, custom) %>% summarise(te_umis = sum(n)) %>% ungroup
	genic_umis = left_join(mcs, genes_tidy) %>% group_by(mc) %>% summarise(gene_umis = sum(n)) %>% ungroup

	mc_te_load = full_join(te_umis, genic_umis) %>% mutate(te_load = te_umis / (te_umis + gene_umis) * 100) %>% 
		select(mc, custom, te_load) %>% spread(mc, te_load) %>% tibble.to.matrix
		
	cols = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'ERV9' = "#4DAF4A", 'ERVL-MaLR/ERVL' = "#004529", 'ERV1/ERVK' = "#FFFFE5", 'SINE' = "#984EA3")
	arranged = c('SINE', 'DNA', 'LINE', 'LTR', 'ERV9', 'ERV1/ERVK', 'ERVL-MaLR/ERVL')
		
	barplot(mc_te_load[arranged, mc_order], col = cols[arranged], las = 2, ylim = c(max(colSums(mc_te_load) + 2), 0), ylab = 'TE load (%)')
}

plot2DProjectionsWrapper = function(mc2d, mc_cc, ct_te, te_load, cell_size, s_score, m_score, s_genes)
{	
	png(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/2D_mc_id_', treatment, '.png'), 
		width = 5000, height = 5000, res = 1500, pointsize = 6)
	par(mar = rep(0,4), bty = 'n', font = 2)
	
	score_value = log2(colMeans(mc_cc@mc_fp[m_genes,]))
		cols = c('darkblue', colorRampPalette(c('darkblue', 'lightblue', 'orange', 'red'))(100), 'darkred')[as.numeric(cut(score_value, breaks = c(-Inf, seq(-.75, 0.75, l = 100), Inf)))]
	plot(mc2d@sc_x, mc2d@sc_y, cex = 1, col = alpha('lightgrey', 0.25), pch = 20, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
	par(new = TRUE)
	plot(mc2d@mc_x, mc2d@mc_y, cex = 8, col = cols, pch = 15, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
	#text(mc2d@mc_x, mc2d@mc_y, labels = names(mc2d@mc_x), cex = 3, col = 'white')
	dev.off()
	
	png(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/2D_te_load_', treatment, '.png'), 
		width = 5000, height = 5000, res = 1500, pointsize = 6)
	par(mar = rep(0,4), bty = 'n')
	plotUMIsOn2D(mc2d, te_load)
	dev.off()

	png(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/2D_cell_size_', treatment, '.png'), 
		width = 5000, height = 5000, res = 1500, pointsize = 6)
	par(mar = rep(0,4), bty = 'n')
	plotUMIsOn2D(mc2d, cell_size)
	dev.off()

	png(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/2D_s_score_', treatment, '.png'), 
			width = 5000, height = 5000, res = 1500, pointsize = 6)
	par(mar = rep(0,4), bty = 'n')
	plotUMIsOn2D(mc2d, s_score)
	dev.off()

	png(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/2D_m_score_', treatment, '.png'), 
			width = 5000, height = 5000, res = 1500, pointsize = 6)
	par(mar = rep(0,4), bty = 'n')
	plotUMIsOn2D(mc2d, m_score)
	dev.off()

	for(i in unique(ct_te))
	{
		umis = colSums(mat_ds[names(ct_te[ct_te == i]), ])

		png(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/2D_te_mod_', treatment, '_', i, '.png'), 
			width = 5000, height = 5000, res = 1500, pointsize = 6)
		par(mar = rep(0,4), bty = 'n')
		plotUMIsOn2D(mc2d, umis)
		dev.off()
	}
	dev.off()
}

calcDiffGeneSig = function(genes_m_filt, mcs)
{
	calcMCCounts = function(matrix, mcs, sum_stat = sum)
	{
		mc_counts = t(apply(matrix, 1, function(x) tapply(x, mcs, sum_stat)))
		return(mc_counts)
	}
	
	shared_barcodes = intersect(colnames(genes_m_filt), names(mcs))
	genes_m_filt = genes_m_filt[, shared_barcodes]
	mcs = mcs[shared_barcodes]
	
	mc_counts = calcMCCounts(genes_m_filt, mcs) %>%
		as.data.frame.table %>%
		rename(gene = Var1, mc = Var2, umis = Freq) %>%
		as_tibble
		
	cont_table_tidy = mc_counts %>% 
		group_by(gene) %>% 
			mutate(umis_outside = sum(umis) - umis) %>% 
		group_by(mc) %>% 
			mutate(others_inside = sum(umis) - umis) %>% 
			ungroup %>% 
		mutate(others_outside = sum(umis) - umis - umis_outside - others_inside)	
		
	doMC::registerDoMC(25)
	
	p_vals = plyr::ddply(cont_table_tidy, plyr::.(mc), function(x) 
	{ 
		x %>%
			rowwise %>%
				do({ broom::tidy(chisq.test(matrix(c(.$umis, .$umis_outside, .$others_inside, .$others_outside), byrow =TRUE, ncol =2)))})
	}, .parallel = TRUE) %>% 
		select(p.value)
		
	stats = cont_table_tidy %>% 
		mutate(p_val = p_vals$p.value, q_val = p.adjust(p_val, 'fdr'))	
	
	return(stats)
}

filterCCGenes = function(mat_genes = genes_m_filt, cc_genes_all, cc_refs, nclust = 20, filter_mods = TRUE, downsamp_n = NULL)
{
	mat_ds_genes = scm_downsamp(mat_genes, downsamp_n)
	mat_trans = as.matrix(log2(1 + 7 * mat_ds_genes))
	#mat_trans = as.matrix(mat_ds_genes)
	cor_m = tgs_cor(t(mat_trans[cc_genes_all,]))

	d = tgs_dist(cor_m)
	#d = as.dist(1 - cor_m)
	hc = hclust(d, 'ward.D')
	ct = cutree(hc, nclust)
	
	annotations = data.frame(marker = as.character(rownames(cor_m) %in% cc_refs), row.names = rownames(cor_m))

	cols = c('black', colorRampPalette(c("darkblue", 'blue', 'white', 'orange', "red"))(100), 'darkred')
	pheatmap(cor_m, breaks = c(-1, seq(-0.5, 0.5, l = 98), 1), clustering_method = 'ward.D', cutree_row = nclust, cutree_col = nclust, 
		annotation_row = annotations, annotation_col = annotations, col = cols, annotation_colors = list('marker' = c('FALSE' = 'lightgrey', 'TRUE' = 'black')))
	
	if(filter_mods == TRUE)
	{
	cc_genes_final = sort(ct[ct %in% unique(ct[cc_refs])]) 
	} else {
	cc_genes_final = sort(ct)
	}
	return(cc_genes_final)
}

plotDomainEnr = function(ct_te, tes_tidy, domain_size = 1.5e7, outdir = 'path')
{
	whole_genome = GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
	ref_genome = whole_genome[seqlevels(whole_genome) %in% paste0('chr', c(1:22, 'X', 'Y', 'MT'))]
	
	ref_genome_bins = tile(ref_genome, width = domain_size) %>% unlist
	
	te_loci = tibble(unique_id = names(ct_te), mod = ct_te) %>% left_join(., tes_tidy %>% select(unique_id, chrom, start, end) %>% distinct) %>%
		mutate(chrom = paste0('chr', chrom)) 
	te_loci_gr = te_loci %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE)
	
	hits = findOverlaps(te_loci_gr, ref_genome_bins)
	te_loci = te_loci %>% mutate(genome_bin = as.data.frame(ref_genome_bins[hits@to]) %>% mutate(bin_id = paste0(seqnames, ': ', round(start/1e6 + 1,1), ' - ', round(end/1e6 + 1, 1), ' Mb')) %>% pull(bin_id))
	
	domain_counts = te_loci %>% count(mod, genome_bin)
	
	cont_table_tidy = domain_counts %>% 
		filter(n > 1) %>%
		group_by(mod) %>% 
			mutate(n_outside = sum(n) - n) %>% 
		group_by(genome_bin) %>% 
			mutate(others_inside = sum(n) - n) %>% 
			ungroup %>% 
		mutate(others_outside = sum(n) - n - n_outside - others_inside)
		
	
	# p_vals = cont_table_tidy %>%
			# rowwise %>%
				# do({ broom::tidy(chisq.test(matrix(c(.$n, .$n_outside, .$others_inside, .$others_outside), byrow =TRUE, ncol =2), simulate.p.value = TRUE)$expected[1,1])})
		
	# pheatmap(log2(cont_table_tidy %>% mutate(expected = p_vals$x, 'obs_vs_exp' = n/expected) %>% select(mod, genome_bin, obs_vs_exp) %>% spread(genome_bin, obs_vs_exp, fill = 1) %>% tibble.to.matrix), cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(-3, seq(-2, 2, l = 98), 3), col = colorRampPalette(c('darkblue', 'white', 'darkred'))(99))
	
	statistics = cont_table_tidy %>%
			rowwise %>%
				do({ broom::tidy(fisher.test(matrix(c(.$n, .$n_outside, .$others_inside, .$others_outside), byrow =TRUE, ncol =2), alternative = 'greater'))}) %>%
			ungroup %>%
			select(p.value) %>%
			mutate(q.value = p.adjust(p.value, 'fdr'),
			       enriched = q.value < 0.05) %>%
			bind_cols(cont_table_tidy, .)
			
	
	m = statistics %>%
			filter(enriched) %>% 
			select(mod, genome_bin, q.value) %>%
			spread(genome_bin, q.value, fill = 1) %>%
			tibble.to.matrix
			
	mod_counts = left_join(domain_counts, statistics %>% select(mod, genome_bin, enriched)) %>% 
		replace_na(list(enriched = FALSE)) %>% 
		group_by(mod, enriched) %>% 
			summarise(n = sum(n)) %>%
			ungroup
	
	# plotting
	annotation = data.frame('chrom' = strsplit2(colnames(m), '\\:', 1) , row.names = colnames(m))
	cols = structure(c("#984EA3", "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 	"#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", 'darkred', 'darkorange2'), names = paste0('chr', c(1:22, 'X', 'Y')))
	
	pdf(paste0(outdir, '/te_domain_enrichment_heatmap_', domain_size, '_.pdf'))
	pheatmap(m[rev(rownames(m)),], 
		cluster_rows = FALSE, 
		cluster_cols = FALSE, 
		breaks = c(0, 0.005, 0.05, 1), 
		col = c('darkred', 'orange', 'white'), 
		border = NA,
		annotation_col = annotation,
		annotation_colors = list('chrom' = cols))
	dev.off()
	
	pdf(paste0(outdir, '/te_domain_enrichment_barplot_', domain_size, '_.pdf'))
	mod_counts %>%
		spread(enriched, n, fill = 0, sep = '_') %>%
		tibble.to.matrix %>%
		t() %>%
		barplot(beside = FALSE, border = NA, las = 2, horiz = TRUE, xlim = c(0, 500), xlab = 'Loci', ylab = 'TE module', col = c('lightgrey', 'red'))
	dev.off()
}

plotModTOR = function(te_mods, te_full_length)
{
	te_loci = tibble(row_id = as.numeric(strsplit2(names(te_mods), '\\|', 3)), mod = te_mods) %>% left_join(., te_full_length %>% select(row_id, chrom, start, end, strand) %>% distinct)
	te_loci_gr = te_loci %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE)
		
	te_loci$repliseq = calcCov(te_loci_gr, repliseq, summary_fun = mean)

	#order mods by median repliscore
	mod_ordering = te_loci %>% group_by(mod) %>% summarise(repliseq_med = median(repliseq)) %>% arrange(repliseq_med) %>% pull(mod)
	
	p = te_loci %>% 
		ggplot(aes(x = factor(mod, levels = mod_ordering), y = repliseq, group = mod, fill = 'blue')) + 
			geom_boxplot(outlier.shape = NA) + 
			coord_flip() +
			scale_fill_manual(values = 'lightblue') +
			theme(legend.position = 'none')

	return(list(p, rev(mod_ordering)))
}

plotDomainComparison = function(mcs, genes_m_filt_ds, mat_ds, mcs1 = c(1:11, 15:20, 22:23), mcs2 = c(12:14, 21), 
								tes_tidy, genes_tidy, startp = 31500000, endp = 46500000)
{
	diff_genes = diff_expr(mcs, genes_m_filt_ds, mcs1, mcs2, min_max_umi = 10, reg = 5)
	diff_genes_anno = left_join(diff_genes %>% rename(final_id = gene), 
		genes_tidy %>% mutate(domain = seqnames == '19' & start >= startp & end <= endp) %>% select(final_id, domain) %>% distinct) %>%
		mutate(type = 'gene')

	diff_tes = diff_expr(mcs, mat_ds, mcs1, mcs2, min_max_umi = 10, reg = 5)
	diff_tes_anno = left_join(diff_tes %>% rename(unique_id = gene), 
		tes_tidy %>% mutate(domain = chrom == '19' & start >= startp & end <= endp) %>% select(unique_id, domain) %>% distinct) %>%
		mutate(type = 'te') %>%
		rename(final_id = unique_id)
		
	combined = bind_rows(diff_genes_anno, diff_tes_anno)
	combined = combined %>% mutate(cols = paste(domain, type, sep = '_'))
	
	p_scatterplot = combined %>% 
		ggplot(aes(tot1, tot2, size = domain, col = cols, alpha = cols)) + 
			geom_point_rast(size = 2) + 
			scale_y_log10() + 
			scale_x_log10() + 
			scale_color_manual(values = c('orange', 'lightblue', 'red', 'blue')) + 
			scale_alpha_discrete(range = c(0.0075, 1, 1, 1)) + 
			scale_size_discrete(range = c(0.5, 1)) 
			#scale_shape_manual(values=c(3,15))

	#boxplot
	p_boxplot = combined %>% 
		mutate(fc = c(tot2 + 5) / c(tot1 + 5)) %>% 
			ggplot(aes(y = fc, x = as.factor(cols), group = as.factor(cols), fill = type)) + 
				coord_cartesian(ylim = c(0,5)) +
				geom_boxplot(outlier.shape = NA) + 
				scale_fill_manual(values = c('orange', 'lightblue')) + 
				geom_abline(intercept = 2, slope = 0, lty = 2, col = 'lightgrey')
				
	p = plot_grid(p_scatterplot, p_boxplot, rel_widths = c(2, 1))
	return(p)
}

plotTEvsIFN = function(genes_tidy, tes_tidy, mc_object, ifn_genes)
{
	# add MC info
	mc_df = mc_object@mc %>% 
		data.frame %>% 
		rownames_to_column %>% 
		rename(CB = rowname, mc = '.')
		
	genes_tidy = inner_join(genes_tidy, mc_df)
	tes_tidy = inner_join(tes_tidy, mc_df)

	ifn_counts = genes_tidy %>% 
		group_by(mc) %>% 
			summarise(genic_umis = sum(n), 
			ifn_umis = sum(n[gene_name %in% ifn_genes])) 

	te_counts = tes_tidy %>% 
		group_by(mc, repname) %>% 
			summarise(te_umis = sum(n)) %>% 
			ungroup

	mc_loads = inner_join(ifn_counts, te_counts)  %>% 
		mutate(ifn_perc = ifn_umis / genic_umis * 100, 
			   te_perc = te_umis  / (te_umis + genic_umis) * 100)
		
	# plotting
	p1 = mc_loads %>%
		group_by(mc) %>%
			summarise(te_perc = sum(te_perc), 
			          ifn_perc = unique(ifn_perc)) %>%
		ggplot(aes(te_perc, ifn_perc, label = mc)) + 
			geom_text_repel() +
			geom_point_rast(size = 2) 
			
	# cor per fam
	p2 = mc_loads %>% 
		group_by(repname) %>% 
		summarise(cors = cor(te_perc, ifn_perc, method = 'pearson'), 
		          total_te_umis = sum(te_umis)) %>%
		left_join(., select(tes_tidy, repname, repclass) %>% distinct) %>%
		ggplot(aes(total_te_umis, cors, col = repclass, label = ifelse(cors > 0.7, repname, ''))) + 
			geom_point_rast(size = 2) + 
			geom_text_repel() +
			scale_y_continuous(limits = c(-1, 1)) +
			scale_x_log10(breaks = c(1, 100, 10000, 1000000), labels = c(1, 100, 10000, 1000000)) +
			scale_color_brewer(palette = 'Set1') +
			theme(legend.position = 'none')

	p_final = plot_grid(p1, p2)
	return(p_final)
}