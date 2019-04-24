source('/net/mraid14/export/data/users/davidbr/general/scripts/utilities.R')

#########################
# FUNCTIONS
#########################
addAnnotation = function(sel_loci, tes_tss, genes_tidy, te_mods, mc_object)
{
	major_tss_bin = tes_tss %>% 
		select(repname, bin_id, row_id) %>% 
		distinct %>% count(repname, bin_id) %>% 
		group_by(repname) %>% 
			top_n(1, n) %>% 
			ungroup %>% 
		select(-n) %>%
		rename(major_tss_bin = bin_id)
	
	tss_gr = genes_tidy %>% 
		filter(tss) %>%  
		select(seqnames, start, end, strand) %>% 
		distinct %>% mutate(seqnames = paste0('chr', seqnames)) %>% 
		makeGRangesFromDataFrame
		
	sel_loci = left_join(sel_loci, major_tss_bin)
	sel_loci$marker = sel_loci$row_id %in% strsplit2(names(te_mods), '\\|', 3)
	sel_loci$dist_tss = mcols(distanceToNearest(makeGRangesFromDataFrame(sel_loci), tss_gr))$distance
	sel_loci$damid = calcCov(sel_loci, damid, summary_fun = mean)
	sel_loci$repliseq = calcCov(sel_loci, repliseq, summary_fun = mean)
	sel_loci = sel_loci %>%
		mutate(repname_custom = ifelse(grepl('Alu', repname), 'SINEs', repname), 
		       repname_custom = ifelse(grepl('LTR12', repname), 'ERV9s', repname_custom),
			   repname_custom = ifelse(grepl('THE1|MST', repname), 'ERVL-MaLRs', repname_custom),
			   repname_custom = ifelse(grepl('LTR10|LTR13', repname), 'ERV1/ERVK', repname_custom))
	sel_loci = addMCFP(tes_tss, sel_loci, mc_object)
	return(sel_loci)
}

msaToLong <- function(msa)
{
  n_loci <- length(msa)
  
  if (is.null(names(msa)) | sum(duplicated(names(msa))) > 0)
  {
    row_labels <- as.character(1:n_loci)
    names(msa) <- row_labels
  } else {
    row_labels <- names(msa)
  }
  
  
  m <- as.matrix(msa)
  colnames(m) <- 1:ncol(m)
  
  #sm <- Matrix::Matrix(m, sparse = TRUE)
  msa_long <- 
    m %>% 
    base::as.data.frame() %>% 
    tibble::rownames_to_column() %>% 
    gather(pos, base, -rowname) %>% 
    dplyr::rename(te_id = rowname) %>% 
    as_tibble %>% 
    #filter(base != '-') %>%
    mutate(pos = as.numeric(pos),
           te_id = factor(te_id, levels = row_labels),
           is_gap = base == '-') %>%
    data.table::as.data.table()
    
	msa_long <- msa_long[order(te_id, pos), ]
  
  # Add column with ungapped pos integer
  msa_long <- 
    msa_long[, ':='(pos_wo_gaps = 1:.N), by = c('te_id', 'is_gap')
    ][, ':='(gap_perc = 1 - sum(base != '-') / n_loci), by = 'pos'] %>%
    mutate(pos_wo_gaps = ifelse(base == '-', NA, pos_wo_gaps)) %>%
    as_tibble
  
  return(msa_long)
} 

addAlignment = function(sel_loci_anno, fam_to_align = c(), type = 'fast')
{
	aligned = sel_loci_anno %>% 
		filter(repname %in% fam_to_align) %>%
			group_by(repname) %>% 
				do(addMSA(., type = type, threads = 18) %>% unnest) %>% 
			ungroup

	res = full_join(aligned, sel_loci_anno)
	return(res)
}

addMCFP = function(tes_tss, sel_loci, mc_object)
{
	mc_compute_fp_custom = function(mc, us)
	{
        f_g_cov = rowSums(us) > 10

        mc_cores = 16
        doMC::registerDoMC(mc_cores)
        all_gs = rownames(us[f_g_cov,])
        n_g = length(all_gs)
        g_splts = split(all_gs, 1+floor(mc_cores*(1:n_g)/(n_g+1)))
        fnc = function(gs) {
                                        .row_stats_by_factor(us[gs,],
                                                                        mc@mc,
                                                                        function(y) {exp(rowMeans(log(1+y)))-1}) }

        clust_geomean = do.call(rbind, mclapply(g_splts, fnc, mc.cores = mc_cores))

        mc_meansize = tapply(colSums(us), mc@mc, mean)
        ideal_cell_size = pmin(1000, median(mc_meansize))
        g_fp = t(ideal_cell_size*t(clust_geomean)/as.vector(mc_meansize))
        #normalize each gene
        fp_reg = 0.1
        #0.1 is defined here because 0.1*mean_num_of_cells_in_cluster
        #is epxected to be 3-7, which means that we regulairze
        #umicount in the cluster by 3-7.
        g_fp_n = (fp_reg+g_fp)/apply(fp_reg+g_fp, 1, median)

        return(g_fp_n)
	}

	#create sparse matrix
	counts_tidy = tes_tss %>% 
		group_by(row_id, CB) %>% 
			summarise(umis = sum(n)) %>% 
			ungroup %>%			
			mutate(row_id = factor(row_id), CB = factor(CB))
	counts_sparse = sparseMatrix(as.integer(counts_tidy$row_id), 
		                         as.integer(counts_tidy$CB), 
								 x = counts_tidy$umis,
								 dimnames = list(levels(counts_tidy$row_id), levels(counts_tidy$CB)))
	######################
	
	shared_barcodes = intersect(colnames(counts_sparse), names(mc_object@mc))
	
	mc_fps = mc_compute_fp_custom(mc_object, counts_sparse[, shared_barcodes])
	colnames(mc_fps) = paste0('mc_', 1:ncol(mc_fps))

	mc_fps = data.frame(mc_fps) %>% 
		rownames_to_column %>% 
		as_tibble %>% 
		rename(row_id = rowname) %>% 
		mutate(row_id = as.numeric(row_id))
	
	res = left_join(sel_loci, mc_fps)
	
	return(res)
}

addMSA = function(intervs, type = 'fast', max_gap_perc = 0.99, threads = 18)
{
	removeGaps = function(alignment, max_gap_perc)
	{
	alignment_bin = as.DNAbin(alignment)
	alignment_bin_wo_gaps = del.colgapsonly(as.matrix(alignment_bin), max_gap_perc)
	
	res = as.character(as.list(alignment_bin_wo_gaps)) %>% 
		lapply(., paste0, collapse = '') %>% 
		unlist %>% 
		DNAStringSet
	return(res)
	}
	
	intervs_gr = makeGRangesFromDataFrame(intervs)
	names(intervs_gr) = intervs$row_id
	
	seqs = getSeq(Hsapiens, intervs_gr)
	intervs$seqs_full = as.character(seqs)
	
	alignment = mafft(seqs, save_alignment = FALSE, type = type, threads = threads)
	alignment_df = alignment %>%	
		data.frame(row_id = as.numeric(names(.)), 'alignment_full' = .) %>%
		as_tibble
		
	alignment_wo_gaps = removeGaps(alignment, max_gap_perc)
	
	alignment_wo_gaps_df = alignment_wo_gaps %>%	
		data.frame(row_id = as.numeric(names(.)), 'alignment_wo_gaps' = .) %>%
		as_tibble
		
	intervs_and_alignment = inner_join(intervs, alignment_df, by = 'row_id')
		
	intervs_and_alignment = inner_join(intervs_and_alignment, alignment_wo_gaps_df, by = 'row_id')
	intervs_and_alignment = intervs_and_alignment %>% mutate(seqs_wo_gaps = gsub('-', '', alignment_wo_gaps))
		
	#res = clustAlignPerMod(intervs_and_alignment, 'indelblock')
	res = intervs_and_alignment
	
	return(res)
}

addTrackValues = function(intervs, track = 'encode.repliseq.wgEncodeUwRepliSeqBg02esWaveSignalRep1', fun = 'global.percentile', sshift = NULL, eshift = NULL)
{
	intervs_transformed = intervs %>% select(chrom, start, end, strand) %>% convertStrand
	gvtrack.create('tor', track, func = fun)
	#gvtrack.create('tor', 'encode.repliseq.wgEncodeUwRepliSeqK562WaveSignalRep1', func = 'global.percentile')
	#gvtrack.create('tor', 'encode.repliseq.wgEncodeUwRepliSeqImr90WaveSignalRep1', func = 'global.percentile')
	if (!is.null(sshift) | !is.null(eshift))
	{
		gvtrack.iterator('tor', sshift = sshift, eshift = eshift)
	}
	
	tor = gextract('tor', intervals = intervs_transformed, iterator = intervs_transformed) %>%
		arrange(intervalID)
		
	intervs$track_value = tor$tor
	return(intervs)	
}

clustAlignment = function(seqs, nclust = NULL, method = 'indelblock')
{
	seqs = DNAStringSet(seqs)
	d = dist.dna(as.DNAbin(seqs), model = method)
	hc = hclust(d, 'ward.D')
	if (!is.null(nclust)) 
	{
		ct = cutree(hc, nclust)
		res = as.numeric(ct)
	} else 
	{
		res = order((1:length(seqs))[hc$order])
	}
	
	return(res)
}

corByDist <- function(comb_genic, comb_tes_tss, min_umis = 10, cell_line_id = 'h1299')
{
  downsamp_n <- ifelse(cell_line_id == 'h1299', 7500, 5000)
  
  # Downsample counts
  genes_ds <- 
    comb_genic %>% 
    filter(cell_line == cell_line_id & treatment == 'dacsb') %>% 
    filter(grepl('protein_coding', gene_type)) %>%
    filter(n_tss == 0 | tss) %>%
    filter(!grepl('RP11', gene_name)) %>%
    filter(seqnames != 'MT') %>%
    dplyr::rename(unique_id = final_id) %>% 
    downsamp(i = 'unique_id', j = 'CB', n = 'n', full = FALSE, ds = downsamp_n) %>%
    filter(tss)
  
  tes_ds <- 
    comb_tes_tss %>% 
    filter(cell_line == cell_line_id & treatment == 'dacsb') %>%
    downsamp(i = 'unique_id', j = 'CB', n = 'n', full = FALSE)
    
  shared_bc <- intersect(genes_ds$CB, tes_ds$CB)
  message(length(shared_bc))
  
  # summarise umis per cell and locus/gene
  te_counts <-
    tes_ds %>%
      group_by(unique_id, CB) %>%
      summarise(umis = sum(n_ds)) %>%
      ungroup
  
  gene_counts <-
    genes_ds %>%
      group_by(unique_id, CB) %>%
      summarise(umis = sum(n_ds)) %>%
      ungroup

  # combine into matrix and calc correlation
  mat_comb <-
    bind_rows(te_counts, gene_counts) %>%
    filter(CB %in% shared_bc) %>%
    tidyToSparse
  
  mat_comb_filt <- mat_comb[rowSums(mat_comb) >= min_umis, ]
  mat_comb_filt <- log2(1 + 7 * mat_comb_filt)
  mat_cor <- tgs_cor(t(as.matrix(mat_comb_filt)), spearman = FALSE)
  
  # only retain TE gene pairs of cor mat
  mat_cor_filt <- mat_cor[rownames(mat_cor) %in% te_counts$unique_id, !colnames(mat_cor) %in% te_counts$unique_id]
  
  # tidy cor matrix
  long_cor <-
    as.data.frame(mat_cor_filt) %>% 
    rownames_to_column(var = 'te_id') %>% 
    gather(-te_id, key = 'gene_id', value = 'cor') %>% 
    as_tibble
   
  # Create gene and TE intervals
  gene_intervals <- 
    genes_ds %>% 
    filter(unique_id %in% long_cor$gene_id) %>%
    select(seqnames, start, end, strand, unique_id) %>% 
    distinct %>% 
    rename(chrom = seqnames) %>% 
    mutate(chrom = paste0('chr', chrom)) %>% 
    filter(chrom != 'chrMT') %>% 
    mutate(start = start - 1) %>%
    dplyr::rename(gene_id = unique_id, gene_strand = strand)
    
  te_intervals <- 
    tes_ds %>% 
    filter(unique_id %in% long_cor$te_id) %>%
    select(chrom, start, end, strand, unique_id) %>% 
    distinct %>% 
    mutate(chrom = paste0('chr', chrom)) %>% 
    mutate(start = start - 1) %>% 
    dplyr::rename(te_id = unique_id, te_strand = strand)
  
  # Get distance for all pairs
  distances <- gintervals.neighbors(gene_intervals %>% select(-gene_strand), te_intervals %>% select(-te_strand), maxdist = 100001, maxneighbors = 1e9)[,c(1,4,8,9)]
  
  # add anno to cor data.frame
  long_cor_dist <- inner_join(long_cor, distances)
  long_cor_dist <- left_join(long_cor_dist, select(gene_intervals, gene_id, gene_strand))
  long_cor_dist <- left_join(long_cor_dist, select(te_intervals, te_id, te_strand))
  long_cor_dist <- long_cor_dist %>% mutate(strandedness = ifelse(gene_strand == te_strand, 'same', 'opposite'))
  long_cor_dist <- long_cor_dist %>% left_join(tes_ds %>% select(unique_id, repclass) %>% dplyr::rename(te_id = unique_id) %>% distinct) # add repclass anno
  
  # plotting 
  p <-
    long_cor_dist %>%
      mutate(dist_bin = cut(dist, breaks = c(1, 1e3, 1e4, 1e5), include.lowest = TRUE),
             label = ifelse(cor >= 0.2, paste(gene_id, te_id), '')) %>%
      filter(!is.na(dist_bin)) %>%
      ggplot(aes(x = dist_bin, y = cor, group = dist_bin, col = strandedness, label = label)) + 
        geom_boxplot(outlier.colour = NA) +
        geom_jitter(size = 1, alpha = 0.05) +
        facet_wrap(~strandedness) +
        coord_cartesian(ylim = c(-0.1, 0.2)) + 
        theme(axis.text.x=element_text(angle=90, hjust=1),
              legend.position = 'none')
  
  return(p)
  # Stats
  # wilcox.test(long_cor_dist %>% filter(strandedness == 'opposite') %>% filter(dist <= 500) %>% pull(cor), 
              # long_cor_dist %>% filter(strandedness == 'opposite') %>% filter(dist >  500) %>% pull(cor))
  # wilcox.test(long_cor_dist %>% filter(strandedness == 'same') %>% filter(dist <= 500) %>% pull(cor), 
              # long_cor_dist %>% filter(strandedness == 'same') %>% filter(dist >  500) %>% pull(cor))
}

plotAlignment = function(df, split_column = NULL)
{
	f = function(alignment)
	{
		consensus_matrix = consensusMatrix(alignment)
		colnames(consensus_matrix) = 1:ncol(consensus_matrix)
		res = consensus_matrix %>% as.data.frame.table(stringsAsFactors = FALSE) %>% rename(base = Var1, pos = Var2, freq = Freq) %>% mutate(pos = as.numeric(pos))
		return(res)
	}	
	consensus_freqs = f(df$alignment_wo_gaps) %>% group_by(pos) %>% mutate(freq = freq / sum(freq)) %>% ungroup %>% rename(consensus_freq = freq)
	indv_freqs = df %>% slice %>% rowwise %>% do(test = f(.$alignment_wo_gaps)) %>% ungroup %>% mutate(row_id = df$row_id) %>% unnest
	
	m = left_join(indv_freqs, consensus_freqs) %>% 
		filter(freq != 0) %>% 
		mutate(divergence = freq - consensus_freq,
		       divergence = ifelse(base == '-', NA, divergence)) %>% 
		select(pos, row_id, divergence) %>% 
		spread(pos, divergence) %>% 
		tibble.to.matrix
		
	m = m[as.character(df$row_id),]
	if (!is.null(split_column))
	{
		rownames(m) =  df %>% pull(split_column)
		n_clusts = length(unique(rownames(m)))
		p_heights = rownames(m) %>% table
		# plotting
		for (index in unique(rownames(m)))
		{
			
			png(paste0('msa_', index, '.png'), width = 1250, height = p_heights[index])
			par(mar = c(0,0,0,0))
			to_plot_mat = m[rownames(m) %in% index, ]
			image(y = 0:nrow(to_plot_mat), x = 0:ncol(to_plot_mat), z = t(to_plot_mat), 
				col = colorRampPalette(c('lightgrey', 'orange', 'darkred', 'black'))(256), 
				zlim = c(0,1), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
			dev.off()
		}
	} else {
	#par(mar = rep(0, 4))
	image(y = 0:nrow(m), x = 0:ncol(m), z = t(m), 
		col = colorRampPalette(c('lightgrey', 'orange', 'darkred', 'black'))(256), 
		zlim = c(0,1), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	#axis(1, at = c(seq(0, ncol(m), 500)), labels = seq(0, ncol(m), by = 500), las = 2, lwd = 0, lwd.ticks = 1)
	#axis(3, at = unique(df$tss_start_on_alignment), labels = 'TSS')
	}
}

plotDistTSS = function(sel_loci_anno)
{
	# plot dist TSS
	sig_levels = sel_loci_anno %>% 
		group_by(repname_custom) %>% 
			summarise(p_value = wilcox.test(dist_tss[expressed == 'expressed'], dist_tss[expressed == 'repressed'])$p.value) %>% 
		mutate(q_value = p.adjust(p_value, 'fdr')) %>% 
		mutate(annotation = ifelse(q_value < 0.05, '*', 'ns'), 
		       annotation = ifelse(q_value < 0.005, '**', annotation), 
			   annotation = ifelse(q_value < 0.0005, '***', annotation)) %>% 
		select(repname_custom, annotation) %>% 
		tibble.to.matrix

	p = sel_loci_anno %>% 
		filter(expressed != 'unclear') %>%
		ggplot(aes(x = repname_custom, y = dist_tss, fill = expressed)) + 
			geom_boxplot(outlier.shape = NA, coef = 0.5) + 
			coord_cartesian(ylim = c(0, 2e6)) + 
			#scale_y_continuous(limits = c(0, 1.75e6), breaks = c(0, 5e5, 1e6, 1.5e6), labels = c(0, 5e5, 1e6, 1.5e6))  + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			scale_fill_manual(values = c('lightblue', 'orange')) + 
			geom_signif(annotations = as.character(sig_levels), xmin = c(1:nrow(sig_levels) - 0.2), xmax = c(1:nrow(sig_levels) + 0.2), y_position  = 1.75e6, tip_length = 0.005)
			
	return(p)
}

plotLamin = function(sel_loci_anno)
{
	# plot dist TSS
	sig_levels = sel_loci_anno %>% 
		group_by(repname_custom) %>% 
			summarise(p_value = wilcox.test(damid[expressed == 'expressed'], damid[expressed == 'repressed'])$p.value) %>% 
		mutate(q_value = p.adjust(p_value, 'fdr')) %>% 
		mutate(annotation = ifelse(q_value < 0.05, '*', 'ns'), 
		       annotation = ifelse(q_value < 0.005, '**', annotation), 
			   annotation = ifelse(q_value < 0.0005, '***', annotation)) %>% 
		select(repname_custom, annotation) %>% 
		tibble.to.matrix

	p = sel_loci_anno %>% 
		filter(expressed != 'unclear') %>%
		ggplot(aes(x = repname_custom, y = damid, fill = expressed)) + 
			geom_boxplot(outlier.shape = NA, coef = 0.5) + 
			scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 1), labels = c(0, 1))  + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			scale_fill_manual(values = c('lightblue', 'orange')) + 
			geom_signif(annotations = as.character(sig_levels), xmin = c(1:nrow(sig_levels) - 0.2), xmax = c(1:nrow(sig_levels) + 0.2), y_position  = 1.15, tip_length = 0.05)
			
	return(p)
}

plotMeth = function(sel_loci_anno)
{

	sel_loci_anno$meth_dmso = calcCov(sel_loci_anno, dmso_meth, summary_fun = mean)
	sel_loci_anno$meth_dacsb = calcCov(sel_loci_anno, dacsb_meth, summary_fun = mean)
	#sel_loci_anno = sel_loci_anno %>% mutate(meth_delta = meth_dacsb - meth_dmso)

	to_plot = sel_loci_anno %>% 
		gather(treatment, 'meth', meth_dmso, meth_dacsb)
		
	p = to_plot %>%
		filter(expressed != 'unclear') %>%
		mutate(treatment = factor(treatment, levels = c('meth_dmso', 'meth_dacsb')),
		       condition = factor(paste(expressed, treatment, sep = '_'), levels = c('expressed_meth_dmso', 'expressed_meth_dacsb', 'repressed_meth_dmso', 'repressed_meth_dacsb'))) %>%
		ggplot(aes(x = condition, y = meth, group = condition, fill = expressed, alpha = treatment)) +
			geom_boxplot(outlier.shape = NA) + 
			facet_wrap(~repname_custom) +
			scale_fill_manual(values = rev(c('darkorange', 'lightblue'))) +
			theme(axis.text.x = element_blank())
			
	return(p)
}

compareGroups <- function(fam_loci, tes_tss, cell_line = NULL)
{
  if (cell_line == 'h1299')
  {
    mc_cc = get(load('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/metacell/genes/cell_cycle/mc.h1299_genes_m_filt_filt_cc_mc.Rda'))
    mc_ordering = c(8,6,7,9,10,3,4,5,1,2)
    mcs <- vectToTib(mc_cc@mc) %>% rename(CB = id) %>% mutate(mc = factor(value, levels = mc_ordering))
  } else {
 		mc_cc = get(load('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/metacell/genes/cell_cycle/mc.hct116_genes_m_filt_cc_mc.Rda'))
		mc_ordering = c(12,8,11,7,6,10,3,4,5,2,1,9)
    mcs <- vectToTib(mc_cc@mc) %>% rename(CB = id) %>% mutate(mc = factor(value, levels = mc_ordering))
    mcs <- mcs %>% mutate(cycling = value %in% c(3,4,5,2,1,9))
  }
  
  print(cell_line)
    
  # Preparing and downsampling TEs
  tes_tss_sum <- # summarise umis per row_id
    tes_tss %>% 
      select(row_id, CB, n) %>% 
      group_by(row_id, CB) %>% 
      summarise(n = sum(n)) %>% 
      ungroup
  
  te_downsamp_n <- tes_tss_sum %>% group_by(CB) %>% summarise(n = sum(n)) %>% pull(n) %>% median
  print(te_downsamp_n)
        
  tes_tss_ds <-
    tes_tss_sum %>%
      downsampTidy(., te_downsamp_n)

  ggdata <-
    inner_join(fam_loci %>% select(seq_clust, row_id) %>% distinct, tes_tss_ds) %>%
    group_by(seq_clust, CB) %>%
    summarise(umis = sum(n)) %>%
    ungroup %>%
    spread(seq_clust, umis, sep = '_', fill = 0) %>%
    inner_join(., mcs)
  
  p <-  
    ggdata %>%
    ggplot(aes(seq_clust_1, seq_clust_2, col = mc)) +
      geom_point_rast(size = 2, alpha = 0.5) +
      scale_color_manual(values = rev(colorRampPalette(c('darkblue', 'blue', 'orange', 'red', 'darkred'))(length(mc_ordering))))
  
  return(p)
}

plotGroupCors <- function(fam_loci, tes_tss, genes_tidy, min_expr = 25, downsamp_n = 7500, cor_delta = 0.2)
{
  print(downsamp_n)
  # Calc gene TE dist
  genes_gr <- genes_tidy %>% select(seqnames, start, end, final_id) %>% distinct %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  tes_gr <- fam_loci %>% filter(expressed == 'expressed') %>% mutate(chrom = gsub('chr', '', chrom)) %>% makeGRangesFromDataFrame
  
  dists <- distanceToNearest(genes_gr, tes_gr)
  gene_te_dists <- tibble(final_id = genes_gr[dists@from]$final_id, TE_dist = mcols(dists)$distance)
  
  # Preparing and downsampling TEs
  tes_tss_sum <- # summarise umis per row_id
    tes_tss %>% 
      select(row_id, CB, n) %>% 
      group_by(row_id, CB) %>% 
      summarise(n = sum(n)) %>% 
      ungroup
  
  te_downsamp_n <- tes_tss_sum %>% group_by(CB) %>% summarise(n = sum(n)) %>% pull(n) %>% median
  
  tes_tss_ds <-
    tes_tss_sum %>%
      downsampTidy(., te_downsamp_n)
      
  clust_counts_ds <-
    inner_join(fam_loci %>% select(row_id, seq_clust), tes_tss_ds) %>%
      group_by(seq_clust, CB) %>%
      summarise(n = sum(n)) %>%
      ungroup %>% 
      spread(seq_clust, n, fill = 0, sep = '_')    
      
  # Preparing and downsamplig genes
  tf_symbols <- tfs %>% filter(tf) %>% pull(gene_symbol) %>% unique
  
  genes_tidy_filt <- genes_tidy %>% 
		filter(grepl('protein_coding', gene_type)) %>% 
		filter(n_tss == 0 | tss) %>%
		filter(!grepl('RP11', gene_name)) %>%
		filter(seqnames != 'MT') 
  
  genes_tidy_filt_ds <- downsampTidy(select(genes_tidy_filt, final_id, CB, n), downsamp_n) %>% # filter lowly expressed genes
    group_by(final_id) %>% mutate(total = sum(n)) %>% ungroup %>% filter(total >= min_expr) %>% # select TFs only
    mutate(tf = strsplit2(final_id, '\\|', 1) %in% tf_symbols) %>% filter(tf)
  
  
  # Combine and calc cor
  clust_counts_ds_all <- inner_join(clust_counts_ds, genes_tidy_filt_ds %>% select(CB, final_id, n)) %>%
    replace_na(list(seq_clust_1 = 0, seq_clust_2 = 0, n = 0))
  
  cors <-
    clust_counts_ds_all %>%
    group_by(final_id) %>%
    summarise(cor1 = cor(n, seq_clust_1, method = 'spearman'), 
              cor2 = cor(n, seq_clust_2, method = 'spearman')) %>%
    mutate(delta = abs(cor1 - cor2))
    
  cors_anno <-
    left_join(cors, gene_te_dists) %>%
    replace_na(list(TE_dist = 0))
  
  p <-
    cors_anno %>%
    mutate(dist_bin = as.character(as.numeric(cut(TE_dist, breaks = c(0, 1000 ,100000, 1e8), include.lowest = TRUE)))) %>%
    ggplot(aes(cor1, cor2, label = ifelse(delta >= cor_delta, as.character(final_id), ''))) + 
      geom_point_rast(size = 2) +
      geom_text_repel(aes(col = dist_bin)) +
      scale_color_manual(values = c('black', 'orange', 'red')) +
      theme(legend.position = 'none')
    
  return(p)
}

plotMCFPs = function(df, split_column = 'comb_clust')
{
	m = df %>% select(row_id, grep('mc_', colnames(df))) %>% tibble.to.matrix %>% log2
	m[is.na(m)] = 1000
	
	cols = c('darkblue', rev(colorRampPalette(c(brewer.pal(11, 'BrBG')))(256)), 'lightgrey')
	
	if (!is.null(split_column))
	{
		#m = m[nrow(m):1,]
		rownames(m) <- df %>% pull(split_column)
		n_clusts = length(unique(rownames(m)))
		p_heights = rownames(m) %>% table
		#p_heights[p_heights <= 10] = p_heights[p_heights <= 10] * 5 # rescale
		# plotting
		for (index in unique(rownames(m)))
		{
			png(paste0('mc_fp_', index, '.png'), width = 300, height = p_heights[index])
			par(mar = c(0,0,0,0))
			to_plot_mat = m[rownames(m) %in% index, ]
			image(y = 0:nrow(to_plot_mat), x = 0:ncol(to_plot_mat), z = t(to_plot_mat), 
				xaxt = 'n', 
				yaxt = 'n',
				breaks = c(c(-999, seq(-1.15, 1.15, l = 256), 999), 1000), 
				col = cols,
				xlab = '',
				ylab = '')
			dev.off()
		}
	}
}

plotFamSeqHighlight = function(sel_loci_anno, fam = 'LTR12C', nclust = 2, outdir = 'path')
{
	fam_loci = sel_loci_anno %>%	
		filter(repname == fam) 
		
	fam_loci = fam_loci %>% 
		mutate(tor_early = repliseq > 0.5)
		
	# cluster sequences
	fam_loci = fam_loci %>% 
		mutate(seq_clust = clustAlignment(alignment_wo_gaps, nclust = nclust, method = 'indelblock')) %>%
		arrange(seq_clust, tor_early, expressed)
		
	print(fam_loci %>% group_by(seq_clust) %>% summarise(pval = wilcox.test(tss_TRUE[tor_early], tss_TRUE[!tor_early])$p.value) %>% mutate(qval = p.adjust(pval)))
		
	n_clusters = fam_loci %>% count(seq_clust, tor_early, expressed)
	print(n_clusters)
	fam_loci$comb_clust = rep(1:nrow(n_clusters), times = n_clusters$n)
	
  clust_ordered_df <- tibble(comb_clust = 1:12, comb_clust_ordered = c(6,4,5,3,1,2,12,10,11,9,7,8))
  fam_loci <- left_join(clust_ordered_df, fam_loci) %>% mutate(comb_clust = comb_clust_ordered) %>% select(-comb_clust_ordered)

	
	mc_fps = fam_loci %>% select(row_id, grep('mc_', colnames(fam_loci))) %>% tibble.to.matrix %>% log2
	mc_fps[mc_fps > 1.15] = 1.16
	mc_fps[mc_fps < -1.15] = -1.16
	
	# cluster loci within clusts based on mc_fp
	row_ordering = as.character(unlist(tapply(rownames(mc_fps), fam_loci$comb_clust , function(x)
	{
		mc_fps_sel = mc_fps[x,] %>% 
			as.data.frame %>% 
			rownames_to_column
			
		print(sum(is.na(mc_fps_sel)))
			
		if (sum(is.na(mc_fps_sel)) == 0)
		{	
			# d = dist(mc_fps_sel)
			# hc = hclust(d, 'ward.D')
			# res = hc$label[hc$order]
			res = names(sort(tglkmeans::TGL_kmeans(mc_fps_sel, k = 6, seed = 19, reorder_func = max)$cluster))
		} else 
		{
			res = x
		}
			return(res)
	})))
	
	fam_loci_sort = tibble(row_id = as.numeric(row_ordering)) %>% left_join(., fam_loci)
	
	#plotting
	oldwd = getwd()
	setwd(outdir)
	
	plotAlignment(fam_loci_sort, split_column = 'comb_clust')
	plotMCFPs(fam_loci_sort, split_column = 'comb_clust')
	
	setwd(oldwd)
}

plotCovOnMSA <- function(fam_loci, max_gap_freq = 0.99, track = 'encode.chipseq.k562.nfya')
{
  msa <- fam_loci$alignment_full %>% DNAStringSet
  names(msa) <- fam_loci$row_id
  #msa <- Repguide:::.rmGaps(msa, max_gap_freq)
  
  msa_long <- msaToLong(msa)
  
  binned_intvs <- fam_loci %>% makeGRangesFromDataFrame %>% tile(., width = 1)
  n_elements = elementNROWS(binned_intvs)
  element_pos <- unlist(sapply(binned_intvs, function(x) 1:length(x)))
	binned_intvs = unlist(binned_intvs) %>% as_tibble %>% rename(chrom = seqnames)
	
	# ADD ANNO #
  binned_intvs$te_id <- as.character(rep(fam_loci$row_id, n_elements))
  binned_intvs$pos_wo_gaps <- element_pos
  scores <- addTrackValues(binned_intvs %>% mutate(start = start -1), track = track, fun = 'global.percentile')
  
  msa_guide_cov <- left_join(msa_long, scores) %>% left_join(., fam_loci %>% select(row_id, comb_clust) %>% rename(te_id = row_id) %>% mutate(te_id = as.character(te_id)))
  
  # PLOTTING
   # exclude alignment pos based on max gapped threshold
  msa_guide_cov_filt <-
    msa_guide_cov %>%
    filter(gap_perc <= max_gap_freq)

  # long to wide format
  m <-
	  msa_guide_cov_filt %>%
	  select(te_id, pos, track_value) %>%
	  spread(pos, track_value, fill = NA) %>%
    left_join(fam_loci %>% arrange(comb_clust) %>% select(row_id) %>% transmute(te_id = as.character(row_id)), .) %>%
    tibble.to.matrix
    
  p <-
    msa_guide_cov_filt %>% 
      group_by(comb_clust, te_id) %>% 
      summarise(mean = mean(-log2(1 - track_value), na.rm = TRUE)) %>% 
      ungroup %>% 
      ggplot(aes(as.factor(comb_clust), mean, group = comb_clust)) + 
        geom_boxplot()
  ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/', track,'_boxplot.pdf'), p, device = 'pdf')
    
  m <- -log2(1 - m)
  m[m > 15] <- 15
  
  # d <- tgs_dist(m)
  # hc <- hclust(d, 'ward.D2')
  
  # Plotting
  par(mar = rep(0, 4))
  image(y = 0:nrow(m), x = 0:ncol(m), z = t(m), 
        col =  c("#F5F5F5", colorRampPalette(c(brewer.pal(9, 'YlOrRd')))(99)),
        breaks = c(0, seq(5, 15, l = 100)),
        xaxt = 'n',
        yaxt = 'n',
        xlab = '',
        ylab = '',
        axes = FALSE)
  abline(v = c(1900, 2025), lwd = 0.5, col = 'black')
  
  aff_mean <- rowMeans(m, na.rm = TRUE)
  seq_clusts <- fam_loci %>% count(seq_clust) %>% slice(1) %>%  pull(n)
  print(wilcox.test(aff_mean[1:seq_clusts], aff_mean[(seq_clusts+1):length(aff_mean)]))
        
  #axis(1, at = c(seq(0, ncol(m), 100)), labels = seq(0, ncol(m), by = 100), las = 2, lwd = 0, lwd.ticks = 1)
}

plotSeqDivBarplot <- function(fam_loci, smoo_window = 5)
{
  a = fam_loci %>% filter(seq_clust == 1 & expressed == 'expressed') %>% pull(alignment_wo_gaps) %>% DNAStringSet %>% consensusMatrix(as.prob = TRUE) %>% as.data.frame.table %>% filter(Var1 %in% c('A', 'C', 'G', 'T', '-')) %>% mutate(Var2 = rep(1:(nrow(.)/5), each = 5)) %>% as_tibble
  b = fam_loci %>% filter(seq_clust == 2 & expressed == 'expressed') %>% pull(alignment_wo_gaps) %>% DNAStringSet %>% consensusMatrix(as.prob = TRUE) %>% as.data.frame.table %>% filter(Var1 %in% c('A', 'C', 'G', 'T', '-')) %>% mutate(Var2 = rep(1:(nrow(.)/5), each = 5)) %>% as_tibble

  comb = full_join(a,b, by = c('Var1', 'Var2')) %>% rename(base = Var1, pos = Var2, clustA = Freq.x, clustB = Freq.y)

  p = comb %>% group_by(pos) %>% summarise(divergence = max(abs(clustA - clustB))) %>% 
    mutate(div_smoo = zoo::rollmean(divergence, k = smoo_window, fill = 'extend')) %>% 
    ggplot(aes(pos, div_smoo, fill = div_smoo)) + 
      geom_bar(stat = 'identity') + 
      scale_fill_gradientn(colours = colorRampPalette(c('lightgrey', 'orange', 'darkred', 'black'))(256)) +
      scale_x_continuous(breaks = seq(0, nrow(comb)/5, 250), expand = c(0,0)) + ylim(0,1) + 
      geom_rect(xmin = 2094, xmax = 2115, ymin = 0, ymax = 1, alpha = 0.01, fill = 'lightgrey')
  
  return(p)
}

plotConsensusSeqs <- function(fam_loci, from = 1900, to = 2025)
{
  seqs1 <- fam_loci %>% filter(seq_clust == '1' & expressed == 'expressed') %>% pull(alignment_wo_gaps) %>% DNAStringSet %>% 
    subseq(., from, to) %>% consensusMatrix
  seqs2<- fam_loci %>% filter(seq_clust == '2' & expressed == 'expressed') %>% pull(alignment_wo_gaps) %>% DNAStringSet %>% 
    subseq(., from, to) %>% consensusMatrix
  
  p1 <- ggseqlogo(seqs1, seq_type = 'dna', method = 'prob') + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) + theme(axis.text.x=element_blank())
  p2 <- ggseqlogo(seqs2, seq_type = 'dna', method = 'prob') + ylim(0,1) + theme(axis.text.x=element_blank())
  
  p <- plot_grid(p2, p1, ncol = 1)
  return(p)
}


selectLoci = function(tes_all, te_full_length, sel_families = NULL, min_expression = 1, max_repressed_loci = 500, s_seed = 19)
{
	sampleLoci = function(df, nloci = 500)
	{
		if (nrow(df) > nloci & unique(df$expressed == 'repressed'))
		{
			res = sample_n(df, nloci)
		} else {
			res = df
		}
	return(res)
	}
	
	# UMI counts per TE locus
	loci_umi_counts = tes_all %>% 
		group_by(row_id, tss) %>% 
			summarise(umis = sum(n)) %>% 
			ungroup %>%
		spread(tss, umis, fill = 0, sep = '_')
		
	# annotate loci with UMI counts inside and outside putative TSS
	te_full_length_anno = full_join(te_full_length, loci_umi_counts) %>% 
		replace_na(list(tss_FALSE = 0, tss_TRUE = 0)) %>%
		mutate(expressed = ifelse(tss_TRUE >= min_expression, 'expressed', 'unclear'), 
               expressed = ifelse(tss_FALSE == 0 & tss_TRUE == 0, 'repressed', expressed))
		
	# filter based on chrom and min size and selected families and expressed/repressed
	te_full_length_filt = te_full_length_anno %>%
		filter(chrom != 'chrY') %>%
		#filter(end - start + 1 >= 100) %>%
		filter(repname %in% sel_families) 
	
	# sample repressed loci per fam
	set.seed(s_seed)
	te_full_length_filt_reduced = te_full_length_filt %>% 
		group_by(repname, expressed) %>%
			do (sampleLoci(., max_repressed_loci) %>% unnest) %>%
			ungroup 
			
	res = te_full_length_filt_reduced %>%
		arrange(row_id)
			
	return(res)
}

sqEuDist = function(a, b)
{
	res = sum(((a - b)^2))
	return(res)
}


		   
###################################
# LOAD
###################################
damid = import.bw('~/hg38/rawdata/damid/4DN/H1-hESC.bw')
damid$score = percent_rank(damid$score)
repliseq = import.bw('~/hg38/rawdata/repliseq/4DN/hct116_early_S_comb.bigWig')
repliseq$score = percent_rank(repliseq$score)
dmso_meth = import.bw('~/hg38/rawdata/wgbs/brocks_2017/h1299_dmso.bigWig')
dacsb_meth = import.bw('~/hg38/rawdata/wgbs/brocks_2017/h1299_dacsb.bigWig')

genes = import.gff('~/davidbr/general/data/gencode/gencode.v28.annotation.gtf')
te_full_length = readRDS('~/davidbr/general/data/repeatmasker/hg38/full_length_te_coords.RDS')

cell_line = 'hct116'

genes_tidy = readRDS(paste0('~/davidbr/proj/epitherapy/data/', cell_line, '/10x/dacsb/counts/gene_counts_tidy.RDS'))
tes_tss = readRDS(paste0('~/davidbr/proj/epitherapy/data/', cell_line, '/10x/dacsb/counts/te_tss_bin_counts_tidy.RDS'))
tes_tss_sm = readRDS(paste0('~/davidbr/proj/epitherapy/data/', cell_line, '/10x/dacsb/counts/te_tss_bin_counts_sm.RDS'))
tes_all = readRDS(paste0('~/davidbr/proj/epitherapy/data/', cell_line, '/10x/dacsb/counts/te_bin_counts.RDS')) %>% 
	mutate(tss = unique_id %in% tes_tss$unique_id)
	
mc_te = readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/dacsb/metacell/tes/output/metacell_object.RDS'))
te_mods = readRDS(paste0('~/davidbr/proj/epitherapy/workspaces/', cell_line, '_te_modules.RDS'))

# load('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/metacell/genes/cell_cycle/gset.cc_genes_filt.Rda')
# cc_umis = colSums(genes_m_ds[names(object@gene_set),])paste0('~/davidbr/proj/epitherapy/workspaces/', cell_line, '_te_anno_df_for_cis_figure.RDS')
####################################
# COMPUTATION
####################################
sel_fams = data.frame(te_mods) %>% 
	rownames_to_column %>% 
	mutate(repname = strsplit2(rowname, '\\|', 1)) %>% 
	count(repname) %>% 
	filter(n >= 5) %>% pull(repname)

sel_loci = selectLoci(tes_all, te_full_length, sel_families = sel_fams, min_expression = 11, max_repressed_loci = 500, s_seed = 19)
sel_loci_anno = addAnnotation(sel_loci, tes_tss, genes_tidy, te_mods, mc_object = mc_te)
sel_loci_anno = addAlignment(sel_loci_anno, fam_to_align = c('LTR12C'), type = 'fast')
saveRDS(sel_loci_anno, paste0('~/davidbr/proj/epitherapy/workspaces/', cell_line, '_te_anno_df_for_cis_figure.RDS'))

####################################
# PLOTTING
####################################
p_dist_tss = plotDistTSS(sel_loci_anno)
p_lamin = plotLamin(sel_loci_anno)
p_meth = plotMeth(sel_loci_anno)

plotFamSeqHighlight(sel_loci_anno, fam = 'LTR12C', nclust = 2, 
  outdir = paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/', cell_line, '/'))
  
p_div_barplot <- plotSeqDivBarplot(fam_loci, smoo_window = 1)
p_seqs <- plotConsensusSeqs(fam_loci, from = 1900, to = 2025)
p_cluster1_vs_2_plot <- compareGroups(fam_loci, tes_tss, cell_line)
p_tf_cors <- plotGroupCors(fam_loci, tes_tss, genes_tidy, min_expr = 25, downsamp_n = 7500, cor_delta = ifelse(cell_line == 'h1299', 0.2, 0.15)) #fam_loci can be h1299 and hct116

png('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/nfya_gm12878_chip_msa_cov.png', width = 7500, height = 5000, res = 900)
plotCovOnMSA(fam_loci, max_gap_freq = 0.99, track = 'encode.chipseq.gm12878.nfya')
dev.off()
  
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/', cell_line, '/', 'p_clust1_vs_2_scatterplot.pdf'), p_cluster1_vs_2_plot, device= 'pdf')  
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/', cell_line, '/', 'p_dist_tss.pdf'), p_dist_tss, device = 'pdf', 
  width = 10, height = 5)
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/', cell_line, '/', 'p_lamin.pdf'), p_lamin, device = 'pdf', 
  width = 10, height = 5)
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/', cell_line, '/', 'p_meth.pdf'), p_meth, device = 'pdf', 
  width = 10, height = 10)
p_comb = plot_grid(p_dist_tss, p_lamin, p_meth, nrow = 1, align ='hv', axis = 'tblrb', rel_widths = c(1,1,2))
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/', cell_line, '/', 'comb_fig.pdf'), p_comb, device = 'pdf', 
  width = 15, height = 5)
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/', cell_line, '/', 'tf_cor_scatterplot.pdf'), p_tf_cors, device = 'pdf')
ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/seq_divergence_barplot.pdf', p_div_barplot, device = 'pdf', width = 10, height = 2)  
ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/consensus_seqs.pdf', p_seqs, device = 'pdf', width = 10, height = 2)   

#################################
# Localized effect correlation
#################################
p_cor_dist <- corByDist(comb_genic, comb_tes_tss, cell_line_id = cell_line, min_umis = 10)
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/cis_regulation/', cell_line, '/', 'cor_by_dist.pdf'), p_cor_dist, device = 'pdf')


##################################
# NFY CELL LINE COMPARISON
##################################
genome_bins <- Repguide:::binGenome(Hsapiens, bin_width = 100)

k562 <- addTrackValues(genome_bins %>% as_tibble %>% mutate(start = start -1) %>% dplyr::rename(chrom = seqnames), track = 'encode.chipseq.k562.nfya', fun = 'max')
gm12 <- addTrackValues(genome_bins %>% as_tibble %>% mutate(start = start -1) %>% dplyr::rename(chrom = seqnames), track = 'encode.chipseq.gm12878.nfya', fun = 'max')

set.seed(19)
p <-
  tibble(k562 = k562$track_value, gm12 = gm12$track_value) %>% 
    filter(!is.na(k562)) %>%
    filter(!is.na(gm12)) %>%
    mutate(k562_enr = k562 > 20) %>%
    group_by(k562_enr) %>%
    sample_n(sum(.$k562_enr)) %>%
    ungroup %>%
    ggplot(aes(x = k562_enr, y = gm12, group = k562_enr)) + 
      geom_boxplot()
      
ggsave('~/davidbr/proj/epitherapy/manuscript/figures/supplements/cis_supps/nfy_enr_comparison.pdf', p, device = 'pdf', width = 2, height = 5)

###################################
