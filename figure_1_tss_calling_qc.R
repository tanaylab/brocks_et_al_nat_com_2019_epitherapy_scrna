source('/net/mraid14/export/data/users/davidbr/general/scripts/utilities.R')

#########################
# FUNCTIONS
#########################
plotCellSize = function(comb_genic, comb_tes)
{
	genic_cell_size = comb_genic %>% 
		group_by(cell_line, treatment, CB) %>% 
			summarise(cell_size = sum(n)) %>% 
			ungroup %>%
		mutate(type = 'genic')
	
	te_cell_size = comb_tes %>% 
		group_by(cell_line, treatment, CB) %>% 
			summarise(cell_size = sum(n)) %>% 
			ungroup %>%
		mutate(type = 'te')

	comb_cell_size = bind_rows(genic_cell_size, te_cell_size) %>%
		mutate(treatment = factor(treatment, levels = c('dmso', 'dacsb')))
	
	# plot
	dodge <- position_dodge(width = 0.25)
	
	p = comb_cell_size %>%
		ggplot(aes(x = type, y = cell_size, fill = treatment, alpha = type)) +
			geom_split_violin() +
			geom_boxplot(width=.1, outlier.colour=NA, position = dodge) + 
			coord_cartesian(ylim = c(1, 1e5)) + 
			scale_y_log10(breaks = c(1, 10, 100, 1000, 1e4, 1e5)) +
			facet_wrap(~cell_line) +
			scale_alpha_discrete(range = c(1, 0.25)) + 
			scale_fill_brewer(palette = 'Pastel1') 

	return(p)
}

plotTELoad = function(comb_genic, comb_tes)
{
	genic_cell_size = comb_genic %>% 
		group_by(cell_line, treatment, CB) %>% 
			summarise(cell_size_g = sum(n)) %>% 
			ungroup
	
	te_cell_size = comb_tes %>% 
		group_by(cell_line, treatment, CB) %>% 
			summarise(cell_size_te = sum(n)) %>% 
			ungroup

	comb_te_load = full_join(genic_cell_size, te_cell_size, by = c('cell_line', 'treatment', 'CB')) %>%
		mutate(treatment = factor(treatment, levels = c('dmso', 'dacsb'))) %>%
		replace_na(list(cell_size_g = 0, cell_size_te = 0)) %>%
		mutate(te_load = cell_size_te / (cell_size_te + cell_size_g))
		
	# plot
	p = comb_te_load %>%
		ggplot(aes(x = te_load, alpha = treatment, fill = treatment)) + 
			geom_histogram(alpha = 0.5, position = 'identity', bins = 100) + 
			facet_wrap(~cell_line) + 
			coord_cartesian(xlim = c(0.0001, 1), ylim = c(0, 300)) + 
			scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.01, 0.1, 1)) + 
			scale_alpha_discrete(range = c(0.15, 1)) +
			scale_fill_brewer(palette = 'Pastel1') 
	
	return(p)
}
plotTELoadPerClass = function(comb_genic, comb_tes)
{
	genic_cell_size = comb_genic %>% 
		group_by(cell_line, treatment, CB) %>% 
			summarise(cell_size_g = sum(n)) %>% 
			ungroup
	
	te_cell_size = comb_tes %>% 
		group_by(cell_line, treatment, CB, repclass) %>% 
			summarise(cell_size_te = sum(n)) %>% 
			ungroup

	comb_te_load_class = inner_join(genic_cell_size, te_cell_size, by = c('cell_line', 'treatment', 'CB')) %>%
		mutate(treatment = factor(treatment, levels = c('dmso', 'dacsb'))) %>%
		replace_na(list(cell_size_g = 0, cell_size_te = 0)) %>%
		mutate(te_load = cell_size_te / (cell_size_te + cell_size_g))
		
	# plot
	dodge <- position_dodge(width = 0.25)
	
	p = comb_te_load_class %>%
		ggplot(aes(x = repclass, y = te_load, alpha = treatment, fill = repclass)) + 
			geom_split_violin() + 
			geom_boxplot(width = .1, outlier.colour = NA, position = dodge) + 
			facet_wrap(~cell_line) + 
			coord_cartesian(ylim = c(1e-5, 1)) + 
			scale_y_log10(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.01, 0.1, 1), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.01, 0.1, 1)) +
			scale_alpha_discrete(range = c(0.15, 1)) +
			scale_fill_brewer(palette = 'Set1')
	
	return(p)
}
plotNTSS = function(comb_genic)
{
	ntss = comb_genic %>%
		group_by(cell_line, treatment) %>% 
			summarise(total_tss = length(unique(final_id[tss]))) %>% 
			ungroup %>%
		mutate(treatment = factor(treatment, levels = c('dmso', 'dacsb')))
		
	p = ntss %>%
		ggplot(aes(x = treatment, y = total_tss, fill = treatment)) +
			geom_bar(stat = 'identity') +
			facet_wrap(~cell_line) +
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			scale_fill_brewer(palette = 'Pastel1')
			
	return(p)
}
plotTSSComp = function(comb_genic)
{
	tss_h1299 = comb_genic %>% 	
		filter(tss & cell_line == 'h1299') %>% 
		select(gene_name, start, treatment) %>% 
		distinct %>%
		rename(start_h1299 = start)
		
	tss_hct116 = comb_genic %>% 	
		filter(tss & cell_line == 'hct116') %>% 
		select(gene_name, start, treatment) %>% 
		distinct %>%
		rename(start_hct116 = start)
		
	tss_dist_counts = inner_join(tss_h1299, tss_hct116, by = c('gene_name', 'treatment')) %>%
		group_by(gene_name, treatment) %>%
			summarise(min_dist = min(abs(start_h1299 - start_hct116))) %>%
			ungroup %>% 
		mutate(dist_bin = cut(min_dist, breaks = c(0, 20, 200, 2000, Inf), include.lowest = TRUE) %>% as.numeric) %>%
		count(treatment, dist_bin)
		
	p = tss_dist_counts %>%
		mutate(treatment = factor(treatment, levels = c('dmso', 'dacsb'))) %>%
		ggplot(aes(x = dist_bin, y = n, fill = treatment)) + 
			geom_bar(stat = 'identity') +
			facet_wrap(~treatment) +
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			scale_fill_brewer(palette = 'Pastel1')
			
	return(p)
}
plotTSSEnr = function(comb_genic)
{
	tss_enr = comb_genic %>% 
		group_by(cell_line, treatment, CB) %>% 
			summarise(tss_enr = sum(n[tss]) / sum(n)) %>% 
			ungroup

	p = tss_enr %>%
		mutate(treatment = factor(treatment, levels = c('dmso', 'dacsb'))) %>%
		ggplot(aes(x = tss_enr, fill = treatment)) +
			geom_histogram(alpha = 0.5, position = 'identity', bins = 100) + 
			facet_wrap(~cell_line) + 
			scale_x_continuous(limits = c(0.25, 1), breaks = c(0.25, 0.5, 0.75, 1)) + 
			scale_alpha_discrete(range = c(0.15, 1)) + 
			scale_fill_brewer(palette = 'Pastel1')
			
	return(p)
}

plotTSSPos <- function(comg_genic, comb_genic_tss, cell_line_id = 'h1299', reg = 1)
{
  tes_filt <-
    comb_tes_tss %>% 
    filter(cell_line == cell_line_id)

  genes_filt <-
    comb_genic %>%
    filter(tss) %>%
    filter(cell_line == cell_line_id) %>%
    select(seqnames, start, end, strand, final_id, gene_name, treatment) %>%
    distinct
  
  loci_derep <-  
    tes_filt %>%
    group_by(treatment, row_id) %>% 
    summarise(umis = sum(n)) %>% 
    ungroup %>% 
    spread(treatment, umis) %>% 
    replace_na(list(dmso = 0)) %>% mutate(fc = log2((dacsb + reg) / (dmso + reg))) %>% 
    filter(fc >= 2) %>% 
    pull(row_id)
    
  loci_derep_coords <-
    tes_filt %>%
    filter(row_id %in% loci_derep) %>%
    select(chrom, start, end, strand) %>%
    distinct %>%
    makeGRangesFromDataFrame
  
  gene_tss_coords <-
    genes_filt %>%
    mutate(tss_center = round((start + end) / 2)) %>%
    group_by(gene_name) %>%
    summarize(tss_dist = min(abs(tss_center[treatment == 'dmso'] - tss_center[treatment == 'dacsb'])),
              seqnames = seqnames[1],
              start = list(unique(tss_center)),
              end = list(unlist(start),
              strand = strand[1]) %>%
    filter(tss_dist != Inf) %>% # exclude genes that only occur in one treatment
    unnest(start, end) %>%
    mutate(te_dist = Inf) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
              
 
  hits <- distanceToNearest(gene_tss_coords, loci_derep_coords)
  gene_tss_coords[hits@from]$te_dist <- mcols(hits)$distance

  ggdata <-
    gene_tss_coords %>%
    as_tibble %>%
    group_by(gene_name) %>%
    summarise(tss_dist = unique(tss_dist),
              te_dist = min(te_dist)) %>%
    mutate(dist_bin = cut(te_dist, breaks = c(1, 1000, 1e4, 1e5), include.lowest = TRUE)) %>%
    filter(!is.na(dist_bin))
    
  ggdata %>%
    ggplot(aes(tss_dist, col = dist_bin)) +
      stat_ecdf(geom = "step") +
      scale_x_log10()
   
  # p value
  wilcox.test(ggdata$tss_dist[ggdata$te_dist <= 1000], ggdata$tss_dist[ggdata$te_dist > 1000])

}

#########################
# LOAD
#########################
gene_anno = import.gff('~/davidbr/general/data/gencode/gencode.v28.annotation.gtf')

# Genic Counts
h1299_genes_dmso_tidy = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/counts/gene_counts_tidy.RDS') %>% 
	mutate(cell_line = 'h1299', treatment = 'dmso', type = 'genic')
hct116_genes_dmso_tidy = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dmso/counts/gene_counts_tidy.RDS') %>% 
	mutate(cell_line = 'hct116', treatment = 'dmso', type = 'genic')
h1299_genes_dacsb_tidy = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/counts/gene_counts_tidy.RDS') %>% 
	mutate(cell_line = 'h1299', treatment = 'dacsb', type = 'genic')
hct116_genes_dacsb_tidy = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/counts/gene_counts_tidy.RDS') %>% 
	mutate(cell_line = 'hct116', treatment = 'dacsb', type = 'genic')
	
comb_genic = bind_rows(h1299_genes_dmso_tidy, hct116_genes_dmso_tidy, h1299_genes_dacsb_tidy, hct116_genes_dacsb_tidy)	

# TE Counts
h1299_tes_dmso = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/counts/te_bin_counts.RDS') %>% 
	mutate(cell_line = 'h1299', treatment = 'dmso', type = 'te')
h1299_tes_dacsb = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/counts/te_bin_counts.RDS') %>% 
	mutate(cell_line = 'h1299', treatment = 'dacsb', type = 'te')
hct116_tes_dmso = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dmso/counts/te_bin_counts.RDS') %>% 
	mutate(cell_line = 'hct116', treatment = 'dmso', type = 'te')
hct116_tes_dacsb = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/counts/te_bin_counts.RDS')	%>% 
	mutate(cell_line = 'hct116', treatment = 'dacsb', type = 'te')
	
comb_tes = bind_rows(h1299_tes_dmso, h1299_tes_dacsb, hct116_tes_dmso, hct116_tes_dacsb) %>%
	filter(grepl('LTR|SINE|DNA|LINE', repclass))
  
# TE counts TSS
h1299_tes_dmso = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/counts/te_tss_bin_counts_tidy.RDS') %>% 
	mutate(cell_line = 'h1299', treatment = 'dmso', type = 'te')
h1299_tes_dacsb = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/counts/te_tss_bin_counts_tidy.RDS') %>% 
	mutate(cell_line = 'h1299', treatment = 'dacsb', type = 'te')
hct116_tes_dmso = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dmso/counts/te_tss_bin_counts_tidy.RDS') %>% 
	mutate(cell_line = 'hct116', treatment = 'dmso', type = 'te')
hct116_tes_dacsb = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/counts/te_tss_bin_counts_tidy.RDS')	%>% 
	mutate(cell_line = 'hct116', treatment = 'dacsb', type = 'te')
  
comb_tes_tss = bind_rows(h1299_tes_dmso, h1299_tes_dacsb, hct116_tes_dmso, hct116_tes_dacsb) %>%
	filter(grepl('LTR|SINE|DNA|LINE', repclass))

###############################################################
####### PLOTTING #############
###############################################################
p_cell_size_violinplot = plotCellSize(comb_genic, comb_tes)
p_te_load_hist = plotTELoad(comb_genic, comb_tes)
p_te_load_repclass_violin = plotTELoadPerClass(comb_genic, comb_tes)
p_te_load_repclass_tss_violin = plotTELoadPerClass(comb_genic, comb_tes_tss)
p_ntss_barplot = plotNTSS(comb_genic)
p_tss_dist_comp_barplot = plotTSSComp(comb_genic)
p_tss_enr_hist = plotTSSEnr(comb_genic)



#############################################################
####### EXPORT ##########
#############################################################
ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/genic_tss_calling/genic_te_cell_size_violin.pdf', p_cell_size_violinplot, device = 'pdf', width = 8, height = 4)

ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/genic_tss_calling/te_load_hist.pdf', p_te_load_hist, device = 'pdf', width = 8, height = 4)

ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/genic_tss_calling/te_load_repclass_violin.pdf', p_te_load_repclass_violin, device = 'pdf', width = 8, height = 4)

ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/genic_tss_calling/te_load_repclass_tss_violin.pdf', p_te_load_repclass_tss_violin, device = 'pdf', width = 8, height = 4)

ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/genic_tss_calling/ntss_all_conditions_barplot.pdf', p_ntss_barplot, device = 'pdf', width = 4, height = 4)

ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/genic_tss_calling/tss_cell_line_overlap_dist.pdf', p_tss_dist_comp_barplot, device = 'pdf', width = 4, height = 4)

ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/genic_tss_calling/tss_enr_sc_hist.pdf', p_tss_enr_hist, device = 'pdf', width = 8, height = 4)
