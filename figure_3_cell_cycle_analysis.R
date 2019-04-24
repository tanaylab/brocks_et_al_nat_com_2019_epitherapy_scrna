source('/net/mraid14/export/data/users/davidbr/general/scripts/utilities.R')
source('/net/mraid14/export/data/users/davidbr/proj/epitherapy/scripts/figure_4_biological_proces_link_functions.R')


#####################################
set.seed(seed = 19)
register_param(param = 'scm_n_downsamp_gstat', default_value = 5000, package = 'metacell')
###################################
# LOAD
###################################
cell_line = 'h1299'
treatment = 'dacsb'

te_full_length = readRDS('~/davidbr/general/data/repeatmasker/hg38/full_length_te_coords.RDS')
aims = fread('/net/mraid14/export/data/users/davidbr/proj/epitherapy/data/aim_genes/ubiquitous_AIMS.txt')
aims_all = unique(unlist(c(read.delim('~/davidbr/proj/epitherapy/data/aim_genes/all_AIMS.txt', header = F, row.names = 1))))
ctas = fread('~/davidbr/proj/epitherapy/data/cancer_testis_antigens/ctas.txt')
ctas_curated = na.omit(unique(c(strsplit2(ctas[,2], '\\/', 1), strsplit2(ctas[,2], '\\/', 2))))
aims_and_ctas = unique(c(aims[,2], aims_all, ctas_curated))
gene_anno = import.gff('~/davidbr/general/data/gencode/gencode.v28.annotation.gtf')
tfs = fread('/net/mraid14/export/data/users/davidbr/general/data/human_tfs/tfs.txt')[,1]
ifns_up = fread(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/ifn_genes_selected.txt'), header = FALSE)[,1]
repliseq = import.bw('~/hg38/rawdata/repliseq/4DN/hct116_early_S_comb.bigWig')
repliseq$score = percent_rank(repliseq$score)

genes_tidy = readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/',treatment,'/counts/gene_counts_tidy.RDS'))
genes_m = readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/',treatment,'/counts/gene_counts_sm.RDS'))
genes_m_filt = filterGeneMatrix(genes_m, genes_tidy)
addSCMat(genes_m_filt, 'genes_m_filt')

tes_tidy = readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/',treatment,'/counts/te_tss_bin_counts_tidy.RDS'))
tes_m = readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/',treatment,'/counts/te_tss_bin_counts_sm.RDS'))

mc_te = readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/dacsb/metacell/tes/output/metacell_object.RDS'))
mcs_te = mc_te@mc
mat_all = readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/dacsb/metacell/tes/output/umi_matrix_all.RDS'))
mat_ds = readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/dacsb/metacell/tes/output/umi_matrix_marker_ds.RDS'))
mat_ds_trans = log2(1 + 7 * as.matrix(mat_ds))
mat = readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/dacsb/metacell/tes/output/umi_matrix_marker_raw.RDS'))
##########################

if (cell_line == 'h1299')
{
	if (treatment == 'dacsb')
	{
		anchors = c('PCNA', 'E2F1' , 'TOP2A|1' , 'MKI67', 'UBE2S|1', 'ATAD2|1', 'ARL6IP1|1', 'AURKA|1', 'RRM2')
		refs = c('ATAD2|1', 'CDC20|1', 'CDKN3|1', 'CDCA8|1', 'RRM2', 'CCNB2|1', 'TOP2A|1', 'INCENP|1', 'CDKN2D|1', 'HELLS|1', 'MCM8', 'PCNA', 'DEK|1', 'CENPK|1')
		cor_threshold = 0.25
		downsamp_n = 7500
		gene_anti = c('PDCD5|1', 'TROAP|1', 'TMEM60|1')
		ignore_cells = NULL
		s_genes = c('CHAF1A', 'ATAD2|1', 'RRM2', 'HIST1H4C|1', 'CDT1', 'H2AFX|1')
		m_genes = c('ARL6IP1|1', 'UBE2C|1', 'AURKA|1', 'CDC20|1', 'MKI67')
		mc2d = get(load('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/metacell/genes/cell_cycle/mc2d.h1299_genes_m_filt_filt_cc_mc2d.Rda'))
		mc_cc = get(load('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/metacell/genes/cell_cycle/mc.h1299_genes_m_filt_filt_cc_mc.Rda'))
		#mc_ordering = c(8,6,3,4,5,1,2,10,9,7)
		mc_ordering = c(8,6,7,9,10,3,4,5,1,2)
		set_param(param = 'scm_n_downsamp_gstat', value = downsamp_n, package = 'metacell')
		mc_s = 3
		mc_m = 1
		mod_names = c('Alu I', 'Alu II', 'ERV9 I', 'ERVL-MaLR/ERVL', 'ERV9 II', 'ERV9 III', 'ERV9 IV', 'ERV9 V', 'ERV9 VII', 'ERV9 VIII')
	} 
	if (treatment == 'dmso')
	{
		anchors = c('PCNA|1', 'E2F1' , 'TOP2A' , 'MKI67', 'UBE2S|1', 'ATAD2|1', 'ARL6IP1|1', 'AURKA|1', 'RRM2')
		refs = c("ATAD2|1", "E2F1", "PCNA|1", "MKI67", "ASPM", "CENPA|1", "CKS1B|1", "CDC20|1", "RRM2", "SGO1|1", "KIF20B", "RAD21", "KIF23",  "CDCA8|1", 'SGO2|1')
		cor_threshold = 0.25
		downsamp_n = 10000
		gene_anti = NULL
		ignore_cells = NULL
		s_genes = c('ATAD2|1', 'CLSPN|1', 'E2F1', 'PCNA|1')
		mat_ds = tes_m[intersect(rownames(tes_m), names(ct_te)),]
		set_param(param = 'scm_n_downsamp_gstat', value = downsamp_n, package = 'metacell')
		mc2d = get(load('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/metacell/genes/cell_cycle/mc2d.h1299_genes_m_filt_cc_mc2d.Rda'))
		mc_cc = get(load('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/metacell/genes/cell_cycle/mc.h1299_genes_m_filt_cc_mc.Rda'))
		mc_m = 6
		mc_s = 3
		mc_ordering = c(8,9,10,4,1,2,3,5,6,7)
	}
}

if (cell_line == 'hct116')
{
	if (treatment == 'dacsb')
	{
		anchors = c('PCNA|1', 'E2F1' , 'TOP2A|1' , 'MKI67', 'UBE2S|1', 'ATAD2|1', 'ARL6IP1|1', 'AURKA|1', 'RRM2')
		refs = c('ARL6IP1|1', 'AURKA|1', 'UBE2T|1', 'ATAD2|1', 'CENPN|1', 'MKI67', 'CCNB2|1', 'CENPE', 'CDKN3|1', 'PCNA|1', 'HIST1H4C|1')
		cor_threshold = 0.14
		downsamp_n = 5000
		gene_anti = c('KPNB1', 'NEU1|1')
		ignore_cells = names(mcs_te[mcs_te == 13])
		mcell_mat_ignore_cells('genes_m_filt', 'genes_m_filt', ig_cells = names(mcs_te[mcs_te == 13]), reverse = FALSE)
		s_genes = c('PCNA|1', 'HIST1H4C|1', 'ATAD2|1', 'RRM2')
		m_genes = c('ARL6IP1|1', 'UBE2C|1', 'AURKA|1', 'CDC20|1', 'MKI67')
		mc2d = get(load('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/metacell/genes/cell_cycle/mc2d.hct116_genes_m_filt_cc_mc2d.Rda'))
		mc_cc = get(load('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/metacell/genes/cell_cycle/mc.hct116_genes_m_filt_cc_mc.Rda'))
		#mc_ordering = c(12,11,7,6,10,5,3,4,2,1,9,8)
		mc_ordering = c(12,8,11,7,6,10,3,4,5,2,1,9)
		mc_s = 5
		mc_m = 1
		set_param(param = 'scm_n_downsamp_gstat', value = downsamp_n, package = 'metacell')
		mod_names = c('Alu I', 'Alu II', 'Alu III', 'ERV9 I', 'LINE/DNA', 'ERV1/ERVK', 'ERV9 II', 'ERV9 III', 'ERV9 IV', 'ERVL-MaLR/ERVL')
	}
	if (treatment == 'dmso')
	{
		anchors = c('PCNA|1', 'E2F1' , 'TOP2A|1' , 'MKI67', 'UBE2S|1', 'ATAD2|1', 'ARL6IP1|1', 'AURKA|1', 'RRM2')
		refs = c('HIST1H4C|1', 'ATAD2|1', 'E2F1', 'PCNA|1', 'INCENP|1', 'CDCA8|1', 'SGO2|1', 'CDCA2|1', 'KIFC1|1', 'RAD21', 'FOXM1', 'MKI67', 'CKS1B|1', 'ASPM|1', 'CDC20|1', 'TUBB4B|1', 'CDKN3|1')
		cor_threshold = 0.15
		downsamp_n = 7500
		gene_anti = NULL
		ignore_cells = NULL
		s_genes = c('HIST1H1E|1', 'HIST1H4C|1', 'ATAD2|1', 'RRM2', 'E2F1', 'PCNA|1')
		set_param(param = 'scm_n_downsamp_gstat', value = downsamp_n, package = 'metacell')
		mat_ds = tes_m[intersect(rownames(tes_m), names(ct_te)),]
		mc2d = get(load('~/davidbr/proj/epitherapy/data/hct116/10x/dmso/metacell/genes/cell_cycle/mc2d.hct116_genes_m_filt_cc_mc2d.Rda'))
		mc_cc = get(load('~/davidbr/proj/epitherapy/data/hct116/10x/dmso/metacell/genes/cell_cycle/mc.hct116_genes_m_filt_cc_mc.Rda'))
		mc_s = 5
		mc_m = 1
		mc_ordering = c(4,3,10,9,5,6,8,7,1,2)
	}
}

genes_m_filt_ds = scm_downsamp(genes_m_filt, downsamp_n)
mat_cor = tgs_cor(t(as.matrix(mat_ds_trans)), spearman = FALSE)
ct_te = clusterCorMatrix(mat_cor, nclust = 10) #ct_height = 4)
##########################
# CALC SOME STATS
##########################

# number of diff genes_ds
statistics = calcDiffGeneSig(genes_m_filt, mcs = mcs_te)
statistics %>% filter(q_val < 0.01) %>% pull(gene) %>% unique %>% length

# TE mod Gene correlation
cors = calcGeneTECor(genes_m_filt_ds, te_vect = calcUMIsPerMod(mat_ds[, shared_barcodes], ct_te, sum)[5,])


####################
# CREATE GENIC/TE META AND PLOT MARKER HEATMAP
####################
createMeta(sparse_matrix = genes_m_filt, id = 'genes_te_meta', mcs_te)

set_param("scm_mc_mark_k_per_clust", value = 6, package = 'metacell') #6
set_param("scm_mc_mark_min_gene_fold", 0.6, 'metacell') #0.6
set_param("scm_mc_mark_min_gene_cov", 0.2, 'metacell') #0.2

# plot Marker Heatmap MC_FP
pdf(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/', 'mc_fp_markers.pdf'))
gene_clustering = plotMarkerGeneHeatmap(genes_tidy, te_full_length, nclust = 8)
dev.off()

# plot TE loci Heatmap MC_FP
pdf(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/', 'mc_fp_te_loci.pdf'))
plotTEFPs(mc_te, col_order = gene_clustering, te_modules = ct_te[ct_te > 2])
dev.off()

# plot TE load Barplot
pdf(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/', 'mc_te_load.pdf'))
plotMCTEload(tes_tidy, genes_tidy, mcs_te = mcs_te, mc_order = gene_clustering)
dev.off()

# plot IFN TE cor
#ifns_up = selectIFNsUp(dmso_tidy, dacsb_tidy, min_max_umi = 25, reg = 5)
p = plotTEvsIFN(genes_tidy, tes_tidy, mc_object = mc_te, ifn_genes = ifns_up)
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/supplements/trans_supps/', cell_line, '_ifn_vs_te_load.pdf'), p, device = 'pdf', width = 8, height = 4)

# plot TE mod Domain Enrichment
plotDomainEnr(ct_te, tes_tidy, domain_size = 1.5e7, outdir = paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/'))

# plot chromosome counts by module
p_chrom_counts <-
  tibble(unique_id = names(ct_te), mod = ct_te) %>% 
    left_join(., tes_tidy) %>% 
    select(mod, chrom, unique_id) %>% 
    distinct %>% 
    count(mod, chrom) %>% 
    mutate(n = cut(n, breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128), include.lowest = TRUE)) %>%
    ggplot(aes(chrom, as.factor(mod), fill = n)) +
      geom_tile() +
      scale_fill_brewer(palette = 'Spectral', direction = -1)
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/chrom_counts_by_mod.pdf'), p_chrom_counts, device = 'pdf', height = 4, width = 8)


# plot Domain comparison between MCs for chr19 region
if (cell_line == 'h1299' & treatment == 'dacsb')
{
	p = plotDomainComparison(mcs = mcs_te, genes_m_filt_ds, mat_ds, mcs1 = c(1:11, 15:20, 22:23), mcs2 = c(12:14, 21), 
									tes_tidy, genes_tidy, startp = 31500000, endp = 46500000)
	ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/domain_mc_count_comparison.pdf'), p, device = 'pdf', width = 10, height = 5)								
} else { # MT metacell 
	diff_expression = diff_expr(mcs_te, filterGeneMatrix(genes_m, genes_tidy, exclude_mt = FALSE), mcs1 = 13, mcs2 = c(1:12, 14:25), reg = 5, min_max_umi = 25) 

	p = diff_expression %>% 
		ggplot(aes(tot2, tot1, col = grepl('MT-', gene, fixed = TRUE), label = ifelse(enr >= 3.75, gene, ''))) + # labels top 10 non-MT most enriched genes
			geom_point_rast(size = 2) + 
			scale_y_log10() + 
			scale_x_log10() + 
			geom_text_repel() +
			theme(legend.position = 'none') + 
			scale_color_manual(values = c('black', 'red'))
			
	ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/apoptotic_MC_scatterplot.pdf'), p, 
		device = 'pdf', width = 5, height = 5)	
}

#######################################
# DIFF EXPRESSION COMPARISON
#######################################

res = plotDiffExpr(mcs1 = c(3,4,12,13), mcs2 = c(15:18, 6))
p = res[[1]]
diff_genes = res[[2]]
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/meta_comparison.pdf'), p, device = 'pdf')


####################
# GET CELL CYCLE GENES AND META
######################
cc_genes = getCCGenes(mat_id = 'genes_m_filt', gene_anchors = anchors, gene_anti = gene_anti, cor_threshold = cor_threshold, 
                     downsamp_n = downsamp_n, ignore_cells = ignore_cells)

pdf(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/cc_gene_selection_', treatment, '.pdf'), width = 20, height = 20)
cc_genes_final = filterCCGenes(mat_genes = genes_m_filt, cc_genes_all = cc_genes, cc_refs = refs, nclust = 20, downsamp_n = downsamp_n)
dev.off()
cc_gset = tgGeneSets(cc_genes_final, "cell cycle")
scdb_add_gset("cc_genes_final", cc_gset)

runMeta(mat_id = 'genes_m_filt', gset_id = 'cc_genes_final')

pdf(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/cc_mods_', treatment, '.pdf'), width = 15, height = 15)
cc_mods = filterCCGenes(mat_genes = genes_m_filt, cc_genes_all = names(cc_genes_final), cc_refs = refs, nclust = 8, filter_mods = FALSE, downsamp_n = downsamp_n)
dev.off()	

fwrite(list(cc_mods), paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/cc_mods_', treatment, '.txt', sep = '\t'))
#################			  

s_genes = names(tail(sort(mc_cc@mc_fp[,mc_s]), 10))
m_genes = names(tail(sort(mc_cc@mc_fp[,mc_m]), 10))

s_score = colSums(genes_m_filt_ds[s_genes,])
m_score = colSums(genes_m_filt_ds[m_genes,])
cell_size = calcCellSize(genes_tidy)
te_load = calcTELoad(tes_tidy, genes_tidy)

plot2DProjectionsWrapper(mc2d, mc_cc, ct_te, te_load, cell_size, s_score, m_score, s_genes)

# plot Cor Gene TE mods and load
shared_barcodes = intersect(colnames(genes_m_filt_ds), colnames(mat_ds))
pdf(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/gene_te_load_cor_', treatment, '.pdf'), width = 6.5, height = 5)
plotCorGeneTE(genes_m_filt_ds, te_load[colnames(genes_m_filt_ds)], top_n = 15)
dev.off()

pdf(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/gene_te_mod_5_cor_', treatment, '.pdf'), width = 6.5, height = 5)
plotCorGeneTE(log2(1 + 7 * genes_m_filt_ds[, shared_barcodes]), log2(1 + 7 * colSums(mat_ds[names(ct_te[ct_te == 5]), shared_barcodes])), top_n = 15)
dev.off()

pdf(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/gene_te_mod_6_cor_', treatment, '.pdf'), width = 6.5, height = 5)
plotCorGeneTE(log2(1 + 7 * genes_m_filt_ds[, shared_barcodes]), log2(1 + 7 * colSums(mat_ds[names(ct_te[ct_te == 6]), shared_barcodes])), top_n = 15)
dev.off()



res = plotModTOR(te_mods = ct_te, te_full_length)
p = res[[1]]
mod_ordering = res[[2]]
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/tor_boxplot.pdf'), p, device = 'pdf', width = 4, height = 10, units = 'cm')

p = plotTEmodMC(tes_ds = mat_ds, te_mods = ct_te, mc_object = mc_cc, mc_ordering = mc_ordering, mod_ordering, cell_size, s_score, m_score, 
	smoothed = T, mcs1 = 1:5, mcs2 = 6:10)
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/main/biological_process/', cell_line, '/te_mods_against_cc_metas_', treatment, '.pdf'), p, device = 'pdf', width = 8, height = 30, units = 'cm')


