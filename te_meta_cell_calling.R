library("devtools")
load_all("/home/davidbr/src/metacell/")
scdb_init("/home/davidbr/tmp/", force_reinit = TRUE)
scfigs_init("/home/davidbr/tmp/figs/")
require(tgconfig)
tgconfig::override_params("~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/metacell/tes/h1299_dacsb.yaml", package="metacell")
 
#####################
# FUNCTION
#####################

metaPipeline = function(mat_id = 'mat_id', T_vm = 0.1, T_tot = 25, T_top3 = 2, n_resamp = 500)
{
	mat_id = mat_id
	graph_id = paste0(mat_id, '_graph')
	gset_id = paste0(mat_id, '_gset')
	gstat_id = paste0(mat_id, '_gstat')
	coc_id = paste0(mat_id, '_coc')
	mc_id = paste0(mat_id, '_mc')
	mc2d_id = paste0(mat_id, '_2d')
	
	#####################
	# MARKER SELECTION
	#####################
	mcell_plot_umis_per_cell(mat_id = mat_id, min_umis_cutoff = NA, bin_for_cutoff = 50)
	mcell_add_gene_stat(mat_id = mat_id, gstat_id = gstat_id, force = TRUE)
	mcell_gset_filter_varmean(gstat_id = gstat_id, gset_id = gset_id, T_vm = T_vm, force_new = T)
	mcell_gset_filter_cov(gstat_id = gstat_id, gset_id = gset_id, T_tot = T_tot, T_top3 = T_top3)
	mcell_plot_gstats(gstat_id = gstat_id, gset_id = gset_id)
	print(paste0('GSET size: ', length(scdb_gset(gset_id)@gene_set)))
	mcell_gset_split_by_dsmat(gset_id = gset_id, mat_id = mat_id, K = 8)
	#mcell_plot_gset_cor_mats(gset_id = gset_id, scmat_id = mat_id, downsamp = T)
	
	#####################
	# CREATE GRAPH
	#####################						  
	mcell_add_cgraph_from_mat_bknn(graph_id=graph_id, mat_id=mat_id, gset_id = gset_id, K = 100, dsamp = TRUE)

	#####################
	# CREATE METACELLS
	#####################
	mcell_coclust_from_graph_resamp(coc_id = coc_id, graph_id = graph_id, min_mc_size = 25, p_resamp = 0.75, n_resamp = n_resamp) 

	####################
	# METACELL COVER
	####################
	mcell_mc_from_coclust_balanced(mc_id = mc_id, coc_id = coc_id, mat_id = mat_id, K = 30, min_mc_size= 25, alpha = 2) # 20
	
	####################
	# EXPORT CLUST_FP
	####################
	mcell_mc_export_tab(mc_id = mc_id, gstat_id = gstat_id, mat_id = mat_id, T_gene_tot = 25, T_fold = 1, metadata_fields = NULL) 

	####################
	# PLOTTING
	####################
	mcell_gset_from_mc_markers(gset_id = 'markers', mc_id = mc_id)
	mcell_mc_plot_marks(mc_id = mc_id, gset_id = 'markers', mat_id = mat_id, plot_cells = TRUE)
	mcell_mc2d_force_knn(mc2d_id = mc2d_id, mc_id = mc_id, graph_id = graph_id)
	mcell_mc2d_plot(mc2d_id = mc2d_id)
	mcell_mc_plot_subheats(mc_id = mc_id, mat_id = mat_id)
	
	#SAVE as RDS
	# mc = scdb_mc(mc_id)
	# mat_all = scdb_mat(mat_id)@mat
	# mat_ds = gset_get_feat_mat(gset_id = gset_id, mat_id = mat_id, downsamp = TRUE)
	# mat = gset_get_feat_mat(gset_id = gset_id, mat_id = mat_id, downsamp = FALSE)
	
	# saveRDS(mc, 'metacell_object.RDS')
	# saveRDS(mat_all, 'umi_matrix_all.RDS')
	# saveRDS(mat_ds, 'umi_matrix_marker_ds.RDS')
	# saveRDS(mat, 'umi_matrix_marker_raw.RDS')
}


set_param(param = 'scm_n_downsamp_gstat', value = NULL, package = 'metacell')
mcell_import_scmat_10x_custom('hct116_tes_dacsb', "~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/counts/te_tss_bin_counts_sm.RDS", force = TRUE)
mcell_import_scmat_10x_custom('h1299_tes_dacsb', "~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/counts/te_tss_bin_counts_sm.RDS", force = TRUE)

mat_id = 'hct116_tes_dacsb'
T_vm = 0.1
T_tot = 9
T_top3 = 2
n_resamp = 500

metaPipeline(mat_id = 'hct116_genes_dacsb', T_vm = 0.1, T_tot = 9, T_top3 = 2, n_resamp = 150) FOR TEs





		  




