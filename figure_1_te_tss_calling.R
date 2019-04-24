source('/net/mraid14/export/data/users/davidbr/general/scripts/utilities.R')

#########################
# FUNCTIONS
#########################
plotFamUMIsConditionComparison = function(condition1, condition2, tss_union)
{
	condition1_tss = left_join(tss_union, condition1)
	condition2_tss = left_join(tss_union, condition2)
	
	fam_umi_counts_condition = full_join(condition1_tss %>% group_by(repname) %>% summarise(umis = sum(as.numeric(n))),
					 condition2_tss %>% group_by(repname) %>% summarise(umis = sum(as.numeric(n))), by = 'repname') %>%
		replace(., is.na(.), 0) %>%
		mutate(fc = umis.y / umis.x)
		
	plot(fam_umi_counts_condition$umis.x + 1, fam_umi_counts_condition$umis.y + 1, log = 'xy', pch = 20)
	return(fam_umi_counts_condition)
}

compareTSSbetweenConditions = function(condition1, condition2)
{
	shared_families = intersect(condition1$repname, condition2$repname)

	con1 = condition1 %>% 
		group_by(repname, bin_id) %>%
		summarise(umis = sum(n), repclass = unique(repclass)) %>%
		ungroup
		
	con2= condition2 %>% 
		group_by(repname, bin_id) %>%
		summarise(umis = sum(n), repclass = unique(repclass)) %>%
		ungroup
		
	tss_bin_stats = full_join(con1, con2, by = c('repname', 'bin_id')) %>%
		filter(repname %in% shared_families) %>%
		rowwise %>%
		mutate(repclass = unique(na.omit(c(repclass.x, repclass.y)))) %>%
		select(-repclass.x, -repclass.y) %>%
		group_by(repname, bin_id) %>%
				mutate(shared = sum(is.na(umis.x), is.na(umis.y)) == 0,
					   max_umis = max(umis.x, umis.y, na.rm = TRUE)) %>%
				ungroup %>% 
		group_by(repclass) %>%
			mutate(expr_bin = ntile(max_umis, n = 4)) %>%
			ungroup %>%
		group_by(repclass, expr_bin) %>%
			mutate(med_expr = median(max_umis)) %>%
			ungroup %>%
		arrange(repname, bin_id)
		
	tss_bin_stats = tss_bin_stats %>% mutate(shared = ifelse(is.na(umis.x), 'h1299', shared), shared = ifelse(is.na(umis.y), 'hct116', shared), shared = ifelse(shared == TRUE, 'both', shared))
	
	p = count(tss_bin_stats, repclass, expr_bin, med_expr, shared) %>% 
			ggplot(aes(fill = repclass, y = n, x = as.character(expr_bin), alpha = shared)) + 
				geom_bar(stat = 'identity', position = 'fill') + 
				scale_alpha_discrete(range = c(1, 0.33)) + 
				scale_fill_brewer(palette = 'Set1') +
				facet_wrap(~ repclass, ncol = 4)
			
	return(p)
}

getTSSBinUnion = function(tss_condition1, tss_condition2)
{
	tss_combined = full_join(tss_condition1 %>% select(repname, bin_id) %>% distinct, 
							 tss_condition2 %>% select(repname, bin_id) %>% distinct, by = c('repname', 'bin_id'))
	return(tss_combined)
}

calcReadHamming <- function(bam_file, tss_bins, fam = 'LTR12C', uniquely = TRUE, max_reads_per_locus = 100)
{
  library(GenomicAlignments)
  
  valid_barcodes <- tss_bins$CB %>% unique
  which_regions <-  tss_bins %>% filter(repname == fam)  %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% unique
  
  reads = readGAlignments(bam_file,
                          index = gsub('.bam$', '.bai', fam_file),
                          use.names=TRUE, 
                          param = ScanBamParam(which = which_regions,
                                             what = 'seq',
                                             tag = c("NH", "CB"),
                                             flag = scanBamFlag(isPaired=TRUE, isProperPair=TRUE, isFirstMateRead=TRUE)))
                                             
  reads <- granges(reads, use.mcols = TRUE, use.names = TRUE)
  reads <- reads[!is.na(reads@elementMetadata$CB)]
  if (!is.null(valid_barcodes)) { reads <- reads[reads@elementMetadata$CB %in% valid_barcodes] }
  
  if (uniquely) { reads <- reads[reads@elementMetadata$NH == 1] }                                            
  
  reads_plus <- reads[strand(reads) == '+']
  end(reads_plus) <- start(reads_plus) + 1
  reads_plus_unique <- reads_plus[!duplicated(names(reads_plus))] # duplicated entries because first mate is considered (2nd mate is distinct)
  
  reads_plus_unique$row_id <- NA
  hits <- findOverlaps(reads_plus_unique, which_regions)
  reads_plus_unique[hits@from]$row_id <- which_regions[hits@to]$row_id
  
  reads_plus_unique <- reads_plus_unique[!is.na(reads_plus_unique$row_id)]
  
  # Select reads mapping to max covered bp per loci
  set.seed(19)
  reads_plus_unique <-
    reads_plus_unique %>% 
      as_tibble %>% 
      dplyr::count(row_id, start) %>% 
      group_by(row_id) %>% 
      top_n(1, n) %>% 
      group_by(row_id) %>%
      sample_n(1) %>% 
      ungroup %>% 
      left_join(., 
                as_tibble(reads_plus_unique))
  
  # Sample reads
  reads_plus_unique_ds <-
    reads_plus_unique %>% 
      group_by(row_id) %>% 
      #mutate(n = length(row_id)) %>% 
      #filter(n >= reads_per_locus) %>%
      sample_n(size = max_reads_per_locus, replace = TRUE) %>%
      ungroup %>%
      distinct
      
  reads_plus_unique_ds <-  
    reads_plus_unique_ds %>%
    as.data.frame %>%
    rownames_to_column %>%
    as_tibble 
  
  # Calc all vs all dist
  all_vs_all_dists <- 
    stringdist::stringdistmatrix(reads_plus_unique_ds$seq, reads_plus_unique_ds$seq, method = 'hamming', nthread = 50) 
    
  colnames(all_vs_all_dists) <- reads_plus_unique_ds$rowname
  rownames(all_vs_all_dists) <- reads_plus_unique_ds$rowname
  
  #all_vs_all_dists <- as(all_vs_all_dists, 'sparseMatrix')
  
  res <-
    all_vs_all_dists %>% 
    as.data.frame %>%
    rownames_to_column %>%
    as_tibble %>%
    gather('read_id2', 'dist', -rowname) %>%
    filter(!is.na(dist)) %>%
    rename(read_id1 = rowname)
  
  # Add locus info
  res2 <- left_join(res, tibble(read_id1 = reads_plus_unique_ds$rowname, row_id1 = reads_plus_unique_ds$row_id), 'read_id1')
  res3 <- left_join(res2, tibble(read_id2 = reads_plus_unique_ds$rowname, row_id2 = reads_plus_unique_ds$row_id), by = 'read_id2')
  res3 <- res3 %>% mutate(same_locus = row_id1 == row_id2)

  # Plotting
  m_dist <- res3 %>% select(read_id1, read_id2, dist) %>% tidyToSparse %>% as.matrix
  m_loci <- res3 %>% select(read_id1, read_id2, same_locus) %>% tidyToSparse %>% as.matrix
  
  d  <- tgs_dist(m_loci * 1)
  hc <- hclust(d, 'ward.D2')
  
  m_dist_ord <- m_dist[hc$order, hc$order]
  m_dist_ord[lower.tri(m_dist_ord)] = NA
  
  m_loci_ord <- m_loci[hc$order, hc$order]
  m_loci_ord[upper.tri(m_loci_ord)] = NA
  
  par(mar = rep(0.1, 4))
  image(m_dist_ord, 
        col = c('purple3', 'darkred', 'orange', 'white'),
        breaks = c(0, 0.9, 1.9, 3.9, 128),
        xaxt = 'n',
        yaxt = 'n')
        
  image(m_loci_ord, 
        col = c('white', 'black'),
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE)
  
  
  return(res3)
}

readPairMappingComp <- function(read_pairs, tes_full_length)
{
  te_coords <-
    te_full_length %>% 
    filter(repclass %in% c('LTR', 'SINE', 'DNA', 'LINE')) %>% 
    mutate(seqnames = substring(chrom, 4, 99)) %>% 
    as_granges() %>%
    filter(countOverlaps(.) == 1)
  
  first_mates <- first_mates <- GenomicAlignments::first(read_pairs) %>% as_granges %>% anchor_5p %>% mutate(width = 1)
  second_mates <- GenomicAlignments::second(read_pairs, real.strand = TRUE) %>% as_granges %>% anchor_3p %>% mutate(width = 1)
  
  hits_first <- 
    findOverlaps(first_mates, te_coords) %>% 
    as_tibble %>% 
    rename(mate_index = queryHits, te_index_first = subjectHits) %>% 
    mutate(repclass = te_coords[te_index_first]$repclass,
           repname = te_coords[te_index_first]$repname)
  hits_second <- findOverlaps(second_mates, te_coords) %>% as_tibble %>% rename(mate_index = queryHits, te_index_second = subjectHits)
  
  ggdata <-
    left_join(hits_first, 
              hits_second, 
              by = 'mate_index') %>% 
    mutate(same_te = te_index_first == te_index_second) %>% 
    replace_na(list(same_te = FALSE)) %>% 
    group_by(repname) %>%
    summarise(repclass = repclass[1],
              total_reads = length(mate_index),
              same_te_reads = sum(same_te),
              diff_te_reads = sum(!same_te),
              prop_same = same_te_reads / (same_te_reads + diff_te_reads))
    count(repclass, repname, same_te)
    
  p <-
  ggplot(ggdata, aes(x = repclass, y = n, fill = repclass, alpha = same_te)) +
    geom_bar(stat = 'identity') +
    scale_fill_brewer(palette = 'Set1')

  return(p)
}


plotTSSComparison <- function(untreated = h1299_tes_tss_dmso, 
                              treated = h1299_tes_tss_dacsb,
                              min_umis = 150)
{
  a <- 
    untreated %>%
    group_by(repname, bin_id) %>%
      mutate(bin_umis = sum(n)) %>%
    group_by(repname) %>% 
      summarise(bins = list(1:nbins_fam[1]), tss_bin_dmso = list(bins %in% bin_id), max_umis = max(bin_umis[tss_bin_dmso])) %>% 
    unnest %>%
    filter(max_umis >= min_umis)
    
  b <- 
    h1299_tes_tss_dacsb %>% 
    group_by(repname, bin_id) %>%
      mutate(bin_umis = sum(n)) %>%
    group_by(repname) %>% 
      summarise(bins = list(1:nbins_fam[1]), tss_bin_dacsb = list(bins %in% bin_id), max_umis = max(bin_umis[tss_bin_dacsb])) %>% 
    unnest %>%
    filter(max_umis >= min_umis)
    
  ggdata <-
    inner_join(a,b, c('repname', 'bins')) %>% 
    mutate(type = ifelse(tss_bin_dmso + tss_bin_dacsb == 2, 'shared', 'none'),
           type = ifelse(tss_bin_dmso & !tss_bin_dacsb, 'dmso', type),
           type = ifelse(!tss_bin_dmso & tss_bin_dacsb, 'dacsb', type)) %>%
    left_join(., untreated %>% select(repname, repclass) %>% distinct)

  p1 <-
    ggdata %>% 
    ggplot(aes(bins, repname, fill = type)) + 
    geom_raster() + 
    scale_fill_manual(values = c('red', 'blue', 'lightgrey', 'darkgreen')) +
    facet_wrap(~repclass, scale = 'free')
    
  p2 <- 
    ggdata %>% 
    filter(type != 'none') %>% 
    count(repclass, type) %>% 
    ggplot(aes(x = repclass, y = n, fill = type)) + 
    geom_bar(stat = 'identity', position = 'fill') +
    coord_flip()
    
  p <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(5, 2), labels = 'auto')
  return(p1)  
}

#########################
# LOAD
#########################
library(Gviz)
gene_anno = import.gff('~/davidbr/general/data/gencode/gencode.v28.annotation.gtf')
te_full_length = readRDS('~/davidbr/general/data/repeatmasker/hg38/full_length_te_coords.RDS')
read_pairs <- GenomicAlignments::readGAlignmentPairs('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/aligned/chunked_bams/possorted_deduplicated_chr1.bam')


h1299_genes_dmso = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/counts/gene_counts_tidy.RDS')
hct116_genes_dmso = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dmso/counts/gene_counts_tidy.RDS')
h1299_genes_dacsb = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/counts/gene_counts_tidy.RDS')
hct116_genes_dacsb = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/counts/gene_counts_tidy.RDS')


h1299_tes_dmso = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/counts/te_bin_counts.RDS')
h1299_tes_tss_dmso = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/counts/te_tss_bin_counts_tidy.RDS')

h1299_tes_dacsb = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/counts/te_bin_counts.RDS')
h1299_tes_tss_dacsb = readRDS('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/counts/te_tss_bin_counts_tidy.RDS')

hct116_tes_dmso = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dmso/counts/te_bin_counts.RDS')
hct116_tes_tss_dmso = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dmso/counts/te_tss_bin_counts_tidy.RDS')

hct116_tes_dacsb = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/counts/te_bin_counts.RDS')
hct116_tes_tss_dacsb = readRDS('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/counts/te_tss_bin_counts_tidy.RDS')
#########################

tss_union = getTSSBinUnion(hct116_tes_tss_dmso, hct116_tes_tss_dacsb)
test = plotFamUMIsConditionComparison(hct116_tes_dmso, hct116_tes_dacsb, tss_union)

p = compareTSSbetweenConditions(h1299_tes_tss_dacsb, hct116_tes_tss_dacsb)
ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/te_tss_calling/tss_comparison.pdf', p, device = 'pdf')

png('~/davidbr/proj/epitherapy/manuscript/figures/main/te_tss_calling/mate1_hamming_LTR12C_mapping.png', width = 5000, height = 5000)
calcReadHamming(bam_file = "~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/aligned/chunked_bams/possorted_deduplicated_chr1.bam", 
                tss_bins = h1299_tes_tss_dacsb, 
                fam = 'LTR12C', 
                uniquely = TRUE, 
                max_reads_per_locus = 100)
dev.off()

p_mate_mapping <- readPairMappingComp(read_pairs, tes_full_length)
ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/te_tss_calling/mate_mapping_same_te_locus_h1299_chr1_mates.pdf', p_mate_mapping, device = 'pdf')

p_tss_comp <- plotTSSComparison(h1299_tes_tss_dmso, h1299_tes_tss_dacsb, min_umis = 150)
ggsave('~/davidbr/proj/epitherapy/manuscript/figures/main/te_tss_calling/tss_comparison.pdf', p_tss_comp, device = 'pdf')