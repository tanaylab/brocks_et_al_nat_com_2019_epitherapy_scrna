cell_line <- 'h1299'
treatment <- 'dacsb'

liste <- list()
i <- 1
for (cell_line in c('h1299', 'hct116'))
{
for (treatment in c('dmso', 'dacsb'))
{

if (cell_line == 'h1299')
{
  if (treatment == 'dacsb')
  {
    mc_cc = get(load('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/metacell/genes/cell_cycle/mc.h1299_genes_m_filt_filt_cc_mc.Rda'))
    cycling_ids <- c(3,4,5,1,2)
  } else {
    mc_cc = get(load('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/metacell/genes/cell_cycle/mc.h1299_genes_m_filt_cc_mc.Rda'))
    cycling_ids <- c(1:7)
  }
}

if (cell_line == 'hct116')
{
	if (treatment == 'dacsb')
  {
    mc_cc = get(load('~/davidbr/proj/epitherapy/data/hct116/10x/dacsb/metacell/genes/cell_cycle/mc.hct116_genes_m_filt_cc_mc.Rda'))
    cycling_ids <- c(3,4,5,2,1,9,12)
  } else 
    mc_cc = get(load('~/davidbr/proj/epitherapy/data/h1299/10x/dmso/metacell/genes/cell_cycle/mc.h1299_genes_m_filt_cc_mc.Rda'))
    cycling_ids <- c(1:8)
}
mc_cc <- tibble(CB = names(mc_cc@mc), mc = mc_cc@mc)
tes <- readRDS(paste0('~/davidbr/proj/epitherapy/data/',cell_line,'/10x/',treatment,'/counts/te_tss_bin_counts_tidy.RDS'))

liste[[i]] <-
  inner_join(mc_cc, tes) %>%
  mutate(cycling = mc %in% cycling_ids) %>%
  group_by(repname, cycling) %>%
  summarise(umis = sum(n)) %>%
  ungroup %>%
  tidyr::complete(repname, cycling, fill = list(umis = 0)) %>%
  group_by(cycling) %>%
  mutate(umis_tot = sum(umis)) %>%
  ungroup %>%
  mutate(umis_norm = umis / umis_tot * min(umis_tot)) %>%
  select(-umis, -umis_tot) %>%
  left_join(., tes %>% select(repname, repclass) %>% distinct) %>%
  spread(cycling, umis_norm, sep = '_') %>%
  mutate(fc = log2((cycling_TRUE + 10)/ (cycling_FALSE + 10)),
         cell_line = cell_line,
         treatment = treatment)
  
  i <- i + 1
}
}

p <-   
  do.call(rbind, liste) %>% 
    mutate(treatment = factor(treatment, levels = c('dmso', 'dacsb'))) %>%
    ggplot(aes(x = repclass, y = fc, group = repclass, col = repclass)) + 
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter(alpha = .5, size = 0.5) +
      scale_color_brewer(palette = 'Set1') +
      theme(legend.position = 'none') +
      ylim(-2, 2) +
      facet_grid(treatment ~ cell_line)
      #geom_text_repel()
      
ggsave(paste0('~/davidbr/proj/epitherapy/manuscript/figures/supplements/cell_cycle_supps/cc_fam_umi_counts_fc.pdf'), p, device = 'pdf', height = 6, width = 6)