picardqc <- Sys.glob('qc/data/*/*_picard_qc.txt') %>%
 	setNames(.,basename(.)) %>%
 	map(. %>%
 		readLines(8) %>%
 	.[7:8] %>%
 	str_split('\t') %>%
 	simplify2array %>%
 	t %>%
 	{set_colnames(.[2,,drop=F],.[1,])} %>%
 	as_data_frame) %>%
 	bind_rows(.id='file') %>%
 	as_tibble
