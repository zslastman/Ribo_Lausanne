picardqc <- Sys.glob('pipeline/qc/data/*/*_picard_qc.txt') %>% setNames(.,basename(.)) %>%
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

select <- dplyr::select


picardqc%>%select(file,PCT_RIBOSOMAL_BASES,PCT_UTR_BASES,PCT_CODING_BASES,PCT_INTERGENIC_BASES)%>%as.data.frame



picardqc%>%
	mutate(city = file%>%str_extract('(?<=ctrl\\d)[L|B]'))%>%
	# filter(!is.na(city))%>%
	arrange(is.na(city),city,file)%>%
	select(file,PCT_RIBOSOMAL_BASES,PCT_UTR_BASES,PCT_CODING_BASES,PCT_INTERGENIC_BASES,matches('BIAS'))%>%as.data.frame


