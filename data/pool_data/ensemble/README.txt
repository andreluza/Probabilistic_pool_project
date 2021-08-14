# README for chelsa worldclim cru data

This data are the results of a species distribution modeling exercise comparing the performance of temporal vs. spatial aggregations of different climate datasets.

# The scripts that produced the data are located here:
https://gitlab.wsl.ch/karger/ch_wc_cru

# The data is sorted into 4 different folders.

- metadata
	- species_names

		[group].csv

		A file containing the names and ids for the species of each group ‘amphibians’, ‘reptiles’, or ‘mammals’. The file has 3 columns. ‘id_ch_wc_cru’ gives the id that is used in the raster stacks or the data frames. ‘id_iucn’ is the respective IUCN id. ‘species_name’ is the name of the species based on the IUCN taxonomy.

- statistics: 
	- [group]
		
		df_[metric]_[group].RDS	
		
		A data frame containing the metric ‘tss’ or ‘auc’ for the group ‘amphibians’, ‘reptiles’, or ‘mammals’.

		Hint: There might be more entries in this data table then in the stacks. If a sdm would not work or crash, there will be still one entry. Better to use the stacks as reference for the species that have been successfully modeled.

- boxplots:
		[group]_[sdm]_[metric].pdf

		A pdf showing a boxplot comparing the differnent climate models [ch, , cru] and their performance based on a metric (tss, auc) for different sdms (rf, glm, gam).

	
- stacks:
	- regional:
		- [group]
			[metric].stk[climate]_[aggregation]_[sdm].RDS
	
			R data files containing raster stacks of a given metric (auc, 	pres_abs, prob, tss) for all species of a given group. 

		
		
	
