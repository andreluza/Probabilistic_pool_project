
# ------------------------
#   Description of files
# ------------------------

hypervolume_local_null_deltaFD* : null deltaFD for local communities, for three different data sets (varying their sampling effort)
hypervolume_obs_deltaFD*: observed deltaFD for local communities, for three different effort datasets
hypervolumePOOL*: functional hypervolume of the local probabilistic species pool

mean_hypervolume_pool_*: results of hypervolume for 100 random samples of composition, per community
- Observation: I had to run this hypervolume analysis in 6 different steps, as this analysis is very
		demanding in terms of computation memory and time.

FD_pool_final: global functional diversity of non-volant small mammal communities, at cells of 2x2 degrees
- Observation: FD (hypervolume) was calculated for cells having >=5 species
SR_pool_final: global species richness of non-volant small mammal communities, at cells of 2x2 degrees

# ------------------------
#   Description of the folder
# ------------------------

figures: folder with raw figures from R (in format 'pdf'), SVG files with edited figures,
	and png format.