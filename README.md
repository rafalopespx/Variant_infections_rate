README
================

# Combining genomic data and infection estimates to characterize the complex dynamics of SARS-CoV-2 Omicron variants in the United States

Rafael Lopes1,#, Kien Pham1, Fayette Klaassen2, Melanie H. Chitwood1,
Anne M. Hahn1, Seth Redmond1, Nicole A. Swartwood2, Joshua A. Salomon3,
Nicolas A. Menzies2, Ted Cohen1,*,#, Nathan D. Grubaugh1,4,*,#

1 Department of Epidemiology of Microbial Diseases and Public Health
Modeling Unit, Yale School of Public Health, New Haven, CT, USA 2
Department of Global Health and Population, Harvard T. H. Chan School of
Public Health, Boston, MA, USA 3 Department of Health Policy, Stanford
University School of Medicine, Stanford, CA, USA 4 Department of Ecology
and Evolutionary Biology, Yale University, New Haven, CT, USA \*
Co-senior authors \# Corresponding authors: <rafael.lopes@yale.edu>,
<theodore.cohen@yale.edu>, <nathan.grubaugh@yale.edu>

## Data Availability

The findings of this study are based on metadata associated with
3,103,250 sequences available on GISAID from September 1st, 2021 up to
April 22nd, 2023, and accessible at
<https://doi.org/10.55876/gis8.231023hd> (GISAID Identifier:
EPI_SET_231023hd). All genome sequences and associated metadata in this
dataset are published in GISAID’s EpiCoV database. To view the
contributors of each individual sequence with details such as accession
number, Virus name, Collection date, Originating Lab and Submitting Lab
and the list of Authors, visit https//doi.org/10.55876/gis8.231023hd

## Pipeline Running order

All the codes to reproduce the paper analysis are at Scripts/ folder. At
2023-10-30, the pipeline running order is:

- **manuscript_figures.R**, make all the manuscript figures.
- **manuscript_table.R**, make all the manuscript tables.
- **01_metadata_cleaning.R**, clean metadata from GISAID and set variant
  categories, count and frequencies.
- (Optional) **02_plot_metadata.R**, plot figures with variant counts
  and frequencies.
- **03_infections_per_variant_estimates.R**, estimates infections per
  variants.
- (Optional) **04_plot_infections_estimates.R**, generate plots of
  infections.
- **05_variant_rt_estimates_daily.R**, estimate Rt per variant per
  state.
- **06_rt_ratios.R**, estimates rt ratio per pairs of variants.
- **07_attack_rate_svi.R**, attack rate vs SVI correlation and figure4
  of the manuscript.
