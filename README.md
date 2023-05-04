README
================

## Pipeline runnig order status

At 2023-05-03, the pipeline running order is:

- 01_metadata_cleaning.R, to clean the metadata from GISAID and set
  variant categories, count and frequencies
- (Optional) 02_plot_metadata.R, to plot the figures with variant counts
  and frequencies
- 05a_infections_per_variant_estimates.R, script to stretch the time
  series to daily basis and estimates infections per variants
- 05b_plot_infections_estimates.R, script to generate the plot of
  infections over all states
- 06a_variant_rt_estimates_daily.R, script to estimate Rt per variant
  per state
- 07a_plot_rt_estimates_daily.R, script to plot the Rt estimates per
  variant per state
- 09_plot_infections_with_rt.R, script to plot Rt estimates and Rt,
  altogether

## To-do list

- [ ] Organize the repo to the scripts have a more logical order
- [ ] Calculate Rt ratios per states per variant
- [ ] Pairs of time of emergence to peak comparison ratio:
  - [ ] BA.1 vs. BA.2
  - [ ] BA.2 vs. XBB
  - [ ] BA.4 vs. BA.5
    - [ ] Their timing to peak with BA.2
      - [ ] BA.5 vs. BA.2
      - [ ] BA.4 vs. BA.2
