# ms-peptide-kd
Code and data for analyzing results from PEDAL (Parallel Equilibrium Dialysis for Affinity Learning) experiments by Zhao and Grigoryan, et al.

# Table of contents
`ms_allData3.xlsx` -- mass spectrometry data (including peptide sequences and peak intencities) obtained as part of our study of PDZ-peptide interactions using PEDAL. The file contains five independent datasets (different peptide libraries or different experimental settings). This is the default dataset read by `analyzeMS.ms`.

`Kd_FP_sum.xlsx` -- peptide affinity validation data obtained using using Fluorescence Polarization (FP) in low throughput.

`analyzeMS.ms` -- the main script for analyzing PEDAL data. It reads a data file (by default, `ms_allData3.xlsx`) and output a CSV file containing Kd, estimated error, observed apparent error, number of occurrences, alpha and cross-correlation value of each unique peptide. The script also reads in datafile `Kd_FP_sum.xlsx` and uses `plotFpMsKd.m` to plot the affinities from PEDAL assay against affinities from FP assay for a group of peptides
	
`best.txt` -- the best estimation of Kd (i.e., the one with the lowest predicted error) for each unique sequence identified from all four peptide libraries.   

`cluster_seqs.m` -- this script reads `best.txt` and perform a t-SNE analysis of all the peptides included.
