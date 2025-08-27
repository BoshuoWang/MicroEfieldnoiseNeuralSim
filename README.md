# MicroEfieldnoiseNeuralSim
Code and data for simulation of neural response to E-field with microscopic noise. Cite as: Wang B, Dataset: MicroEfieldnoiseNeuralSim, https:doi.org/10.5281/zenodo.16949817

Associated publications:
[1]  Wang B, Hussain M A, Worbs T, Thielscher A, Grill W M and Peterchev A V 2025 Single pointwise samples of electric field on a neuron model cannot predict activation threshold by brain stimulation 2025.07.31.667968. https://doi.org/10.1101/2025.07.31.667968 
[2]  Wang B, Aberra A S, Grill W M and Peterchev A V 2018 Modified cable equation incorporating transverse polarization of neuronal membranes for accurate coupling of electric fields J. Neural Eng. 15 026003. https://doi.org/10.1088/1741-2552/aa8b7c

Field distribution folder: files that analyze and generate E-field distributions.
	distribution_analysis.m: analyze probablity distribution of E-field ratio from Weise et al. 2025 (from digitized figure, data in Weise2025.csv). Creates pSE_distribution.mat
	fit_dist: generates 10 E-field profiles with probablity distribution that best match Weise et al. 2025. Creates E-field_profiles.mat in root folder.
	plot_dist: Visualization and further analysis of the 10 E-field profiles.
	process_coaxial_field: Processes E-fields from spherical electrode simulation in Wang and Aberra 2025. Reads homogeneous_y118z81.txt and micro3D_y118z81.txt and creates PointSource.mat in root folder.

Neural simulations: adapted from [2]

main_ES.m: main function running simulations.
		
		mod_prmtr: structure specifying model parameter:
			model_name                Model: 'UF', 'PE'
			model_axon                Axon:  'HH', 'MRG'
			model_size                Axon size: '3um', '0.3um'
			id                        Parameter ID of test case: 1-700 for UF and 1-280 for PE
		
		out_ctrl: structure specifying output control:
			if_save_data   	Whether to save results in a .mat file for each simulation (logical or 0/1)
			if_write_log  	Whether to write threshold finding process in a .txt log (logical or 0/1)
			if_plot        	Whether to plot threshold finding process (logical or 0/1)

specify_model_ES_HH.m & specify_model_ES_RMG.m: defining models and simulation setup

simulate_cable_HH.m & simulate_cable_RMG.m: cable solver

Folder "Shared functions and data": auxillary functions for running simulations

compile_all_datasets.m: collecting all individual simulation data into datasets

plot_results.m: visualize results

Other subfolders: simulation results and logs, and compiled datasets
