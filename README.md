# MicroEfieldnoiseNeuralSim
Code and data for simulation of neural response to E-field with microscopic noise

Associated publications:
[1]  Wang B, Hussain M A, Worbs T, Thielscher A, Grill W M and Peterchev A V 2025 Single pointwise samples of electric field on a neuron model cannot predict activation threshold by brain stimulation 2025.07.31.667968. https://doi.org/10.1101/2025.07.31.667968 
[2]  Wang B, Aberra A S, Grill W M and Peterchev A V 2018 Modified cable equation incorporating transverse polarization of neuronal membranes for accurate coupling of electric fields J. Neural Eng. 15 026003. https://doi.org/10.1088/1741-2552/aa8b7c

Field distribution folder: files that analyze and generate E-field distributions.
	distribution_analysis.m: analyze probablity distribution of E-field ratio from Weise et al. 2025 (from digitized figure, data in Weise2025.csv). Creates pSE_distribution.mat
	fit_dist: generates 10 E-field profiles with probablity distribution that best match Weise et al. 2025. Creates E-field_profiles.mat in root folder.
	plot_dist: Visualization and further analysis of the 10 E-field profiles.
	process_coaxial_field: Processes E-fields from spherical electrode simulation in Wang and Aberra 2025. Reads homogeneous_y118z81.txt and micro3D_y118z81.txt and creates PointSource.mat in root folder.


