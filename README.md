# GrainBoundaries-DC
 Code used to create the paper on conduction through grain boundaries

 All files written and maintained by Eva Benford and Iddo Amit
 Current file version written on MacOS Sonoma 14.5 running on Apple Silicon M3

 Content:
    
   Alpha - The Python code to create the alpha matrix, and the results.
 	
   Distributions - including the Python file to calculate the charge density / doping distribution and the presented results
 	
   Look-up_netlists:
 		LUTs_builder - The Python file to create the LUTs
 		LUT_net_writer# - A Python file that write LTSpice netlists. The # represents the number of grains (0 is 10)
 		Analyser - The Python file that runs and analyses the results. Requires LTSpice package for Python as well as a local installation of LTSpice
 		Sample_LUTs - A library of LUTs, where n_D,1 = 0.9e16 1/cm^3 and n_D,2 ranges from 0.9 to 1.1e16 1/cm^3. Only a sample is provided due to space restrictions. All files can be rebuilt with LUTs_builder.py
 		Sample_Netlists - A library of the first 50 (out of 1,000) netlists used for the 10G:9B system presented in the paper
 	
   Potential_map - A repository of files required to create the voltage map presented in the paper
 
   Weibull_Statistics:
 		Weibull - A Python library with Weibull commands
 		overall_statistics - Calculation of statistics presented in the paper
 		Distributions - Results of charge density distributions for grain radius ranging from 100 to 500 nm
