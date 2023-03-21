------------------------
This README Accompanies:
------------------------


--------------------
Corresponding Author:
---------------------

Conor McMeniman
cmcmeni1@jhu.edu

------------------------------
Matlab scripts in this repository:
------------------------------

Script to process corrected trajectory field from ivTrace "full.txt". This script outputs a "trace_data.mat" file with information such as landings, platform occupancy, experiment type, date

process_tractories.m 

Script to concatenate the "trace_data.mat" files from each experiment type (i.e. CO2, 1-human, 2-human) onto a single cell to be further analyzed:

concatenate_data_cell.m

Script to analyze landing and weather data from each experiment type:

Trajectories_landing_analysis.m

Script to load compound abundance matrix and run statistical analysis:

compound_abundance_analysis.m

------------------------------
Data files in this repository:
------------------------------

data_lab.mat contains a cell array with landing information for laboratory experiments. To be analyzed with "Trajectories_landing_analysis.m"


data_field_2020.mat contains a cell array with landing information for semi-field experiments carried out in 2020 (i.e. CO2, 1-human, 2-human). To be analyzed with "Trajectories_landing_analysis.m"


data_6-human.mat contains a cell array with landing information for semi-field experiments carried out in 2022 (i.e. 6-human). To be analyzed with "Trajectories_landing_analysis.m"


weather_2020.mat contains weather information (i.e. temperature, humidity, wind speed, wind direction, wind gust) during experimental time (8pm-4am) for experiments carried out in 2020 (i.e. CO2, 1-human, 2-human). To be analyzed with "Trajectories_landing_analysis.m"


weather_2022.mat contains weather information (i.e. temperature, humidity, wind speed, wind direction, wind gust) during experimental time (8pm-1am) for experiment carried out in 2026 (i.e. 6-human). To be analyzed with "Trajectories_landing_analysis.m"

compound_abundance.mat contains a table with the relative abundances of compounds measured for 6-human experiments








