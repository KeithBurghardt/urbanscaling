# Urban Scaling Code
Code for "Analyzing urban scaling laws in the United States over 115 years"

Citation:
Keith Burghardt, Johannes H. Uhl, Kristina Lerman, Stefan Leyk. (2022). Analyzing urban scaling laws in the United States over 115 years. Environment and Planning B: Urban Analytics and City Science (in press)



To reproduce results, after extracting all data referenced in the manuscript:
- Reproduce CIUPs from urban_scaling_create_urban_boundaries.py where auxilary data is from https://github.com/johannesuhl/USRoadNetworkEvolution/tree/main/auxiliary_data
- Extract urban properties from files shown in https://github.com/johannesuhl/USRoadNetworkEvolution/tree/main/CBSA_statistics
- Compile files into one larger file and add indoor area (BUI) statistics with data_analysis/FullStats.py
- Extract CIUP statistics and incorporate BUFA statistics with extract_CIUP_BUFA.py
- Finally extract all scaling statistics using the Global Human Settlement Layer with WorldScaling.py. This file namely utilizes statistics extracted by Boeing, 2020*
- Create generalized urban areas for a range of parameters, in order to assess the sensitivity of scaling laws to the choice of urban boundary modeling parameters (urban_boundary_modeling_sensitivity.py)
- Fitting scaling laws are in data_analysis/ScalingStats.py, data_analysis/ScalingStats_robust, and nl_fit.r
Data:
- Auxilary data: https://github.com/johannesuhl/USRoadNetworkEvolution/tree/main/auxiliary_data
- BUA: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/J6CYUJ
- BUI: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1WB9E4 
- BUFA: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/HXQWNJ
- Road lengths: https://figshare.com/projects/USRoadNetworkEvolution/137044


* Boeing, G. (2022), Street Network Models and Indicators for Every Urban Area in the World. Geogr Anal, 54: 519-535. https://doi.org/10.1111/gean.12281
