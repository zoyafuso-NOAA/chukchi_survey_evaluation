Introduction
============

This accession contains physical, chemical, and biological data collected during research cruises for the Arctic Nearshore Impact Monitoring in Development Area (ANIMIDA) III Study.  The study occurred at 107 stations in the Beaufort Sea area in July and August of 2014 and 2015, and involved measurements and sampling of water, sediment, and biota conducted by investigators from several universities and research organizations.  The study also includes a spring river survey conducted in May, 2015.  The dataset includes more than 50,000 data values across 300+ variables, 500+ taxonomic names, 20+ collection methods, and 40+ lab analysis methods.  Seawater samples include variables such as salinity, total suspended solids, and dissolved elements.  Samples from the benthic zone include taxonomic counts and biomass. Sediment samples include variables such as hydrocarbons, metals, grain size distribution, and stable isotope ratios. All samples include numerical data value, date of sample collection, location in latitude and longitude coordinates, station identifier, sampling method, and the researcher responsible for the measurement. Furthermore, these data were harmonized according to a community driven set of controlled vocabularies for variable names, chemical speciations, units of measure, sample mediums, and sample types maintained by the Consortium of Universities for the Advancement of Hydrologic Science, Inc.  The data are accompanied by an entity relationship diagram and data dictionary explaining the data structure. CTD, ADCP, and thermosalinograph sensor data are also provided following NCEI's netCDF templates for such data. All data are provided under the Creative Commons Attribution 3.0 United States (CC BY 3.0 US) license which states that the data are available for all to use provided that proper credit is given.

Except for sensor data (CTD and TSG), these files are stored in Esri file geodatabase, comma separated values (CSV), and Excel format. Each format contains the same sample results as other formats. Multiple formats are provided so that the user can choose the format most compatible with the user's software.  See "ERD and Data Dictionary.pdf" for an entity relationship diagrams of the tables included with this product, and a data dictionary describing each table.  Original sample data files prior to loading into the common database structure are also provided. Original data files are typically structured for optimal viewing of the data by the original investigator, though they may not include the same standard data descriptors as described in the data dictionary.

Sensor data include results of conductivity temperature depth (CTD) data sondes and thermosalinographs (TSG). Sensor data are stored in netCDF format using NODC templates described at https://www.nodc.noaa.gov/data/formats/netcdf/v2.0/. 

For more information about the project, see:
http://arcticstudies.org/animida_iii/index.html

Citation
========

Except for sensor data, sampling results are stored as numerical data values in the DataValues table (or file, or worksheet, depending on which format you are using). Each row representing a data value will also have a number in the SourceID column.  This SourceID points to a record in the Sources table. The Sources table includes a Citation column that indicates how the data value should be cited.  In other words, every data value is linked to its proper citation. You'll find that a single citation works for all data values in a dataset, such as sediment metals data from a particular investigator.

If you need to cite the entire project database, please use:
Kasper, Jeremy, et al. (2015). Sampling Results from 2014 and 2015 field cruises for ANIMIDA III. 

License
=======

These data are licensed under the Creative Commons Attribution 3.0 United States (CC BY 3.0 US) license. Please provide attribution.
https://creativecommons.org/licenses/by/3.0/us/

Disclaimer
==========

We attempt to maintain the highest accuracy of content in our data. Any errors or omissions can be reported for investigation. We make no claims, promises, or guarantees about the absolute accuracy, completeness, or contents of this data and expressly disclaim liability for errors and omissions.

Contact
=======

Dr. Tim Whiteaker
The University of Texas at Austin
whiteaker@utexas.edu
512-471-0570

Version
=======

Version 1.0
2017-04-27