Introduction
============

These files represent sampling results for the 2014 and 2015 ANIMIDA III cruises. 

These files are stored in comma separated values (CSV) format, a text format which can be read by many software applications.  These files duplicate results also stored as Esri file geodatabase. These files are for users who prefer working with CSV or who do not have ArcGIS.

DataValues.csv contains the results, while the other files provide additional metadata about the results. See "ERD and Data Dictionary.pdf" for a description of each column in each file.

Issues when adding data to ArcGIS
=================================

If you decide to start with these CSV files and then add them to ArcGIS, keep the following considerations in mind.

Some site names are interpreted by ArcGIS as coordinate values, so it adds SiteName_x and SiteName_y fields and populates them with the coordinate values that it thinks are represent by the site name. These two fields are meaningless and can be deleted.

Synonyms, Method, LabMethodDescription, and SourceDescription columns are longer than 255 characters and may not display correctly in ArcMap. We recommend creating a copy of the data file and removing those columns from the copy, and then adding the copy to ArcMap.  You will still have those columns in the original data file if you need to look them up.

