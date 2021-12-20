Introduction
============

This file represents sampling results for the 2014 and 2015 ANIMIDA III cruises. 

This file is stored in Excel format. This file duplicates results also stored as Esri file geodatabase and CSV format. This file is for users who prefer working with Excel.

Each worksheet represents a table or view from the project database. The DataValues worksheet contains the results, while the other worksheets provide additional metadata about the results. See "ERD and Data Dictionary.pdf" for a description of each column in each worksheet.

Tips for viewing in Excel
=========================

Changing column width
---------------------

You may want to widen columns to view all contents of a cell. For example, date columns might show up as ###### if the column is not wide enough.

Try these key strokes for widening columns.

CTRL+a to select all cells.
ALT, h, o, i to autofit column contents.
ALT, h, o, w to enter your own column width.

Looking up values
-----------------

The DataValues worksheet is intended to include enough descriptive columns to use the data for most cases. If more information about a particular item, such as collection method, is required, then you can use the MethodID column in the DataValues worksheet to map to methods in the Methods worksheet.

Using methods as an example, if you want to append method description as a new column in the DataValues worksheet, use an Excel function called VLOOKUP.

1. Insert a column after MethodID in the DataValues worksheet. Let's assume MethodID is column S, so the new column would be column T. Call the new column MethodDescription. 

2. Enter this formula in cell T2:
=VLOOKUP(S2,Methods!$A$2:$B$37,2,FALSE)

3. Copy cell T2 down the rest of the column.

In a nutshell, the VLOOKUP formula uses the MethodID in cell S2 to find a matching method from the table on the Methods worksheet defined by $A$2:$B$37. We use the dollar sign so that Excel doesn't unintentionally update the row numbers defining the table. Otherwise, the table would copy down as A3:B38, A4:B39, and so on. The next function argument, 2, tells VLOOKUP to use the value in the second column of the table as the value for cell T2. Finally, FALSE means use an exact match, not an approximate match.

Filtering values
----------------

You can filter data in Excel from the Home tab, in the Editing panel (typically on the far right). In that panel, click Sort & Filter, and then click Filter. You'll then see drop down arrows next to column names. As an example, suppose you want to only see arsenic data:

1. While viewing the DataValues worksheet, enable filtering (as described above).
2. Click the drop down arrow in the VariableName column.
3. Type arsenic in the Search box, and then press ENTER or click OK.

Issues when adding data to ArcGIS
=================================

If you decide to start with the XLSX file and then add it to ArcGIS, keep the following considerations in mind.

Some site names are interpreted by ArcGIS as coordinate values, so it adds SiteName_x and SiteName_y fields and populates them with the coordinate values that it thinks are represent by the site name. These two fields are meaningless and can be deleted.

Synonyms, Method, LabMethodDescription, and SourceDescription columns are longer than 255 characters and may not display correctly in ArcMap. We recommend creating a copy of the data file and removing those columns from the copy, and then adding the copy to ArcMap.  You will still have those columns in the original data file if you need to look them up.

