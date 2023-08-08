


create table  hauls_nbs 
  as select a.CRUISEJOIN, a.HAULJOIN, a.REGION, a.VESSEL, a.CRUISE,
  a.HAUL, a.HAUL_TYPE, a.PERFORMANCE, a.START_TIME, a.DURATION,
  a.DISTANCE_FISHED, a.NET_WIDTH, a.NET_MEASURED, a.NET_HEIGHT,
  a.STRATUM, a.START_LATITUDE, a.END_LATITUDE, a.START_LONGITUDE,
  a.END_LONGITUDE, a.STATIONID, a.GEAR_DEPTH, a.BOTTOM_DEPTH,
  a.BOTTOM_TYPE, a.SURFACE_TEMPERATURE, a.GEAR_TEMPERATURE,
  a.WIRE_LENGTH, a.GEAR, a.ACCESSORIES, a.SUBSAMPLE, a.AUDITJOIN,
  floor(a.cruise/100) year
  from racebase.haul a
  where ((cruise = 198202 and vessel = 21 and haul_type in (0,3) and (stratum in (0,81,70,71) or stratum is NULL) and start_latitude > 60)
      or (cruise = 198502 and vessel = 60 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 198808 and vessel = 21 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL) and start_latitude > 59)
      or (cruise = 199101 and vessel = 37 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 199101 and vessel = 78 and haul != 154 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 199102 and vessel = 37 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 199102 and vessel = 78 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 201002 and vessel = 89 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 201002 and vessel = 94 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 201002 and vessel = 162 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 201702 and vessel = 94 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 201702 and vessel = 162 and haul_type = 3 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 201801 and haul_type = 13 and (stratum in (0,81,70,71) or stratum is NULL))
      or (cruise = 201902 and haul_type = 3 and (stratum in (81,70,71)))
      or (cruise = 202102 and haul_type = 3 and (stratum in (81,70,71) )))
    and region = 'BS'
    and performance >= 0
    and not(cruise = 198808 and vessel = 21 and haul in (9,14,105,106,107,108,113,114,115))
    and not(cruise = 198202 and vessel = 21 and haul in (6,7,10,40))
    and not(cruise = 198502 and vessel = 60 and haul in (253));
  
  drop table hauls_ebs_nbs;

  create table hauls_ebs_nbs as  
  select CRUISEJOIN, HAULJOIN, REGION, VESSEL, CRUISE,
  HAUL, HAUL_TYPE, PERFORMANCE, START_TIME, DURATION,
  DISTANCE_FISHED, NET_WIDTH, NET_MEASURED, NET_HEIGHT,
  STRATUM, START_LATITUDE, END_LATITUDE, START_LONGITUDE,
  END_LONGITUDE, STATIONID, GEAR_DEPTH, BOTTOM_DEPTH,
  BOTTOM_TYPE, SURFACE_TEMPERATURE, GEAR_TEMPERATURE,
  WIRE_LENGTH, GEAR, ACCESSORIES, SUBSAMPLE, AUDITJOIN,
  YEAR from hauls_survey
  UNION ALL SELECT * FROM hauls_nbs;
  
    
    
    
    