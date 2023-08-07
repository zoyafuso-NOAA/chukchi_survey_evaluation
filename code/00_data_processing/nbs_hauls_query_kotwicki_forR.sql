select * from racebase.haul
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
      or (cruise = 202102 and haul_type = 3 and (stratum in (81,70,71)))
      or (cruise = 202202 and haul_type = 3 and (stratum in (81,70,71)) or stratum is NULL)
      or (cruise = 202302 and haul_type = 3 and (stratum in (81,70,71)) or stratum is NULL))
    and region = 'BS'
    and performance >= 0
    and not(cruise = 198808 and vessel = 21 and haul in (9,14,105,106,107,108,113,114,115))
    and not(cruise = 198202 and vessel = 21 and haul in (6,7,10,40))
    and not(cruise = 198502 and vessel = 60 and haul in (253));