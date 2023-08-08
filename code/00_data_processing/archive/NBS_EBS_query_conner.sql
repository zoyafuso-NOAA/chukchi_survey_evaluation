select * from racebase.haul
where cruisejoin in (select cruisejoin
                    from safe.survey
                    where survey in ('NBS_SHELF')
                    )
and distance_fished is not null
and net_width is not null;