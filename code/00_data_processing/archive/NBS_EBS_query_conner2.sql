select *
from racebase.haul
where region = 'BS'
and cruisejoin not in (select cruisejoin
                    from safe.survey
                    where survey in ('EBS_SLOPE')
                    )
and cruise >= 198201
and distance_fished is not null
and net_width is not null
and ( gear between 30 and 44
      or
      gear in (47, 178, 180, 181)
);

