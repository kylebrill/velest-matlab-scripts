#
#awk script
#extract events from velest single-event-mode .loc output
#
#command is awk -f extract_singleloc.awk [filename]

$2=="DATE" && $3=="ORIGIN" { getline; lon=substr(-$4,1,8);                 lat=substr($3,1,7); z=substr($5-2.0,1,5);
    print lon, lat, z >> "/local/kabrill_res/velest-tight/singlemods/m_mean/locations_raw.dat"}

$1=="AMX" && $2=="PRX" { 
	fnd=1
	if (fnd=1 && !NF) {fnd=0}
	while (fnd=1) {
		getline
		if (fnd=1 && !NF) {fnd=0;next}
		print $6 >> "/local/kabrill_res/velest-tight/singlemods/m_mean/weights.dat"
		}
	}
