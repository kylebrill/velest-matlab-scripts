#
#awk script
#extract final hypocenters from velest simultaneous-mode.OUT output
#
#

$1=="~~~" && $2=="output" && $3=="final" {fnd=1; first=1
next }

{
if (fnd==1 && first==1 && $1=="date" && $2=="origin")
while (fnd=1) {
    getline

    lon=substr(-$6,1,8);lat=substr($5,1,7); z=substr($7-2.0,1,5);
    print lon, lat, z >> "locations_raw.dat";

if (!NF)
break;

}}


