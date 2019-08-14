#
# awk script
#
#
#Find one specific earthquake location #515

BEGIN { fnd=0; lon=0 }

$1=="~~~" && $2=="output" && $3=="final"  {  fnd=1; first=1
next }

{
if (fnd==1 && first==1 && $1=="315" ) {
   lon=substr(-$6,1,8); lat=substr($5,1,7); z=substr($7-1.5,1,5);
   print lon, lat, z }
}


