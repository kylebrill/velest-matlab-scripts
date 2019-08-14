#
# awk-script
#
# Erstellung des Vp-Final-Modelles als Input fuer GMT
# aus dem VELEST-output

BEGIN { fnd=0; vp=0 }

$1=="nlay"   {  fnd=1; first=1
		next }

($1=="" || ($1>2 && $2=="0.00..."))  { if (fnd==1 && vp!=0) {
	   print "80.0",vp
	   fnd=0; vp=0 }
          else next
	}

  { 
    if (fnd==1 && first==0) {
        print substr($2,1,length($2)-3),vp
	vp=$5; if (vp=="km/s") vp=$4
        z=substr($2,1,length($2)-3)
	print z,vp }
    if (fnd==1 && first==1) {
	print substr($2,1,length($2)-3),$5
	vp=$5; z=substr($2,1,length($2)-3)
	first=0 }
  }

