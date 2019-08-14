#
# awk-script
#
# Erstellung des Vp-Start-Modelles als Input fuer GMT
# aus dem VELEST-output

BEGIN { fnd=0 }

$1=="layer" && $2=="vel" && $3=="depth" {  fnd=1;next }
		
fnd==1  {	vp=$2; z=$3						
                print $3,$2
		fnd=2
		next }

$1==""  { if (fnd==2) {
	  print "80.00",vp 
	  exit }
	}

  { if (fnd==2) { 
	print $3,vp
	vp=$2; z=$3
	print $3,$2 }
  }

