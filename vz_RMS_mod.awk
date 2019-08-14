#
# awk script
#
#
#Find minimum RMS in vel.out file

$7=="RMS" && $8=="RESIDUAL=" { rms=$9 }
	
END { print rms }
