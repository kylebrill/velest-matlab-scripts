#
# awk script
#
#
#Find station with delays

$1=="sta" && $2=="phase" {flag=1;next}
/WARNING/ {next}
$1=="~~~" {flag=0}

flag {
delay=substr($NF,length($NF)-7,length($NF)); 
print $1
}
