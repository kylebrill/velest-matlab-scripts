#
# awk script
#
#
#Find station delays

$1=="sta" && $2=="phase" {flag=1;next}
/WARNING/ {next}
$1=="~~~" {flag=0}

flag {
delay=substr($NF,length($NF)-7,length($NF)); 
print delay
}
