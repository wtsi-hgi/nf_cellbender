#!/usr/bin/env bash
set -e
set -o pipefail
# script to check that farm home dir has enough disk space for Nextflow.
# By default Nextflow writes a small amount of files into home dir...

####### CLI ARGUMENTS:
export MIN_QUOTA="${1:-500M}"
echo first script argument MIN_QUOTA is "$MIN_QUOTA"
# $1 is minimum quota required in farm home dir to run Nextflow. Use human readable format. Can be a fraction.
# example: 0.5T 
# example: 400GiB 
# example: 2T

####### SCRIPT:

# function to convert human readable disk space (like 2Tb or 500G) to machine readable numbers:
dehumanise() {
  for v in "${@:-$(</dev/stdin)}"
  do  
    echo $v | awk \
      'BEGIN{IGNORECASE = 1}
       function printpower(n,b,p) {printf "%u\n", n*b^p; next}
       /[0-9]$/{print $1;next};
       /K(iB)?$/{printpower($1,  2, 10)};
       /M(iB)?$/{printpower($1,  2, 20)};
       /G(iB)?$/{printpower($1,  2, 30)};
       /T(iB)?$/{printpower($1,  2, 40)};
       /KB$/{    printpower($1, 10,  3)};
       /MB$/{    printpower($1, 10,  6)};
       /GB$/{    printpower($1, 10,  9)};
       /TB$/{    printpower($1, 10, 12)}'
  done
} 

MACHINE_MIN_QUOTA=$(dehumanise "$MIN_QUOTA")
echo MACHINE_MIN_QUOTA is "$MACHINE_MIN_QUOTA"

quota -s 
DISK_USED=$(quota -s | head -n 4 | tail -n 1 | xargs echo | cut -f1 -d" ")
echo DISK_USED is "$DISK_USED"
MACHINE_DISK_USED=$(dehumanise "$DISK_USED")
echo MACHINE_DISK_USED is "$MACHINE_DISK_USED"

DISK_TOTAL=$(quota -s | head -n 4 | tail -n 1 | xargs echo | cut -f3 -d" ")
echo DISK_TOTAL is "$DISK_TOTAL"
MACHINE_DISK_TOTAL=$(dehumanise "$DISK_TOTAL")
echo MACHINE_DISK_TOTAL is "$MACHINE_DISK_TOTAL"

MACHINE_DISK_LEFT="$(($MACHINE_DISK_TOTAL-$MACHINE_DISK_USED))"
echo MACHINE_DISK_LEFT is "$MACHINE_DISK_LEFT"
DISK_LEFT=$(printf %s\\n "$MACHINE_DISK_LEFT" | numfmt --to=iec-i)
echo DISK_LEFT is "$DISK_LEFT"

if (( "$MACHINE_DISK_LEFT" < "$MACHINE_MIN_QUOTA" )); then
    echo -e "\nnot enough disk quota in home dir (according to \"quota -s\")"
    echo "quota disk left $DISK_LEFT in home dir < min quota required $MIN_QUOTA"
    echo "will now exit with exitcode 1"
    exit 1
fi
echo -e "\nenough disk quota in home dir (according to \"quota -s\")"
echo "quota disk left $DISK_LEFT for Unix group $UNIX_GROUP > min quota required $MIN_QUOTA"
