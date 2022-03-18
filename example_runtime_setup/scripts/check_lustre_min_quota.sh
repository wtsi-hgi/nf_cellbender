#!/usr/bin/env bash
set -e
set -o pipefail
# script to check that Nextflow BASE_DIR has enough disk space quota available to run pipeline.

[[ -z "$1" ]] && { echo "Error: parameter 1 is empty" ; exit 1; }
[[ -z "$2" ]] && { echo "Error: Parameter 2 is empty" ; exit 1; }
####### CLI ARGUMENTS:
export BASE_DIR="$1" 
echo first script argument BASE_DIR is "$BASE_DIR"
# $1 must be an absolute /lustre path, which is part of the humgen Unix groups.
# example: /lustre/scratch119/humgen/projects/team227_bioinfo/pipelines/dgrna/inputs_crams
export MIN_QUOTA="$2"
echo second script argument MIN_QUOTA is "$MIN_QUOTA"
# $2 is minimum quota required to run the pipeline. Use human readable format. Can be a fraction.
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
[[ -z "$MACHINE_MIN_QUOTA" ]] && { echo "Error: Min quota not in correct format, check script argument." ; exit 1; }

UNIX_GROUP=$(ls -ld "$BASE_DIR" | cut -f4 -d' ')
echo UNIX_GROUP is "$UNIX_GROUP"
SCRATCH=$(echo "$BASE_DIR" | grep -o ".lustre.scratch[0-9]*")
echo SCRATCH is "$SCRATCH"
[[ -z "$SCRATCH" ]] && { echo "Error: scratch area not determined, check script argument." ; exit 1; }
lfs quota -gh "$UNIX_GROUP" "$SCRATCH"
DISK_USED=$(lfs quota -gh "$UNIX_GROUP" "$SCRATCH" | tail -n 1 | xargs echo | cut -f1 -d" ")
echo DISK_USED is "$DISK_USED"
DISK_TOTAL=$(lfs quota -gh "$UNIX_GROUP" "$SCRATCH" | tail -n 1 | xargs echo | cut -f3 -d" ")
echo DISK_TOTAL is "$DISK_TOTAL"
MACHINE_DISK_USED=$(dehumanise "$DISK_USED")
echo MACHINE_DISK_USED is "$MACHINE_DISK_USED"
MACHINE_DISK_TOTAL=$(dehumanise "$DISK_TOTAL")
echo MACHINE_DISK_TOTAL is "$MACHINE_DISK_TOTAL"

MACHINE_DISK_LEFT="$(($MACHINE_DISK_TOTAL-$MACHINE_DISK_USED))"
echo MACHINE_DISK_LEFT is "$MACHINE_DISK_LEFT"
DISK_LEFT=$(printf %s\\n "$MACHINE_DISK_LEFT" | numfmt --to=iec-i)
echo DISK_LEFT is "$DISK_LEFT"

if (( "$MACHINE_DISK_LEFT" < "$MACHINE_MIN_QUOTA" )); then
    echo -e "\nnot enough disk quota at BASE_DIR $BASE_DIR (according to \"lfs quota -hg $UNIX_GROUP $SCRATCH\")"
    echo "quota disk left $DISK_LEFT for Unix group $UNIX_GROUP < min quota required $MIN_QUOTA"
    echo "will now exit with exitcode 1"
    exit 1
fi
echo -e "\nenough disk quota at BASE_DIR $BASE_DIR (according to \"lfs quota -hg $UNIX_GROUP $SCRATCH\")"
echo "quota disk left $DISK_LEFT for Unix group $UNIX_GROUP > min quota required $MIN_QUOTA"
