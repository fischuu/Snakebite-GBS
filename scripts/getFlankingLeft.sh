input=$1
output=$2
cut -f1,2,3,5,6 $input | awk 'BEGIN { OFS = "\t" } {a=0; b=0; for (i=1;i<=NF;i++) if ($i < a || i == 1) a = $i; else if($i > b|| i == 1)b = $i; print $1,a-6, a}' > $output