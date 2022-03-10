awk -v P=$2 '{print P"/"$2".count.txt"}' $1
