awk 'BEGIN { OFS = "\t" } { $5 = $3 - $2 ; $6 = $4 / $5} 1' $1