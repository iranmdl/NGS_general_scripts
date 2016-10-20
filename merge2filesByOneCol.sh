awk 'BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  hash[ $1 ] = $2
  next
  }
FNR < NR {
    if ( $1 in hash ) { 
        printf "%s\t%s\n", $0, ( hash[ $1 ] ) 
    }
    else
    	printf "%s\n", $0
}' gene_names.txt wt80Tvswt80.txt > wt80Tvswt80_renamed.txt
