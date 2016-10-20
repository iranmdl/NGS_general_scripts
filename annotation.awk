### Run: awk -f annotation.awk File1.txt File2.txt | column -t
## Script para añadir el functional annotation column a un fichero que contiene "Loci - fold2change - pvalue -padj"
# 

##################################
# Los ficheros NO tienen header  #
##################################

## Process first file of arguments. Save 'id' as key and 'No' as value
## of a hash.
FNR == NR {

    hash[ $1 ] = $2
    next
}

## Process second file of arguments.
FNR < NR {
    if ( $1 in hash ) { 
        printf "%s %s\n", $0, ( hash[ $1 ] ) 
    }
    else
    	printf "%s\n", $0
}


## Functional annotation
# awk -f ../annotation.awk ../ITAG2.4_NAMES_parsed_uniq_forAWK.txt sign_up.txt | column -t > FA_up_BI.txt
# awk -f ../annotation.awk ../ITAG2.4_NAMES_parsed_uniq_forAWK.txt sign_down.txt | column -t > FA_down_BI.txt

# sed 's/ \+ /\t/g' FA_up_BI.txt > FA_up_BI_1.txt 	#Replace more than 1 whitespace with tab
# sed 's/;/ /g' FA_up_BI_1.txt > FA_up_BI_2.txt 	#Replace ; with one whitespace

