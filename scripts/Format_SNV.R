library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

get_Info = function( mut ){
	mutType = sapply( mut , function(x){ unlist( strsplit( x , "|" , fixed = TRUE ) )[2] } )
	gene = sapply( mut , function(x){ unlist( strsplit( x , "|" , fixed = TRUE ) )[4] } )
	names( mutType ) = names( gene ) = NA
	
	info = cbind( mutType , gene )

	info
}

get_Type = function( ref , alt ){
	ifelse( nchar(ref) != nchar(alt) , "INDEL", "SNV" )

}

################################################
################################################

files = list.files(file.path(input_dir, "annot_vcf"))
data = NULL
for(i in 1:length(files)){
	lines = as.numeric( unlist( strsplit( str_trim(system( paste( "wc -l", file.path(input_dir, "annot_vcf", files[i])) , intern = TRUE ), side=("both")) , " " , fixed=TRUE) )[1] )
	if( lines > 0 ){
		vcf = read.csv( paste( file.path(input_dir, "annot_vcf") , files[i] , sep = "/" ) , sep = "\t" , header = FALSE , skip = 5 )
		data = rbind( data , cbind( unlist( strsplit( files[i] , "." , fixed = TRUE ))[1] , vcf ) )
	}
}


data[ , 5 ] = ifelse( data[ , 5 ] %in% "-" , "" , data[ , 5 ] )
data[ , 6 ] = ifelse( data[ , 6 ] %in% "-" , "" , data[ , 6 ] )


variant = c( "^AU$" , "^CU$" , "^GU$" , "^TU$" , "^TIR$" )
names(variant) = c( "A" , "C" , "G" , "T" , "I" )

snv_data = NULL
for( i in 1:nrow(data) ){

	print( paste( i , "in" , nrow(data) , '(' , round( i/nrow(data) * 100 ) , '%)' ) )

	dat = data[i , ]

	annot = as.character( dat[ "V9" ] )
	dp = grep( "^DP$", unlist( strsplit( annot , ":" , fixed = TRUE )) )
	alt = grep( variant[ ifelse( nchar( as.character( dat[ "V5" ] ) ) > 1 , "I" , as.character( dat[ "V5" ] ) )  ], unlist( strsplit( annot , ":" , fixed = TRUE )) )
	
	dp = as.numeric( as.character( sapply( dat[ "V11" ] , function(x){ unlist( strsplit( x , ":" , fixed = TRUE ) )[dp] } ) ) )
	alt = as.numeric( as.character( sapply( dat[ "V11" ] , function(x){ unlist( strsplit(  unlist( strsplit( x , ":" , fixed = TRUE ) )[alt] , "," , fixed = "TRUE" ))[1] } ) ) )

	vaf = alt/dp *100
	
	if( !is.na(vaf) & dp >= 20 & vaf > 10 ){
		snv_data = rbind( snv_data , 
							cbind( data[ i , c( 1:3 , 5:6 ) ] , get_Info( mut = data[ i , 9 ] ) , get_Type( ref = data[ i , 5 ] , alt = data[ i , 6 ] ) ) )
	}
}

colnames( snv_data ) = c( "Sample" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "Gene" , "MutType" )

snv_data$Effect[ snv_data$Effect %in% "downstream_gene_variant" ] = "3'Flank" 
snv_data$Effect[ snv_data$Effect %in% "upstream_gene_variant" ] = "5'Flank" 
snv_data$Effect[ snv_data$Effect %in% "3_prime_UTR_variant" ] = "3'UTR" 
snv_data$Effect[ snv_data$Effect %in% c( "5_prime_UTR_premature_start_codon_gain_variant" , "5_prime_UTR_variant" ) ] = "5'UTR" 
snv_data$Effect[ snv_data$Effect %in% c( "synonymous_variant" , "initiator_codon_variant" ) ] = "Silent" 
snv_data$Effect[ snv_data$Effect %in% c( "missense_variant" , "missense_variant&splice_region_variant" ) ] = "Missense_Mutation" 
snv_data$Effect[ snv_data$Effect %in% c( "intergenic_region" , "intron_variant" , "non_coding_transcript_exon_variant" ) ] = "Intron" 
snv_data$Effect[ snv_data$Effect %in% c( "splice_acceptor_variant&intron_variant" , "splice_donor_variant&intron_variant" , "splice_region_variant" , "splice_region_variant&intron_variant" , "splice_region_variant&non_coding_transcript_exon_variant" , "splice_region_variant&synonymous_variant" ) ] = "Splice_Site" 
snv_data$Effect[ snv_data$Effect %in% c( "missense_variant" , "missense_variant&splice_region_variant" ) ] = "Missense_Mutation" 

snv_data$Effect[ snv_data$Effect %in% c( "frameshift_variant" , "frameshift_variant&splice_region_variant" , "frameshift_variant&start_lost" , "frameshift_variant&stop_gained" ) ] = "Frame_Shift_Del" 
snv_data$Effect[ snv_data$Effect %in% c( "conservative_inframe_insertion" , "disruptive_inframe_insertion" , "disruptive_inframe_insertion&splice_region_variant" ) ] = "In_Frame_Ins" 

snv_data$Effect[ snv_data$Effect %in% c( "stop_gained" , "stop_gained&conservative_inframe_insertion" , "stop_gained&splice_region_variant" , "start_lost" , "start_retained_variant" ) ] = "Nonsense_Mutation" 


################################################
################################################
data = snv_data

rna = read.csv( file.path(input_dir, "annot_WES.txt"), stringsAsFactors=FALSE , sep="\t"  )
rownames(rna)=rna[,3]

patient = intersect( unique( sort( data$Sample ) ) , rna$patient )

data = data[ data$Sample %in% patient , ]
data$Sample = rna[ data$Sample , ]$run

data = data[ , c( "Sample" , "Gene" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" ) ]

write.table( data , file=file.path(output_dir, "SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

