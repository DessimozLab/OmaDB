#' Map the Protein Sequence Function
#' 
#' The function to identify a sequence.
#'
#' @param query the sequence to be searched, it can be either a string or an AAString object from the Biostrings package
#' @param search argument to choose search strategy. Can be set to 'exact', 'approximate' or 'mixed'. Defaults to 'mixed', meaning first tries to find exact match. If no target can be found, uses approximate search strategy to identify query sequence in database.
#' @param full_length a boolean indicating whether or not for exact matches, the query sequence must be matching the full target sequence. By default, a partial exact match is also reported as exact match.
#' @return a data.frame containing the information of matches for the query sequence
#' @export
#' @examples
#' mapSequence(query="MNDPSLLGYPNVGPQQQQQQQQQQHAGLLGKGTPNALQQQLHMNQLTGIPPPGLMNNSDVHTSSNNNSRQLLDQLANGNANMLNMNMDNNNNNNNNNNNNNNNGGGSGVMMNASTAAVNSIGMVPTVGTPVNINVNASNPLLHPHLDDPSLLNNPIWKLQLHLAAVSAQSLGQPNIYARQNAMKKYLATQQAQQAQQQAQQQAQQQVPGPFGPGPQAAPPALQPTDFQQSHIAEASKSLVDCTKQALMEMADTLTDSKTAKKQQPTGDSTPSGTATNSAVSTPLTPKIELFANGKDEANQALLQHKKLSQYSIDEDDDIENRMVMPKDSKYDDQLWHALDLSNLQIFNISANIFKYDFLTRLYLNGNSLTELPAEIKNLSNLRVLDLSHNRLTSLPAELGSCFQLKYFYFFDNMVTTLPWEFGNLCNLQFLGVEGNPLEKQFLKILTEKSVTGLIFYLRDNRPEIPLPHERRFIEINTDGEPQREYDSLQQSTEHLATDLAKRTFTVLSYNTLCQHYATPKMYRYTPSWALSWDYRRNKLKEQILSYDSDLLCLQEVESKTFEEYWVPLLDKHGYTGIFHAKARAKTMHSKDSKKVDGCCIFFKRDQFKLITKDAMDFSGAWMKHKKFQRTEDYLNRAMNKDNVALFLKLQHIPSGDTIWAVTTHLHWDPKFNDVKTFQVGVLLDHLETLLKEETSHNFRQDIKKFPVLICGDFNSYINSAVYELINTGRVQIHQEGNGRDFGYMSEKNFSHNLALKSSYNCIGELPFTNFTPSFTDVIDYIWFSTHALRVRGLLGEVDPEYVSKFIGFPNDKFPSDHIPLLARFEFMKTNTGSKKV")
#' mapSequence(search="mixed",query="NKLLQPTDFQQSHIAEASKSLVDCTKQALMEMADTLTDSKTAKKQQPTGDSTPSGTATNSAVSTPLTPKIELFANGKDEANQALLQHKKLSQYSIDEDDDIENRMVMPKDSKYDDQLWHALDLSNLQIFNISANIFKYDFLTRLYLNGNSLTELPAEIKNLSNLRVLDLSHNRLTSLPAELGSCFQLKYFYFFDNMVTTLPWEFGNLCNLQFLGVEGNPLEKQFLKILTEKSVTGLIFYLRDNRPEIPLPHERRFIEINTDGEPQREYDSLQQSTEHLATDLAKRTFTVLSYNTLCQHYATPKMYRYTPSWALSWDYRRNKLKEQILSYDSDLLCLQEVESKTFEEYWVPLLDKHGYTGIFHAKARAKTMHSKDSKKVDGCCIFFKRDQFKLITKDAMDFSGAWMKHKKFQRTEDYLNRAMNKDNVALFLKLQHIPSGDTIWAVTTHLHWDPKFNDVKTFQVGVLLDHLETLLKEETSHNFRQDIKKFPVLICGDFNSYINSAVYELINTGRVQIHQEGNGRDFGYMSEKNFSHNLALKSSYNCIGELPFTNFTPSFTDVIDYIWFSTHALRVRGLLGEVDPEYVSKFIGFPNDKFPSDHIPLLARFEFMKTNTGSKKV")



mapSequence <- function(query,search,full_length=FALSE){
	
	if(missing(query)){
		stop("You must provide a sequence to query.")
	}

	if(class(query)=="AAString"){
		query = as.character(query)
	}
	
	if(missing(search)){
		url = urlGenerator(type="sequence",query_param1="query",query_param1_value=query,
			query_param2="full_length",query_param2_value=as.character(full_length))
	}

	else{

		search_values = c('approximate','exact','mixed')
		
		if(!search %in% search_values ){
			stop("search parameter invalid. Must be one of 'approximate', 'exact', 'mixed'")
		}
		else{
			url = urlGenerator(type="sequence",query_param1="query",query_param1_value=query,query_param2="full_length",
				query_param2_value=as.character(full_length),query_param3="search",query_param3_value=search)

		}
	}
	

	return(requestFactory(url))
}