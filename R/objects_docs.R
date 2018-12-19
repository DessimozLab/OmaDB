#' An example HOG object. 
#'
#' An object containing information for the HOG:0273533.1b.
#'
#' @format An S3 object with 8 variables:
#' \describe{
#'   \item{hog_id}{hog identifier}
#'   \item{level}{the taxonomic level of this hog}
#'   \item{levels_url}{url pointer to the hog information at a given level}
#'   \item{members_url}{url pointer to the list of gene members for this hog}
#'   \item{alternative_members}{a dataframe object containing the rest of the taxonomic levels in this hog}
#'   \item{roothog_id}{the root taxonomic level of this hog}
#'   \item{parent_hogs}{a dataframe containing information on the parent hogs to the current hogs}
#'   \item{children_hogs}{a dataframe containing information on the children hogs to the current hogs}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/hog/HOG:0273533.1b/}
"hog"

#' An example OMA group object. 
#'
#' An object containing information for the OMA group number 737636.
#'
#' @format An S3 object with 4 variables:
#' \describe{
#'   \item{group_nr}{group number, not stable across releases}
#'   \item{fingerprint}{fingerprint of the oma group, stable across releases}
#'   \item{related_groups}{url to the endpoint containing the list of oma groups that share some of the orthologs with this oma group}
#'   \item{members}{list of protein members of this oma group}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/group/YEAST58/}
"group"

#' An example orthologs object. 
#'
#' A dataframe containing information for the orthologs of protein YEAST00058.
#'
#' @format A dataframe object with 15 variables:
#' \describe{
#'   \item{entry_nr}{entry number of the ortholog}
#'   \item{omaid}{oma identifier of the ortholog}
#'   \item{canonicalid}{canonicalid of the ortholog}
#'   \item{sequence_md5}{sequence_md5 of the ortholog}
#'   \item{oma_group}{oma_group of the ortholog}
#'   \item{oma_hog_id}{hog id of the ortholog}
#'   \item{chromosome}{chromosomal location of the ortholog}
#'   \item{locus.start}{start locus of the ortholog}
#'   \item{locus.end}{end locus of the ortholog}
#'   \item{locus.strand}{locus strand of the ortholog}
#'   \item{is_main_isoform}{true/false}
#'   \item{rel_type}{relationship type of the ortholog to the gene}
#'   \item{distance}{ortholog distance}
#'   \item{score}{ortholog score}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/protein/YEAST00058/orthologs}
"orthologs"

#' An example genome alignment object. 
#'
#' A dataframe containing information for the whole genome aligment of YEAST and ASHGO.
#'
#' @format A dataframe object with 12 variables for each member of the pair, as well some 3 additional variables:
#' \describe{
#'   \item{entry_nr}{entry number of the ortholog}
#'   \item{omaid}{oma identifier of the ortholog}
#'   \item{canonicalid}{canonicalid of the ortholog}
#'   \item{sequence_md5}{sequence_md5 of the ortholog}
#'   \item{oma_group}{oma_group of the ortholog}
#'   \item{oma_hog_id}{hog id of the ortholog}
#'   \item{chromosome}{chromosomal location of the ortholog}
#'   \item{locus.start}{start locus of the ortholog}
#'   \item{locus.end}{end locus of the ortholog}
#'   \item{locus.strand}{locus strand of the ortholog}
#'   \item{is_main_isoform}{true/false}
#'   \item{rel_type}{relationship type of the ortholog to the gene}
#'   \item{distance}{ortholog distance}
#'   \item{score}{ortholog score}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/pairs/YEAST/ASHGO/}
"pairs"

#' An example protein object. 
#'
#' An object containing information for the YEAST00058 protein.
#'
#' @format A S3 object with 23 variables:
#' \describe{
#'   \item{entry_nr}{entry number of the protein}
#'   \item{entry_url}{url pointer to the protein}
#'   \item{omaid}{oma identifier of the protein}
#'   \item{canonicalid}{canonicalid of the protein}
#'   \item{sequence_md5}{sequence_md5 of the protein}
#'   \item{oma_group}{oma_group of the protein}
#'   \item{oma_hog_id}{hog id of the protein}
#'   \item{chromosome}{chromosomal location of the protein}
#'   \item{locus}{GRanges object with the locus information for the protein}
#'   \item{is_main_isoform}{true/false}
#'   \item{roothog_id}{root taxonomic level of the relevant hog}
#'   \item{roothog_id}{taxonomic levels of the hog in which the protein is present}
#'   \item{sequence_length}{length of the protein sequence}
#'   \item{sequence}{AAString of the protein sequence}
#'   \item{cdna}{DNAString of the protein sequence}
#'   \item{domains}{url pointer to the list of protein domains}
#'   \item{xref}{url pointer to the list of protein cross references}
#'   \item{orthologs}{url pointer to the list of protein orthologs}
#'   \item{homeologs}{url pointer to the list of protein homeologs}
#'   \item{gene_ontology}{url pointer to the list of protein GO ontologies}
#'   \item{oma_group_url}{url pointer to the protein oma group}
#'   \item{oma_hog_members}{url pointer to the protein hog members}
#'   \item{alternative_isoforms_urls}{list of url pointers to the protein isoforms}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/protein/6633022/}
"protein"

#' An example newick format taxonomy object. 
#'
#'
#' @format An S3 with 2 variables:
#' \describe{
#'   \item{root_taxon}{sequence that was queried}
#'   \item{newick}{taxonomy newick}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/taxonomy/Alveolata/?type=newick}
"taxonomy"


#' An example dataframe containing proteins identified from a given sequence. 
#'
#'
#' @format A dataframe with 3 variables:
#' \describe{
#'   \item{query}{sequence that was queried}
#'   \item{identified_by}{type of identification}
#'   \item{targets}{list of protein targets identified}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/sequences/?query=MNDPSLLGYPNVGPQQQQQQQQQQHAGLLGKGTPNALQQQLHMNQLTGIPPPGLMNNSDVHTSSNNNSRQLLDQLANGNANMLNMNMDNNNNNNNNNNNNNNNGGGSGVMMNASTAAVNSIGMVPTVGTPVNINVNASNPLLHPHLDDPSLLNNPIWKLQLHLAAVSAQSLGQPNIYARQNAMKKYLATQQAQQAQQQAQQQAQQQVPGPFGPGPQAAPPALQPTDFQQSHIAEASKSLVDCTKQALMEMADTLTDSKTAKKQQPTGDSTPSGTATNSAVSTPLTPKIELFANGKDEANQALLQHKKLSQYSIDEDDDIENRMVMPKDSKYDDQLWHALDLSNLQIFNISANIFKYDFLTRLYLNGNSLTELPAEIKNLSNLRVLDLSHNRLTSLPAELGSCFQLKYFYFFDNMVTTLPWEFGNLCNLQFLGVEGNPLEKQFLKILTEKSVTGLIFYLRDNRPEIPLPHER}
"sequence_map"

#' An example xref object. 
#'
#'
#' @format A dataframe with 8 variables:
#' \describe{
#'   \item{xref}{cross reference}
#'   \item{source}{source of the cross reference}
#'   \item{entry_nr}{oma database entry number}
#'   \item{oma_id}{oma id of the cross reference}
#'   \item{genome.code}{genome_id of the cross reference}
#'   \item{genome.taxon_id}{taxon_id of the cross reference}
#'   \item{genome.species}{species of the cross reference}
#'   \item{genome.genome_url}{genome url pointer of the cross reference}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/xref/?search=MAL}
"xref"

#' An example dataframe containing GO annotations identified from a given sequence. 
#'
#'
#' @format A dataframe with 13 variables:
#' \describe{
#'   \item{Qualifier}{qualifier of the annotation}
#'   \item{GO_ID}{GO term for the annotation}
#'   \item{With}{GO term for the annotation}
#'   \item{Evidence}{evidence for the annotation}
#'   \item{Date}{date}
#'   \item{DB_Object_Type}{identified object type}
#'   \item{DB_Object_Name}{identified object name}
#'   \item{Aspect}{aspect}
#'   \item{Assigned_By}{assignment of the annotation}
#'   \item{GO_name}{GO term name}
#'   \item{DB}{database}
#'   \item{DB.Reference}{database reference}
#'   \item{Synonym}{synonym}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/function/?query=MNDPSLLGYPNVGPQQQQQQQQQQHAGLLGKGTPNALQQQLHMNQLTGIPPPGLMNNSDVHTSSNNNSRQLLDQLANGNANMLNMNMDNNNNNNNNNNNNNNNGGGSGVMMNASTAAVNSIGMVPTVGTPVNINVNASNPLLHPHLDDPSLLNNPIWKLQLHLAAVSAQSLGQPNIYARQNAMKKYLATQQAQQAQQQAQQQAQQQVPGPFGPGPQAAPPALQPTDFQQSHIAEASKSLVDCTKQALMEMADTLTDSKTAKKQQPTGDSTPSGTATNSAVSTPLTPKIELFANGKDEANQALLQHKKLSQYSIDEDDDIENRMVMPKDSKYDDQLWHALDLSNLQIFNISANIFKYDFLTRLYLNGNSLTELPAEIKNLSNLRVLDLSHNRLTSLPAELGSCFQLKYFYFFDNMVTTLPWEFGNLCNLQFLGVEGNPLEKQFLKILTEKSVTGLIFYLRDNRPEIPLPHER}
"sequence_annotation"
