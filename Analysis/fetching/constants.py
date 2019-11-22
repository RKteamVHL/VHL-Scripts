import os
LINE_DELIMITER = '\n'
ROW_DELIMITER= '\t'
LIST_DELIMITER= ';'


FILES_DIRECTORY = os.path.join('ApiScripts','files')

CLINVAR_NAME = 'ClinVar'
CLINVAR_FILENAME = "ClinVar.txt"
CLINVAR_HREF = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_HEADER = "AlleleID{0}Type{0}Name{0}GeneID{0}GeneSymbol{0}HGNC_ID{0}ClinicalSignificance{0}ClinSigSimple{0}LastEvaluated{0}RS# (dbSNP){0}nsv/esv (dbVar){0}RCVaccession{0}PhenotypeIDS{0}PhenotypeList{0}Origin{0}OriginSimple{0}Assembly{0}ChromosomeAccession{0}Chromosome{0}Start{0}Stop{0}ReferenceAllele{0}AlternateAllele{0}Cytogenetic{0}ReviewStatus{0}NumberSubmitters{0}Guidelines{0}TestedInGTR{0}OtherIDs{0}SubmitterCategories{0}VariationID\n".format(
	ROW_DELIMITER)


STUDENTS_NAME = 'KimStudents2019'
STUDENTS_FILENAME = 'CiVIC Extracted Data (Summer 2019-Present).txt'
STUDENTS_SHEETIDS = ['1826130937', '0', '763995397', '511845134', '301841530']
STUDENTS_HREF = [f'https://docs.google.com/spreadsheets/d/1YYXHLM9dOtGyC1l8edjQnfcAf0IM1EuTrz8ZKbQkQGU/export?format=tsv&gid={name}' for name in STUDENTS_SHEETIDS ]
STUDENTS_HEADER = "PMID{0}cDNA_Position{0}Multiple Mutants in Case{0}Mutation Event c.DNA.{0}Predicted Consequence Protein Change{0}variant_name{0}Mutation Type{0}Kindred Case{0}Familial/Non-familial{0}Phenotype{0}Reference{0}Age{0}Notes{0}Evidence Statements\n".format(
	ROW_DELIMITER)


GNOMAD_POST_DATA = {
	"query": '''{
		transcript(transcript_id: "ENST00000256474", reference_genome: GRCh37){
			variants(dataset: gnomad_r2_1) {
				consequence
				flags
				hgvs
				hgvsc
				hgvsp
				lof
				lof_filter
				lof_flags
				pos
				rsid
				variantId
				exome {
					ac
					ac_hemi
					ac_hom
					an
					af
					filters
					populations {
						id
						ac
						an
						ac_hemi
						ac_hom
					}
				}
				genome {
					ac
					ac_hemi
					ac_hom
					an
					af
					filters
					populations {
						id
						ac
						an
						ac_hemi
						ac_hom
					}
				}
			}
		}
	}'''
}

GNOMAD_NAME = 'gnomAD_v2_1'
GNOMAD_FILENAME = 'gnomAD_v2_1.txt'
GNOMAD_HREF = 'http://gnomad.broadinstitute.org/api'

# 58 corresponds to the VHL gene 
CIVIC_VHL_GENE_ID =  58
CIVIC_NAME = 'CIViC'
CIVIC_FILENAME = "CIViC.txt"