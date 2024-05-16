check_pkgs_path <- file.path(dirname(sys.frame(1)$ofile), 
                             "../../utils/check_and_install_packages.R")
source(check_pkgs_path)

# loading dependencies
# Define the list of required packages
required_packages <- c("hash", "jsonlite")

# Check and install the required packages
check_and_install_packages(required_packages)

get_max_kmer_count_peptides_df <- function(list_of_sign_peptide_ids, kmer_db_path, key_table){
  # helper functions
  get_kmer_hash <- function(list_of_sign_peptide_ids, key_table, kmer_db_path) {  
    # sort first
    pep_order <- key_table[key_table$target_id %in% list_of_sign_peptide_ids, c('target_inter_id', 'target_id')]
    pep_order <- pep_order[order(pep_order$target_inter_id),]
    peptides_to_find <- pep_order$target_id

    # creating selected peptides hash:
    peptides_kmers <- hash()

    # going over the kmer db to extract kmers for our peptides:
    file_connection <- file(kmer_db_path, open='r')

    pep_idx <- 1
    while(length(line <- readLines(file_connection, n = 1, warn = F)) > 0 & 
          pep_idx <= length(peptides_to_find)) {
      tryCatch({
        json_data <- fromJSON(line)
        k <- json_data$peptide_id
        if(k == peptides_to_find[pep_idx]) {
          peptides_kmers[[k]] <- json_data$kmers
          pep_idx <- pep_idx + 1
        }
      }, error = function(e) {
        cat("the length ")
        cat("Error parsing JSON on line:", line, "\n")
        cat("Error message: ", e$message, "\n")
      })
    }
    close(file_connection)

    if (length(keys(peptides_kmers)) != length(list_of_sign_peptide_ids)) {
      print("Error! something went wrong: # of found peptides does not match # of provided peptide ids")
      return(NA)
    }

    return(peptides_kmers)
  }
  get_kmer_tally <- function(curr_kmer_hash){
    kmer_tally <- hash()

    peptide_ids <- keys(curr_kmer_hash)
    for(pep_id in peptide_ids) {
      tryCatch({
        lapply(curr_kmer_hash[[pep_id]], function(x) {
          if(is.null(kmer_tally[[x]])) {
            kmer_tally[[x]] <<- 1
          } else {
            kmer_tally[[x]] <<- kmer_tally[[x]] + 1
          }
          NULL
        })
      }, error = function(e) {
        cat(pep_id, "\n")
        cat("Error message: ", e$message, "\n")
      })
    }

    return(kmer_tally)
  }
  get_pep_max_kmer_df <- function(curr_kmer_hash, kmer_tally){
    peptied_max_kmer_count_df <- data.frame(target_id=keys(curr_kmer_hash), max_kmer_count=NA)

    for(row_idx in 1:nrow(peptied_max_kmer_count_df)) {
      pep_id <- peptied_max_kmer_count_df[row_idx,1]
      all_kmers <- curr_kmer_hash[[pep_id]]
      max_count <- 0
      lapply(all_kmers, function(k) {
        if(kmer_tally[[k]] > max_count) max_count <<- kmer_tally[[k]]
        NULL
      })
      peptied_max_kmer_count_df[row_idx,2] <- max_count
    }
    return(peptied_max_kmer_count_df)
  }

  # extracting kmers for each peptide of interest from pre-made kmer db jsonl
  curr_kmer_hash <- get_kmer_hash(list_of_sign_peptide_ids, key_table, kmer_db_path)
  
  # calculating tally for all kmers present in our data
  kmer_tally <- get_kmer_tally(curr_kmer_hash)
  
  # calculating max kmer count for each peptide
  my_df_pep_kmer_count <- get_pep_max_kmer_df(curr_kmer_hash, kmer_tally)
  
  # returning table of peptides of interest and max kmer count for each of them
  return(my_df_pep_kmer_count)
}

# Usage Instructions
# This function takes as an argument:
#   List of unique peptide ids of interest (list_of_sign_peptide_ids)

#   Kmer database .jsonl file path (kmer_db_path)for all peptides in the library, 
#   generated separatley by a python script located in:
#   bin/reusable_small_scripts/phipseq-key-kmer-parser/phipseq-key-kmer-parser.py

#   Library key table for current phipseq library. It must contain
#   target_inter_id and target_id columns. target_id column must contain all
#   peptides of interest(list_of_sign_peptide_ids). This file also should be used
#   with above mentioned python script to generate .jsonl kmer database file.

# It returns a data frame with two columns:
#   target_id - contains peptides of interest (list_of_sign_peptide_ids)
#   max_kmer_count - integer value, max count of one of the peptide kmers
#     (kmer database contains only these kmers that are unique for a given 
#     peptide)


# This function can be sourced and used in your R scripts or R Markdown documents.

# Example usage:
#
# source("../../../../../bin/R/phipseq/max_kmer_count/v0.0.0.R")
# list_of_sign_peptide_ids <- c("peptide1", "peptide2", "peptide3")
# kmer_db_path <- "path/to/your/kmer_db.jsonl"
# key_table <- read.csv("path/to/your/key_table.csv")
# result <- get_max_kmer_count_peptides_df(list_of_sign_peptide_ids, kmer_db_path, key_table)

# P.S. Developed and tested in:
# /Users/mac_studio_np/Documents/data_analysis/phipseq/2024-04-30_tprC_manual_analysis_exploration/analysis/scripts/2024-06-09_single_rabbit_deseq.qmd
