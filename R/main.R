
#' Encode DNA Sequences into a Feature Matrix
#'
#' This function converts a set of DNA sequences into a
#' nucleotide encoding matrix, where each column represents
#' a nucleotide position in the input sequences.
#'
#' @param dna_strings A character vector containing DNA sequences
#' of equal length.
#'
#' @return A data frame in which each column represents a
#' nucleotide position and contains factor values with levels
#' A, T, C, and G.
#'
#' @examples
#' dna_encoding(c("ATGCA", "GGTAC"))
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Sites for Multiple Samples
#'
#' This function predicts m6A probabilities and prediction
#' statuses for multiple RNA samples using a trained
#' random forest model.
#'
#' @param ml_fit A trained random forest model.
#' @param feature_df A data frame containing the input
#' features required for m6A prediction.
#' @param positive_threshold A numeric value between 0 and 1
#' used as the classification threshold for positive prediction.
#'
#' @return A data frame containing the original input features,
#' predicted m6A probabilities, and predicted m6A statuses.
#'
#' @examples
#' # Load the example random forest model
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'
#' # Load example input data
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#'
#' prediction_multiple(ml_fit = ml_fit, feature_df = example_df, positive_threshold = 0.5)
#'
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}


#' Single Sample Prediction for m6A Sites
#'
#' This function takes individual feature values and uses
#' a pre-trained machine learning model to predict the m6A
#' probability and status for a single sample.
#'
#' @param ml_fit A trained random forest model.
#' @param gc_content A numeric value representing GC content.
#' @param RNA_type A character string specifying the RNA type.
#' @param RNA_region A character string specifying the RNA region.
#' @param exon_length A numeric value representing exon length.
#' @param distance_to_junction A numeric value representing
#' the distance to the nearest exon junction.
#' @param evolutionary_conservation A numeric value representing
#' evolutionary conservation.
#' @param DNA_5mer A character string representing the
#' five-nucleotide sequence surrounding the candidate site.
#' @param positive_threshold A numeric value used as the
#' classification threshold for positive prediction.
#'
#' @return A named character vector containing the predicted
#' m6A probability and prediction status.
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'
#' prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
