getTCR <- function(input_dir,samples) {
  #This function will return a dataframe with unique cell barcodes, all TRA/TRB chain clonotype information and MAIT/iNKT evidence in a single row.
  ##Load filtered_contig_annotations.csv files for all samples
  df_l <- lapply(samples,function(x) {read.csv(paste0(input_dir,"/",x,"/outs/filtered_contig_annotations.csv"))})
  ##Load clonotypes.csv files for all samples
  cl_l <- lapply(samples,function(x) {read.csv(paste0(input_dir,"/",x,"/outs/clonotypes.csv"))})
  ##Parse filtered_contig_annotations.csv files for each sample
  for (i in seq_along(df_l)){
    df <- df_l[[i]]
    df <- df[(df$chain %in% c("TRA","TRB")) & (df$productive %in% c("true","TRUE","True")),]
    ##For every chain per cell barcode, keep the row with the most UMIs (remove multiple rows per chain)
    ###Sorting by UMIs rather than reads.
    df <- df %>% group_by(barcode,chain) %>% arrange(desc(umis)) %>% slice_head(n=1)
    ##For every cell barcode, keep the two chains with the most UMIs (remove multiple chains per barcode)
    ###Using slice_head instead of top_n to handle ties in the number of UMIs per chain.
    df <- df %>% group_by(barcode) %>% arrange(desc(umis)) %>% slice_head(n=2)
    ##Add columns with TRA/TRB chain gene-level clonotype
    df <- df %>% 
      mutate(TRAct = ifelse(chain == "TRA", paste(interaction(v_gene,  j_gene, c_gene)), NA)) %>%
      mutate(TRBct = ifelse(chain == "TRB", paste(interaction(v_gene,  j_gene, c_gene)), NA))
    ##Create dataframe with unique cell barcodes in rows and TRA/TRB chain gene-, amino acid- and nucleotide-level clonotype information 
    out <- parseTCR(df)
    ##Add MAIT/iNKT evidence per cell
    cl <- cl_l[[i]]
    out$MAIT <- cl$mait_evidence[match(as.character(out$clonotypeID), as.character(cl$clonotype_id))]
    out$iNKT <- cl$inkt_evidence[match(as.character(out$clonotypeID), as.character(cl$clonotype_id))]
    out$clonotype_freq <- cl$proportion[match(as.character(out$clonotypeID), as.character(cl$clonotype_id))]
    df_l[[i]] <- out
  }
  names(df_l) <- samples
  return(df_l)
}


parseTCR <- function(df){
  #This function will output a dataframe of unique cell barcodes and their TRA & TRB chain
  #clonotype information in a single row
  ##Create vector of unique cell barcodes in dataframe
  unique_cb <- unique(df$barcode)
  ##Create output dataframe where each cell barcode is matched to heavy and light chain clonotype information in a single row
  out <- data.frame(matrix(nrow=length(unique_cb),ncol=10))
  colnames(out) <- c("barcode","TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene","TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene","clonotypeID")
  out$barcode <- unique_cb
  ##For every unique cell barcode:
  for (i in seq_along(unique_cb)){
    barcode.i <- out$barcode[i]
    location.i <- which(df$barcode == barcode.i)
    ##Write in its clonotype ID
    if(length(location.i)>1){
      out[i,"clonotypeID"] <- df[location.i[1],"raw_clonotype_id"]
    } else {
      out[i,"clonotypeID"] <- df[location.i,"raw_clonotype_id"]
    }
    ##If there are two rows corresponding to it in the input data
    if (length(location.i) == 2){
      ##and the first row is not a TRA chain row
      if (is.na(df[location.i[1],c("TRAct")])) {
        ##then, write in the TRA chain clonotype info from the second row
        out[i,c("TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene")]<-df[location.i[2], c("TRAct","cdr3","cdr3_nt","v_gene")]
        ##and write in the TRB chain clonotype info from the first row
        out[i,c("TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene")]<-df[location.i[1], c("TRBct","cdr3","cdr3_nt","v_gene")]
      } else { 
        ##otherwise, write in the TRA chain clonotype info from the first row
        out[i,c("TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene")]<-df[location.i[1], c("TRAct","cdr3","cdr3_nt","v_gene")]
        ##and write in the TRB chain clonotype info from the second row
        out[i,c("TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene")]<-df[location.i[2], c("TRBct","cdr3","cdr3_nt","v_gene")]
      }
      ##otherwise, if there is only one row for this cell barcode
    } else if (length(location.i) == 1) {
      chain.i <- df$chain[location.i]
      ##and it corresponds to a TRA chain
      if (chain.i == "TRA"){
        ##write in the TRA chain clonotype information from that row
        out[i,c("TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene")]<-df[location.i, c("TRAct","cdr3","cdr3_nt","v_gene")]
        ##otherwise
      } else {
        ##write in the TRB chain clonotype information from that row
        out[i,c("TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene")]<-df[location.i, c("TRBct","cdr3","cdr3_nt","v_gene")]
      }
    }
  }
  return(out)
}
