# call_introns_Sorghum

Sbicolor_730_v5.1_primary_exons <- fread("Sbicolor_730_v5.1_primary_exons.csv", data.table = F)
# change names
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "five_prime_UTR")] <- "5'UTR"
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "three_prime_UTR")] <- "3'UTR"


# Load required packages
library(dplyr)
library(IRanges)

# Your input data frame (make sure it's already loaded as `Sbicolor_730_v5.1_primary_exons`)

# Function to find introns for a single gene
find_introns <- function(df) {
  # Get mRNA annotation
  mrna_row <- df %>% filter(type == "mRNA")
  if (nrow(mrna_row) != 1) return(NULL)  # skip genes with missing or multiple mRNA annotations
  
  # Define mRNA range
  mrna_range <- IRanges(start = mrna_row$start, end = mrna_row$end)
  
  # Define all annotated exon-like regions (UTRs + CDS)
  exon_parts <- df %>% filter(type %in% c("5'UTR", "CDS", "3'UTR"))
  if (nrow(exon_parts) == 0) return(NULL)  # no exon annotations
  
  exon_ranges <- IRanges(start = exon_parts$start, end = exon_parts$end)
  
  # Calculate introns = mRNA minus exon parts
  introns <- setdiff(mrna_range, reduce(exon_ranges))
  
  # Format result
  if (length(introns) == 0) return(NULL)
  
  data.frame(
    gene = mrna_row$gene,
    type = "intron",
    chr = mrna_row$chr,
    start = start(introns),
    end = end(introns),
    string = mrna_row$string,
    stringsAsFactors = FALSE
  )
}

# Apply function to each gene
intron_df <- Sbicolor_730_v5.1_primary_exons %>%
  split(.$gene) %>%
  lapply(find_introns) %>%
  bind_rows()

# Combine introns with original data
combined <- bind_rows(Sbicolor_730_v5.1_primary_exons, intron_df) %>%
  arrange(gene, start)

# View results
head(combined)

fwrite(combined, "Sbicolor_730_v5.1_primary_exons_introns.csv")
