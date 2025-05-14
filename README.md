# Call Introns from Sorghum Gene Annotations

This script processes gene annotation data for *Sorghum bicolor* to identify **intronic regions** between annotated UTRs and CDS regions. It uses the IRanges package to calculate introns by subtracting exon-like regions from full mRNA annotations.

## 📦 Requirements

Make sure you have the following R packages installed:

```r
install.packages("data.table")
install.packages("dplyr")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IRanges")
```

## 📥 Input

The script expects a CSV file named:

- `Sbicolor_730_v5.1_primary_exons.csv`

This file should include a gene annotation table with at least these columns:
- `gene`
- `type` (e.g., `"mRNA"`, `"CDS"`, `"five_prime_UTR"`, `"three_prime_UTR"`)
- `chr`
- `start`
- `end`
- `string` (strand)

---

## 🚀 Step-by-Step Instructions

### 1. Load the annotation data

```r
library(data.table)

Sbicolor_730_v5.1_primary_exons <- fread("Sbicolor_730_v5.1_primary_exons.csv", data.table = FALSE)
```

### 2. Clean up the type names for clarity

We rename UTR types to match common biological conventions:

```r
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "five_prime_UTR")] <- "5'UTR"
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "three_prime_UTR")] <- "3'UTR"
```

---

### 3. Load required libraries

```r
library(dplyr)
library(IRanges)
```

---

### 4. Define the function to find introns for a single gene

This function:
- Filters for the mRNA record
- Collects all exon-like parts (UTRs and CDS)
- Subtracts them from the mRNA range to get introns

```r
find_introns <- function(df) {
  mrna_row <- df %>% filter(type == "mRNA")
  if (nrow(mrna_row) != 1) return(NULL)
  
  mrna_range <- IRanges(start = mrna_row$start, end = mrna_row$end)
  
  exon_parts <- df %>% filter(type %in% c("5'UTR", "CDS", "3'UTR"))
  if (nrow(exon_parts) == 0) return(NULL)
  
  exon_ranges <- IRanges(start = exon_parts$start, end = exon_parts$end)
  introns <- setdiff(mrna_range, reduce(exon_ranges))
  
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
```

---

### 5. Apply the function to each gene and combine the results

```r
intron_df <- Sbicolor_730_v5.1_primary_exons %>%
  split(.$gene) %>%
  lapply(find_introns) %>%
  bind_rows()
```

---

### 6. Merge introns with the original annotations

This keeps the original UTR/CDS/mRNA annotations and adds introns:

```r
combined <- bind_rows(Sbicolor_730_v5.1_primary_exons, intron_df) %>%
  arrange(gene, start)
```

---

### 7. View and export results

You can inspect the first few rows:

```r
head(combined)
```

Then export the complete data frame:

```r
fwrite(combined, "Sbicolor_730_v5.1_primary_exons_introns.csv")
```

---

## 📤 Output

- **`Sbicolor_730_v5.1_primary_exons_introns.csv`**: a CSV file that contains the original gene annotations **plus** inferred introns.

---
