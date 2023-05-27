# Workflow for peptide and chemical analysis of LC-MS/MS queries (spectra)
This workflow is for classifying and extracting spectral informaton from LC-MS/MS queries from Mascot report.
In particular, those queries which do not result in an accepted hit are of interest.

![Data processing workflow](Figure_S2_Data_processing_workflow_Github.pdf)

Figure: Overview data processing workflow


## Step 1: Export "significant matches" and "Display non-significant matches" of Mascot Summary Reports for MS/MS

### 1.A) Export Protein-families with only significant matches
1. Pick FDR 1% recalculate, note significance threshold at homology
2. Go to Export Search Results
3. Choose Export format .csv
4. Check that significance threshold at homology is the same to achieve FDR 1%
5. Pick FDR type distinct PSMs
6. Leave box unticked “Display non-significant matches”
…
7. Do not “include sub-set protein hits” by typing in zero
8. Tick “Group protein-families”
9. “Protein Hit Information” pick Score, Description, Mass (Da), Number of queries matched, emPAI
10. “Peptide Match Information” tick all available boxes, i.e. “Unassigned queries”
11. “Query Level Information” tick all available boxes

### 1.B) Export Protein-families with display non-significant matches (Anchored)
re-do export same as above except in 6) tick box for “Display non-significant matches”

## Step 2) Mascot report.csv file sorting and filtering

### 2.A) Sort and filter Protein-families with only significant matches
1. Open Mascot file .csv which contain significant matches export in table calculation tool (e.g. Libre calc, Excel)
2. Remove duplicate pep_query leaving only unique queries
3. Apply filter across headers of section “Protein hits”
4. Untick protein-families to contaminants such as keratin, human/primate-like, bovine, trypsin, ovalbumin (BLASTp hit peptide sequence)
5. Copy query number into new sheet with leading column entry for "Assigned_QC" (The pep_seq entries are used for further analyses e.g. in Unipept)
6. In this step, select contaminants instead and copy and paste query number with leading column entry for "Assigned_contams"
7. Apply filter across headers of section "Peptide matches not assigned to protein hits" - further down in Mascot report
8. Untick "empty" in pep_seq and copy the query numbers with leading column entry for "Unassigned" - these are proto-unassigned and anchored will need to be removed in later step in R
9. Select "empty" in pep_seq and copy the query numbers with leading colums entry for "Unmatched" 
10. In new sheet add new columns pep_query, moverz, charge, StringTitle, Retention time range, TotalIonsIntensity, NumVals, StringIons1
11. Entries for new information for each query can be retrieved from section "Queries" further below in Mascot report. Empty rows in Queries are removed. E.g. VLOOKUP may be used using query as unique ID
12. Move to part 2.B

### 2.B) Sort and filter Protein-families with display non-significant matches (Anchored)
1. Open Mascot file .csv which contain display non-significant matches export in table calculation tool
2. Make note of Significance threshold (p_signif) value
3. Remove duplicate pep_query leaving only unique queries
4. Apply filter across headers of section “Protein hits”
5. Untick protein-families to contaminants such as keratin, human/primate-like, bovine, trypsin, ovalbumin (BLASTp hit peptide sequence)
6. Filter for queries with pep_expect > significance threshold
7. Copy query number into new sheet with leading column entry for "Anchored_QC"
8. In this step, select protein family contaminants instead 
9. Filter for queries with pep_expect > significance threshold
10. Copy and paste query number with leading column entry for "Anchored_contams"
11. Loop up information on queries as in step 11 above
12. Save file as .csv

## Step 3) Remove Anchored IDs from Proto-Unassigned in R script "Remove.R"
1. Execute "Remove.R" for each datafile and save as x_prep.csv prepared files
The file now contains Assigned_QC, Assigned_contams, Anchored_QC, Anchored_contams, Unassigned, Unmatched
Summary statistics are done in R script "Summary_stats.R"

## Step 4) Extract MS2 fragments
1. Execute "ExtractMS2.R" for each datafile and save as .Rda files

## Step 5) Prepare MS2 ions as m/z list for Mass defect (Modulo 1) analysis
1. Execute "PrepMS2modulo.R" for each datafile and save .csv files

## Step 6) Mass defect analysis (Modulo 1)
1. Execute scripts forked to site as detailed in manuscript
