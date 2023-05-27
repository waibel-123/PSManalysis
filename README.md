# Export Mascot Summary Reports for MS/MS

## 1) Export Protein-families with only significant matches
1. Pick FDR 1% recalculate, note significance threshold at homology
2. Go to Export Search Results
3. Choose Export format .csv
4. Check that significance threshold at homology is the same to achieve FDR 1%
5. Pick FDR type distinct PSMs
6. Leave box unticked “Display non-significant matches”
…
7. Do not “include sub-set protein hits” by setting to zero
8. Tick “Group protein-families”
9. “Protein Hit Information” pick Score, Description, Mass (Da), Number of queries matched, emPAI
10. “Peptide Match Information” tick all available boxes, i.e. “Unassigned queries”
11. “Query Level Information” tick all available boxes

## 2) Export Protein-families with display non-significant matches (Anchored)
re-do export same as above except 6) Tick “Display non-significant matches”

