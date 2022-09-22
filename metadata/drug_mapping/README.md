# Normalization of drug names to ChEMBL

## Columns 
- UNNAMED: meaningless 
- Index: meaningless
- Name: from original data
- Dataset: from original data
- Found Names: List (of list) of names suggested by normalization software  
- Found Name IDs: List (of list) of IDs associated to each name in `Found Names` 
- Match Type: type of match reported by normalization software. Can be: `exact`, `substring`, `partial`, `none`
- Edit Distance: Measure of dis-similarity between `Name` and `Found Names` 
- Query Time: well... 
- ChEMBL: **Final gold ChEMBL ID** to be assigned to `Name` 
- Manual normalization: comment if normalization was done manually, e.g. `typo` for names that probably have typos 

## Files 

- `drugnames_raw.csv` : original data

### TODO:

- `missing_normalized_drug_names_NONE`: no match whatsoever could be found

### READY:

Divided by type of match of normalization software:
- `normalized_drug_names_EXACT.csv` : assign ChEMBL ID directly 
- `normalized_drug_names_SUBSTRING.csv` : these entries were in the form `NAME (CODE_NAME_1, CODE_NAME_2, ...)`. With a little processing I could assign a ChEMBL ID 
- `manual_normalized_drug_names_PARTIAL.csv` : partial match, due to typo, add manually by searching amended name in ChEMBL web page
- `manual_normalized_drug_names_NONE.csv` : no match found, missing entry in gold file: add manually by searching amended name in ChEMBL web page


