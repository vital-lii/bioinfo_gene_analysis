## 1. Data Processing

* Use procd.py to process raw data (such as GSE189792_Processed_data.xlsx), retain specified rows, filter zero-value rows
* Manually remove column names containing '1' (these are results from original experiment corrections/processing, not original measurement data) cleaned_data.xlsx

## 2. Secondary Data Processing for Group Mean Comparison

* cal_mean.py calculates inter-group means and generates five CSVs (clea.csv down up specific_gf/scfa). Can be specified using command line arguments -f and -p. Currently defaults to keeping only log2fc=1 and p≤0.05. Output files can be directly used for mapping, sequence number 5.

## 3. Gene Pathway Mapping

Use clear.csv directly with calc_means.py, outputs five files including up and down-regulated genes, but without differential genes filtering and deduplication, resulting in cle4.csv (already includes P-value filtered results).

* (Direct mapping) Perform pathway enrichment analysis using enrich.py library. Specify -i input file -o output file -k database -l view how many libraries, library names -p specify p-value -n specify top n rows

## 4. Extracting Differential Genes

* Extract rows with reject_null=TRUE from cle4.csv (using diff_genes.py), resulting in diff_genes.csv and other related files

## 5. Intersection with SAGE Data Differential Genes

* Take intersection of gene names from both datasets, remove duplicates, resulting in non_common.csv, temporarily stored in results/jiaoji folder. Subsequent pathway mapping and interactions.

## 5. Remaining Single Group Data Processing

- Differential genes from both datasets can undergo pathway mapping and interactions independently. Logic same as above.

## 6. Keyword Search

* autophage_search.py for keyword searching, utils is a required folder
  Multiple results presented, see code for search terms.

# Pathway Hypotheses

Hypothesis 1: SCFAs-HDAC-M2 Polarization-Acupuncture Regulation in Asthma
Scientific Basis:
- SCFAs (especially butyrate) are recognized HDAC inhibitors
- HDAC inhibition promotes macrophage polarization to M2 (anti-inflammatory) type
- M2 macrophages play crucial role in asthma inflammation control
- Acupuncture proven to regulate immune cell polarization
Research Value:
- Reveal molecular mechanism of acupuncture regulation through SCFAs-HDAC-M2 axis
- Connect traditional acupuncture with modern epigenetics research
- Provide theoretical foundation for combined acupuncture and sodium butyrate treatment

Hypothesis 2: SCFA-T cells-Acupuncture-Asthma
Scientific Basis:
- SCFAs shown to influence T cell differentiation and function
- SCFAs promote Treg cell production and inhibit Th2, Th17 responses
- Asthma pathology core relates to Th1/Th2 imbalance and Treg dysfunction
- Acupuncture can regulate T cell subgroup balance and function
Research Value:
- Elucidate mechanism of acupuncture regulating lung T cell response via gut SCFAs
- Explore acupuncture-SCFAs-T cell axis role in asthma pathogenesis and treatment
- Provide clues for new immunoregulatory targets

Hypothesis 3: SCFAs-Acupuncture-Pulmonary FFAR Expression
Scientific Basis:
- FFARs (especially FFAR2/GPR43 and FFAR3/GPR41) are key SCFA receptors
- Lungs express these receptors, particularly in epithelial and immune cells
- SCFAs influence cell function and inflammation via FFAR signaling
- G protein-coupled receptor signaling pathways are known acupuncture targets
Research Value:
- Investigate if/how acupuncture regulates pulmonary FFAR expression
- Reveal SCFAs-FFAR signaling role in acupuncture treatment of asthma
- Establish molecular connection between "gut-lung axis" and acupuncture efficacy

# Changes in P-value Calculation in cal_mean.py

* For one group with expression and one without: Higher expression levels (larger gf_mean/np.max(gf_vals)); lower variation (larger 1-gf_cv); results in smaller P-value; final P-value restricted to [1e-10, 0.05]
* For both groups with expression: First calculate relative difference (relative_diff) and average coefficient of variation (avg_cv); if both relative difference and variation < 10%, groups considered not significantly different, p_val = 1.0; otherwise use Welch's t-test (no assumption of equal variances)
* Finally, apply multiple testing correction to all P-values
* Avoids zero variance issues in statistical testing. More reliable significance. —Different calculation methods for specific expression and lung-specific expression, multiple testing correction controls false positive rate

# Handling Gene Name Duplicates

   To avoid gene name duplication in subsequent visualization, use non_duplicate_genes.py, output to nondu.csv. Gene count reduced from 30000+ to 10000+, file usable for heatmap and volcano plot visualization (processed by gene_analysis.py)

# Gene Pathway Mapping

* Use gseapy for enrichment analysis in enrich.py library. Specify -i input file -o output file -k database -l view how many libraries, library names -p specify p-value -n specify top n rows

