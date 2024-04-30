# sporeCostsVer2
Calculation of energetic costs of spore formation based on protein abundance and gene expression data. 



The flat files are formatted as follows:
"Description \t data1,data2,data3.."

To parse the file you first split by tabs, then split the second object by commas to get the data array.


#### consumer_resource_spore_model.tsv
- Line 1 = Hours
- Line 2 = Resource concentration
- Line 3 = Cell concentration
- Line 4 = Spore concentration


#### cr_spores_efficiency.png 
- Line 1 = Range of possible spore/cell ATP ratios
- Line 2 = Sporulation efficiencies




#### evo_ratio.tsv
- Line 1 = Range of deletion sizes
- Line 2 = Ratio of rates for lower Ne
- Line 3 = Ratio of rates for higher Ne

The neutral expectation (horizontal line) is 0.2668621700879765


#### evo_ratio_conditional.tsv
- Line 1 = Range of deletion sizes
- Line 2 = Ratio of rates for lower Ne
- Line 3 = Ratio of rates for higher Ne
- Line 4 = Neutral expectation of ratio

