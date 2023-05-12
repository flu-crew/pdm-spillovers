## A(H1N1)pdm09 human/swine spillover and evolution analysis 2009-2021 ##
Here we provide scripts and describe how to reproduce the computational analysis in <br/>
*Markin, A., Zanella, G.C., Arendsee, Z.W., Zhang, J., Krueger, K.M., Gauger, P.C., Vincent Baker, A.L.
and Anderson, T.K., 2022. **Reverse-zoonoses of 2009 H1N1 pandemic influenza A viruses and evolution in United States swine results in viruses with zoonotic potential**. bioRxiv 2022.12.15.520479. https://doi.org/10.1101/2022.12.15.520479*

### Project structure ###
See the respective README.txt files within the folders for details/instructions.
- [detection-stats](detection-stats/): R code and data for analyzing swine/human pdm09 detection statistics between 2010 and 2021.
- [HI-data](HI-data/): Hemagglutinin Inhibition assays results and R code for Figure 7.
- [blast](blast/): Whole-genome sequences analysis of putative (swine-to-human) variants from 2020 and 2021.
- [aadiff](aadiff/): Code and HA sequences for Supplemental Table S2.
- [spillovers-by-season](spillovers-by-season/): R code for persistence of spillovers over time analysis (Fg 3A).
- [usda-h1-surveillance](usda-h1-surveillance/): R code for swine H1pdm09s to H1-all ratio visulization by season (Fig 3B).
- [spillovers-geo-plot](spillovers-geo-plot/): Python code for US map plotting, color-coded by frequency of human-to-swine spillovers (Fig 3C).
- [host-parsimony](host-parsimony/): Maximum parsimony analysis of host transitions (spillovers) with Python.


Other folders and code are specific to the genetic analysis of HA pdm09 sequences (see below).

### Genetic and phylodynamic analyses of HA pdm09 sequences ###

#### Swine IAV sequence data from the ISU Veterinary Diagnostic Laboratory ####
Unpublished sequences from the Iowa State University Veterinary Diagnostic Laboratory were submitted to NCBI Genbank, Accesssions OQ179020 - OQ179505.

9 HA genes did not meet quality control criteria for NCBI Genbank and are available at AgDataCommons, (dataset) Markin, Alexey; Ciacci Zanella, Giovana; Arendsee, Zebulun W.; Zhang, Jianqiang; Krueger, Karen M.; Gauger, Phillip C.; Vincent Baker, Amy L.; Anderson, Tavis K. (2023). Data from: Reverse-zoonoses of 2009 H1N1 pandemic influenza A viruses and evolution in United States swine results in viruses with zoonotic potential. Ag Data Commons. https://doi.org/10.15482/USDA.ADC/1528393.

#### Environment setup ####
Executing the scripts requires Python 3.7 or higher. For the full list of the required python packages see [requirements.txt](requirements.txt). To install these packages we recommend using virtualenv and run `pip install -r requirements.txt` within the virtual environment. If you never used virtualenv, below we provide instructions on how to set up a virtual environment.

Additionally, we assume that IQ-Tree is installed and can be executed by calling `iqtree`.


#### Download human HA genes ####
We cannot publicly upload the genetic sequences from GISAID, but we provide a list of accessions. See [pdm09-data/README.txt](pdm09-data/README.txt) for the instructions on how to obtain the sequences.

Then run the following commands:
```
cd human-swine-analysis/
python ./build_working_dataset.py  # Align sequences and construct a similarity matrix.
cd ../
```

#### Treetime correction ####
Treetime v0.8.4 may fail on large trees due to the insufficient recursion depth allowed in Python3 by default.
To fix this problem, please run `which treetime` after installing phylo-treetime==0.8.4 and modify the treetime executable file as follows:
```
<...>
if __name__ == '__main__':
    sys.setrecursionlimit(15000)

    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
```
Note the `sys.setrecursionlimit(15000)` line, which makes sure that treetime will process our dataset with around 14,000 strains.

#### Run phylodynamic anlyses ####
We run the phylogenetic inference + phylodynamic analyses 20 times independently. Therefore, it is best to be run on a server in parallel, e.g., as follows:
```
for i in {1..20}
do
	bash run-treetime_multi.sh $i > log_pdm09_treetime_run${i}.log &
done
wait
```

#### Analyze spillover dynamics ####
After executing the phylogenetic pipeline, run the [spillover_analysis.py](human-swine-analysis/spillover_analysis.py) python code to compute spillover statistics (across all 20 replicates) as well as replicate the machine learning analysis of the 2020-21 swine pdm09 HA sequences. The machine learning analysis will generate Figure 5.

#### Aggregate spillover statistics across 20 replicates ####
Please run the R code in [log-combine.R](trees/log_combine.R) to reproduce the spillover statistics presented in the manuscript. Note that the previous steps are not required to run this code, since all necessarily log files are already included in this project.

### Setting up a virtual environment ###
To set up a virtual environment in Python, first make sure `virtualenv` is installed:
```
pip install virtualenv
```
Then run `virtualenv pdm-venv` to create a virtual envirnoment and activate it using `source pdm-venv/bin/activate` (assuming you are on Mac/Linux). Once you are done using the virtual environment run `deactivate`.
