# CassetteCrafter
CassetteCrafter automates the creation of DNA sequence libraries containing all combinations of desired mutations, formatted for Golden Gate Assembly workflows.  


This software generates a library of DNA sequence cassettes that contain all possible combinations of specified mutations, with the necessary recognition sequences for Golden Gate Assembly. A user provided starting sequence is used to generate a library of mutant sequences. These sequences are then split into cassettes with optimized overhangs to fit oligo length requirements. The user chooses which Type IIs enzyme will be used for plasmid construction, and recognition sites are added onto the ends of each cassette. If neccessary, additional bases are added to meet length requirements.

For more information on cloning terminology and techniques, see [Traditional Cloning Basics from ThermoFisher](https://www.thermofisher.com/us/en/home/life-science/cloning/cloning-learning-center/invitrogen-school-of-molecular-biology/molecular-cloning/cloning/traditional-cloning-basics.html) 

## Installation and Setup

### Requirements

- Python 3.x
- Flask 
- BioPython
- Pandas
- Optional: Pytest
  
### Setup

1. Clone the repository:  
   `git clone https://github.com/jennastanislaw/cassettecrafter.git`

2. In the cassettecrafter directory, create the conda environment containing packages needed to run the software:   
 `conda env create -f environment.yml`

3. Activate the environment
   
4. Run the App:  
   From the cassettecrafter directory, run
   `python3 GUI/app_script.py`

5. Open Web App:  
   Open a web browser and go to http://127.0.0.1:5000

### For command line usage, see [Command Line Instructions:](./commandlineinstructions.md)

## Demo Instructions
1. Upload the gene sequence file, demo file is located at `<cassette/crafter/parent/directory>/cassettecrafter/test_data/test_seq_single.csv`  
2. Upload the mutation file, demo file is located at `<cassette/crafter/parent/directory>/cassettecrafter/test_data/demo_mutation_list.csv`
3. Choose your enzyme, for the demo we use BbsI.
4. Specify your minimum oligo size, for the demo we use 25.
5. Specify your maximum oligo size, for the demo we use 100.  
- Note: If your minimum and maximum oligo sizes are such that a large number of split sites (>7) are required, determining the locations of these split sites can take an extended ammount of time. 

## Golden Gate Assembly 

Golden Gate assembly is a cloning technique used to do scar-less DNA sequence assembly using Type IIs restriction enzymes. These enzymes have cut sites outside of their recognition sequences, allowing for a final construct that does not contain enzyme recognition sites. 

![Type IIP vs IIs](./Restriction_enzyme_cuts.png)  
![Golden Gate Assembly](./Golden_Gate_multi-insert_diagram.png) 
(Images via [Snapgene](https://www.snapgene.com/guides/golden-gate-assembly))

More information on Golden Gate assembly can be found here:  
[Golden Gate Assembly](https://www.snapgene.com/guides/golden-gate-assembly)  
[Plasmids 101: Golden Gate Cloning](https://blog.addgene.org/plasmids-101-golden-gate-cloning) 







