# StarvationTraits: Microbial Starvation Traits

***Sub-Project of:***
**Dimensions: Collaborative Research: Microbial seed banks: processes and patterns of dormancy-driven biodiversity**

This repository contains open-source code, data, & text files for an REU project related to a National Science Foundation Dimensions of Biodiversity grant (#1442246) awarded to Dr.'s Jay Lennon and Ken Locey of Indiana University (Bloomington) and to Dr. Stuart Jones of the University of Notre Dame (#...).

For information regarding the main project, please visit: 

1. *NSF*: http://www.nsf.gov/awardsearch/showAward?AWD_ID=1442246&HistoricalAwards=false
2. *Dimensions Repo*: https://github.com/LennonLab/Dimensions

## Project Goals

* **Aim 1.)** 

* **Aim 2.)** 

* **Aim 3.)** 

### Repo Contents

* **analyses:**

* **bin:** 
	* *PreSensInteractiveRegression.R*: An R script written by Mario Muscarella (Indiana University) containing functions used in the analysis of PreSens oxygen respiration data.
	* *ReadSynergy.R*: An R script written by Mario Muscarella (Indiana University) containing functions used to import data from BioTek Synergy MX plate reader.
	* *ReadAB1.py*: A Python script written by Mario Muscarella (Indiana University) used to read raw ABI sequencer files and export *.fasta and *.qual files.
	* *TrimMovingAverage.py*: A Python script written by Mario Muscarella (Indiana University) used to trim raw sequence files (*.fasta) based on Phred quality scores.
	* *MergeSeqs.py*: A Python script written by Mario Muscarella (Indiana University) used to merge individual *.fasta files into a multi fasta file.

* **data:**
	* Sequences: Contains raw and processed Sanger sequencing files used for colony identifications.
	* Respiration: Contains raw output files (.txt) from PreSens SensorDish Oxygen system
	* GrowthCurves: Contains raw output files (.txt) from BioTek Synergy MX plate reader. 
	* DeathCurves: Contains plate count data from long term starvation experiments

* **figures:**






## Contributors

Rachel Ferrill: REU Research, Undergraduate, Department of Biology, Transylvania University. 

[Mario Muscarella](http://mmuscarella.github.io/): Ph.D. candidate in the [Lennon Lab](http://www.indiana.edu/~microbes/people.php)

[Dr. Jay Lennon](http://www.indiana.edu/~microbes/people.php): Principle Investigator, Associate Professor, Department of Biology, Indiana University, Bloomington. Head of the [Lennon Lab](http://www.indiana.edu/~microbes/people.php).
