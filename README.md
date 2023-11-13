# Fungal bioActive compound Resistant Target Seeker (FunARTS) Overview

FunARTS is a specific and efficient mining tool for the identification of fungal bioactive compounds with interesting and novel targets.
The aims of FunARTS are to (i) automate the process of target-directed (also called resistance-guided) genome mining in fungi, 
(ii) screen for potential novel bioactive compound targets, and (iii) prioritize putative Biosynthetic Gene Clusters (BGCs) for their subsequent characterization.

FunARTS can be installed locally, or you can use the free public webserver located at https://funarts.ziemertlab.com

See https://github.com/ziemertlab/funartswebapp for a guide on installing the webserver independently.


# Installation of FunARTS

There are three options for installing FunARTS:

- Using Docker Images 
- Using Anaconda/Miniconda
- Manual Installation for Linux/Ubuntu

## 1- Using Docker Image:

- Firstly, if you don't have Docker, you should install the Docker engine on your computer. Please check out the latest version of Docker 
[on the official website.](https://docs.docker.com/get-docker/)

- To run FunARTS Image, you should download the "docker_run_funarts.py" file from the command line or from the repository using a web browser.
```bash
    mkdir FunARTSdocker && cd FunARTSdocker
    wget https://github.com/ziemertlab/funarts/raw/master/docker_run_funarts.py
```
**Note:** Python 3.x is needed to run "docker_run_funarts.py".

- FunARTS Image include only Ascomycota reference set. If you need other reference sets, please download (~450MB) and unzip all of them.
```bash
    mkdir FunARTSdocker && cd FunARTSdocker
    wget https://funarts.ziemertlab.com/static/zip_refsets/reference.zip
    unzip reference.zip 
```
- Enter the required arguments and run the script
```bash
    python docker_run_funarts.py [-h] [input] [resultdir] [-optional_arguments]
```
- You can see the other details [on Docker Hub](https://hub.docker.com/r/ziemertlab/funarts)

## 2- Using Anaconda/Miniconda:
We recommend [Anaconda3/Miniconda3](https://docs.anaconda.com/free/anaconda/install/index.html) (with python >=3.8) and 
it is necessery for the [conda](https://docs.conda.io/en/latest/index.html) package manager.

- Clone/Download the repository (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/funarts
```
- Enter the funarts folder:
```bash
    cd funarts
```
- download (~450MB) and unzip reference sets:
```bash
    wget https://funarts.ziemertlab.com/static/zip_refsets/reference.zip
    unzip reference.zip 
```
- Create a new environment and install all the packages using the environment.yml file with conda:
```bash
    conda env create -f environment.yml
```
- Activate funarts environment and run FunARTS (See [Usage](https://github.com/ZiemertLab/FunARTS/blob/master/README.md#usage) for more):
```bash
    conda activate funarts
```
```bash
    python funartspipeline.py [-h] [input] [refdir] [-optional_arguments]
```

## 3- Manual Installation for Linux/Ubuntu:
The analysis server will start a local antiSMASH job if cluster annotation is not already provided as input. We recommend antiSMASH version >= 6.0.1.
See [antiSMASH](https://docs.antismash.secondarymetabolites.org/install/) for installation instructions.

**Note:** Python version 3.8 or higher is recommended.

- Clone/Download the repository (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/funarts
```
- Enter the funarts folder:
```bash
    cd funarts
```
- download (~450MB) and unzip reference sets:
```bash
    wget https://funarts.ziemertlab.com/static/zip_refsets/reference.zip
    unzip reference.zip 
```
- Install required libraries and applications (root / sudo required):
```bash
    apt-get update
    apt-get install -y hmmer2 hmmer diamond-aligner fasttree prodigal ncbi-blast+ muscle
    pip install -r requirements.txt
```
- Run FunARTS (See [Usage](https://github.com/ZiemertLab/FunARTS/blob/master/README.md#usage) for more):
```bash
    python funartspipeline.py [-h] [input] [refdir] [-optional_arguments]
```


## Comparing the results of multi-genome analysis:

The BiG-SCAPE algorithm is used to compare the results of multi-genome analysis. 
All clustered BGCs from antiSMASH results are analyzed to determine BGC similarity. 
The BiG-SCAPE algorithm generates sequence similarity networks of BGCs and classifies them into gene cluster families (GCFs).

To install [the BiG-SCAPE](https://bigscape-corason.secondarymetabolites.org/index.html), please see https://github.com/medema-group/BiG-SCAPE/wiki/installation

**Note:** Make sure that the Pfam database is in the same folder as bigscape.py

# Running FunARTS
FunARTS uses a webserver to queue jobs to the analysis pipeline. Details on webserver usage can be found at: https://funarts.ziemertlab.com/help 

Alternatively jobs can be run directly using the funartspipeline.py script (see -h for options).

````
usage: funartspipeline.py [-h] [-hmms HMMDBLIST] [-khmms KNOWNHMMS] [-duf DUFHMMS] [-cchmms CUSTCOREHMMS] [-chmms CUSTOMHMMS]
                          [-rhmm RNAHMMDB] [-t THRESH] [-td TEMPDIR] [-rd RESULTDIR] [-cpu MULTICPU] [-opt OPTIONS]
                          [-org ORGNAME] [-ras] [-asp ANTISMASHPATH] [-bcp BIGSCAPEPATH] [-rbsc]
                          input refdir

Start from genbank file and compare with pre-computed reference for Duplication and Transfers

positional arguments:
  input                 genbank file to start query - 
  refdir                Directory of precomputed reference files

optional arguments:
  -h, --help            show this help message and exit
  -hmms HMMDBLIST, --hmmdblist HMMDBLIST
                        hmm file, directory, or list of hmm models for core gene id
  -khmms KNOWNHMMS, --knownhmms KNOWNHMMS
                        Resistance models hmm file
  -duf DUFHMMS, --dufhmms DUFHMMS
                        Domains of unknown function hmm file
  -cchmms CUSTCOREHMMS, --custcorehmms CUSTCOREHMMS
                        User supplied core models. hmm file
  -chmms CUSTOMHMMS, --customhmms CUSTOMHMMS
                        User supplied resistance models. hmm file
  -rhmm RNAHMMDB, --rnahmmdb RNAHMMDB
                        RNA hmm models to run (default: None)
  -t THRESH, --thresh THRESH
                        Hmm reporting threshold. Use global bitscore value or Model specific options: gathering= GA, trusted=TC, noise= NC(default: none)
  -td TEMPDIR, --tempdir TEMPDIR
                        Directory to create unique results folder
  -rd RESULTDIR, --resultdir RESULTDIR
                        Directory to store results
  -cpu MULTICPU, --multicpu MULTICPU
                        Turn on Multi processing set # Cpus (default: Off, 1)
  -opt OPTIONS, --options OPTIONS
                        Analysis to run. expert=Exploration mode, kres=Known resistance, duf=Domain of unknown function (default:None)
  -org ORGNAME, --orgname ORGNAME
                        Explicitly specify organism name
  -ras, --runantismash  Run input file through antismash first
  -asp ANTISMASHPATH, --antismashpath ANTISMASHPATH
                        Location of the executable file of antismash or location of antismash 'run_antismash.py' script
  -bcp BIGSCAPEPATH, --bigscapepath BIGSCAPEPATH
                        Location of bigscape 'bigscape.py' script
  -rbsc, --runbigscape  Run antismash results through bigscape
````

# Usage 

- For basic run with positional arguments;
````
    python funartspipeline.py /PATH/input_genome.gbk /PATH/funarts/reference/ascomaycota
````

- To save all output data files: `-rd`, `--resultdir`
````
    python funartspipeline.py /PATH/input_genome.gbk /PATH/funarts/reference/ascomaycota -rd /PATH/result_folder
````

- To use antiSMASH: `-asp`, `--antismashpath` and to run antiSMASH: `-ras`, `--runantismash`
````
    python funartspipeline.py /PATH/input_genome.gbk /PATH/funarts/reference/ascomaycota -asp /PATH/antismash -ras -rd /PATH/result_folder
````

- If there is an exsiting antiSMASH job, .json files of antiSMASH results are available fo FunARTS: `-asp`, `--antismashpath`
````
    python funartspipeline.py /PATH/antismash_result.json /PATH/funarts/reference/ascomaycota -asp /PATH/antismash -rd /PATH/result_folder
````

- To run FunARTS with exploration mode, please use `-opt`, `--options` parameter;
````
    python funartspipeline.py /PATH/input_genome.gbk /PATH/funarts/reference/ascomaycota -asp /PATH/antismash -ras -opt 'expert' 
````

- To identify known resistance, please use `-khmms`, `--knownhmms` and `-opt`, `--options` parameters;
````
    python funartspipeline.py /PATH/input_genome.gbk /PATH/funarts/reference/ascomaycota -asp /PATH/antismash -ras -khmms /PATH/funarts/reference/knownresistance.hmm -opt 'kres'
````

- To identify domain of unknown function(DUF), please use `-duf`, `--dufhmms` and `-opt`, `--options` parameters;
````
    python funartspipeline.py /PATH/input_genome.gbk /PATH/funarts/reference/ascomaycota -asp /PATH/antismash -ras -khmms /PATH/funarts/reference/dufmodels.hmm -opt 'duf'
````

- For multi-genome input, it is enough to put commas without any space between the paths of genome files;
````
    python funartspipeline.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk,/PATH/input_genome3.gbk /PATH/funarts/reference/ascomaycota -rd /PATH/result_folder
````

- To run the BiG-SCAPE algorithms, please use `-bcp`, `--bigscapepath` and `-rbsc`, `--runbigscape`
````
    python funartspipeline.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk /PATH/funarts/reference/ascomaycota -bcp /PATH/BiG-SCAPE_1.1.5/bigscape.py -rbsc -rd /PATH/result_folder
````

# Support
If you have any issues please feel free to contact us at arts-support@ziemertlab.com

# Licence
This software is licenced under the GPLv3. See LICENCE.txt for details.

# Publication
If you found FunARTS to be helpful, please [cite us](https://doi.org/10.1093/nar/gkad386):

Yılmaz, T. M., Mungan, M. D., Berasategui, A., & Ziemert, N. (2023). FunARTS, the Fungal bioActive compound Resistant Target Seeker, an exploration engine for target-directed genome mining in fungi. Nucleic Acids Research
