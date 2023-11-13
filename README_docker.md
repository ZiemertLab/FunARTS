# Fungal bioActive compound Resistant Target Seeker (FunARTS) Overview

FunARTS is a specific and efficient mining tool for the identification of fungal bioactive compounds with interesting and novel targets.
The aims of FunARTS are to (i) automate the process of target-directed (also called resistance-guided) genome mining in fungi, 
(ii) screen for potential novel bioactive compound targets, and (iii) prioritize putative Biosynthetic Gene Clusters (BGCs) for their subsequent characterization.

FunARTS can be installed locally, or you can use the free public webserver located at https://funarts.ziemertlab.com

See https://github.com/ziemertlab/funarts for a guide on installing the main analysis pipeline independently.

See https://github.com/ziemertlab/funartswebapp for a guide on installing the webserver independently.


## Quick start with docker:

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
## Arguments of docker_run_funarts.py
````
usage: docker_run_funarts.py [-h] [-refdir REFDIR] [-hmms HMMDBLIST] [-cchmms CUSTCOREHMMS] [-chmms CUSTOMHMMS]
                             [-rhmm RNAHMMDB] [-t THRESH] [-cpu MULTICPU] [-opt OPTIONS] [-org ORGNAME] [-ras] [-rbsc]
                             input resultdir

Start from genbank file and compare with pre-computed reference for Duplication and Transfers

positional arguments:
  input                 gbk file to start query. For multi-genome input, put commas without any space between the
                        paths.
  resultdir             Directory to store results

optional arguments:
  -h, --help            show this help message and exit
  -refdir REFDIR        Directory of precomputed reference files (default="~/ascomycota/")
  -hmms HMMDBLIST, --hmmdblist HMMDBLIST
                        hmm file, directory, or list of hmm models for core gene id
  -cchmms CUSTCOREHMMS, --custcorehmms CUSTCOREHMMS
                        User supplied core models. hmm file
  -chmms CUSTOMHMMS, --customhmms CUSTOMHMMS
                        User supplied resistance models. hmm file
  -rhmm RNAHMMDB, --rnahmmdb RNAHMMDB
                        RNA hmm models to run (default: None)
  -t THRESH, --thresh THRESH
                        Hmm reporting threshold. Use global bitscore value or Model specific options: gathering= GA,
                        trusted= TC, noise= NC(default: none)
  -cpu MULTICPU, --multicpu MULTICPU
                        Turn on Multi processing set # Cpus (default: Off, 1)
  -opt OPTIONS, --options OPTIONS
                        Analysis to run. expert=Exploration mode, kres=Known resistance, duf=Domain of unknown function (default: None)
  -org ORGNAME, --orgname ORGNAME
                        Explicitly specify organism name
  -ras, --runantismash  Run input file through antismash first
  -rbsc, --runbigscape  Run antismash results through bigscape
````
## Usage 

- For basic run with positional arguments;
````
    python docker_run_funarts.py /PATH/input_genome.gbk /PATH/result_folder
````

- To use antiSMASH, please use: `-ras`, `--runantismash`;
````
    python docker_run_funarts.py /PATH/input_genome.gbk /PATH/result_folder -ras
````
- To use other reference sets (default: ascomycota);
````
    python docker_run_funarts.py /PATH/input_genome.gbk /PATH/result_folder -refdir /PATH/reference_set_file
````
- For multi-genome input, it is enough to put commas without any space between the paths of genome files;
````
    python docker_run_funarts.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk,/PATH/input_genome3.gbk /PATH/result_folder
````

- It is possible to run FunARTS if there is a .json files of antiSMASH 6.1.1 results;
````
    python docker_run_funarts.py /PATH/antismash_result.json /PATH/result_folder
````

- To run FunARTS analysis options, please use `-opt`, `--options` parameter (default: None);
````
    python docker_run_funarts.py /PATH/input_genome.gbk /PATH/result_folder -ras -opt kres,duf,expert 
````

- To run the BiG-SCAPE algorithms, please use `-rbsc`, `--runbigscape`;
````
    python docker_run_funarts.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk /PATH/result_folder -ras -rbsc
````

# Support
If you have any issues please feel free to contact us at arts-support@ziemertlab.com

# Licence
This software is licenced under the GPLv3. See LICENCE.txt for details.

# Publication
If you found FunARTS to be helpful, please [cite us](https://doi.org/10.1093/nar/gkad386):

Yılmaz, T. M., Mungan, M. D., Berasategui, A., & Ziemert, N. (2023). FunARTS, the Fungal bioActive compound Resistant Target Seeker, an exploration engine for target-directed genome mining in fungi. Nucleic Acids Research
