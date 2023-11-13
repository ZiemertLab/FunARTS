# FunARTS Web Server Overview

This is a sub repository for the Fungal bioActive compound Resistant Target Seeker (FunARTS).

This can be used to view results generated from the public server at https://funarts.ziemertlab.com
or using output from the main analysis pipeline at https://github.com/ziemertlab/funarts

For usage of the web server see https://funarts.ziemertlab.com/help


# Installation of FunARTS Web Server

There are three options for installing FunARTS Web Server:

- Using Docker Images 
- Using Anaconda/Miniconda
- Manual Installation for Linux/Ubuntu

## 1- Using Docker Image:

- Firstly, if you don't have Docker, you should install the Docker engine on your computer. Please check out the latest version of Docker 
[on the official website.](https://docs.docker.com/get-docker/)

- Edit desired paths to run Docker Images:
````
I- Desired Result Directory: /my/path/to/result_folders:/results/
II- Desired Log Directory: /my/path/to/log_file:/run/
````
- Enter the required paths and run the command:

````bash
docker run -it -v /my/path/to/results:/result_folders/ -v /my/path/to/log_file:/run/ -p 5000:5000 funartswebbapp:latest
````

## 2- Using Anaconda/Miniconda:
We recommend [Anaconda3/Miniconda3](https://docs.anaconda.com/free/anaconda/install/index.html) (with python >=3.8) and 
it is necessery for the [conda](https://docs.conda.io/en/latest/index.html) package manager.

- Clone/Download the repository (~15MB) (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/funartswebapp
```
- Enter the funartswebapp folder:
```bash
    cd funartswebapp
```
- Create a new environment and install all the packages using the environment.yml file with conda:
```bash
    conda env create -f environment.yml
```
- Activate funartswebapp environment:
```bash
    conda activate funartswebapp
```
- Edit desired folders in configs (config/funartsapp_default.conf and config/uwsgi.conf)
  (See [Confugiration](https://github.com/ZiemertLab/FunARTS/tree/master/webapp#configuration-of-funarts-web-server) for more):
- Run server (from funartswebapp folder)
  (See [Usage](https://github.com/ZiemertLab/FunARTS/tree/master/webapp#usage-of-funarts-web-server) for more):
```bash
    uwsgi --ini config/uwsgi.conf
```


## 3- Manual Installation for Linux/Ubuntu:

**Note:** Python version 3.8 or higher is recommended.

- Clone/Download the repository (~15MB) (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/funartswebapp
```
- Enter the funartswebapp folder:
```bash
    cd funartswebapp
```
- Install required libraries and applications (root / sudo required):
```bash
    pip install -r requirements.txt
```
- Edit desired folders in configs (config/funartsapp_default.conf and config/uwsgi.conf)
  (See [Confugiration](https://github.com/ZiemertLab/FunARTS/tree/master/webapp#configuration-of-funarts-web-server) for more):
- Run server (from funartswebapp folder)
  (See [Usage](https://github.com/ZiemertLab/FunARTS/tree/master/webapp#usage-of-funarts-web-server) for more):
```bash
    uwsgi --ini config/uwsgi.conf
```


# Configuration of FunARTS Web Server
- Edit desired folders in configs (config/funartsapp_default.conf and config/uwsgi.conf) and write your working directories instead of "~PATH":
````
EXAMPLE:
config/funartsapp_default.conf:
        ...
        UPLOAD_FOLDER = "~PATH/uploads"
        RESULTS_FOLDER = "~PATH/results"
        ARCHIVE_FOLDER = "~PATH/archive"
        ...
config/uwsgi.conf:
        ...
        logto = ~PATH/funartswebapp.log
        stats = ~PATH/uwsgi.stats.sock
        touch-reload = ~PATH/uwsgi.reload
        pidfile = ~PATH/uwsgi.pid
        ...
````

# Usage of FunARTS Web Server

- Run server (from funartswebapp folder):
```bash
    uwsgi --ini config/uwsgi.conf
```
**Note:** It may need to run "redis-server" on the terminal.

- Click your local server on your browser:

```
    http://127.0.0.1:5000/
```

**Note:** The link may differ according to your configuration.

- Click result page to view FunARTS results:
```
    http://127.0.0.1:5000/results
```
- Enter the result file name and click "View Report":

**Note:** Make sure your result folders are in the specified results folder path of funartswebapp.conf.
````
config/funartsapp_default.conf:
        ...
        RESULTS_FOLDER = "~PATH/results"
        ...
````
- You can view your FunARTS results now!

### Optional:
 If you do not have a FunARTS result file yet, you can download the sample result from the link below (~85MB). 
After extracting the zip file, follow the relevant steps.

````bash
  wget https://funarts.ziemertlab.com/archive/GCF_001890705.1.zip
  tar -C ~PATH/results ~PATH/-xzvf GCF_001890705.1.zip
````

# Optional - Submission a job using local webserver:
To start the analaysis on local webserver, please see https://github.com/ziemertlab/funarts and install FunARTS.

- Firstly, run server (from funartswebapp folder)::
```bash
    uwsgi --ini config/uwsgi.conf
```
**Note:** It may need to run "redis-server" on the terminal.

- To submit an input file, run "runjobs.py" on the terminal:
```bash
    cd funarts
    python runjobs.py run -pid /tmp/runjobs.pid
```
- Local webserver is ready to analyse!!



# Support
If you have any issues please feel free to contact us at arts-support@ziemertlab.com

# Licence
This software is licenced under the GPLv3. See LICENCE.txt for details.

# Publication
If you found FunARTS to be helpful, please [cite us](https://doi.org/10.1093/nar/gkad386):

Yılmaz, T. M., Mungan, M. D., Berasategui, A., & Ziemert, N. (2023). FunARTS, the Fungal bioActive compound Resistant Target Seeker, an exploration engine for target-directed genome mining in fungi. Nucleic Acids Research
