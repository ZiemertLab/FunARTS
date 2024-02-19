#!/usr/bin/env python
# Copyright (C) 2023 Turgut Mesut YILMAZ
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Lab of Nadine Ziemert, Div. of Microbiology/Biotechnology
# Funding by the German Centre for Infection Research (DZIF)
#
# This file is part of FunARTS
# FunARTS is free software. you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version
#
# License: You should have received a copy of the GNU General Public License v3 with FunARTS
# A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>.

import argparse, os

def call_startquery(args):
    reference= "/funarts/reference/ascomycota/"
    container_volume_list = []
    container_argument_list = []
    docker_run_cmd = "docker run -it --rm "
###Input volume - path in container
    input_paths = []
    if "," in args.input:
        input_list = args.input.split(",")
        for i in input_list:
            path,file = os.path.split(i)
            input_paths.append("/funarts/uploads/%s" % file)
            container_volume_list.append("-v %s:/funarts/uploads/" % path)
    else:
        path, file = os.path.split(args.input)
        input_paths.append("/funarts/uploads/%s" % file)
        container_volume_list.append("-v %s:/funarts/uploads/" % path)

 ##Reference volume - path
    if args.refdir != reference:
        path, file = os.path.split(args.refdir)
        container_ref_arg = "/custom_reference/%s" % file
        container_volume_list.append("-v %s:/custom_reference/" % path)
    else:
        container_ref_arg = reference

### Directory to store results
    path, file = os.path.split(args.resultdir)
    container_volume_list.append("-v %s:/funarts/results/" % path)
    container_argument_list.append("-rd /funarts/results/%s" % file)

### Antismash path
    if args.runantismash is True:
        container_argument_list.append("-ras -asp /usr/local/bin/antismash")
    else:
        container_argument_list.append("-asp /usr/local/bin/antismash")

### Custom hmm models list
    if args.hmmdblist is not None:
        path, file = os.path.split(args.hmmdblist)
        container_argument_list.append("-hmms /funarts/custom_files/%s"%file)
        container_volume_list.append("-v %s:/funarts/custom_files/" % path)

### Known Resistance models hmms
    if "kres" in args.options:
        container_argument_list.append("-khmms /funarts/reference/knownresistance.hmm")

### Domains of unknown function hmm file
    if "duf" in args.options:
        container_argument_list.append("-duf /funarts/reference/dufmodels.hmm")

### User supplied core models
    if args.custcorehmms is True:
        path, file = os.path.split(args.custcorehmms)
        container_argument_list.append("-cchmms /funarts/custom_files/%s"%file)
        container_volume_list.append("-v %s:/funarts/custom_files/" % path)

### User supplied resistance models
    if args.customhmms is True:
        path, file = os.path.split(args.customhmms)
        container_argument_list.append("-chmms /funarts/custom_files/%s"%file)
        container_volume_list.append("-v %s:/funarts/custom_files/" % path)

### RNA hmm models to run
    if args.rnahmmdb is True:
        path, file = os.path.split(args.rnahmmdb)
        container_argument_list.append("-rhmm /funarts/custom_files/%s"%file)
        container_volume_list.append("-v %s:/funarts/custom_files/" % path)

### Hmm reporting threshold
    if args.thresh is not None:
        container_argument_list.append("-t %s"%args.thresh)

### Multi processing set - CPU
    if args.multicpu > 1:
        container_argument_list.append("-cpu %s" % args.multicpu)
        docker_run_cmd += "--cpus=%s"%args.multicpu + " "

### Opions for Analysis to run
    container_argument_list.append("-opt %s"%args.options)

### Explicitly specify organism name
    if args.orgname is not None:
        container_argument_list.append("-org %s"%args.orgname)

### Run BiG-Scape Path:
    if args.runbigscape is True:
        container_argument_list.append("-rbsc -bcp /funarts/BiG-SCAPE/bigscape.py")


## Docker Run FunARTS
    docker_run_cmd += " ".join(list(set(container_volume_list)))
    docker_run_cmd += " ziemertlab/funarts:latest "
    docker_run_cmd += ",".join(input_paths) + " " + container_ref_arg + " "
    docker_run_cmd += " ".join(container_argument_list)
    print(docker_run_cmd)
    os.system(docker_run_cmd)


# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Start from genbank file and compare with pre-computed reference for Duplication and Transfers""")
    parser.add_argument("input", help="gbk file to start query. For multi-genome input, put commas without any space between the paths.")
    parser.add_argument("resultdir", help="Directory to store results")
    parser.add_argument("-refdir", help="Directory of precomputed reference files (default='~/ascomycota/')",default="/funarts/reference/ascomycota/")
    parser.add_argument("-hmms","--hmmdblist", help="hmm file, directory, or list of hmm models for core gene id",default=None)
    parser.add_argument("-cchmms","--custcorehmms", help="User supplied core models. hmm file",default=False)
    parser.add_argument("-chmms","--customhmms", help="User supplied resistance models. hmm file",default=False)
    parser.add_argument("-rhmm","--rnahmmdb", help="RNA hmm models to run (default: None)",default=None)
    parser.add_argument("-t","--thresh", help="Hmm reporting threshold. Use global bitscore value or Model specific options: gathering= GA, trusted= TC, noise= NC(default: None)",default=None)
    parser.add_argument("-cpu", "--multicpu", help="Turn on Multi processing set # Cpus (default: Off, 1)", type=int, default=1)
    parser.add_argument("-opt", "--options", help="Analysis to run. expert=Exploration mode, kres=Known resistance, duf=Domain of unknown function (default: None)",default=None)
    parser.add_argument("-org", "--orgname", help="Explicitly specify organism name", default=None)
    parser.add_argument("-ras", "--runantismash", help="Run input file through antismash first", action='store_true', default=False)
    parser.add_argument("-rbsc", "--runbigscape",help="Run antismash results through bigscape", action='store_true', default=False )
    args = parser.parse_args()
    call_startquery(args)
