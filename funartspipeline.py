#!/usr/bin/env python
# Copyright (C) 2015,2016 Mohammad Alanjary
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Lab of Nadine Ziemert, Div. of Microbiology/Biotechnology
# Funding by the German Centre for Infection Research (DZIF)
#
# This file is part of ARTS
# ARTS is free software. you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version
#
# License: You should have received a copy of the GNU General Public License v3 with ARTS
# A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>.

import argparse, tempfile, os, shutil, pickle, sys, re, glob, json, subprocess, shlex, time, multiprocessing as mp
import parsegbk, makeseqsql, makehmmsql, seqsql2fa, extractdbgenes, getrnagenes, setlog
from combine_results import *
from subprocess import Popen, PIPE
from threading import Timer
from distutils.dir_util import copy_tree


def makequerydb(infasta,tdir,fname,orgname,idprfx=""):
    makeseqsql.runlist([tdir+infasta],tdir+fname+".db",transonly=True,orgname=orgname)
    seqsql2fa.writefasta(tdir+fname+".db",tdir+fname+".faa",idprfx=idprfx)
    seqsql2fa.writefasta(tdir+fname+".db",tdir+fname+".fna",nuc=True,idprfx=idprfx)

def getcutvals(fname):
    cut_vals={}
    with open(fname,"r") as fil:
        non_decimal = re.compile(r'[^\d.]+')
        for line in fil:
            if line.startswith("NAME"):
                modelname=line.split()[1]
                if modelname not in cut_vals:
                    cut_vals[modelname]={}
            if line.startswith("GA ") or line.startswith("TC ") or line.startswith("NC "):
                cv=line.split()
                cut_vals[modelname][cv[0]]=[float(non_decimal.sub('',cv[1])),float(non_decimal.sub('',cv[2]))]
    return cut_vals

def converthmm(hmmfile,outfile):
    if hmmfile and os.path.exists(hmmfile):
        with open(outfile,"w") as ofil:
            status = subprocess.call(['hmmconvert',str(hmmfile)],stdout=ofil)
        if not status:
            return outfile
    return False

def validatehmms(hmmfile):
    if hmmfile and os.path.exists(hmmfile):
        with open(os.devnull,"w") as devnull:
            status = subprocess.call(['hmmstat',hmmfile],stderr=devnull,stdout=devnull)
        if status==0:
            return True
    return False

def renamehmms(hmmdbs,tdir=None):
    if os.path.exists(hmmdbs):
        tf = tempfile.NamedTemporaryFile(prefix="hmmdb_",suffix=".hmm",dir=tdir,delete=False)

def concathmms(hmmdbs,tdir=None):
    if type(hmmdbs) is list:
        flist=[x for x in hmmdbs if os.path.isfile(x)]
    elif type(hmmdbs) is str and "," in hmmdbs:
        flist=[x for x in hmmdbs.split(",") if os.path.isfile(x)]
    else:
        return False
    if len(flist):
        tf = tempfile.NamedTemporaryFile(prefix="hmmdb_",suffix=".hmm",dir=tdir,delete=False)
        for fname in flist:
            if validatehmms(fname):
                with open(fname,"rb") as fil:
                    shutil.copyfileobj(fil, tf)
        tf.close()
        return tf.name
    else:
        return False

def parsegmatrix(fname,avgs=True):
    temp={}
    with open(fname,"r") as fil:
        for line in fil:
            if line.startswith("#Gene"):
                temp["_orgs"]=line.strip().split("\t")[5:]
            if line.startswith("#Singles"):
                if "," in line:
                    temp["_singles"]=set(line.strip().split("\t")[1].split(","))
            elif not line.startswith("#"):
                x = line.strip().split("\t")
                if avgs:
                    temp[x[0]]=[float(v) for v in x[1:5]]
                    # temp[x[0]].append(float(x[4])/len(x[5:]))
                    temp[x[0]].extend([float(v) for v in x[5:]])
                else:
                    temp[x[0]]=[float(x[1])]
    return temp

def getdupgenes(refgm,genematrix,maxcount=3,minsr=0.2,minubiq=0.2):
    glist = list(set(refgm.keys())&set(genematrix.keys()))
    duplist = []
    for k in glist:
        if k.startswith("_"):
            continue
        # Get counts that are > Avg + STDev of ref genes
        # ubiq = float(refgm[k][-2])/refgm[k][-1] #ubiquity numorgs with gene / total orgs
        genematrix[k] = list(genematrix[k])
        genematrix[k].extend(refgm[k])
        # if genematrix[k][0] > refgm[k][0]+refgm[k][1] and refgm[k][0] <= maxcount and refgm[k][2] <= maxrsd and refgm[k][3] >= minubiq:
        if genematrix[k][0] > refgm[k][0]+refgm[k][1] and refgm[k][0] <= maxcount and refgm[k][2] >= minsr:
            duplist.append(k)
    return duplist

def runhmmer(fname,hmmdb,tdir,cut=None,mcpu=1):
    tf = tempfile.NamedTemporaryFile(prefix="domrslt_",suffix=".domhr",dir=tdir,delete=False)
    tf.close()
    # model_data = getmodeldata(fname=os.path.join(refdir, "model_metadata.json"))
    # print(model_data)
    with open(tf.name+".log","w") as logfil:
        cmd=["hmmsearch", "--domtblout", tf.name, "--noali", "--notextw", hmmdb, fname]
        if mcpu>1:
            cmd[1:1] = ["--cpu", str(mcpu)]
        if cut and cut.lower() in ("ga","tc","nc"):
            cmd[1:1] = ["--cut_%s"%cut.lower()]
        else:
            cmd[1:1] = ["-E","0.01"]
        subprocess.call(cmd, stdout=logfil, stderr=logfil)
    return tf.name

def checkbgcprox(clusters,gtypes=[],glist=[],modeldata={}):
    bgchits = {"cluster-%s"%bgc[0]:{"row":["cluster-%s"%bgc[0],bgc[1],bgc[2],"%s - %s"%(bgc[3],bgc[4])],"hits":[]} for bgc in clusters}
    corebgchits = {}
    corehitnum = 0
    for i,ghits in enumerate(glist):
        if "seqs" in ghits.keys():
            for seqid,gene in ghits['seqs'].items():
                #check all start end overlap gene[2]=scaffold gene[3]=start gene[4]=end bgc[2..3..4]
                for bgc in clusters:
                    clustid = "cluster-%s"%bgc[0]
                    desc = func = ""
                    if gene[2]==bgc[2] and not ((int(gene[4]) < int(bgc[3])) or (int(bgc[4]) < int(gene[3]))):
                        if gene[1] in modeldata.keys():
                            desc = modeldata[gene[1]][0]+": "+modeldata[gene[1]][1]
                            func = modeldata[gene[1]][2]
                        if gene[1] not in corebgchits.keys():
                            corebgchits[gene[1]] = {}
                        if gtypes[i] == "Core":
                            corehitnum+=1
                            corebgchits[gene[1]][seqid] = [str(x) for x in bgc[:5]]+[gene[3],gene[4]]+["Core"]

                        bgchits[clustid]["hits"].append([seqid,gene[1],gene[3],gene[4],gtypes[i],desc,func])
                        break
    return bgchits, corebgchits, corehitnum

#use this if local antismash is not run, turns gbk file into AS results webpage
# def makeantismashresults(aspath,infile,outdir):
#     ashelpdir = os.path.join(os.path.dirname(os.path.realpath(__file__)),"antismash_helper")
#     #This hard code needs to be changed for users that will locally install funarts, maybe make another argument with different antismash path?
#     aspath = "~/anaconda3/envs/funarts/lib/python3.8/site-packages/antismash" ###updated on 21.06.2023
#     # aspath = "/home/turgutyilmaz/funarts_d/antismash-6.0.1/run_antismash.py"
#     if not aspath:
#         log.warning("No antismash installation found")
#     else:
#         # added on 21.06.2023 - start
#         if aspath.endswith("run_antismash.py"):
#             aspath = os.path.split(aspath)[0]
#         # added on 21.06.2023 - end
#         sys.path.insert(0, aspath)
#         try:
#             from Bio import SeqIO
#             from antismash.config import load_config
#             from antismash.config import set_config
#             from antismash.output_modules import html
#             import straight.plugin
#
#             shutil.copytree(os.path.join(ashelpdir,"output_template"),outdir)
#             try:
#                 seqrecs = list(SeqIO.parse(infile,'genbank'))
#             except ValueError:
#                 log.error("Biopython could not parse genbank. please ensure no external sequence references are present")
#                 log.exception("exception")
#             with open(os.path.join(ashelpdir,"defaultopts.pkl"),"r") as fil:
#                 opts = pickle.load(fil)
#
#             opts.full_outputfolder_path = opts.outputfoldername = os.path.realpath(outdir)
#             opts.sep='/'
#             # If Antismash Version is <=4, user needs to choose "antismash.specific_modules" instead of "antismash.modules".
#             # opts.plugins=list(straight.plugin.load('antismash.specific_modules'))
#             opts.plugins = list(straight.plugin.load('antismash.modules')) ###updated on 21.06.2023
#             opts.smcogs=False
#
#             load_config(opts)
#             set_config(opts)
#             html.generator.generate_webpage(seqrecs,opts)
#         except ImportError:
#             log.error("Error importing antismash library, check antismash configuration")

def makeantismashresults_json(aspath,infile,outdir,mcpu):
    # added on 22.06.2023 = if json file is created by antismash v7, pipeline need to use antismash v7 for reuse-results
    # with open(infile, "r") as f:
    #     js = json.load(f)
    # if str(js["version"]).split(".")[0] == "7":
    #     aspath = "/home/turgutyilmaz/anaconda3/envs/as7/bin/antismash"
    #     cmd = ["/home/turgutyilmaz/anaconda3/envs/as7/bin/python3", aspath, "--reuse-results", infile, "--output-dir", outdir, "--cpus", str(mcpu)]
    # else:
    #     cmd = ["python3", aspath, "--reuse-results", infile, "--output-dir", outdir,"--cpus",str(mcpu)]
    # added on 22.06.2023 = if json file is created by antismash v7, pipeline need to use antismash v7 for reuse-results
    cmd = ["python3", aspath, "--reuse-results", infile, "--output-dir", outdir,"--cpus",str(mcpu)]
    if os.path.exists(outdir) != True:
        os.mkdir(outdir)
    if os.path.exists(os.path.join(outdir, "knownclusterblast")) != True:
        os.mkdir(os.path.join(outdir, "knownclusterblast"))
    result_gbk_file = ""
    try:
        aslog_dir= os.path.abspath(os.path.join(outdir, os.pardir))
        with open(aslog_dir+"/aslog.txt",'w') as aslog:
            cmd = " ".join(cmd)
            subprocess.call(cmd,stderr=aslog,stdout=aslog, shell=True)
        result_gbk_file = os.path.join(outdir, os.path.basename(infile)[0:-5]+".gbk")
        #log.info(result_gbk_file)
        # added on 11.04.2023 = json file error by antismash version - antismash 7.0.0beta2
        if os.path.isfile(result_gbk_file):
            log.info(result_gbk_file)
        else:
            log.error("Error running Antismash reuse-results. The error may be caused by the Antismash version. It is recommended to use Antismash v6.")
            return "json_error"
        # added on 11.04.2023 = json file error by antismash version - antismash 7.0.0beta2
    except:
        log.error("Error running Antismash reuse-results")
    return result_gbk_file

def runantismash(aspath,infile,outdir,mcpu):
    asdir = os.path.join(outdir,"antismash")
    if os.path.exists(aspath) and os.path.exists(infile):
        ##hard-coded quick solution for new version##
        python_version = "python3"
        output_command = "--output-dir"
        # with open("/home/ubuntu/disk/antismash/run_antismash.py", "r") as run_antismash_file:
        #     run_antismash_lines = run_antismash_file.readlines()

        #find antismash version
        # if "python3" in run_antismash_lines[0]:
        #     python_version = "python3"
        #     output_command = "--output-dir"
        # else:
        #     python_version = "python"
        #     output_command = "--outputfolder"
        cmd = python_version+ " " + aspath + " " +infile + " --taxon fungi --minimal " +output_command + " " + asdir + " -v --cpus " + str(mcpu)
        # cmd2 = python_version+ " " + aspath + " " +infile + " --taxon fungi --minimal " +output_command + " " + asdir + " -v --cpus " + str(mcpu) + " --genefinding-tool prodigal"  ### changed 25.11.2022
        cmd3 = python_version + " " + aspath + " " + infile + " --taxon fungi --minimal " + output_command + " " + asdir + " -v --cpus " + str(mcpu) + " --genefinding-tool glimmerhmm" ## added 07.11.2022
        with open(outdir+"/aslog.txt",'w') as aslog:
            proc = Popen(shlex.split(cmd) ,stderr=aslog,stdout=aslog)
            #kill antismash after 3 hours
            timeout = 10800
            timer = Timer(timeout, proc.kill)
            try:
                start = time.time()
                timer.start()
                stdout, stderr = proc.communicate()
            finally:
                end = time.time()
                if end - start >= timeout:
                    log.error("Antismash taking too long, please use antismash servers for job and try again with json file")
                timer.cancel()
            log.info("first as done")
        for line in open(outdir+"/aslog.txt", "r"):
            #antismash 5 doesnt accept WGS record without sequences so download them
            if "whole genome shotgun records are not supported" in line:
                log.info("wgs")
                wgs_dir = os.path.join(outdir,"wgs_record")
                os.mkdir(wgs_dir)
                wgs_gbk_path = os.path.join(wgs_dir, os.path.basename(infile))
                with open(outdir + "/wgs_dl_log.txt", 'w') as dllog:
                    #another package with pip (one of Kai's useful tools), but this hard code needs to be changed as well
                    cmd_ncbi = ["/home/ubuntu/.local/bin/ncbi-acc-download --recursive " + os.path.splitext(os.path.split(infile)[-1])[0] + " -o " +  wgs_gbk_path]
                    p = Popen(cmd_ncbi, stdout=dllog, stderr=dllog, shell=True)
                    out, err = p.communicate()
                log.info(os.listdir(wgs_dir))
                wgs_cmd = python_version+ " " + aspath + " " +wgs_gbk_path + " --taxon fungi --minimal " +output_command + " " + asdir + " -v --cpus " + str(mcpu)
                # wgs_cmd2 = python_version+ " " + aspath + " " +wgs_gbk_path + " --taxon fungi --minimal " +output_command + " " + asdir + " -v --cpus " + str(mcpu) + " --genefinding-tool prodigal"
                wgs_cmd2 = python_version+ " " + aspath + " " +wgs_gbk_path + " --taxon fungi --minimal " +output_command + " " + asdir + " -v --cpus " + str(mcpu) + " --genefinding-tool glimmerhmm" ## added 25.11.2022
                with open(outdir+"/aslog.txt",'w') as aslog:
                    proc = Popen(shlex.split(wgs_cmd2) ,stderr=aslog,stdout=aslog)
                    timer = Timer(timeout, proc.kill)
                    try:
                        start = time.time()
                        timer.start()
                        stdout, stderr = proc.communicate()
                    finally:
                        end = time.time()
                        if end - start >= timeout:
                            log.error("Antismash taking too long, please use antismash servers for job and try again with json file")
                        timer.cancel()
                        break
            #If antismash cant find any genes it gives an error
            if "ValueError: Called find_genes, but genefinding disabled" in line:
                with open(outdir+"/aslog.txt",'w') as aslog:
                    # proc = Popen(shlex.split(cmd2), stderr=aslog, stdout=aslog)
                    proc = Popen(shlex.split(cmd3), stderr=aslog, stdout=aslog) ### changed 25.11.2022
                    timer = Timer(timeout, proc.kill)
                    try:
                        start = time.time()
                        timer.start()
                        stdout, stderr = proc.communicate()
                    finally:
                        end = time.time()
                        if end - start >= timeout:
                            log.error("Antismash taking too long, please use antismash servers for job and try again with json file")
                        timer.cancel()
                        break
            ## added 07.11.2022 - start
            #If antismash gives an error for no genes and no genefinding
            if "contains no genes and no genefinding tool specified" in line:
                with open(outdir+"/aslog.txt",'w') as aslog:
                    proc = Popen(shlex.split(cmd3) ,stderr=aslog,stdout=aslog)
                    timer = Timer(timeout, proc.kill)
                    try:
                        start = time.time()
                        timer.start()
                        stdout, stderr = proc.communicate()
                    finally:
                        end = time.time()
                        if end - start >= timeout:
                            log.error("Antismash taking too long, please use antismash servers for job and try again with json file")
                        timer.cancel()
                        break
            ## added 07.11.2022 - end
        if python_version == "python3":
            fil = [x for x in os.listdir(asdir) if x.endswith(".gbk") and "region" not in x and "cluster" not in x]
        else:
            fil = [x for x in os.listdir(asdir) if x.endswith(".final.gbk")]
        if len(fil):
            return os.path.join(asdir,fil[0])
        else:
            log.error("Could not find Antismash output")
            return infile
    else:
        log.error("Could not find Antismash executable and/or input file")
        return infile

def writecoretable(corelist,gm,krlist,modeldata,outfile):
    summary = {"data":[],"seqs":corelist["seqs"]}
    tablestr = "#Core_gene\tDescription\tFunction\tDuplication\tBGC_Proximity\tPhylogeny\tKnown_target\t[Hits_listed]\n"
    funcstats = {}
    if "matrix" in gm.keys():
        for k in gm["matrix"].keys():
            if k not in modeldata:
                log.debug("%s has no model metadata!"%k)
            metadata = modeldata.get(k,["N/A","N/A","N/A","N/A","N/A","N/A","N/A"])
            if metadata[2]:
                if metadata[2] not in funcstats:
                    funcstats[metadata[2]] = 0
                funcstats[metadata[2]] += 1
            temp = {"coregene":k,"TC":metadata[3],"dNdS":metadata[4],"SC":metadata[5],"Ubiq":metadata[6],"description":"%s: %s"%(metadata[0],metadata[1]),"func":metadata[2],"duplicate":"N/A","proximity":"N/A",
                    "phylogeny":"N/A","hits":corelist["core"].get(k,{}).get("seqs",[]),"proxhits":[],"known_hit":"N/A"}
            if "duplicates" in gm.keys():
                if k in gm["duplicates"]:
                    temp["duplicate"] = "Yes"
                else:
                    temp["duplicate"] = "No"
            if "proximity" in gm.keys():
                if k in gm["proximity"].keys():
                    temp["proximity"] = "Yes"
                    temp["proxhits"] = gm["proximity"][k]
                else:
                    temp["proximity"] = "No"
            if "phylogeny" in gm.keys():
                if k in gm["phylogeny"].keys():
                    temp["phylogeny"] = "Yes"
                    temp["phylhits"] = gm["phylogeny"][k]
                else:
                    temp["phylogeny"] = "No"
            if "seqs" in krlist:
                temp["known_hit"] = "No"
                for seqid in temp["hits"]:
                    if seqid in krlist["seqs"].keys():
                        temp["known_hit"] = "Yes"
            summary["data"].append(temp)
            temp["allhits"] = "[" + "; ".join([str(corelist["seqs"][seqid][-2]).replace(";","") for seqid in temp["hits"]]) + "]"
            tablestr += "{coregene}\t{description}\t{func}\t{duplicate}\t{proximity}\t{phylogeny}\t{known_hit}\t{allhits}\n".format(**temp)
        summary["funcstats"]=funcstats
    if outfile:
        with open(outfile+".json","w") as jfil, open(outfile+".tsv","w") as tfil:
            json.dump(summary,jfil,indent=2)
            tfil.write(tablestr)
    else:
        return summary

#reformat duplicate table and write to file
def writeduptable(matrix,dups,corelist,modeldata,outfile):
    # metadata = {k:"%s: %s"%tuple(modeldata.get(k,["n/a","n/a"])[:2]) for k in dups}
    hitlist = {k:[str(corelist["seqs"].get(seqid,["n/a","n/a"])[-2]).replace(";","") for seqid in corelist["core"].get(k,{}).get("seqs",[])] for k in dups}
    dupmat = [[x]+matrix.get(x,["n/a"])[:5]+["["+"; ".join(hitlist.get(x,""))+"]"]+["%s: %s"%tuple(modeldata.get(x,["n/a","n/a"])[:2])] for x in dups]
    with open(outfile+".json","w") as jfil, open(outfile+".tsv","w") as tfil:
        temp = {"data":dupmat[:-1],"hits":hitlist}
        json.dump(temp,jfil,indent=2)
        tfil.write("\t".join(["#Core_gene","Count","Ref_median","Ref_stdev","Ref_RSD","Ref_ubiquity","[Hits_listed]","Description"])+"\n")
        for x in dupmat:
            tfil.write("\t".join([str(y) for y in x])+"\n")

#reformat bgc table and write
def writebgctable(temp,outfile):
    data=[]
    for key,x in temp.items():
        if key.startswith("cluster-"):
            corehits = len([hit for hit in x["hits"] if hit[4] == "Core"])
            data.append(x["row"]+[str(corehits),str(len(x["hits"])-corehits),x["hits"]])
    temp["data"] = data
    with open(outfile+".json","w") as jfil, open(outfile+".tsv","w") as tfil:
        json.dump(temp,jfil,indent=2)
        tfil.write("\t".join(["#Cluster","Type","Source","Location","Core hits","Other hits","Genelist"])+"\n")
        for x in temp["data"]:
            tfil.write("\t".join([str(y) for y in x])+"\n")

#get model metadata from json, tsv or hmm file
def getmodeldata(fname="",hmmname=""):
    modeldata = {}
    log.debug("Getting model metadata...")
    if ".json" in fname.lower() and os.path.exists(fname):
        with open(fname,"r") as fil:
            modeldata = json.load(fil)
    elif ".tsv" in fname.lower() and os.path.exists(fname):
        with open(fname,"r") as fil:
            for line in fil:
                if not line.startswith("#"):
                    x = line.strip().split("\t")
                    if x[0] not in modeldata:
                        modeldata[x[0]] = x[1:]
    elif ".hmm" in hmmname.lower() and os.path.exists(hmmname):
        with open(hmmname,"r") as fil:
            x = ["N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A"]
            for line in fil:
                if line.startswith("ACC"):
                    x[0] = " ".join(line.split()[1:])
                if line.startswith("NAME"):
                    x[1] = " ".join(line.split()[1:])
                if line.startswith("DESC"):
                    x[2] = " ".join(line.split()[1:])
                if line.startswith("TC"):
                    x[4] = line.split()[1]
                if line.startswith("HMM "):
                    if x[0] not in modeldata:
                        modeldata[x[0]] = x[1:]
                    x = ["N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A"]
    else:
        log.error("Could not get model metadata, ensure file is .json, .tsv, or .hmm")
        return False
    return modeldata

#simple sub strings from location text
def parsesourcelocation(x):
    x=x.strip()
    source = ""
    loc = ["","",""]

    if "|source|" in x:
        idx1 = x.index("|source|")
        if "|" in x[idx1+8:]:
            idx2 = x.index("|",idx1+8)
            source = x[idx1+8:idx2]
        else:
            source = x[idx1+8:]
    if "|loc|" in x:
        idx1 = x.index("|loc|")
        if "|" in x[idx1+5:]:
            idx2 = x.index("|",idx1+5)
            loc = x[idx1:idx2].split("_")
        else:
            loc = x[idx1+5:].split("_")
    return source,loc

def startquery(infile=None, refdir=None, td=None, rd=None, hmmdbs=None, rnahmm=None, cut=None,
                toconsole=False, mcpu=1, asrun=False, knownhmms=False, dufhmms=False,
                custcorehmms=False, custhmms=False, aspath=None, options="kres,duf", custorgname=None, bcp=None, run_bsc=False):
    try:
        #Set Working directory
        if type(rd) is str and not rd.endswith("/"):
            rd+="/"
        if type(rd) is str and os.path.exists(rd):
            tdir = os.path.realpath(rd)
        else:
            tdir=tempfile.mkdtemp(prefix="funarts-query-",dir=td)
        if not tdir.endswith("/"):
            tdir+="/"
        if not refdir.endswith("/"):
            refdir+="/"
        try:
            mcpu = int(mcpu)
        except ValueError:
            mcpu = mp.cpu_count()
        oldfname = os.path.splitext(os.path.split(infile)[-1])[0]
        #Start log
        global log
        log = setlog.init(tdir+"funarts-query.log",toconsole=toconsole)

        #Set folders if not created already
        if not os.path.exists(tdir+"coregenes"):
            os.mkdir(tdir+"coregenes")
            os.mkdir(tdir+"tables")

        log.info("Init folders complete")
        log.debug("Starting query: dir=%s; pid=%s; args=%s"%(tdir,os.getpid(),locals()))
    except Exception as e:
        log.error("Problem with folder creation / disk error (%s)"%e)
        log.exception("exception")
        raise e

    #Start Antismash run if needed
    if asrun:
        try:
            os.mkdir(tdir+"antismash")
            log.info("Starting Antismash job...")
            infile = runantismash(aspath,os.path.realpath(infile),tdir,mcpu)
            log.info("Finished Antismash")
        except Exception as e:
            log.error("Antismash failed to execute (%s)"%e)
            log.exception("exception")
    else:
        # if infile.endswith(".gbk") or infile.endswith(".gbff"):
        #     makeantismashresults(aspath, os.path.realpath(infile), os.path.realpath(tdir + "antismash"))
        # elif infile.endswith(".json"):
        #     infile = makeantismashresults_json(aspath, os.path.realpath(infile), os.path.realpath(tdir + "antismash"), mcpu)
        # added on 27.06.2023 - start : The situation that local antismash is not run and the solution offers
        if infile.endswith(".json"):
            infile = makeantismashresults_json(aspath, os.path.realpath(infile), os.path.realpath(tdir + "antismash"),mcpu)
        else:
            log.warning("Local antismash is not run, check antismash configuration or check your arguments and be sure '-ras/--runantismash' is exist")
            log.warning("Local antismash is not run, check antismash_result.json file. It is recommended to use Antismash v6 or v7")
        # added on 27.06.2023 - end : The situation that local antismash is not run and the solution offers

    try:
        ### PREP SQL
        if not os.path.isfile(infile):
            log.error("Error: No such input path %s"%infile)
            raise IOError
        try:
            queryorg, querygenus, clusters, locus_names = parsegbk.convertgenes(infile,tdir,plasmid=True,userecnum=True,clust=True,rename=oldfname)
            with open(os.path.join(tdir,"locus_to_region.pickle"), "wb") as locus_file:
                pickle.dump(locus_names, locus_file)
            if custorgname:
                #remove any bad chars
                #if the orgname has characters like :/# then the astral-rangerdtl part of phylogeny also fails so maybe we can change the line below?
                # queryorg = ''.join(ch for ch in custorgname.replace(".","").replace("#","").replace(":","").replace("/","").replace(",","").replace("(","_").replace(")","_").replace("-","_").replace(" ","_") if ch.isalnum() or "_")
                queryorg = ''.join(ch for ch in custorgname.replace(".","").replace(",","").replace("(","_").replace(")","_").replace("-","_").replace(" ","_") if ch.isalnum() or "_")
                querygenus = queryorg.replace("_"," ").split()[0]

        except ValueError:
            log.error("Biopython could not parse genbank. please ensure no external sequence references are present")
            log.exception("exception")

        infasta = oldfname+".fna"
        makequerydb(infasta,tdir,"queryseqs",orgname=queryorg)
        log.info("query: org=%s"%queryorg)
        log.info("query: genus=%s"%querygenus)

        #use list of hmmmodels if directory is specified
        if type(hmmdbs) is list and all([validatehmms(x) for x in hmmdbs]):
            hmmdb = concathmms(hmmdbs,tdir)
            modeldata = getmodeldata(hmmname=hmmdb)
        else:
            hmmdb = os.path.join(refdir,"coremodels.hmm")
            if "expert" in options and os.path.exists(os.path.join(refdir,"coremodels_exp.hmm")):
                hmmdb = os.path.join(refdir,"coremodels_exp.hmm")
            modeldata = getmodeldata(fname=os.path.join(refdir,"model_metadata.json"),hmmname=hmmdb)
        if custcorehmms and validatehmms(custcorehmms):
            hmmdb = concathmms([converthmm(hmmdb,os.path.join(tdir,"core.hmm")),converthmm(custcorehmms,os.path.join(tdir,"custcore.hmm"))],tdir)
            modeldata.update(getmodeldata(hmmname=custcorehmms))
            # with open(tdir+"model_metadata.json","w") as expfil:
            #     json.dump(modeldata,expfil,indent=2)
    except Exception as e:
        log.error("Problem init database")
        log.exception("exception")
        raise e

    #defaultcut="tc" #changed on 08.11.2023
    defaultcut = None #added on 08.11.2023

    try:
        #Set cutoff thresholds for core search
        cut_vals = getcutvals(hmmdb)
        log.debug("Setting thresholds...")
        if cut and cut.upper() in ("GA","NC"):
            log.debug("Setting %s cutoff..."%cut)
            cut = {k:x.get(cut.upper(),[0])[0] for k,x in cut_vals.items()}
            defaultcut = None
        elif cut and cut.upper()[0] in "E" and cut.upper()[1] in "123":
            try:
                factors = {"1":0.90,"2":0.75,"3":0.50,"4":0.30}
                f = factors[cut[1]]
                cut = {k:float(x.get("NC",[0])[0])*f for k,x in cut_vals.items()}
                defaultcut = None
                log.debug("Setting %s NC cutoff..."%f)
            except ValueError:
                cut = None
                log.warning("Warning: invalid value for threshold, using defaults")
        else:
            cut = None
    except Exception as e:
        log.warning("Could not set thresholds using default trusted cutoffs")
        cut = None

    try:
        ### RUN RNA Search and add to sql
        if rnahmm and os.path.exists(rnahmm):
            log.info("Starting RNA hmmsearch...")
            hmmdomrslt=runhmmer(tdir+"queryseqs.fna",rnahmm,tdir,cut="tc",mcpu=mcpu)
            log.info("Adding RNA hmmresults...")
            makehmmsql.run(hmmdomrslt,tdir+"queryseqs.db",ev=0.1,rna=True)

            ### Extract RNA hits
            getrnagenes.runfile(tdir+"queryseqs.db",outdir=tdir+"coregenes")
    except Exception as e:
        log.warning("Could not get RNA hits")
        log.exception("exception")

    try:
        #Run Known resistance models and custom supplied models
        knownhits={}
        if knownhmms and os.path.exists(knownhmms) and "kres" in options.lower():
            #Confirm and add custom hmm models
            log.info("Checking customhmms %s..."%custhmms)
            if validatehmms(custhmms):
                knownhmms = concathmms([converthmm(knownhmms,os.path.join(tdir,"kres.hmm")),converthmm(custhmms,os.path.join(tdir,"custres.hmm"))],tdir)
                log.info("Combined knownhmms = %s"%knownhmms)
            modeldata.update(getmodeldata(hmmname=knownhmms))
            log.info("Start known resistance search...")
            hmmdomrslt=runhmmer(tdir+"queryseqs.faa",knownhmms,tdir,cut="tc",mcpu=mcpu)
            os.rename(hmmdomrslt,tdir+"knownhits.domhr")
            #write results to json file
            knownhits = {"data":[],"seqs":{}}
            with open(os.path.join(tdir,"knownhits.domhr"),"r") as fil, open(os.path.join(tdir,"tables","knownhits.json"),"w") as ofil, open(os.path.join(tdir,"tables","knownhits.tsv"),"w") as tfil:
                tfil.write("#Model\tDescription\tSequence id\tevalue\tbitscore\tSequence description\n")
                for line in fil:
                    if not line.startswith("#"):
                        x = line.strip().split()
                        x[-1].replace("'","")
                        source,loc = parsesourcelocation(line)
                        seqid = x[0].split("|")[-1]
                        if seqid not in knownhits["seqs"]:
                            knownhits["data"].append([x[4],x[3],seqid,x[6],x[7],x[-1]])
                            tfil.write("\t".join([x[4],x[3],seqid,x[6],x[7],x[-1]])+"\n")
                            # [Seqid, modelname, scaffold_source, loc_start, loc_end, strand, evalue, bitscore]
                            knownhits["seqs"][seqid] = [seqid,x[4],source]+loc+[x[6],x[7],x[-1]]
                        elif float(knownhits["seqs"][seqid][5]) < float(x[7]):
                            knownhits["seqs"][seqid] = [seqid,x[4],source]+loc+[x[6],x[7],x[-1]]
                        if x[4] not in knownhits:
                            knownhits[x[4]] = {}
                        knownhits[x[4]][seqid] = [x[4],x[3],seqid,x[6],x[7],x[-1]]
                json.dump(knownhits,ofil,indent=2)
                log.info("Known Resistance Hits: %s"%len(knownhits["seqs"].keys()))
    except Exception as e:
        log.warning("Could not get known resistance hits")
        log.exception("exception")

    try:
        #Run DUF models
        dufhits={}
        if dufhmms and os.path.exists(dufhmms) and "duf" in options.lower():
            log.info("Start DUF search...")
            hmmdomrslt=runhmmer(tdir+"queryseqs.faa",dufhmms,tdir,cut="tc",mcpu=mcpu)
            os.rename(hmmdomrslt,tdir+"dufhits.domhr")
            #write results to json file
            dufhits = {"data":[],"seqs":{}}
            # with open(os.path.join(tdir,"dufhits.domhr"),"r") as fil, open(os.path.join(tdir,"tables","dufhits.json"),"w") as ofil:
            with open(os.path.join(tdir, "dufhits.domhr"), "r") as fil, open(os.path.join(tdir, "tables", "dufhits.json"), "w") as ofil, open(os.path.join(tdir,"tables","dufhits.tsv"),"w") as tfil:
                tfil.write("#Model\tDescription\tSequence id\tevalue\tbitscore\tSequence description\n")
                for line in fil:
                    if not line.startswith("#"):
                        # x = line.split()[:8]
                        x = line.strip().split()
                        x[-1].replace("'", "")
                        source,loc = parsesourcelocation(line)
                        seqid = x[0].split("|")[-1]
                        if seqid not in dufhits["seqs"]:
                            dufhits["data"].append([x[4],x[3],seqid,x[6],x[7]])
                            tfil.write("\t".join([x[4], x[3], seqid, x[6], x[7], x[-1]]) + "\n")
                            # [Seqid, modelname, scaffold_source, loc_start, loc_end, strand, evalue, bitscore]
                            dufhits["seqs"][seqid] = [seqid,x[4],source]+loc+[x[6],x[7]]
                        elif float(dufhits["seqs"][seqid][5]) < float(x[7]):
                            dufhits["seqs"][seqid] = [seqid,x[4],source]+loc+[x[6],x[7]]
                        if x[4] not in dufhits:
                            dufhits[x[4]] = {}
                        dufhits[x[4]][seqid] = [x[4],x[3],seqid,x[6],x[7]]
                json.dump(dufhits,ofil,indent=2)
                log.info("DUF Hits: %s"%len(dufhits["seqs"].keys()))
    except Exception as e:
        log.warning("Could not get DUF hits")
        log.exception("exception")

    log.info("Milestone_1_complete")

    ### RUN CORE Search and add to sql
    if hmmdb:
        log.info("Starting Core Gene hmmsearch... (%s)"%(hmmdb))
        hmmdomrslt=runhmmer(tdir+"queryseqs.faa",hmmdb,tdir,cut=defaultcut,mcpu=mcpu)
        log.info("Adding hmmresults...")
        # makehmmsql.run(hmmdomrslt, tdir + "queryseqs.db", ev=0.1) ### orginale = 29.08.2022
        makehmmsql.run(hmmdomrslt, tdir + "queryseqs.db", ev=0.1, filt2=os.path.join(refdir,'scores_cutoff')) ##v2
    else:
        log.error("No hmm models found")
        raise IOError

    ### EXTRACT Core and Count duplicates and check proximity
    log.info("Extracting core genes...")
    # corelist = getcoregenes.writeall(tdir+"coregenes/",tdir+"queryseqs.db")
    rslt = {}
    corelist = {}
    try:
        # corelist = getcoregenes.writeall(tdir+"coregenes/",tdir+"queryseqs.db",filt2=cut)
        corelist = extractdbgenes.writeall(tdir+"coregenes/",tdir+"queryseqs.db",filt2=cut)
        refgm = parsegmatrix(refdir + "genematrix.txt")
        # log.debug("refgm::::%s"%refgm)
        genematrix = parsegmatrix(tdir+"coregenes/"+"genematrix.txt",avgs=False)
        # log.debug("refgm::::%s"%genematrix)
        rslt["singles"] = list(refgm.get("_singles",set())&genematrix.get("_singles",set()))
        log.info("Found %s single copy genes for MLST: %s"%(len(rslt["singles"]),rslt["singles"]))
        rslt["duplicates"] = list(getdupgenes(refgm,genematrix))
        rslt["orgs"] = list(genematrix["_orgs"])+list(refgm["_orgs"])
        log.info("Found %s duplicate genes"%len(rslt["duplicates"]))
        if "_singles" in genematrix:
            del genematrix["_singles"]
        if "_orgs" in genematrix:
            del genematrix["_orgs"]
        rslt["matrix"] = genematrix
        #check core gene and known resistance factor proximity to BGC clusters
        log.info("Writing duplicates table...")

        #write out display tables
        writeduptable(rslt["matrix"],rslt["duplicates"],corelist,modeldata,os.path.join(tdir,"tables","duptable"))

    except Exception as e:
        log.error("Could not extract coregenes")
        log.exception("exception")

    try:
        log.info("Writing bgc and core tables...")
        bgchits,rslt["proximity"],numhits = checkbgcprox(clusters,gtypes=["Core","ResModel","DUF"],glist=[corelist,knownhits,dufhits],modeldata=modeldata)
        log.info("Proximity hits found: %s"%numhits)
        writebgctable(bgchits,os.path.join(tdir,"tables","bgctable"))
        writecoretable(corelist,rslt,knownhits,modeldata,os.path.join(tdir,"tables","coretable"))
        log.info("Milestone_2_complete")
    except Exception as e:
        log.error("Could not write output tables")
        log.exception("exception")

    log.info("Exporting results...")
    try:
        for x in [f for f in os.listdir(tdir) if f.startswith("queryseqs") or f.startswith("domrslt") or f.endswith(".hmm")]:
            os.remove(os.path.join(tdir,x)) #Save disk space
    except Exception as e:
        log.error("Problem with export")
        log.exception("exception")

    #count multi criteria hits
    twoplus = set()
    if "proximity" in rslt and "duplicates" in rslt:
        twoplus |= set(rslt["proximity"])&set(rslt["duplicates"])

    # log.info("Hits with two or more criteria: %s : %s"%(len(twoplus),twoplus))
    log.info("Hits with Duplication and Proximity criteria: %s : %s" % (len(twoplus), twoplus))

    #update table
    writecoretable(corelist,rslt,knownhits,modeldata,os.path.join(tdir,"tables","coretable"))

    log.info("SUCCESS! job finished")
    return True

def parse_input_orgs(input):
    c = 0
    final_input_list = []
    if "," in input:
        input_list = input.split(",")
        for i in input_list:
            if i in final_input_list:
                c += 1
                i_final = i + "_" + str(c)
                final_input_list.append(i_final)
            else:
                final_input_list.append(i)
        return final_input_list
    else:
        final_input_list = [input]
    return final_input_list

def run_bigscape(antismash_result_directories, bspath, maindir):
    combined_log.info("bigscape run start")
    all_antismash_results = os.path.join(maindir , "all_antismash")
    result_directories = antismash_result_directories
    regions_to_filename_dict = {}
    if os.path.exists(all_antismash_results) != True:
        os.mkdir(all_antismash_results)
    for directory in result_directories:
        antismash_results_path = os.path.join(directory, "antismash")
        copy_tree(antismash_results_path, os.path.join(all_antismash_results, os.path.basename(directory)))
        for file in os.listdir(antismash_results_path):
            if (".region" in file or ".cluster" in file) and ".gbk" in file:
                full_file_path = os.path.join(os.path.join(directory, "antismash"), file)
                regions_to_filename_dict[full_file_path] = result_directories[directory]
    bigscape_result_folder = os.path.join(maindir , "antismash_bigscape_result")
    # cmd = "python %s -i %s -o %s --mix" % (bspath, all_antismash_results, bigscape_result_folder)
    bslog = os.path.join(maindir , "bslog.txt")
    bserror = os.path.join(maindir, "bserror.txt")
    cmd = "python3 %s -i %s -o %s --mix > %s 2> %s" % (bspath, all_antismash_results, bigscape_result_folder, bslog, bserror)
    try:
        subprocess.call([cmd], shell=True)
    except subprocess.CalledProcessError as e:
        for k, v in e._dict_.items():
            print(k)
            combined_log.error(k)
            if isinstance(v, bytes):
                err = v.decode("utf-8")
            else:
                err = v
            if k == "stderr":
                combined_log.error(err)
    combined_log.info("bigscape run saved to: %s"%bigscape_result_folder)
    return bigscape_result_folder, regions_to_filename_dict


# def create_symlinks(result_dict, maindir):
#     for result_path in result_dict:
#         os.symlink(result_path, os.path.join(os.path.abspath(os.path.join(maindir, os.pardir)), os.path.basename(result_path)))

def write_paths_to_file(result_dict,maindir):
    path_file = open(os.path.join(maindir,"all_paths"), "w")
    path_file.write("#Result Directory\tInput File\n")
    for result_path in result_dict:
        input = result_dict[result_path]
        path_file.write(result_path + "\t" + input + "\n")
    path_file.close()

def call_startquery(args):
    #parse input function expects a string of comma seperated organisms as in "org1,org2,org3" or if its just one organism "org1"
    input_list = parse_input_orgs(args.input)
    tdir_increment = 0

    #initialize folders
    if args.resultdir != None:
        main_dir = args.resultdir
    else:
        main_dir = tempfile.mkdtemp(prefix="funarts-query-", dir=args.tempdir)

    if os.path.exists(main_dir) != True:
        os.mkdir(main_dir)

    start_query_count = 0
    result_directories = []
    result_dict = {}

    global combined_log
    combined_log = setlog.init(os.path.join(main_dir, "combined.log",), toconsole=True, logname="combinedlog")

    if len(input_list) == 1:
        if os.path.exists(os.path.join(main_dir, "combined.log")) ==True:
            os.remove(os.path.join(main_dir, "combined.log"))
        #if only one input call artspipeline once
        #run_bsc and bigscape path must be added to front end part since it will give argument errors otherwise
        startquery(infile=args.input, refdir=args.refdir, td=args.tempdir, rd=args.resultdir, hmmdbs=args.hmmdblist,
                   rnahmm=args.rnahmmdb, cut=args.thresh,
                   toconsole=True, mcpu=args.multicpu, asrun=args.runantismash,
                   knownhmms=args.knownhmms, dufhmms=args.dufhmms,
                   custcorehmms=args.custcorehmms, custhmms=args.customhmms, aspath=args.antismashpath,
                   options=args.options, custorgname=args.orgname,
                   bcp = args.bigscapepath, run_bsc = args.runbigscape)
    else:
        for input in input_list:
            start_query_count += 1
            resultdir = os.path.join(main_dir, os.path.basename(args.resultdir) + "_" + str(tdir_increment))
            if os.path.exists(resultdir) != True:
                os.mkdir(resultdir)
                #make symlink immediately so progress can be viewed during processing
                os.symlink(resultdir,os.path.join(os.path.dirname(args.resultdir),os.path.basename(resultdir)))
            #keep result directory paths and corresponding input organisms for further usage from frontend
            result_directories.append(resultdir)
            result_dict[resultdir] = input
            tdir_increment += 1
            combined_log.info("funartspipeline number " + str(start_query_count) + " start")
            combined_log.info("file: " + str(list(os.path.split(input))[-1]))
            startquery(infile=input, refdir=args.refdir, td=args.tempdir, rd=resultdir, hmmdbs=args.hmmdblist,
                       rnahmm=args.rnahmmdb, cut=args.thresh,
                       toconsole=True, mcpu=args.multicpu, asrun=args.runantismash,
                       knownhmms=args.knownhmms, dufhmms=args.dufhmms,
                       custcorehmms=args.custcorehmms, custhmms=args.customhmms, aspath=args.antismashpath,
                       options=args.options, custorgname=args.orgname,
                       bcp=args.bigscapepath, run_bsc = args.runbigscape)
            combined_log.info("funartspipeline number " + str(start_query_count) + " end")
        run_bsc = args.runbigscape
        combined_log.info("Run bigscape: %s %s"%(run_bsc,args.bigscapepath))
        if run_bsc:
            #bspath must be added to front end like aspath
            bspath = args.bigscapepath
            bigscape_result_folder, regions_to_filename_dict = run_bigscape(result_dict, bspath, main_dir)
        #File containing all the paths for individual ARTS results and inputs
        write_paths_to_file(result_dict,main_dir)
        generate_summary(result_dict, main_dir)
        #combine all the tables from arts results
        combine_core_results(result_dict, main_dir)
        combine_known_results(result_dict, main_dir)
        combine_dup_results(result_dict, main_dir)
        if run_bsc:
            combine_bigscape_results(bigscape_result_folder,main_dir, regions_to_filename_dict)
        try:
            parse_json(result_dict, main_dir, regions_to_filename_dict)
        except:
            combined_log.error("Problem in parse_json")
        #Create the symlinks on results folder, should only be available on server
        #create_symlinks(result_dict, main_dir)
        generate_plots(main_dir)



# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Start from genbank file and compare with pre-computed reference for Duplication and Transfers""")
    parser.add_argument("input", help="gbk file to start query")
    parser.add_argument("refdir", help="Directory of precomputed reference files")
    parser.add_argument("-hmms","--hmmdblist", help="hmm file, directory, or list of hmm models for core gene id",default=None)
    parser.add_argument("-khmms","--knownhmms", help="Resistance models hmm file",default=False)
    parser.add_argument("-duf","--dufhmms", help="Domains of unknown function hmm file",default=False)
    parser.add_argument("-cchmms","--custcorehmms", help="User supplied core models. hmm file",default=False)
    parser.add_argument("-chmms","--customhmms", help="User supplied resistance models. hmm file",default=False)
    parser.add_argument("-rhmm","--rnahmmdb", help="RNA hmm models to run (default: None)",default=None)
    parser.add_argument("-t","--thresh", help="Hmm reporting threshold. Use global bitscore value or Model specific options: gathering= GA, trusted= TC, noise= NC(default: None)",default=None)
    parser.add_argument("-td", "--tempdir", help="Directory to create unique results folder", default=None)
    parser.add_argument("-rd", "--resultdir", help="Directory to store results", default=None)
    parser.add_argument("-cpu", "--multicpu", help="Turn on Multi processing set # Cpus (default: Off, 1)", type=int, default=1)
    parser.add_argument("-opt", "--options",help="Analysis to run. expert=Exploration mode, kres=Known resistance, duf=Domain of unknown function (default: None)",default=None)
    parser.add_argument("-org", "--orgname", help="Explicitly specify organism name", default=None)
    parser.add_argument("-ras", "--runantismash", help="Run input file through antismash first", action='store_true', default=False)
    parser.add_argument("-asp", "--antismashpath", help="Location of the executable file of antismash or location of antismash 'run_antismash.py' script", default=False)
    parser.add_argument("-bcp", "--bigscapepath", help="Location of bigscape 'bigscape.py' script", default=False)
    parser.add_argument("-rbsc", "--runbigscape",help="Run antismash results through bigscape", action='store_true', default=False )
    args = parser.parse_args()
    call_startquery(args)
    # startquery(infile=args.input,refdir=args.refdir,td=args.tempdir,rd=args.resultdir,hmmdbs=args.hmmdblist,rnahmm=args.rnahmmdb,cut=args.thresh,
    #            astjar=args.astral,toconsole=True,mcpu=args.multicpu,asrun=args.runantismash,knownhmms=args.knownhmms,dufhmms=args.dufhmms,
    #            custcorehmms=args.custcorehmms,custhmms=args.customhmms,aspath=args.antismashpath,options=args.options,custorgname=args.orgname,prebuilttrees=args.prebuilttrees)

