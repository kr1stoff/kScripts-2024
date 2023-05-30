#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/23 15:00
import os
import yaml
import sys
import click

PATH = os.path.dirname(os.path.abspath(__file__))  
DIR=os.path.dirname(PATH)

sys.path.append(os.path.dirname(PATH))
import lib.common as common
import lib.pipe as pipe

# read YAML
p_cfg = os.path.join(DIR, "modules/cfg.yaml")
f_config = yaml.safe_load(open(p_cfg , mode='r',encoding='utf-8_sig').read())
b_base = f"{f_config['b_base']}"
th = f"{f_config['threads']}"
p_bin = os.path.join(DIR, "bin") 

class SeqQc():
    def __init__(self,**kwargs):
        self.id = kwargs['id']
        self.out = os.path.abspath(kwargs['out'])
        self.p_sh = os.path.abspath(kwargs['sh'])
        self.fq1 = os.path.abspath(kwargs['fq1'])
        self.fq2 = '' if kwargs['fq2'] == 'None' else os.path.abspath(kwargs['fq2'])


    def make_dir(self):
        self.p_fastqc = os.path.join(self.out,'fastqc')
        self.p_filter_fastqc = os.path.join(self.out,'filter_fastqc')
        self.p_fastp = os.path.join(self.out,'fastp')
        os.makedirs(self.p_fastqc ,exist_ok=True)
        os.makedirs(self.p_fastp ,exist_ok=True)
        os.makedirs(self.p_sh ,exist_ok=True)
        os.makedirs(self.p_filter_fastqc ,exist_ok=True)


    def fastqc(self):
        """ run fastp """
        sh = os.path.join(self.p_sh,"1.1fastqc.sh")
        cmd = []
        
        if self.fq2: #PE
            l_fq1 ,l_fq2 = pipe.pe_fq_link(self.id,self.fq1,self.fq2,self.p_fastqc) 
            cmd.append(f"{b_base}/fastqc -q -t {th} {l_fq1} {l_fq2} -o {self.p_fastqc}")
            cmd.append(f"unzip -qn {self.p_fastqc}/*1_fastqc.zip -d {self.p_fastqc}")
            cmd.append(f"unzip -qn {self.p_fastqc}/*2_fastqc.zip -d {self.p_fastqc}")
        else: #SE
            l_fq = pipe.se_fq_link(self.id,self.fq1,self.p_fastqc)
            cmd.append(f"{b_base}/fastqc -q -t {th} {l_fq}  -o {self.p_fastqc}")
            cmd.append(f"unzip -qn {self.p_fastqc}/*_fastqc.zip -d {self.p_fastqc}")

        common.cmd2shell(cmd,sh)
        return sh


    def fastp(self):
        """ run fastp """
        sh = os.path.join(self.p_sh,"1.2fastp.sh")
        cmd=[]
        if self.fq2: #PE
            cmd.append(f"{b_base}/fastp -n 1 -l 50 -y -w 16 -i {self.fq1} -I {self.fq2} -j {self.p_fastp}/{self.id}.clean.fq.stat.json -o {self.p_fastp}/{self.id}_1.clean.fq.gz -O {self.p_fastp}/{self.id}_2.clean.fq.gz -h {self.p_fastp}/{self.id}.clean.html ")
            cmd.append(f"{b_base}/python3 {p_bin}/1.qc/parse_fastp.json.py -id {self.id} -i {self.p_fastp}/{self.id}.clean.fq.stat.json -m PE ")

        else: #SE
            cmd.append(f"{b_base}/fastp -n 1 -l 50 -y -w 16 -i {self.fq1} -j {self.p_fastp}/{self.id}.clean.fq.stat.json -h {self.p_fastp}/{self.id}.clean.html -o {self.p_fastp}/{self.id}_1.clean.fq.gz ")
            cmd.append(f"{b_base}/python3 {p_bin}/1.qc/parse_fastp.json.py -id {self.id} -i {self.p_fastp}/{self.id}.clean.fq.stat.json -m SE ")

        common.cmd2shell(cmd,sh)
        return sh

    def filter_fastqc(self):
        """ run fastqc """
        sh = os.path.join(self.p_sh,"1.3filter_fastqc.sh")
        cmd = []
        
        if self.fq2: #PE
            cmd.append(f"{b_base}/fastqc -q -t {th} {self.p_fastp}/{self.id}_1.clean.fq.gz {self.p_fastp}/{self.id}_2.clean.fq.gz -o {self.p_filter_fastqc}")
            cmd.append(f"unzip -qn {self.p_filter_fastqc}/*_1.clean_fastqc.zip -d {self.p_filter_fastqc}")
            cmd.append(f"unzip -qn {self.p_filter_fastqc}/*_2.clean_fastqc.zip -d {self.p_filter_fastqc}")
        else: #SE
            cmd.append(f"{b_base}/fastqc -q -t {th}{self.p_fastp}/{self.id}_1.clean.fq.gz -o {self.p_filter_fastqc}")
            cmd.append(f"unzip -qn {self.p_filter_fastqc}/*_fastqc.zip -d {self.p_filter_fastqc}")

        common.cmd2shell(cmd,sh)
        return sh




# Option
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-id',required=True,type=click.STRING,help="Sample id")
@click.option('-o', '--out',required=True,type=click.Path(),default='.',show_default=True,help="Output dir")
@click.option('-sh',required=True,type=click.Path(),help="Write shell dir")
@click.option('-fq1',required=True,type=click.Path(),help="Input SE_fq / fq1")
@click.option('-fq2',required=False,type=click.Path(),default='None',help="Input fq2 if seq_mode=PE")


def main(id,out,fq1,fq2,sh):
    #log
    logfile=f"{out}/log"
    os.makedirs(logfile,exist_ok=True)

    project = SeqQc(
        id = id,
        out =out,
        fq1 = fq1,
        fq2 = fq2,
        sh = sh
    )

    project.make_dir()
    sh1 = project.fastqc()
    sh2 = project.fastp()
    shlist=[sh1,sh2]
    common.mul_pool(shlist, logfile)
    sh3 = project.filter_fastqc()
    common.prun((sh3, logfile))

if __name__ == "__main__":
    main()
