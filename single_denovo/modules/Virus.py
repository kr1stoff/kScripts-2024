#!/usr/bin/env python3
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

# read YAML
p_cfg = os.path.join(DIR, "modules/cfg.yaml")
p_bin = os.path.join(DIR, "bin")

f_config = yaml.safe_load(open(p_cfg , mode='r',encoding='utf-8_sig').read())
b_base = f"{f_config['b_base']}"
b_denovo = f"{f_config['b_denovo']}"
b_act = f"{f_config['ACTIVATE']}"
p_db = f"{f_config['p_db']}"
th = f"{f_config['threads']}" 


class VIRUS():

    def __init__(self,**kwargs):
        self.id = kwargs['id']
        self.out = os.path.abspath(kwargs['out'])
        self.ref = os.path.abspath(kwargs['ref'])
        self.p_sh = os.path.abspath(kwargs['sh'])
        self.fq1 = os.path.abspath(kwargs['fq1'])
        self.fq2 = '' if (kwargs['fq2'] == 'None') else os.path.abspath(kwargs['fq2'])


    def make_dir(self):
        os.makedirs(self.p_sh ,exist_ok=True)
        self.p_ass = os.path.join(self.out,'align')
        os.makedirs(self.p_ass ,exist_ok=True)


    def align(self):
        sh  = os.path.join(self.p_sh,"2.ass.sh")
        cmd=[]
        cmd.append(f"""#Virus consensus
cd {self.p_ass}
source {b_act} denovo
{b_base}/bwa index -a bwtsw {self.ref}
""")
        if self.fq2:
            cmd.append(f"{b_base}/bwa mem -t 80 -Y {self.ref} {self.fq1} {self.fq2} |{b_base}/samtools view -@ 24 -hS -bF 12 |{b_base}/samtools sort -@ 24 -o {self.p_ass}/{self.id}_bf12_sort.bam")
        else:
            cmd.append(f"{b_base}/bwa mem -t 80 -Y {self.ref} {self.fq1} |{b_base}/samtools view -@ 24 -hS -bF 12 |{b_base}/samtools sort -@ 24 -o {self.p_ass}/{self.id}_bf12_sort.bam")

        cmd.append(f"""
{b_base}/samtools depth -a {self.p_ass}/{self.id}_bf12_sort.bam >{self.p_ass}/bf12_sort.depth
{b_base}/samtools mpileup -aa -A -d 0 -Q 0 {self.p_ass}/{self.id}_bf12_sort.bam | ivar consensus -p {self.id}_consensus
{b_base}/samtools mpileup -aa -A -d 0 -B -Q 0 --reference {self.ref} {self.p_ass}/{self.id}_bf12_sort.bam |ivar variants -p variants -r {self.ref}

#Depth dic
python3 {p_bin}/2.ass_qc/depth_pic.py {self.p_ass}/bf12_sort.depth {self.p_ass} {self.ref}
python3 {p_bin}/2.ass_qc/ivar_var.py variants.tsv {self.p_ass}
""")
        common.cmd2shell(cmd,sh)
        return sh



# Option
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-id',required=True,type=click.STRING,help="Sample id")
@click.option('-o', '--out',required=True,type=click.Path(),default='.',show_default=True,help="Output dir")
@click.option('-sh',required=True,type=click.Path(),default='.',show_default=True,help="Write shell dir")
@click.option('-fq1',required=True,type=click.Path(),help="Input SE_fq / fq1")
@click.option('-fq2',required=False,type=click.Path(),default='None',help="Input fq2 if seq_mode=PE")
@click.option('-ref',required=True,type=click.Path(),default='None',help="Virus reference")

def main(id,out,fq1,fq2,sh,ref):
    #log
    logfile=f"{out}/log"
    os.makedirs(logfile,exist_ok=True)

    project = VIRUS(
        id = id,
        out =out,
        fq1 = fq1,
        fq2 = fq2,
        sh = sh,
        ref = ref
    )
    project.make_dir()

    sh1 = project.align()
    common.prun((sh1, logfile))


if __name__ == "__main__":
    main()