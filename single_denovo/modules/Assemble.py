#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/23 15:00
import os
import yaml
import sys
import click
import sys

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
cpu ,threads = common.get_cfg()

class Assemble():
    def __init__(self,**kwargs):
        self.id = kwargs['id']
        self.kind = kwargs['kind']
        self.ass_flag = kwargs['ass']
        self.if_pollute = kwargs['if_pollute']
        self.out = os.path.abspath(kwargs['out'])
        self.p_sh = os.path.abspath(kwargs['sh'])
        self.fq1 = os.path.abspath(kwargs['fq1'])
        self.fq2 = '' if (kwargs['fq2'] == 'None') else os.path.abspath(kwargs['fq2'])

        if (kwargs['cut_reads']) and (self.kind == 'bacteria'):
            self.num = int(kwargs['cut_reads'])*1000000  #M
        else: self.num = ''


    def make_dir(self):
        os.makedirs(self.p_sh ,exist_ok=True)
        if self.num:
            self.p_deal = os.path.join(self.out,'deal_reads')
            os.makedirs(self.p_deal ,exist_ok=True)
        if self.ass_flag:
            self.p_ass = os.path.join(self.out,'spades')
            os.makedirs(self.p_ass ,exist_ok=True)
    


    def Advanced_Option(self):
        """ run cut reads """
        sh = os.path.join(self.p_sh,"2.1Advanced_Option.sh")
        cmd = []

        #Depollute
        if self.if_pollute:
            if self.fq2:  #PE
                cmd.append(f"{b_base}/python3 {p_bin}/1.qc/depollute/depollute_by_kraken2.py -p {self.p_deal}/depollute --config {p_bin}/1.qc/depollute/depollute_by_kraken2.yaml {self.fq1} {self.fq2}")
                self.fq1 = os.path.join(self.p_deal,'depollute_1.fq')
                self.fq2 = os.path.join(self.p_deal,'depollute_2.fq')

            else: #SE
                cmd.append(f"{b_base}/python3 {p_bin}/1.qc/depollute/depollute_by_kraken2.py -p {self.p_deal}/depollute --config {p_bin}/1.qc/depollute/depollute_by_kraken2.yaml {self.fq1}")
                self.fq1 = os.path.join(self.p_deal,'depollute_1.fq')
            cmd.append(f"{b_base}/python3 {p_bin}/1.qc/parse_kraken.report.py -i {self.p_deal}/temp_depollute/depollute.report")


        #Cut reads
        if self.num:    
            if self.fq2:  #PE
                cmd.append(f"{b_denovo}/seqtk sample -s 11 {self.fq1} {self.num} >{self.p_deal}/cut_1.fq")
                cmd.append(f"{b_denovo}/seqtk sample -s 11 {self.fq2} {self.num} >{self.p_deal}/cut_2.fq")
                cmd.append(f"{b_denovo}/iTools Fqtools stat -InFq {self.p_deal}/cut_1.fq -InFq {self.p_deal}/cut_2.fq -OutStat {self.p_deal}/cutfq.stat")
                cmd.append(f"{b_denovo}/python3 {p_bin}/1.qc/parse.itools.stat.py -i {self.p_deal}/cutfq.stat -m PE")
                self.fq1 = os.path.join(self.p_deal,'cut_1.fq')
                self.fq2 = os.path.join(self.p_deal,'cut_2.fq')

            else:   #SE
                cmd.append(f"{b_denovo}/seqtk sample -s 11 {self.fq1} {self.num} >{self.p_deal}/cut_1.fq")
                cmd.append(f"{b_denovo}/iTools Fqtools stat -InFq {self.p_deal}/cut_1.fq -OutStat {self.p_deal}/cutfq.stat")
                cmd.append(f"{b_denovo}/python3 {p_bin}/1.qc/parse.itools.stat.py -i {self.p_deal}/cutfq.stat -m SE")
                self.fq1 = os.path.join(self.p_deal,'cut_1.fq')

        common.cmd2shell(cmd,sh)
        return sh



    def spades(self):
        sh = os.path.join(self.p_sh,"2.2spades.sh")
        cmd=[]
        if self.ass_flag:    #choose assemble
            cmd.append(f"source {b_act} denovo")
            if self.fq2:  #PE
                cmd.append(f"time spades.py -t {int(threads*3/4)} --careful -o {self.p_ass} -1 {self.fq1} -2 {self.fq2} >{self.p_ass}/spades.o ")
            else: #SE
                cmd.append(f"time spades.py -t {int(threads*3/4)} --careful -o {self.p_ass} -s {self.fq1} >{self.p_ass}/spades.o ")
            cmd.append(f"{b_denovo}/python3 {p_bin}/2.ass_qc/ass_fa.filter.py -fa {self.p_ass}/scaffolds.fasta -p {self.id}_scaffolds -l 500")
            cmd.append(f"{b_denovo}/python3 {p_bin}/2.ass_qc/ass_fa_stat.py -fa {self.p_ass}/{self.id}_scaffolds.fasta -p {self.id}_scaffolds")
        common.cmd2shell(cmd,sh)
        return sh




# Option
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-id',required=True,type=click.STRING,help="样本编号")
@click.option('-o', '--out',required=True,type=click.Path(),default='.',show_default=True,help="输出目录")
@click.option('-sh',required=True,type=click.Path(),help="写shell脚本目录")
@click.option('-k','--kind',required=True,type=click.Choice(['bacteria','fungi','virus']),help="选择病原体类型")
@click.option('-fq1',required=True,type=click.Path(),help="过滤后的fq1")
@click.option('-fq2',type=click.Path(),default='None',help="过滤后的fq2，单端只需要输出fq1")
@click.option('--if_pollute', is_flag=True,help="是否选择去污染")
@click.option('--cut_reads',type=click.INT,help="切reads(M),真菌不切")
@click.option('--ass',is_flag=True,help="运行组装流程")

def main(id,out,sh,fq1,fq2,cut_reads,ass,kind,if_pollute):

    #log
    logfile=f"{out}/log"
    os.makedirs(logfile,exist_ok=True)

    project = Assemble(
        id = id,
        out =out,
        fq1 = fq1,
        fq2 = fq2,
        sh = sh,
        cut_reads = cut_reads,
        ass = ass,
        kind = kind,
        if_pollute = if_pollute
    )
    project.make_dir()
    sh1 = project.Advanced_Option()
    common.prun((sh1, logfile))

    sh2 = project.spades()
    common.prun((sh2, logfile))


if __name__ == "__main__":
    main()
