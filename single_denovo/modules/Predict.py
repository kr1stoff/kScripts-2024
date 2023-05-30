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

# read YAML
p_cfg = os.path.join(DIR, "modules/cfg.yaml")
p_bin = os.path.join(DIR, "bin")
f_config = yaml.safe_load(open(p_cfg , mode='r',encoding='utf-8_sig').read())
b_base = f"{f_config['b_base']}"
b_denovo = f"{f_config['b_denovo']}"
b_act = f"{f_config['ACTIVATE']}"
cpu ,threads = common.get_cfg()

class Predict():
    def __init__(self,**kwargs):

        self.id = kwargs['id']
        self.fa = kwargs['fa']
        self.out = os.path.abspath(kwargs['out'])
        self.p_sh = os.path.abspath(kwargs['sh'])
        self.kind = os.path.abspath(kwargs['kind'])
        self.bam = os.path.abspath(kwargs['bam']) if kwargs['bam'] else ''


    def make_dir(self):
        os.makedirs(self.p_sh ,exist_ok=True)
        self.p_predict = os.path.join(self.out,'predict')
        os.makedirs(self.p_predict,exist_ok=True)
        self.p_island = os.path.join(self.out,'island')
        os.makedirs(self.p_island,exist_ok=True)


    
    def bacteria_predict(self):
        cmd = []
        sh = os.path.join(self.p_sh,"4.precict.sh")
        cmd.append(f"""#prokka
source {b_act} denovo
cd {self.p_predict}
time prokka {self.fa} --prefix predict --cpus {int(cpu*3/4)} --outdir {self.p_predict} --kingdom Bacteria --force --addgenes --quiet

if [ -f "{self.p_predict}/predict.gff" ]; then
    {b_base}/python3 {p_bin}/3.pre/chinese_prokka.py -i {self.p_predict}/predict.txt -o {self.p_predict}/predict_kind.txt

    # bed
    grep -v '#' {self.p_predict}/predict.gff|grep 'CDS' |grep 'ID='|cut -f 1,4,5,9|awk -F';' '{{print $1}}'|sed 's/ID=//g' >{self.p_predict}/predict.bed
    bedtools coverage -mean -a {self.p_predict}/predict.bed -b {self.bam} >{self.p_predict}/bed_bam.merge.depth
    bedtools coverage -a {self.p_predict}/predict.bed -b {self.bam} >{self.p_predict}/bed_bam.merge.cov
    paste {self.p_predict}/bed_bam.merge.depth {self.p_predict}/bed_bam.merge.cov |cut -f 4,5,13 >{self.p_predict}/bed_bam.merge_depth_cov.txt

    # gene_length
    {b_denovo}/python3 {p_bin}/2.ass_qc/ass_fa_stat.py -fa {self.p_predict}/predict.ffn -p predict
fi
    # island
    time {b_base}/python3 {p_bin}/3.pre/multi_Island.py -i {self.p_predict}/predict.gbk -o {self.p_island} >{self.p_island}/island.log 2>&1
    #phispy
    time {b_base}/python3 {p_bin}/3.pre/run_phispy.py -i {self.p_predict}/predict.gbk -o {self.out}/phispy
""")

        common.cmd2shell(cmd,sh)
        return sh

    
    def fungi_predict(self):
        cmd = []
        sh = os.path.join(self.p_sh,"4.precict.sh")
        cmd.append(f"""#genmark
source {b_act} denovo
cd {self.p_predict}
time {b_denovo}/gmes_petap.pl --ES --fungus --cores {int(cpu*3/4)} --sequence {self.fa}  >{self.p_predict}/pre.log

if [ -f "{self.p_predict}/genemark.gtf" ]; then
    echo "gtf ==> gff3 ==> cds.fa ==> cds.faa"
    gffread {self.p_predict}/genemark.gtf -o ->{self.p_predict}/predict.gff3
    {b_base}/perl {p_bin}/tools/getGene.pl {self.p_predict}/predict.gff3 {self.fa} -type mrna >{self.p_predict}/predict.fa
    {b_base}/perl {p_bin}/tools/cds2aa.pl {self.p_predict}/predict.fa >{self.p_predict}/predict.faa

    grep 'CDS' {self.p_predict}/predict.gff3|grep 'Parent='|cut -f 1,4,5,9|sed -e 's/Parent=//g' >{self.p_predict}/predict.gff3.bed
    bedtools coverage -mean -a {self.p_predict}/predict.gff3.bed -b {self.bam} >{self.p_predict}/bed_bam.merge.depth
    bedtools coverage -a {self.p_predict}/predict.gff3.bed -b {self.bam} >{self.p_predict}/bed_bam.merge.cov
    paste {self.p_predict}/bed_bam.merge.depth {self.p_predict}/bed_bam.merge.cov |cut -f 4,5,13 >{self.p_predict}/bed_bam.merge_depth_cov.txt

    {b_base}/python3 {p_bin}/3.pre/gene_mark_count.py {self.p_predict}/predict.gtf {self.p_predict}/predict_kind.txt
    {b_denovo}/python3 {p_bin}/2.ass_qc/ass_fa_stat.py -fa {self.p_predict}/predict.fa -p predict
fi
""")

        common.cmd2shell(cmd,sh)
        return sh



# Option
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-id',required=True,type=click.STRING,help="样本编号")
@click.option('-fa',required=True,type=click.Path(),help="组装序列 contigs/scaffolds fasta")
@click.option('-o', '--out',required=True,type=click.Path(),default='.',show_default=True,help="输出目录")
@click.option('-sh',required=True,type=click.Path(),default='.',show_default=True,help="写shell脚本目录")
@click.option('-k','--kind',required=True,type=click.Choice(['bacteria','fungi']),help="选择病原体类型")
@click.option('-b', '--bam',required=True,type=click.Path(),help="输入将组装序列比对过滤序列，经过排序后的bam文件，用于合成bed文件")


def main(id,fa,out,sh,kind,bam):
    #log
    logfile=f"{out}/log"
    os.makedirs(logfile,exist_ok=True)

    project = Predict(
                id = id,
                fa = fa,
                out = out,
                sh = sh,
                kind = kind,
                bam = bam
            )
    project.make_dir()
    sh1 = project.bacteria_predict() if (kind == 'bacteria') else project.fungi_predict()
    common.prun((sh1, logfile))



if __name__ == "__main__":
    main()
