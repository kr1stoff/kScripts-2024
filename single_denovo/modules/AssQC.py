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
picard = f"{f_config['picard']}" 

class AssembleQC():
    def __init__(self,**kwargs):

        self.id = kwargs['id']
        self.out = os.path.abspath(kwargs['out'])
        self.p_sh = os.path.abspath(kwargs['sh'])
        self.fa = os.path.abspath(kwargs['fa'])

        self.quast_flag = kwargs['quast']
        self.checkm_flag = kwargs['checkm']

        self.ref = '' if (kwargs['ref'] == 'None') else os.path.abspath(kwargs['ref'])

        if kwargs['depth_stat'] and kwargs['fq1'] :
            self.depth_flag = kwargs['depth_stat']
            self.fq1 = os.path.abspath(kwargs['fq1']) if (kwargs['fq1']) else ''
            self.fq2 = '' if (kwargs['fq2'] == 'None') else os.path.abspath(kwargs['fq2'])
        else: self.depth_flag = ''

        if kwargs['depth_stat'] and not (kwargs['fq1']):
            print(f"Input error: If choose 'depth_stat' ,must input `fq1`  ")
            exit(-1)



    def make_dir(self):
        os.makedirs(self.p_sh ,exist_ok=True)
        if self.quast_flag:
            self.p_quast = os.path.join(self.out,'quast')
            os.makedirs(self.p_quast ,exist_ok=True)
        if self.checkm_flag:
            self.p_checkm = os.path.join(self.out,'checkm')
            os.makedirs(self.p_checkm ,exist_ok=True)
        if self.depth_flag:
            self.p_depth = os.path.join(self.out,'depth')
            os.makedirs(self.p_depth ,exist_ok=True)
        if self.ref:
            self.p_ref = os.path.join(self.out,'ref')
            os.makedirs(self.p_ref ,exist_ok=True)


    def quast(self):
        sh = os.path.join(self.p_sh,"3.1quast.sh")
        cmd = []
        if self.quast_flag:
            cmd.append(f"source {b_act} denovo ")
            cmd.append(f"python3 {b_denovo}/quast.py {self.fa} --plots-format png -o {self.p_quast} --silent -t {int(threads/2)}")
        common.cmd2shell(cmd,sh)
        return sh



    def checkm(self):
        sh = os.path.join(self.p_sh,"3.2checkm.sh")
        cmd = []
        if self.checkm_flag:
            dir = os.path.dirname(self.fa)  #dir of fasta
            fa = os.path.basename(self.fa)
            cmd.append(f"source {b_act} denovo ")
            cmd.append(f"checkm lineage_wf -t {int(threads*3/4)} -q --tab_table -f {self.p_checkm}/checkm.txt -x {fa} {dir} {self.p_checkm}  >{self.p_checkm}/checkm.o 2>&1 ")
            cmd.append(f"{b_base}/python {p_bin}/2.ass_qc/parse_checkm_result.py -i {self.p_checkm}/storage/bin_stats_ext.tsv -o {self.p_checkm}/check.txt ")

        common.cmd2shell(cmd,sh)
        return sh



    def depth(self):
        sh = os.path.join(self.p_sh,"3.3depth_stat.sh")
        cmd = []
        if self.depth_flag:
            depth_fa = f"{self.p_depth}/{os.path.basename(self.fa)}"
            cmd.append(f"cp {self.fa} {self.p_depth}")
            cmd.append(f"{b_base}/bwa index -a bwtsw {depth_fa}")
            cmd.append(f"{b_base}/samtools faidx {depth_fa}")


            if self.fq2:  #PE
                cmd.append(f"mode=PE")
                cmd.append(f"{b_base}/bwa mem -t {int(threads*3/4)} {depth_fa} {self.fq1} {self.fq2} |{b_base}/samtools view -@ 24 -bS |{b_base}/samtools sort -@ 24 > {self.p_depth}/{self.id}.sort.bam")
            else:   #SE
                cmd.append(f"mode=SE")
                cmd.append(f"{b_base}/bwa mem -t {int(threads*3/4)} {depth_fa} {self.fq1} |{b_base}/samtools view -@ 24 -bS |{b_base}/samtools sort -@ 24 > {self.p_depth}/{self.id}.sort.bam")

            cmd.append(f"""#samtools
{b_base}/samtools depth -a {self.p_depth}/{self.id}.sort.bam >{self.p_depth}/{self.id}.bam.depth

# Insert size
if [ "$mode" == 'PE' ];then
    {picard} CollectInsertSizeMetrics --Histogram_FILE {self.p_depth}/picard_insertsize.pdf -I {self.p_depth}/{self.id}.sort.bam -O {self.p_depth}/insertsize.txt
    {b_base}/Rscript {p_bin}/2.ass_qc/insertsize.R {self.p_depth}/insertsize.txt {self.p_depth}
fi

# GC-depth png
{b_base}/python {p_bin}/2.ass_qc/depth_base_stat.py -l 2000 -g {depth_fa} -d {self.p_depth}/{self.id}.bam.depth -s {self.p_depth}/depth_base.stat
{b_base}/Rscript {p_bin}/2.ass_qc/depth_GC_plot.r -i {self.p_depth}/depth_base.stat -o {self.p_depth}/depth_base.stat.depth_GC
{b_base}/convert {self.p_depth}/depth_base.stat.depth_GC.pdf {self.p_depth}/depth_base.stat.depth_GC.png

# Uniformity
{b_base}/python {p_bin}/2.ass_qc/count_depth.py -i {self.p_depth}/{self.id}.bam.depth -ref {depth_fa} -o {self.p_depth}/uniformity.txt -a {self.p_depth}/uniformity_all.txt
""")

        common.cmd2shell(cmd,sh)
        return sh

    def align_ref(self):
        sh = os.path.join(self.p_sh,"3.4align_ref.sh")
        cmd = []
        if self.ref:
            ref_name = os.path.basename(self.ref)
            cmd.append(f"""
cd {self.p_ref}
cp {self.ref} {self.p_ref}
makeblastdb -in {self.p_ref}/{ref_name} -out ref -dbtype nucl
blastn -num_threads {int(threads/2)} -evalue 1e-10 -query {self.fa} -db {self.p_ref}/ref -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send pident length mismatch gapopen evalue bitscore' -out {self.p_ref}/ref_blastn.txt
{b_base}/python3 {p_bin}/tools/cal_blast.py -i {self.p_ref}/ref_blastn.txt --mismatch_rate 0.1 --identity_rate 0.8
""")
        common.cmd2shell(cmd,sh)
        return sh


# Option
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-id',required=True,type=click.STRING,help="样本编号")
@click.option('-fa',required=True,type=click.Path(),help="组装的序列 contigs/scaffolds fasta")
@click.option('-ref',required=False,type=click.Path(),default='None',help="病原体参考序列")
@click.option('-o', '--out',required=True,type=click.Path(),default='.',show_default=True,help="输出目录")
@click.option('-sh',required=True,type=click.Path(),default='.',show_default=True,help="Write shell dir")
@click.option('--quast',is_flag=True,help="运行quast质控")
@click.option('--checkm',is_flag=True,help="运行checkm质控")
@click.option('--depth_stat',is_flag=True,help="运行深度计算")
@click.option('-fq1',required=False,type=click.Path(),help="过滤后的fq1，用于计算深度")
@click.option('-fq2',required=False,type=click.Path(),default='None',help="过滤后的fq2，用于计算深度，单端则只需要输出fq1")


def main(id,fa,fq1,fq2,ref,out,sh,quast,checkm,depth_stat):
    #log
    logfile=f"{out}/log"
    os.makedirs(logfile,exist_ok=True)

    project = AssembleQC(
                id = id,
                fa = fa,
                fq1 = fq1,
                fq2 = fq2,
                ref = ref,
                out = out,
                sh = sh,
                quast =quast,
                checkm = checkm,
                depth_stat = depth_stat
            )
    project.make_dir()
    sh1 =project.quast()
    sh2 = project.checkm()
    sh3 = project.depth()
    sh4 = project.align_ref()
    shlist=[sh1,sh2,sh3,sh4]
    common.mul_pool(shlist, logfile)



if __name__ == "__main__":
    main()