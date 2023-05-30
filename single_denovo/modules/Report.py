#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/23 15:00
import os
import yaml
import sys
import click

PATH = os.path.dirname(os.path.abspath(__file__))  
DIR = os.path.dirname(PATH)

sys.path.append(os.path.dirname(PATH))
import lib.common as common

# read YAML
p_cfg = os.path.join(DIR, "modules/cfg.yaml")
p_bin = os.path.join(DIR, "bin")
f_config = yaml.safe_load(open(p_cfg , mode='r',encoding='utf-8_sig').read())
b_base = f"{f_config['b_base']}"
p_GDHR = f"{f_config['p_GDHR']}"

class Report():
    def __init__(self,**kwargs):

        self.id = kwargs['id']
        self.out = os.path.abspath(kwargs['out'])
        self.analy_dir = os.path.dirname(self.out)
        self.p_sh = os.path.abspath(kwargs['sh'])
        self.kind = kwargs['kind']
        self.upload = '' if kwargs['upload'] == 'None' else os.path.abspath(kwargs['upload'])


    def make_dir(self):
        os.makedirs(self.p_sh ,exist_ok=True)
        os.makedirs(self.upload,exist_ok=True)
        if self.kind != 'virus':
            for var in ['1.qc','2.ass','3.pre','4.anno']:   #report dir
                tmp = os.path.join(self.out,var)
                os.makedirs(tmp,exist_ok=True)
        else:
            for var in ['1.qc','2.ass']:   #report dir
                tmp = os.path.join(self.out,var)
                os.makedirs(tmp,exist_ok=True)

    def run_report(self):
        sh  = os.path.join(self.p_sh,"5.report.sh")
        cmd = []
        
        cmd.append(f"cd {self.analy_dir}")
        cmd.append(f"""# 1.qc
cp -r {self.analy_dir}/1.qc/fastp/{self.id}.basic.stat.txt {self.out}/1.qc
if [ -f "{self.analy_dir}/1.qc/fastqc/{self.id}_1_fastqc.html" ]; then
    cp -r {self.analy_dir}/1.qc/fastqc/{self.id}_1_fastqc {self.out}/1.qc
    cp -r {self.analy_dir}/1.qc/fastqc/{self.id}_2_fastqc {self.out}/1.qc
    cp -r {self.analy_dir}/1.qc/fastqc/{self.id}_1_fastqc.html {self.out}/1.qc
    cp -r {self.analy_dir}/1.qc/fastqc/{self.id}_2_fastqc.html {self.out}/1.qc

    cp -r {self.analy_dir}/1.qc/filter_fastqc/{self.id}_1.clean_fastqc {self.out}/1.qc
    cp -r {self.analy_dir}/1.qc/filter_fastqc/{self.id}_2.clean_fastqc {self.out}/1.qc
    cp -r {self.analy_dir}/1.qc/filter_fastqc/{self.id}_1.clean_fastqc.html {self.out}/1.qc
    cp -r {self.analy_dir}/1.qc/filter_fastqc/{self.id}_2.clean_fastqc.html {self.out}/1.qc

elif [ -f "{self.analy_dir}/1.qc/fastqc/{self.id}_fastqc.html" ];then 
    cp -r {self.analy_dir}/1.qc/fastqc/{self.id}_fastqc {self.out}/1.qc
    cp -r {self.analy_dir}/1.qc/fastqc/{self.id}_fastqc.html {self.out}/1.qc

    cp -r {self.analy_dir}/1.qc/filter_fastqc/{self.id}_clean_fastqc {self.out}/1.qc
    cp -r {self.analy_dir}/1.qc/filter_fastqc/{self.id}_clean_fastqc.html {self.out}/1.qc
fi
""")
        if self.kind == 'bacteria' or self.kind == 'fungi':
            cmd.append(f"""# 2.ass
cp -r {self.analy_dir}/2.ass/deal_reads/cutfq.stat.txt        {self.out}/1.qc

# depollute
if  [ -f "{self.analy_dir}/2.ass/deal_reads/temp_depollute/depollute.report" ];then 
    cp -r {self.analy_dir}/2.ass/deal_reads/temp_depollute/abundance.txt  {self.out}/1.qc
fi

cp -r {self.analy_dir}/2.ass/spades/{self.id}_scaffolds_fa.stat.txt {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/spades/scaffolds.fasta    {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/spades/{self.id}_scaffolds.fasta    {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/spades/{self.id}_scaffolds.length.png    {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/quast            {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/checkm/check.txt {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/depth/depth_base.stat.depth_GC.png {self.out}/2.ass

if [ -f "{self.analy_dir}/1.qc/fastqc/{self.id}_1_fastqc.html" ]; then
    cp -r {self.analy_dir}/2.ass/depth/insertsize.png {self.out}/2.ass
fi

cp -r {self.analy_dir}/2.ass/depth/depth_base.stat    {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/depth/uniformity.txt     {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/depth/uniformity_all.txt {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/ref/ref_cov_stat.txt {self.out}/2.ass
""")
            cmd.append(f"""# 3.precict
cp -r {self.analy_dir}/3.pre/predict/predict_kind.txt {self.out}/3.pre
cp -r {self.analy_dir}/3.pre/predict/predict.length.png {self.out}/3.pre
""")
        if self.kind == 'bacteria':
            cmd.append(f"""
cp -r {self.analy_dir}/3.pre/island/island.txt {self.out}/3.pre
cp -r {self.analy_dir}/3.pre/phispy/prophage_coordinates.txt {self.out}/3.pre
cp -r {self.analy_dir}/3.pre/phispy/phage.fasta {self.out}/3.pre
""")
        if self.kind == 'bacteria' or self.kind == 'fungi':
            cmd.append(f"""# 4.anno
cp -r {self.analy_dir}/4.anno/VFDB/show_virulence_gene.txt  {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/VFDB/show_virulence_gene.xlsx {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/VFDB/top1_virulence_gene.txt  {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/VFDB/top1_virulence_gene.xlsx {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/VFDB/virulence_gene.faa {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/drug/detail_drug_resistance.txt  {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/drug/detail_drug_resistance.xlsx {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/drug/drug.faa {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/eggmapper/COG/all.COG.bar.png {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/eggmapper/COG/all.COG.class.xls {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/eggmapper/GO/GO_anno_stats.xls {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/eggmapper/GO/GO_anno.xls {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/eggmapper/GO/GO_anno_stats_level2.png {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/eggmapper/KEGG/KEGG_anno.txt {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/eggmapper/KEGG/KEGG_anno.stat.txt {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/eggmapper/KEGG/KEGG_anno_stats.png {self.out}/4.anno
""")

            cmd.append(f"""# 4.anno
# cp -r {self.analy_dir}/4.anno/pfam/pfam.txt {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/CAZy/CAZY.txt {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/CAZy/CAZy_anno_stats.png {self.out}/4.anno
cp -r {self.analy_dir}/4.anno/swissprot/swissprot_result.tsv {self.out}/4.anno
""")

        # 病毒有参拼接 
        elif self.kind == 'virus' and os.path.exists(f"{self.analy_dir}/2.ass/align/{self.id}_consensus.fa"):
            cmd.append(f"""
cp -r {self.analy_dir}/2.ass/align/{self.id}_bf12_sort.bam {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/align/bed.depth.stat.png {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/align/bam.stat.txt  {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/align/variants.xlsx {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/align/variants.txt  {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/align/{self.id}_consensus.fa  {self.out}/2.ass
""")

        # 病毒无参拼接 
        elif self.kind == 'virus' and os.path.exists(f"{self.analy_dir}/2.ass/spades/{self.id}_scaffolds.fasta"):
            cmd.append(f"""
cp -r {self.analy_dir}/2.ass/spades/{self.id}_scaffolds_fa.stat.txt {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/spades/scaffolds.fasta    {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/spades/{self.id}_scaffolds.fasta    {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/spades/{self.id}_scaffolds.length.png    {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/quast            {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/depth/depth_base.stat.depth_GC.png {self.out}/2.ass
if [ -f "{self.analy_dir}/1.qc/fastqc/{self.id}_1_fastqc.html" ]; then
    cp -r {self.analy_dir}/2.ass/depth/insertsize.png {self.out}/2.ass
fi

cp -r {self.analy_dir}/2.ass/depth/depth_base.stat    {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/depth/uniformity.txt     {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/depth/uniformity_all.txt {self.out}/2.ass
cp -r {self.analy_dir}/2.ass/ref/ref_cov_stat.txt {self.out}/2.ass
""")

        # src
        pipe_src = os.path.join(DIR,"lib","report","pipe_src")
        cmd.append(f"cp -r {p_GDHR}/src {self.out}")
        cmd.append(f"cp -r {pipe_src}/{self.kind}.png {self.out}/src/image")
        cmd.append(f"cp -r {pipe_src}/icon {self.out}/src/image")

        # zip report
        report_pl = os.path.join(DIR,"lib","report",f"{self.kind}_report.pl")
        
        cmd.append(f"{b_base}/perl {report_pl} {self.id} {self.analy_dir}")
        cmd.append(f"zip -qr {self.id}.zip {self.id}")

        #log
        cmd.append(f"cat {self.analy_dir}/log/*sh.e >{self.analy_dir}/log.txt")

        # Upload
        if self.upload:
            cmd.append(f"cp -r {self.analy_dir}/{self.id}.zip {self.upload}")
            cmd.append(f"cp -r {self.out} {self.upload}")
            # cmd.append(f"cp -r {self.analy_dir}/log.txt {self.upload}")

        common.cmd2shell(cmd,sh)
        return sh



# Option
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-id',required=True,type=click.STRING,help="样本编号")
@click.option('-kind', required=True,type=click.Choice(["bacteria", "fungi", "virus"]),help="选择病原体类型")
@click.option('-o', '--out',required=True,type=click.Path(),default='.',show_default=True,help="输出目录")
@click.option('-u', '--upload',required=False,type=click.Path(),default='None',help="上传目录")
@click.option('-sh',required=True,type=click.Path(),default='.',show_default=True,help="写shell脚本目录")

def main(id,out,sh,upload,kind):
    #log
    logfile=f"{os.path.dirname(out)}/log"
    os.makedirs(logfile,exist_ok=True)

    project = Report(
                id = id,
                out = out,
                kind = kind,
                sh = sh,
                upload = upload
            )
    project.make_dir()
    sh1 = project.run_report()
    common.prun((sh1, logfile))
    



if __name__ == "__main__":
    main()