#!/home/earthtest/miniconda3/bin/python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/26 16:00
import click
import os
import yaml
import sys

PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(PATH))
import lib.common as common

Module = os.path.join(PATH,'modules')
p_cfg = os.path.join(PATH,'modules','cfg.yaml')
f_config = yaml.safe_load(open(p_cfg, mode='r',encoding='utf_8_sig').read())
python = f"{f_config['python3']}"

class RunDenovo():
    def __init__(self,**kwargs):

        self.id = kwargs['id']
        self.kind = kwargs['kind']
        self.out = os.path.abspath(kwargs['out'])

        self.fq1 = os.path.abspath(kwargs['fq1'])
        self.fq2 = os.path.abspath(kwargs['fq2']) if (kwargs['fq2'] != 'None') else None
        self.upload = os.path.abspath(kwargs['upload']) if (kwargs['upload'] != 'None') else None
        self.ref = os.path.abspath(kwargs['ref']) if (kwargs['ref'] != 'None') else None
        self.if_pollute = '--if_pollute' if (kwargs['if_pollute']) else ''
        self.if_denovo = '--if_denovo' if (kwargs['if_denovo']) else ''



    def make_dir(self):
        self.analy_dir = os.path.join(self.out,self.id)
        self.var_dic = {}
        if self.kind != 'virus':
            for var in [self.id,'0.shell','1.qc','2.ass','3.pre','4.anno']:
                tmp = os.path.join(self.analy_dir,var)
                os.makedirs(tmp,exist_ok=True)
                # creat report file
                if (var != self.id) and (var != '0.shell'):  
                    self.var_dic[str(var.split('.')[1])] = tmp
        else:
            for var in [self.id,'0.shell','1.qc','2.ass']:
                tmp = os.path.join(self.analy_dir,var)
                os.makedirs(tmp,exist_ok=True)
                if (var != self.id) and (var != '0.shell'):  
                    self.var_dic[str(var.split('.')[1])] = tmp

    def run_pipeline(self):
        
        cmd = []
        sh  = os.path.join(self.analy_dir,'0.shell','work.sh')

        clean_fq1 =os.path.join(self.var_dic['qc'],'fastp',f"{self.id}_1.clean.fq.gz")
        clean_fq2 =os.path.join(self.var_dic['qc'],'fastp',f"{self.id}_2.clean.fq.gz") if self.fq2 else None
        
        if self.kind == 'bacteria' or self.kind == 'fungi':
            ass_fa = os.path.join(self.var_dic['ass'],f"spades/{self.id}_scaffolds.fasta")
            bam = os.path.join(self.var_dic['ass'],'depth',f"{self.id}.sort.bam")
            pre_faa = os.path.join(self.var_dic['pre'],'predict/predict.faa')
            merge_file = os.path.join(self.var_dic['pre'],"predict/bed_bam.merge_depth_cov.txt")

        basic_argvs = f"-id {self.id} -sh {self.analy_dir}/0.shell"

        if self.kind == 'bacteria' or self.kind == 'fungi':
            # 1.qc
            cmd.append(f"{python} {Module}/SeqQc.py {basic_argvs} -o {self.var_dic['qc']} -fq1 {self.fq1} -fq2 {self.fq2}")
            # 2.assemble
            cmd.append(f"{python} {Module}/Assemble.py {basic_argvs} {self.if_pollute} --ass --cut_reads 9 -k {self.kind} -o {self.var_dic['ass']} -fq1 {clean_fq1} -fq2 {clean_fq2}")
            # 2.1assemble_qc
            cmd.append(f"{python} {Module}/AssQC.py {basic_argvs} -ref {self.ref} -o {self.var_dic['ass']} --quast --checkm --depth_stat -fq1 {clean_fq1} -fq2 {clean_fq2} -fa {ass_fa}")
            # 3.predict
            cmd.append(f"{python} {Module}/Predict.py {basic_argvs} -fa {ass_fa} -o {self.var_dic['pre']} -k {self.kind} -b {bam}")
            # 4.anno
            cmd.append(f"{python} {Module}/Anno.py {basic_argvs} -fa {pre_faa} -o {self.var_dic['anno']} -k {self.kind} --virulence --drug --eggmapper --cazy --pfam --swiss_prot -m {merge_file}")
            # 5.report
            cmd.append(f"{python} {Module}/Report.py {basic_argvs} -kind {self.kind} -o {self.analy_dir}/{self.id} -u {self.upload}")


        # 病毒有参拼接
        elif self.kind == 'virus' and (not self.if_denovo):
            # 1.qc
            cmd.append(f"{python} {Module}/SeqQc.py {basic_argvs} -o {self.var_dic['qc']} -fq1 {self.fq1} -fq2 {self.fq2}")
            # 2.align
            cmd.append(f"{python} {Module}/Virus.py {basic_argvs} -o {self.var_dic['ass']} -fq1 {clean_fq1} -fq2 {clean_fq2} -ref {self.ref}")
            # 3.report
            cmd.append(f"{python} {Module}/Report.py {basic_argvs} -kind {self.kind} -o {self.analy_dir}/{self.id} -u {self.upload}")


        # 病毒无参组装
        elif self.kind == 'virus' and self.if_denovo:
             # 1.qc
            cmd.append(f"{python} {Module}/SeqQc.py {basic_argvs} -o {self.var_dic['qc']} -fq1 {self.fq1} -fq2 {self.fq2}")
            # 2.assemble
            cmd.append(f"{python} {Module}/Assemble.py {basic_argvs} {self.if_pollute} --ass --cut_reads 9 -k {self.kind} -o {self.var_dic['ass']} -fq1 {clean_fq1} -fq2 {clean_fq2}")
            # 2.1assemble_qc
            ass_fa = os.path.join(self.var_dic['ass'],f"spades/{self.id}_scaffolds.fasta")
            cmd.append(f"{python} {Module}/AssQC.py {basic_argvs} -ref {self.ref} -o {self.var_dic['ass']} --quast --depth_stat -fq1 {clean_fq1} -fq2 {clean_fq2} -fa {ass_fa}")
            # 3.report
            cmd.append(f"{python} {Module}/Report.py {basic_argvs} -kind {self.kind} -o {self.analy_dir}/{self.id} -u {self.upload}")


        common.cmd2shell(cmd,sh)
        return sh





####参数################################################################################################################
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-id', required=True,type=click.STRING,help="样本名称")
@click.option('-fq1', required=True,type=click.Path(exists=True),help="输入fq1")
@click.option('-fq2',required=False,type=click.Path(),help="输出fq2，单端仅需输入fq1")
@click.option('--if_pollute', is_flag=True,help="是否选择去污染")
@click.option('-out', required=True,type=click.Path(),default='.',show_default=True,help="输出目录")
@click.option('-kind', required=True,type=click.Choice(["bacteria", "fungi", "virus"]),help="选择病原体类型")
@click.option('-ref', required=False,type=click.Path(),default='None',help="参考基因组")
@click.option('-upload', required=False,type=click.Path(),default='None',help="上传目录")
@click.option('--if_denovo',is_flag=True,help="病毒是否进行无参组装")



def main(id,fq1,fq2,out,kind,ref,upload,if_pollute,if_denovo):
    """创建单菌流程中，需要的分析目录"""
    #log
    logfile=f"{out}/{id}/log"
    os.makedirs(logfile,exist_ok=True)

    project = RunDenovo(
                id = id,
                fq1 = fq1,
                fq2 = fq2,
                out = out,
                kind = kind,
                ref = ref,
                upload = upload,
                if_pollute = if_pollute,
                if_denovo = if_denovo
            )
    project.make_dir()
    sh = project.run_pipeline()
    common.prun((sh, logfile))

    #cp log
    os.system(f"cp -r {out}/{id}/*/log/* {out}/{id}/log")




if __name__ == "__main__":
    main()