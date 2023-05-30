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
p_db = f"{f_config['p_db']}"
cpu ,threads = common.get_cfg()

class Predict():
    def __init__(self,**kwargs):

        self.id = kwargs['id']
        self.fa = kwargs['fa']
        self.kind = kwargs['kind']
        self.virulence_flag = kwargs['virulence']
        self.drug_flag = kwargs['drug']
        self.eggmapper_flag = kwargs['eggmapper']
        self.CAZy_flag = kwargs['cazy']
        self.pfam_flag = kwargs['pfam']
        self.swiss_prot_flag = kwargs['swiss_prot']
        
        self.out = os.path.abspath(kwargs['out'])
        self.p_sh = os.path.abspath(kwargs['sh'])

        self.merge_file = os.path.abspath(kwargs['merge_file']) if kwargs['merge_file'] else ''


    def make_dir(self):
        os.makedirs(self.p_sh ,exist_ok=True)
        if self.virulence_flag:
            self.p_VFDB = os.path.join(self.out,'VFDB')
            os.makedirs(self.p_VFDB ,exist_ok=True)
        if self.drug_flag:
            self.p_drug = os.path.join(self.out,'drug')
            os.makedirs(self.p_drug ,exist_ok=True)
        if self.eggmapper_flag:
            self.p_eggmapper= os.path.join(self.out,'eggmapper')
            os.makedirs(self.p_eggmapper ,exist_ok=True)
        if self.CAZy_flag:
            self.p_CAZy= os.path.join(self.out,'CAZy')
            os.makedirs(self.p_CAZy ,exist_ok=True)
        if self.pfam_flag:
            self.p_pfam= os.path.join(self.out,'pfam')
            os.makedirs(self.p_pfam ,exist_ok=True)
        if self.swiss_prot_flag:
            self.p_swissprot= os.path.join(self.out,'swissprot')
            os.makedirs(self.p_swissprot ,exist_ok=True)


    def VFDB_anno(self):
        sh  = os.path.join(self.p_sh,"4.1vfdb.sh")
        cmd=[]
        if self.virulence_flag:
            cmd.append(f"""#VFDB
source {b_act} denovo
cd {self.p_VFDB}
time blastp -num_threads {int(threads/2)}  -evalue 1e-10 -query {self.fa} -db {p_db}/VFDB/VFDB_setB_pro/VFDB -outfmt '6 qseqid qlen sseqid sgi slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' -out {self.p_VFDB}/virulence_gene.txt

if [ -f "{self.p_VFDB}/virulence_gene.txt" ]; then
    {b_base}/python3 {p_bin}/4.anno/VFDB_filter_blastp.py {self.p_VFDB}/virulence_gene.txt {p_db}/VFDB/VFDB_setB_pro/clean_VFDB_setB_des.txt {self.merge_file} {self.p_VFDB}
    {b_base}/python3 {p_bin}/tools/get_blast_top.py -i {self.p_VFDB}/show_virulence_gene.txt -n 1 -o {self.p_VFDB}/top1_virulence_gene.txt
    
    # Output fasta
    awk 'NR>1' {self.p_VFDB}/show_virulence_gene.txt |cut -f1 |sort |uniq >{self.p_VFDB}/show_virulence_gene.id
    {b_base}/perl {p_bin}/tools/fishInWinter.pl --bformat table --fformat fasta {self.p_VFDB}/show_virulence_gene.id {self.fa} >{self.p_VFDB}/virulence_gene.faa

fi

""")
        common.cmd2shell(cmd,sh)
        return sh

    def drug_anno(self):
        sh  = os.path.join(self.p_sh,"4.2drug.sh")
        cmd =[] 
        if self.drug_flag:
            cmd.append(f"""
# rgi
source {b_act} rgi
cd  {self.p_drug}
time rgi main -n {int(threads/2)} --input_type protein --clean -i {self.fa} --output_file {self.p_drug}/drug
source {b_act} denovo
cut -f 1 {self.p_drug}/drug.txt |cut -f1 -d ' '|sed 's/ORF_ID/Contig/g' > {self.p_drug}/contig.id
cut -f 8,10,11,15,16,21 {self.p_drug}/drug.txt > {self.p_drug}/cut_drug.txt
paste {self.p_drug}/contig.id {self.p_drug}/cut_drug.txt >{self.p_drug}/cut_id_drug.txt
{b_base}/python {p_bin}/4.anno/parse.rgi.py {self.p_drug}/cut_id_drug.txt {self.merge_file}

# Output fasta
awk 'NR>1' {self.p_drug}/drug.txt |cut -f 1,20  |awk -F'\\t' '{{print ">"$1 "\\n"$2}}' >{self.p_drug}/drug.faa
""")
        common.cmd2shell(cmd,sh)
        return sh
    
    def eggmapper(self):
        sh  = os.path.join(self.p_sh,"4.3eggmapper.sh")
        cmd =[] 
        if self.eggmapper_flag:
            Kind = 'bact' if self.kind == 'bacteria' else 'fungi'
            cmd.append(f"""#Eggmapper
source {b_act} denovo
cd {self.p_eggmapper} 
time emapper.py -m diamond --cpu {int(cpu*2/3)} -d {Kind} --override -i {self.fa} --output {self.id} --output_dir {self.p_eggmapper} --data_dir {p_db}/eggNOG
if [ -f "{self.p_eggmapper}/{self.id}.emapper.annotations" ]; then
    {b_base}/python {p_bin}/4.anno/EggnogParser.py -n {self.p_eggmapper}/{self.id}.emapper.annotations -f {self.fa} -o {self.p_eggmapper}
    # COG
    {b_base}/Rscript {p_bin}/4.anno/draw_COG.R {self.p_eggmapper}/COG/all.COG.class.xls {self.p_eggmapper}/COG
    
    # GO
    python3 {p_bin}/4.anno/GO_anno.py -i {self.p_eggmapper}/{self.id}.emapper.annotations -db {p_db}/GO/go-basic.obo
    {b_base}/Rscript {p_bin}/4.anno/draw_GO.R {self.p_eggmapper}/GO/GO_anno_stats.xls {self.p_eggmapper}/GO

    # KEGG
    python3 {p_bin}/4.anno/KEGG_anno.py -i {self.p_eggmapper}/{self.id}.emapper.annotations -db {p_db}/KEGG/db
    {b_base}/Rscript {p_bin}/4.anno/draw_KEGG.R {self.p_eggmapper}/KEGG/KEGG_anno.stat.txt {self.p_eggmapper}/KEGG

fi
rm -rf {{emappertmp_dmdn_*,Rplots.pdf}}
""")
        common.cmd2shell(cmd,sh)
        return sh
    
    def CAZy_anno(self):
        sh  = os.path.join(self.p_sh,"4.4CAZy.sh")
        p_des = f"{p_db}/CAZY/CAZyDB.07302020.fam-activities-hmm.txt"
        cmd =[] 
        if self.CAZy_flag:
            if not os.path.exists(p_des):
                print(f"Database error: No such file [ {p_des} ]")
                exit(-1)

            cmd.append(f"""#CAZy
source {b_act} denovo
cd {self.p_CAZy}
time run_dbcan {self.fa} protein --db_dir {p_db}/CAZY  --out_dir {self.p_CAZy}
csvtk join -t -L -f 1 {self.p_CAZy}/hmmer.out {p_des} |csvtk round -t -f 'Coverage' -n 2 |csvtk mutate -t -f 1 -p '^(\D+)' -n 'kind' >{self.p_CAZy}/CAZY.txt
csvtk freq -t -f 'kind' {self.p_CAZy}/CAZY.txt > {self.p_CAZy}/CAZY.kind_stat.txt
{b_base}/Rscript {p_bin}/4.anno/draw_CAZy.R {self.p_CAZy}/CAZY.kind_stat.txt {self.p_CAZy}
""")
        common.cmd2shell(cmd,sh)
        return sh


    def pfam_anno(self):
        sh  = os.path.join(self.p_sh,"4.5pfam.sh")
        cmd =[] 
        if self.pfam_flag:
            cmd.append(f"""#pfam
source {b_act} denovo
cd {self.p_pfam}
time pfam_scan.pl -cpu {int(cpu*2/3)} -fasta {self.fa} -dir {p_db}/Pfam -outfile {self.p_pfam}/pfam.out
grep -v '^#' {self.p_pfam}/pfam.out >{self.p_pfam}/pfam.txt
sed -i -e '1i seq_id\\talignment_start\\talignment_end\\tenvelope_start\\tenvelope_end\\thmm_acc\\thmm_name\\ttype\\thmm_start\\thmm_end\\thmm_length\\tbit_score\\tE-value\\tsignificance\\tclan' -e '2d' {self.p_pfam}/pfam.txt
""")
        common.cmd2shell(cmd,sh)
        return sh

    def swiss_prot_anno(self):
        Kind = self.kind.capitalize()
        sh  = os.path.join(self.p_sh,"4.6swiss_prot.sh")
        cmd =[] 
        if self.swiss_prot_flag:
            cmd.append(f"{b_base}/python3 {p_bin}/4.anno/UniProt/swissprot_db.py analysing -i {self.fa} -d {p_db}/UniProt/2020-09-23/swissprot_db --softw_align diamond --db_select {Kind} -o {self.p_swissprot}/swissprot_result.tsv")
        common.cmd2shell(cmd,sh)
        return sh


# Option
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-id',required=True,type=click.STRING,help="样本编号")
@click.option('-fa',required=True,type=click.Path(),help="预测的蛋白质faa序列")
@click.option('-o', '--out',required=True,type=click.Path(),default='.',show_default=True,help="输出目录")
@click.option('-k','--kind',required=True,type=click.Choice(['bacteria','fungi']),help="选择病原体类型")
@click.option('-sh',required=True,type=click.Path(),default='.',show_default=True,help="写shell脚本目录")
@click.option('--virulence',is_flag=True,help="运行毒力因子注释")
@click.option('--drug',is_flag=True,help="运行耐药基因注释")
@click.option('--eggmapper',is_flag=True,help="运行eggmapper注释")
@click.option('--pfam',is_flag=True,help="运行pfam注释")
@click.option('--cazy',is_flag=True,help="运行CAZy注释")
@click.option('--swiss_prot',is_flag=True,help="运行swiss_prot注释")
@click.option('-m','--merge_file',required=False,type=click.Path(),help="深度覆盖度文件merge_depth_cov")



def main(id,fa,out,kind,sh,virulence,drug,eggmapper,cazy,pfam,swiss_prot,merge_file):
    #log
    logfile=f"{out}/log"
    os.makedirs(logfile,exist_ok=True)

    project = Predict(
                id = id,
                fa = fa,
                out = out,
                kind = kind,
                sh = sh,
                virulence = virulence,
                drug =drug,
                merge_file = merge_file,
                eggmapper = eggmapper,
                cazy = cazy,
                pfam = pfam,
                swiss_prot = swiss_prot
            )
    project.make_dir()

    sh1 = project.VFDB_anno()
    sh2 = project.drug_anno()
    sh3 = project.eggmapper()
    sh4 = project.CAZy_anno()
    #sh5 = project.pfam_anno() #pfam注释过久，半个小时以上
    sh6 = project.swiss_prot_anno()

    shlist=[sh1,sh2,sh3,sh4,sh6]
    common.mul_pool(shlist, logfile)




if __name__ == "__main__":
    main()