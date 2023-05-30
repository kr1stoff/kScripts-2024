#!/home/earthtest/miniconda3/bin/python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/5 18:00
import os
import yaml
import sys
from subprocess import run

PATH = os.path.dirname(os.path.abspath(__file__))  
sys.path.append(os.path.dirname(PATH))
import lib.common as common
import lib.pipe as pipe


###获取需要读取的文件路径
DIR=os.path.dirname(PATH)
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
p_conf = os.path.join(DIR, "lib/config.yaml")
p_bin=f"{DIR}/bin"
p_src=f"{DIR}/lib/report"

###读取yaml文件
f_config = yaml.safe_load(open(p_conf , mode='r',encoding='utf-8_sig').read())
b_base   = f"{f_config['b_base']}"
b_denovo = f"{f_config['b_denovo']}"
p_db = f"{f_config['p_db']}"            #数据库(/sdbb/share/database)
b_act = f"{f_config['ACTIVATE']}"
kraken2 = f"{f_config['kraken2']}"
ivar = f"{f_config['ivar']}"

###参数
th = f"{f_config['threads']}" 
seqtk_num = f"{f_config['seqtk_num']}"  #抽取的reads数
minlen=f"{f_config['minlen']}"
mincov=f"{f_config['mincov']}"


#todo：定义类
class Denovo():
    def __init__(self,out,id,kind,upload):

        #变量
        self.id=id
        self.kind=kind
        self.d_analy= os.path.abspath(out)
        self.out   = os.path.abspath(os.path.join(out,id))
        self.upload= os.path.abspath(upload)

        self.p_shell= os.path.join(self.out,"0.shell")
        self.p_qc  = os.path.join(self.out,"1.qc")
        self.p_ass = os.path.join(self.out,"2.ass")
        self.p_pre = os.path.join(self.out,"3.pre")
        self.p_anno= os.path.join(self.out,"4.anno")
        self.p_res = os.path.join(self.out,id)


###########################################################################################################################################
#ps:创建分析目录
    def creat_env(self):
        kind_dic={
                "bacteria":{
                        "1.qc":["fastp","fastqc","kraken2","UMSI"],
                        "2.ass":["spades","quast","bwa","checkm"],
                        "3.pre":["prokka","island"],
                        "4.anno":["rgi","VFDB","emapper"],
                        self.id:["1.qc","2.ass","3.pre","4.anno"]
                    },
                "fungi":{
                        "1.qc":["fastp","fastqc","kraken2"],
                        "2.ass":["spades","quast","bwa","checkm"],
                        "3.pre":["genemark"],
                        "4.anno":["rgi","VFDB","emapper"],
                        self.id:["1.qc","2.ass","3.pre","4.anno"]
                    },
                "virus":{
                        "1.qc":["fastp","fastqc","kraken2"],
                        "2.ass":["ivar"],
                        self.id:["1.qc","2.ass"]
                    }
                }
        common.mkdirs(self.p_shell)    #创建0.shell
        env_dic=kind_dic[self.kind]
        for dirs ,subdirs in env_dic.items():
            for subdir in subdirs:
                f_path=os.path.join(self.out,dirs,subdir)
                common.mkdirs(f_path)



###########################################################################################################################################
#ps:定义部分全局变量
    def global_var(self):
        self.p_qc_res= os.path.join(self.p_res,"1.qc")
        self.p_fastp = os.path.join(self.p_qc,"fastp")

        self.p_ass_res= os.path.join(self.p_res,"2.ass")
        self.p_spades = os.path.join(self.p_ass,"spades")
        self.p_bwa    = os.path.join(self.p_ass,"bwa")
        self.p_anno_res = os.path.join(self.p_res,"4.anno")



###########################################################################################################################################
#ps:Fastqc
    def fastqc(self,fq1,fq2):
        p_fastqc= os.path.join(self.p_qc,"fastqc")
        p_sh = os.path.join(self.p_shell,"1.1fastqc.sh")
        
        l_fq1 ,l_fq2=pipe.fq_link(self.id,fq1,fq2,p_fastqc) #根据前缀创建soft-link
        cmd=[]
        cmd.append(f"""#fastqc
{b_base}/fastqc -q -t {th} {l_fq1} {l_fq2} -o {p_fastqc}
unzip -qn {p_fastqc}/*1_fastqc.zip -d {self.p_qc_res}
unzip -qn {p_fastqc}/*2_fastqc.zip -d {self.p_qc_res}
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh



###########################################################################################################################################
#ps:Fastp
    def fastp(self,fq1,fq2):
        p_sh = os.path.join(self.p_shell,"1.2fastp.sh")
        cmd=[]

        cmd.append(f"""#fastp
{b_base}/fastp -n 1 -l 50 -y -w 16 -i {fq1} -I {fq2} -j {self.p_fastp}/{self.id}.clean.fq.stat.json -h {self.p_fastp}/{self.id}.clean.html -o {self.p_fastp}/{self.id}_1.clean.fq.gz -O {self.p_fastp}/{self.id}_2.clean.fq.gz

{b_base}/python3 {p_bin}/1.qc/parse_fastp.json.py -id {self.id} -i {self.p_fastp}/{self.id}.clean.fq.stat.json
cp {self.p_fastp}/{self.id}.basic.stat.txt  {self.p_qc_res}
cp {self.p_fastp}/{self.id}.detail.stat.txt {self.p_qc_res}
{b_base}/python3 {p_bin}/1.qc/fq_qc_stat.py {self.id} {self.p_fastp}/{self.id}.clean.fq.stat.json {self.p_qc_res}/{self.id}_qc.txt
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh



###########################################################################################################################################
#ps:Seqtk
    def cut_reads(self):
        p_sh = os.path.join(self.p_shell,"1.3seqtk.sh")

        cmd=[]

        #细菌 seqtk抽reads，并统计
        if self.kind == "bacteria" :
            cmd.append(f"""#seqtk
{b_denovo}/seqtk sample -s 11 {self.p_fastp}/{self.id}_1.clean.fq.gz  {seqtk_num} >{self.p_fastp}/{self.id}_1.clean.cut.fq
{b_denovo}/seqtk sample -s 11 {self.p_fastp}/{self.id}_2.clean.fq.gz  {seqtk_num} >{self.p_fastp}/{self.id}_2.clean.cut.fq

{b_denovo}/iTools Fqtools stat -InFq {self.p_fastp}/{self.id}_1.clean.cut.fq -InFq {self.p_fastp}/{self.id}_1.clean.cut.fq -OutStat {self.p_fastp}/cutfq.stat
{b_base}/python3 {p_bin}/1.qc/parse.itools.stat.py -i {self.p_fastp}/cutfq.stat -o {self.p_fastp}/cutfq.stat.txt
cp {self.p_fastp}/cutfq.stat.txt {self.p_qc_res}
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh


###########################################################################################################################################
#ps:Kraken2
    def kraken2(self):
        p_kraken2 = os.path.join(self.p_qc,"kraken2")
        p_sh = os.path.join(self.p_shell,"1.4kraken2.sh")
        cmd=[]
        cmd.append(f"""#kraken2
{kraken2} --db {p_db}/kraken2 --threads {th} --output {p_kraken2}/kraken2.out --report {p_kraken2}/kraken2.report {self.p_fastp}/{self.id}_1.clean.cut.fq {self.p_fastp}/{self.id}_1.clean.cut.fq
cut -f3 {p_kraken2}/kraken2.out |taxonkit reformat -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F -P -I 1 >{p_kraken2}/karnken.out.lineage

""")
        common.cmd2shell(cmd,p_sh)
        return p_sh

###########################################################################################################################################
#ps:UMSI
    def umsi(self):
        p_umsi = os.path.join(self.p_qc,"UMSI")
        p_sh = os.path.join(self.p_shell,"1.5UMSI.sh")
        cmd=[]
        common.cmd2shell(cmd,p_sh)
        return p_sh
        #todo:补充。。。。。。。。。。。。。。。。。。。


###########################################################################################################################################
#ps:Spades
    def spades(self):
        p_sh = os.path.join(self.p_shell,"2.1spades.sh")

        cmd=[]
        cmd.append(f"source {b_act} denovo")
        if self.kind == "bacteria": #细菌抽reads组装
            cmd.append(f"spades.py -t {th} --careful -o {self.p_spades} -1 {self.p_fastp}/{self.id}_1.clean.cut.fq -2 {self.p_fastp}/{self.id}_2.clean.cut.fq >{self.p_spades}/spades.o ")

        elif self.kind == "fungi":
            cmd.append(f"spades.py -t {th} --careful -o {self.p_spades} -1 {self.p_fastp}/{self.id}_1.clean.fq.gz -2 {self.p_fastp}/{self.id}_2.clean.fq.gz >{self.p_spades}/spades.o ")

        #filter
        cmd.append(f"python3 {p_bin}/2.ass_qc/contig.filter.py -fa {self.p_spades}/scaffolds.fasta -l {minlen} -c {mincov} -p scaffolds")
        
        cmd.append(f"cp {self.p_spades}/filter_scaffolds.stat.txt {self.p_ass_res} ")
        cmd.append(f"cp {self.p_spades}/scaffolds.fasta {self.p_ass_res}/{self.id}_scaffolds.fasta")
        common.cmd2shell(cmd,p_sh)
        return p_sh



###########################################################################################################################################
#ps:Quast
    def quast(self):
        p_quast  = os.path.join(self.p_ass,"quast")
        p_sh = os.path.join(self.p_shell,"2.2quast.sh")

        cmd=[]
        cmd.append(f"source {b_act} denovo ")
        cmd.append(f"mv {self.p_spades}/scaffolds.fasta {self.p_spades}/raw_scaffolds.fasta ")
        cmd.append(f"mv {self.p_spades}/filter.scaffolds.fasta {self.p_spades}/scaffolds.fasta ")
        cmd.append(f"python3 {b_denovo}/quast.py {self.p_spades}/scaffolds.fasta --plots-format png -o {p_quast} -t {th}")
        cmd.append(f"cp -r {p_quast} {self.p_ass_res}")
        common.cmd2shell(cmd,p_sh)
        return p_sh



###########################################################################################################################################
#ps:Depth_stat
    def count_depth(self):
        p_sh = os.path.join(self.p_shell,"2.3depth.sh")
        cmd=[]
        cmd.append(f"""#build index
cp {self.p_spades}/scaffolds.fasta {self.p_bwa}
{b_base}/bwa index -a bwtsw {self.p_bwa}/scaffolds.fasta
{b_base}/samtools faidx {self.p_bwa}/scaffolds.fasta

{b_base}/bwa mem -t {th} {self.p_bwa}/scaffolds.fasta {self.p_fastp}/{self.id}_1.clean.cut.fq {self.p_fastp}/{self.id}_2.clean.cut.fq |{b_base}/samtools view -@ 24 -bS |{b_base}/samtools sort -@ 24 > {self.p_bwa}/{self.id}.sort.bam
{b_base}/samtools depth -a {self.p_bwa}/{self.id}.sort.bam > {self.p_bwa}/{self.id}.bam.depth

#Insert size
{b_base}/samtools view -h {self.p_bwa}/{self.id}.sort.bam |grep -v '^@' |awk -F '\\t' '$9>0 {{print $9}}' >{self.p_bwa}/fragment.length.txt

#GC-depth png
{b_base}/python {p_bin}/2.ass_qc/depth_base_stat.py -l 2000 -g {self.p_bwa}/scaffolds.fasta -d {self.p_bwa}/{self.id}.bam.depth -s {self.p_bwa}/depth_base.stat
{b_base}/Rscript {p_bin}/2.ass_qc/depth_GC_plot.r -i {self.p_bwa}/depth_base.stat -o {self.p_bwa}/depth_base.stat.depth_GC
{b_base}/convert {self.p_bwa}/depth_base.stat.depth_GC.pdf {self.p_bwa}/depth_base.stat.depth_GC.png
cp -r {self.p_bwa}/{{depth_base.stat.depth_GC.png,depth_base.stat}} {self.p_ass_res}

#uniformity
{b_base}/python {p_bin}/2.ass_qc/count_depth.py -i {self.p_bwa}/{self.id}.bam.depth -ref {self.p_bwa}/scaffolds.fasta -o {self.p_bwa}/uniformity.txt -a {self.p_bwa}/uniformity_all.txt
cp {self.p_bwa}/uniformity.txt {self.p_bwa}/uniformity_all.txt {self.p_ass_res}
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh



###########################################################################################################################################
#ps:Checkm
    def checkm(self):
        p_checkm = os.path.join(self.p_ass,"checkm")
        p_sh = os.path.join(self.p_shell,"2.4checkm.sh")

        cmd=[]
        cmd.append(f"source {b_act} denovo ")
        cmd.append(f"checkm lineage_wf -t {th} -q --tab_table -x scaffolds.fasta {self.p_spades} {p_checkm}  >{p_checkm}/checkm.o ")
        cmd.append(f"{b_base}/python {p_bin}/2.ass_qc/parse_checkm_result.py -i {p_checkm}/storage/bin_stats_ext.tsv -o {p_checkm}/check.txt ")
        cmd.append(f"cp -r {p_checkm}/check.txt {self.p_ass_res}")
        common.cmd2shell(cmd,p_sh)
        return p_sh



###########################################################################################################################################
#ps:Predict
    def predict(self):
        p_pre_res = os.path.join(self.p_res,"3.pre")
        p_sh  = os.path.join(self.p_shell,"3.pre.sh")
        Kind=self.kind.capitalize()

        cmd=[]
        cmd.append(f"source {b_act} denovo ")
########细菌用prokka预测
        if self.kind == "bacteria":
            self.p_prokka =  os.path.join(self.p_pre,"prokka")
            #注释参数所需字典
            self.kind_dic={
                "bacteria":{"pre_faa":f"{self.p_prokka}/predict_bacteria.faa",
                            "egg":"bact",
                            "merge_depth_cov":f"{self.p_prokka}/bed_bam.merge_depth_cov.txt"}
                        }

            cmd.append(f"""#prokka
cd {self.p_prokka}
prokka {self.p_spades}/scaffolds.fasta --prefix predict_{self.kind} --cpus {th} --outdir {self.p_prokka} --kingdom {Kind} --force --addgenes --quiet

if [ -f "{self.p_prokka}/predict_{self.kind}.gff" ]; then
    {b_base}/python3 {p_bin}/3.pre/chinese_prokka.py -i {self.p_prokka}/predict_{self.kind}.txt -o {self.p_prokka}/predict_kind.txt
    cp {self.p_prokka}/predict_kind.txt {p_pre_res}

    #bed
    grep -v '#' {self.p_prokka}/predict_{self.kind}.gff|grep 'CDS' |grep 'ID='|cut -f 1,4,5,9|awk -F';' '{{print $1}}'|sed 's/ID=//g' >{self.p_prokka}/predict_{self.kind}.bed
    bedtools coverage -mean -a {self.p_prokka}/predict_{self.kind}.bed -b {self.p_bwa}/{self.id}.sort.bam >{self.p_prokka}/bed_bam.merge.depth
    bedtools coverage -a {self.p_prokka}/predict_{self.kind}.bed -b {self.p_bwa}/{self.id}.sort.bam >{self.p_prokka}/bed_bam.merge.cov
    paste {self.p_prokka}/bed_bam.merge.depth {self.p_prokka}/bed_bam.merge.cov |cut -f 4,5,13 >{self.p_prokka}/bed_bam.merge_depth_cov.txt

    #gene_length
    {b_base}/perl {p_bin}/3.pre/gene_length.pl -t predict -k {self.id} {self.p_prokka}/predict_{self.kind}.fna
    cp {self.p_prokka}/predict_{self.kind}.tsv {p_pre_res}
    cp {self.p_prokka}/predict_{self.kind}.txt {p_pre_res}
    cp {self.p_prokka}/{self.id}.predict.length_distribution.png {p_pre_res}

fi
""")

########真菌用genemark预测
        elif self.kind == "fungi":
            self.p_gmark  =  os.path.join(self.p_pre,"genemark")
            #注释参数所需字典
            self.kind_dic={
                "fungi":{"pre_faa":f"{self.p_gmark}/predict_cds.faa",
                        "egg":"fungi",
                        "merge_depth_cov":f"{self.p_gmark}/bed_bam.merge_depth_cov.txt"}
                        }   

            cmd.append(f"""#genmark
cd {self.p_gmark}
{b_denovo}/gmes_petap.pl --ES --fungus --cores 60 --sequence {self.p_bwa}/scaffolds.fasta  >{self.p_gmark}/pre.log

if [ -f "{self.p_gmark}/genemark.gtf" ]; then
    echo "gtf ==> gff3 ==> cds.fa ==> cds.faa"
    gffread {self.p_gmark}/genemark.gtf -o ->{self.p_gmark}/predict.gff3
    mv {self.p_gmark}/genemark.gtf {self.p_gmark}/predict.gtf
    {b_base}/perl {p_bin}/tools/getGene.pl {self.p_gmark}/predict.gff3 {self.p_bwa}/scaffolds.fasta -type mrna >{self.p_gmark}/predict_cds.fa
    {b_base}/perl {p_bin}/tools/cds2aa.pl {self.p_gmark}/predict_cds.fa >{self.p_gmark}/predict_cds.faa

    grep 'CDS' {self.p_gmark}/predict.gff3|grep 'Parent='|cut -f 1,4,5,9|sed -e 's/Parent=//g' >{self.p_gmark}/predict.gff3.bed
    bedtools coverage -mean -a {self.p_gmark}/predict.gff3.bed -b {self.p_bwa}/{self.id}.sort.bam >{self.p_gmark}/bed_bam.merge.depth
    bedtools coverage -a {self.p_gmark}/predict.gff3.bed -b {self.p_bwa}/{self.id}.sort.bam >{self.p_gmark}/bed_bam.merge.cov
    paste {self.p_gmark}/bed_bam.merge.depth {self.p_gmark}/bed_bam.merge.cov |cut -f 4,5,13 >{self.p_gmark}/bed_bam.merge_depth_cov.txt

    {b_base}/python3 {p_bin}/3.pre/gene_mark_count.py {self.p_gmark}/predict.gtf {p_pre_res}/gene_count.txt
    {b_base}/perl {p_bin}/3.pre/gene_length.pl -t predict -k {self.id} {self.p_gmark}/data/training.fna
    cp  {self.p_gmark}/predict.gtf  {p_pre_res}
    cp  {self.p_gmark}/data/{self.id}.predict.length_distribution.png {p_pre_res}
fi
""")    
        common.cmd2shell(cmd,p_sh)
        return p_sh



###########################################################################################################################################
#ps:Island
    def island(self):
        cmd=[]
        p_sh  = os.path.join(self.p_shell,"3.1island.sh")
        if self.kind == 'bacteria': 
            p_pre_res = os.path.join(self.p_res,"3.pre")
            p_island  =  os.path.join(self.p_pre,"island")

            cmd.append(f"""#island_predict
if [ -f "{self.p_prokka}/predict_{self.kind}.gbk" ]; then
    source {b_act} denovo
    cd {p_island}
    {b_base}/python3 {p_bin}/3.pre/multi_Island.py {self.p_prokka}/predict_{self.kind}.gbk {p_island}
    cp -r {p_island}/island.txt {p_pre_res}
fi
""")

        common.cmd2shell(cmd,p_sh)
        return p_sh




###########################################################################################################################################
#ps:Virulence
    def vfdb(self):

        p_VFDB   = os.path.join(self.p_anno,"VFDB")
        p_sh  = os.path.join(self.p_shell,"4.1vfdb.sh")

        cmd=[]
        cmd.append(f"""#VFDB
source {b_act} denovo
cd {p_VFDB}
blastp -num_threads {th}  -evalue 1e-10 -query {self.kind_dic[self.kind]['pre_faa']} -db {p_db}/VFDB/VFDB_setB_pro/VFDB -outfmt '6 qseqid qlen sseqid sgi slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp' -out {p_VFDB}/virulence_gene.txt
if [ -f "{p_VFDB}/virulence_gene.txt" ]; then
    {b_base}/python3 {p_bin}/4.anno/VFDB_filter_blastp.py {p_VFDB}/virulence_gene.txt {p_db}/VFDB/VFDB_setB_pro/clean_VFDB_setB_des.txt {self.kind_dic[self.kind]['merge_depth_cov']} {p_VFDB}
    {b_base}/python3 {p_bin}/tools/get_blast_top.py -i {p_VFDB}/show_virulence_gene.txt -n 1 -o {p_VFDB}/top1_virulence_gene.txt
fi
cp -r {p_VFDB}/show_virulence_gene.txt  {self.p_anno_res}
cp -r {p_VFDB}/top1_virulence_gene.txt  {self.p_anno_res}
cp -r {p_VFDB}/show_virulence_gene.xlsx {self.p_anno_res}
cp -r {p_VFDB}/top1_virulence_gene.xlsx {self.p_anno_res}
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh



###########################################################################################################################################
#ps:Drug
    def drug(self):

        p_rgi = os.path.join(self.p_anno,"rgi")
        p_sh = os.path.join(self.p_shell,"4.2drug.sh")

        cmd=[]
        cmd.append(f"""#rgi
source {b_act} rgi
cd  {p_rgi}
rgi main -n {th} --input_type protein --clean -i {self.kind_dic[self.kind]['pre_faa']} --output_file {p_rgi}/drug
source {b_act} denovo
cut -f 1 {p_rgi}/drug.txt |cut -f1 -d ' '|sed 's/ORF_ID/Contig/g' > {p_rgi}/contig.id
cut -f 8,10,11,15,16,21 {p_rgi}/drug.txt > {p_rgi}/cut_drug.txt
paste {p_rgi}/contig.id {p_rgi}/cut_drug.txt >{p_rgi}/cut_id_drug.txt
{b_base}/python {p_bin}/4.anno/parse.rgi.py {p_rgi}/cut_id_drug.txt {self.kind_dic[self.kind]['merge_depth_cov']}
cp -r {p_rgi}/{{detail_drug_resistance.txt,detail_drug_resistance.xlsx}} {self.p_anno_res}
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh



###########################################################################################################################################
#ps:Eggmapper
    def eggmapper(self):
        p_egg = os.path.join(self.p_anno,"emapper")
        p_sh  = os.path.join(self.p_shell,"4.3eggmapper.sh")

        cmd=[]
        cmd.append(f"""#Eggmapper
source {b_act} denovo
cd {p_egg} 
emapper.py -m diamond --cpu {th} -d {self.kind_dic[self.kind]['egg']} --override -i {self.kind_dic[self.kind]['pre_faa']} --output {self.id} --output_dir {p_egg} --data_dir {p_db}/eggNOG
{b_base}/python {p_bin}/4.anno/EggnogParser.py -n {p_egg}/{self.id}.emapper.annotations -f {self.p_spades}/scaffolds.fasta -o {p_egg}
{b_base}/Rscript {p_bin}/4.anno/draw_COG.R {p_egg}/COG/all.COG.class.xls {p_egg}/COG
cp -r {p_egg}/{self.id}.emapper.annotations {self.p_anno_res}
cp -r {p_egg}/COG/all.COG.bar.png  {self.p_anno_res}
rm -rf {{emappertmp_dmdn_*,Rplots.pdf}}
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh


###########################################################################################################################################
#ps：病毒有参拼接
    def virus_ref(self,ref):
        p_ass_res = os.path.join(self.p_res,"2.ass")
        p_sh  = os.path.join(self.p_shell,"2.ass.sh")
        
        cmd=[]
        cmd.append(f"""#Virus consensus
cd {p_ass_res}
{b_base}/bwa index -a bwtsw {ref}
{b_base}/bwa mem -t 80 -Y {ref} {self.p_fastp}/{self.id}_1.clean.fq.gz {self.p_fastp}/{self.id}_2.clean.fq.gz |{b_base}/samtools view -@ 24 -hS -bF 12 |{b_base}/samtools sort -@ 24 -o {self.p_ass}/{self.id}_bf12_sort.bam
{b_base}/samtools depth -a {self.p_ass}/{self.id}_bf12_sort.bam >{self.p_ass}/bf12_sort.depth
{b_base}/samtools mpileup -aa -A -d 0 -Q 0 {self.p_ass}/{self.id}_bf12_sort.bam | {ivar} consensus -p consensus
{b_base}/samtools mpileup -aa -A -d 0 -B -Q 0 --reference {ref} {self.p_ass}/{self.id}_bf12_sort.bam |{ivar} variants -p variants -r {ref}

#Depth dic
python3 {p_bin}/2.ass_qc/depth_pic.py {self.p_ass}/bf12_sort.depth {self.p_ass}
python3 {p_bin}/2.ass_qc/ivar_var.py variants.tsv {self.p_ass}

cp -r {self.p_ass}/{{bed.depth.stat.png,bed.depth.stat.txt,bam.stat.txt}} {p_ass_res}
cp {self.p_ass}/consensus.fa {p_ass_res}/{self.id}_consensus.fa
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh

###########################################################################################################################################
#ps:Report
    def report(self):
        p_sh  = os.path.join(self.p_shell,"report.sh")
        cmd=[]

        #细菌/真菌
        if self.kind == "bacteria" or self.kind == "fungi":
            cmd.append(f"{b_base}/perl {p_src}/{self.kind}_report.pl {self.id} {self.out}")

        #病毒
        elif self.kind == "virus":
            cmd.append(f"{b_base}/perl {p_src}/virus_report.pl {self.id} {self.p_res}")

        #Upload报告
        cmd.append(f"""#Upload
cd {self.out}
cp -r {p_src}/{self.kind}_src {self.p_res}/src
zip -qr {self.id}.zip {self.id}
cp -r {self.id}.zip {self.upload}
cp -r {self.id} {self.upload}
cat {self.out}/log/*.e > {self.out}/log.txt
""")

        common.cmd2shell(cmd,p_sh)
        return p_sh
