#!/usr/bin/env python

import os
import sys
from subprocess import run
import shutil
import logging
from pathlib import Path
import pdb
sys.path.append(str(Path(__file__).parents[0]))
import common


# 设置运行日志
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)

class EI1010():
    """
    工作流脚本处理类
        indict      任务-样本名字典
        task_id     任务id
        fh          log文件句柄
        cursor      MySQL启用中的游标
    """
    def __init__(self, indict, fh, cursor):
        self.indict = indict
        self.fh = fh
        self.cursor = cursor
        self.task_id = self.indict["task_id"]
        self.dir_analysis = f"/sdbb/Earth/Analysis/{self.task_id}"
        self.dir_results = f"/sdbb/Earth/AnalysisResults/{self.task_id}"
        self.samples = self.indict["samples"].keys() # 样本number列表
        self.get_ref_and_tree_dict()
        self.popen_returncode = 0 # 默认有异常,运行成功后更新

    def prepare_bed(self):
        """准备bed文件"""
        in_excel = self.indict["primer_seq"]
        out_bed = common.excel2bed(in_excel, f"{self.dir_analysis}/primers.bed")
        return out_bed

    def get_ref_and_tree_dict(self):
        """获取参考序列和进化树参考"""
        # 参考序列
        sql_cmd = "SELECT * FROM weiyuan.tb_dict_ref_gene"
        self.cursor.execute(sql_cmd)
        refs = self.cursor.fetchall()
        self.dict_ref = {str(ref["ref_id"]):ref["ref_seq"] for ref in refs}
        # 进化树
        sql_cmd = "SELECT * FROM weiyuan.tb_dict_evolutionary_tree"
        self.cursor.execute(sql_cmd)
        refs = self.cursor.fetchall()
        self.dict_tree = {str(ref["tree_id"]):ref["tree_ref_fa"] for ref in refs}

    def my_get_infastq(self, outdict):
        """获取样本表fastq, 写入outdict字典"""
        for name in self.indict["samples"]:
            outdict["samples"][name] = []
            outdict["samples"][name].append(self.indict["samples"][name]["sequencing_data1"])
            if self.indict["samples"][name]["sequencing_data2"] != "":  # 双端
                outdict["samples"][name].append(self.indict["samples"][name]["sequencing_data2"])

    def my_get_infasta(self, outdict):
        """获取样本表fasta, 写入outdict字典"""
        for name in self.indict["samples"]:
            outdict["samples"][name] = self.indict["samples"][name]["fasta"]

    def my_get_tree_ref(self, outdict):
        """
        获取进化树参考序列, 写入outdict字典. 
        1.进化树参考是文件夹,里面有很多FA; 2.可以是一个文件夹
        """
        tree_ids = self.indict["tree_data"].split(",")
        for tid in tree_ids:
            dir_tree_ref = self.dict_tree[tid]
            for rfile in os.listdir(dir_tree_ref):
                outdict["samples"][os.path.basename(rfile)] = os.path.join(dir_tree_ref, rfile)

    def my_get_ref(self, outdict):
        """获取参考序列,写入outdict"""
        ref_id = self.indict["reference_genome"]
        ref_fa = self.dict_ref[ref_id]
        outdict["reference"] = ref_fa # 通用流程不做核酸蛋白注释

    def my_get_pathogen_type(self, outdict):
        """获取病原类型,写入outdict"""
        dict_pathogen = {
            "病毒": "virus",
            "细菌": "bacteria",
            "真菌": "fungi"
        }
        outdict["kindom"] = dict_pathogen[self.indict["pathogen_type"]]

    def trace_tree_shell(self, tree_method, outdict):
        """给trace_tree准备shell"""
        inyaml = f"{self.dir_analysis}/trace_tree.yaml"
        common.dict2yaml(outdict, inyaml) 
        if tree_method == "WGS":
            shell = f"python3 /sdbb/share/pipeline/TracePathogen/main.py -p wgs -i {inyaml}"
        elif tree_method in ["SNP-FA", "SNP-FQ"]:
            shell = f"python3 /sdbb/share/pipeline/TracePathogen/main.py -p snp -i {inyaml}"
        elif tree_method == "CORE":
            shell = f"python3 /sdbb/share/pipeline/TracePathogen/main.py -p core -i {inyaml}"
        return shell

    def ngs_sars2_wga(self):
        """二代测序新冠全基因组扩增子分析路程"""
        # 准备
        dict_sars_wga = {
            "result_dir": "/sdbb/Earth/Analysis",
            "library": str(self.task_id),
            "samples": {},
            "global_sars2": False,
            "bed": ""
        }
        out_bed = self.prepare_bed()
        if out_bed != "": dict_sars_wga["bed"] = out_bed # bed提交了吗
        if self.indict["sars_cov_2_database"] == "是": dict_sars_wga["global_sars2"] = True # 新冠全球库做不做
        self.my_get_infastq(dict_sars_wga)
        yaml_sars_wga = f"{self.dir_analysis}/ngs_sars2_wga.yaml"
        common.dict2yaml(dict_sars_wga, yaml_sars_wga)
        # 运行
        shell = f"python /sdbb/share/pipeline/Sars_Cov2_Amplicon/main.py -i {yaml_sars_wga}"
        logging.debug(f"shell: {shell}")
        self.popen_returncode = common.wrap_popen(shell, self.task_id, self.cursor, self.fh)
        if self.popen_returncode == 0:
            common.lunx_upload(self.dir_analysis, self.dir_results, task_id=self.task_id)
            common.return2mysql(self.task_id, self.samples, self.cursor, return_fasta=True)
        else:
            print(f"popen_returncode: {self.popen_returncode}")

    def ngs_amplicon_general(self):
        """二代测序通用扩增子分析流程"""
        dict_amplicon_general = {
            "result_dir": "/sdbb/Earth/Analysis",
            "library": str(self.task_id),
            "samples": {},
            "bed": ""
        }
        # 获取ref
        self.my_get_ref(dict_amplicon_general) # 通用流程不做核酸蛋白注释
        # bed
        out_bed = self.prepare_bed()
        if out_bed != "": dict_amplicon_general["bed"] = out_bed # bed提交了吗
        dict_amplicon_general["trace"] = False  # 通用流程不做进化树
        dict_amplicon_general["kindom"] = "virus"
        self.my_get_infastq(dict_amplicon_general)
        yaml_amplicon_general = f"{self.dir_analysis}/ngs_amplicon_general.yaml"
        common.dict2yaml(dict_amplicon_general, yaml_amplicon_general)
        # 运行
        shell = f"python /sdbb/share/pipeline/AmpliconGP/main.py -i {yaml_amplicon_general}"
        logging.debug(f"shell: {shell}")
        self.popen_returncode = common.wrap_popen(shell, self.task_id, self.cursor, self.fh)
        if self.popen_returncode == 0:
            common.lunx_upload(self.dir_analysis, self.dir_results, task_id=self.task_id)
            common.return2mysql(self.task_id, self.samples, self.cursor, return_fasta=True)
        else:
            print(f"popen_returncode: {self.popen_returncode}")

    def trace_tree(self):
        """溯源进化树, 全基因组、核心基因、SNP"""
        """二代测序通用扩增子分析流程"""
        dict_trace_tree = {
            "result_dir": "/sdbb/Earth/Analysis",
            "library": str(self.task_id),
            "samples": {}
        }
        # WGS全长 1.没有进化树参考 2.有进化树参考
        if self.indict["tree_method"] == "全长序列" and self.indict["tree_data"] == "":
            self.my_get_infasta(dict_trace_tree)
            shell = self.trace_tree_shell("WGS", dict_trace_tree)
        elif self.indict["tree_method"] == "全长序列" and self.indict["tree_data"] != "":
            self.my_get_infasta(dict_trace_tree)
            self.my_get_tree_ref(dict_trace_tree)
            shell = self.trace_tree_shell("WGS", dict_trace_tree)
        # 基因片段，默认和全基因组一样
        elif self.indict["sequence_type"] == "基因片段":
            self.my_get_infasta(dict_trace_tree)
            shell = self.trace_tree_shell("WGS", dict_trace_tree)
        # 核心基因 细菌真菌 1.没有进化树参考 2.有进化树参考
        elif self.indict["tree_method"] == "coregene" and self.indict["tree_data"] == "":
            self.my_get_infasta(dict_trace_tree)
            self.my_get_pathogen_type(dict_trace_tree)
            shell = self.trace_tree_shell("CORE", dict_trace_tree)
        elif self.indict["tree_method"] == "coregene" and self.indict["tree_data"] != "":
            self.my_get_infasta(dict_trace_tree)
            self.my_get_pathogen_type(dict_trace_tree)
            self.my_get_tree_ref(dict_trace_tree)
            shell = self.trace_tree_shell("CORE", dict_trace_tree)
        # SNP 1.fasta 2.fastq
        elif self.indict["tree_method"] == "SNP" and self.indict["input_fq_or_fa"] == "fasta":
            self.my_get_infasta(dict_trace_tree)
            self.my_get_pathogen_type(dict_trace_tree)
            self.my_get_ref(dict_trace_tree)
            shell = self.trace_tree_shell("SNP-FA", dict_trace_tree)
        elif self.indict["tree_method"] == "SNP" and self.indict["input_fq_or_fa"] == "fastq":
            self.my_get_infastq(dict_trace_tree)
            self.my_get_pathogen_type(dict_trace_tree)
            self.my_get_ref(dict_trace_tree)
            shell = self.trace_tree_shell("SNP-FQ", dict_trace_tree)
        else:
            logging.error(f"不存在的进化树参数组合!!! sequence_type:{self.indict['sequence_type']}, "
            f"tree_method:{self.indict['tree_method']}, tree_data: {self.indict['tree_data']}")
            return None
        # 运行
        logging.debug(f"shell: {shell}")
        self.popen_returncode = common.wrap_popen(shell, self.task_id, self.cursor, self.fh)
        if self.popen_returncode == 0:
            common.lunx_upload(self.dir_analysis, self.dir_results, task_id=self.task_id, batch=True)
            common.return2mysql_batch(self.task_id, self.samples, self.cursor)
        else:
            print(f"popen_returncode: {self.popen_returncode}")
