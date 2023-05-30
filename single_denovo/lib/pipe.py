#!/home/earthtest/miniconda3/bin/python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/3/11 14:06

import os
import sys
import logging

logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M')
sys.path.append(os.path.dirname(__file__))
import common


#根据fq后缀情况，创建软连接
def pe_fq_link(id,fq1,fq2,out):
    s_fix=str(os.path.basename(fq1))
    if s_fix.endswith(".gz"):
        l_fq1, l_fq2=f"{out}/{id}_1.fastq.gz" ,f"{out}/{id}_2.fastq.gz"
        common.creat_link(fq1,l_fq1)
        common.creat_link(fq2,l_fq2)

    elif s_fix.endswith((".fq",".fastq")):
        l_fq1, l_fq2=f"{out}/{id}_1.fastq" ,f"{out}/{id}_2.fastq"
        common.creat_link(fq1,l_fq1)
        common.creat_link(fq2,l_fq2)

    else:
        logging.debug(f"Fastq's suffix can't read ----({s_fix}) ")
        l_fq1 ,l_fq2=f"{out}/{id}_1.fastq.gz" ,f"{out}/{id}_1.fastq.gz"
        common.creat_link(fq1,l_fq1)
        common.creat_link(fq2,l_fq2)

    return l_fq1 ,l_fq2


def se_fq_link(id,fq,out):
    s_fix=str(os.path.basename(fq))
    if s_fix.endswith(".gz"):
        l_fq=f"{out}/{id}.fastq.gz" 
        common.creat_link(fq,l_fq)


    elif s_fix.endswith((".fq",".fastq")):
        l_fq=f"{out}/{id}.fastq" 
        common.creat_link(fq,l_fq)


    else:
        logging.debug(f"Fastq's suffix can't read ----({s_fix}) ")
        l_fq =f"{out}/{id}.fastq.gz" 
        common.creat_link(fq,l_fq)


    return l_fq