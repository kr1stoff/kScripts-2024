#!/usr/bin/env python
# @CreateTime       : 2023/5/19
# @Author           : mengxf
# @version          : v1.1
# @LastModified     : 2023/5/30
# @Description      : 基于primer3的qPCR设计脚本, 替代Windows版PrimerExpress3.0.1实现自动化引物&探针设计. 
# @From             : 脚本基础来自 ChaoboRen/TracePlatform 

import logging
import sys
from pathlib import Path
import click
from lib.qPCR import Qprimer
from lib.Batch import BatchQPCR


logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.group()
def cli():
    """qPCR引物设计软件合集."""

@cli.command()
@click.option('-f', '--fasta', required=True, type=click.Path(), help='输入基因区域FASTA文件.')
@click.option('-o', '--out', required=True, help="引物设计结果输出目录.")
@click.option('-r', '--size_range', default='70-110', show_default=True, help='设计产物长度范围.')
@click.option('-n', '--dryrun', is_flag=True, help='不执行程序,只生成shell文件.')
def qpcr(out, fasta, size_range, dryrun):
    """单区域qPCR设计流程, 输入GeneRegionFASTA, 输出设计结果."""
    fasta = Path(fasta).absolute()
    d_out = Path(out).absolute()
    project = Qprimer(fasta, d_out, dryrun)
    project.format()
    project.internal_design(size_range=size_range)
    project.finish()


@cli.command()
@click.option('-r', '--reference', required=True, type=click.Path(), help='输入考基因组FASTA文件.')
@click.option('-b', '--bed', required=True, type=click.Path(), help='输入保守区域BED文件.')
@click.option('-o', '--out', required=True, help="批量设计结果输出目录.")
def qpcr_batch(reference, bed, out):
    """批量qPCR设计, 输入参考基因组FASTA和保守区域BED文件."""
    mypath = Path(__file__).resolve() #本脚本路径
    Path(out).mkdir(exist_ok=True, parents=True)
    batch = BatchQPCR(mypath, reference, bed, out)
    batch.get_region_fasta()
    batch.design1()
    batch.design2()
    batch.get_exp_txt()


if __name__ == '__main__':
    cli()
