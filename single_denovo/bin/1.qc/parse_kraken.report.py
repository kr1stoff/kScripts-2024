#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/10/26 16:01

import click
import pandas as pd
import os 

def parse_report(infile):
    out = os.path.dirname(os.path.abspath(infile))

    data = pd.read_csv(infile,
                sep='\t',
                encoding='utf-8',
                usecols=[0,1,3,4,5],
                names=['物种丰度占比(%)','序列数','分类级别','Taxid','物种名称'],
                dtype={"序列数": int ,"物种丰度占比":float})

    df = pd.DataFrame(data)
    df = df[df['分类级别'] == 'S']                  #筛选种
    df['物种名称'] = df['物种名称'].map(str.strip)  #去空格
    df.sort_values(by=["序列数"],ascending=False,inplace=True)

    new_df = df[['Taxid','物种名称','序列数','物种丰度占比(%)']]
    new_df.to_csv(f"{out}/abundance.txt",sep='\t',index=False)



## Options
@click.command(context_settings = dict(help_option_names=['-h', '--help']))
@click.option('-i','--infile',required=True,type=click.Path(),help="path of kraken report")

def main(infile):
    """ Parse kraken2 report abundance by species """
    parse_report(infile)




if __name__ == "__main__":
    main()