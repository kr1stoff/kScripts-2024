#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: MingJia
# @Date:   2020-07-14 11:07:40
# @Last Modified by:   MingJia
# @Last Modified time: 2020-07-14 11:07:40
import logging

import click
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
__version__ = '1.0.0'

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    """
    Deal with excel files
    """
    pass


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input',
              required=True,
              type=click.Path(),
              help="The input excel file")
@click.option('-o', '--output',
              required=True,
              type=click.Path(),
              help="The out put table file")
@click.option('-s', '--separator',
              default="\t",
              show_default=True,
              help="The separator in table file")
@click.option('--header/--no-header',
              default=True,
              show_default=True,
              help="Whether read the first line as header")
def excel2table(input, output, separator, header):
    """
    Convert Excel file to Table file
    """
    statue_header = 1 if header else 0
    logging.info(f"Read the excel file {input}")
    df = pd.read_excel(input, header=statue_header)
    str_to_replace = ['\r\n', '\r', '\n']
    df.replace(to_replace=str_to_replace,
               value=' ',
               regex=True,
               inplace=True)
    logging.info(f"Out put the table file {output}")
    df.to_csv(output, sep=separator, index=False)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input',
              required=True,
              type=click.Path(),
              help="The input excel file")
@click.option('-o', '--output',
              required=True,
              type=click.Path(),
              help="The out put table file")
@click.option('-s', '--separator',
              default="\t",
              show_default=True,
              help="The separator in table file")
@click.option('--header/--no-header',
              default=True,
              show_default=True,
              help="Whether read the first line as header")
def table2excel(input, output, separator, header):
    """
    Convert Table file to Excel file
    """
    statue_header = 1 if header else 0
    logging.info(f"Read the table file {input}")
    df = pd.read_csv(input, sep=separator, header=statue_header)
    logging.info(f"Out put the excel file {output}")
    df.to_excel(output, index=False)


cli.add_command(excel2table)
cli.add_command(table2excel)

if __name__ == "__main__":
    cli()
