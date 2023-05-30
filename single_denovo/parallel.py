#!/usr/bin/env python
import os
import logging
import yaml
import click
import lib.common as common
import sys

logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)
PATH = os.path.dirname(os.path.abspath(__file__))

##YAML
p_cfg = os.path.join(PATH,'modules','cfg.yaml')
p_src = os.path.join(PATH, "lib/report/Tree_Typing/src")
f_config = yaml.safe_load(open(p_cfg, mode='r',encoding='utf_8_sig').read())

###读取配置文件的变量
b_base=f"{f_config['b_base']}"
jobs=int(f"{f_config['jobs']}")

###软件 & 脚本
TREE = f"{f_config['p_tree']}"
TYPING = f"{f_config['p_typing']}"
python = f"{f_config['python3']}"
perl = f"{f_config['perl']}"
p_GDHR = f"{f_config['p_GDHR']}"

class MultiRun():
    def __init__(self,out,upload,kind):
        self.kind = kind
        self.out = os.path.abspath(out)
        self.upload = os.path.abspath(upload)
        os.makedirs(f"{self.out}/0.shell",exist_ok=True)


    def single_denovo(self,table,ref):
        path_sh_list = []
        _yaml = yaml.safe_load(open(table, 'r').read())
        self.id_list = list(_yaml['samples'].keys())
        ref = 'None' if ref == 'None' else os.path.abspath(ref)

        #todo: 去污染暂定字段。。。。。？？
        ifpollute = "--if_pollute" if _yaml['sars_cov_2_database'] == "是" else ''
        #todo: 病毒无参流程字段。。。。？？
        ifdenovo = "--if_denovo" if _yaml['sars_cov_2_database'] == "是" else ''

        
        mode = _yaml['pe_or_se']    #单端/双端
        for id in self.id_list:
            cmd = []
            p_sh = f"{self.out}/0.shell/run_{id}.sh"
            path_sh_list.append(p_sh)   #记录每个样本的shell脚本路径
            fq1 = _yaml['samples'][id]['sequencing_data1']
            fq2 = _yaml['samples'][id]['sequencing_data2'] if mode == '双端' else None
            
            if self.kind=='bacteria':
                self.tree_mode='core'
                cmd.append(f"time {python} {PATH}/main.py -id {id} -kind bacteria -fq1 {fq1} -fq2 {fq2} -out {self.out} -ref {ref} -upload {self.upload} {ifpollute}")
            elif self.kind=='fungi':
                self.tree_mode='core'
                cmd.append(f"time {python} {PATH}/main.py -id {id} -kind fungi -fq1 {fq1} -fq2 {fq2} -out {self.out} -ref {ref} -upload {self.upload} {ifpollute}")
            elif self.kind=='virus':    
                self.tree_mode='wgs'
                cmd.append(f"time {python} {PATH}/main.py -id {id} -kind virus -fq1 {fq1} -fq2 {fq2} -out {self.out} -ref {ref} -upload {self.upload}  {ifdenovo}")

            common.cmd2shell(cmd,p_sh)
            logging.info(f"测序模式={mode} ,fq1={fq1} ,fq2={fq2} ,病原体类型={self.kind} ,去污染={ifpollute} ,参考序列={ref}，进化树模式={self.tree_mode}")
        return path_sh_list


    def prepare_fa(self):
        p_sh = f"{self.out}/0.shell/prepare_fa.sh"
        self.p_fa = f"{self.out}/fa"
        cmd = []

        cmd.append(f"mkdir -p {self.p_fa}")
        for id in self.id_list:
            if self.kind != 'virus':
                cmd.append(f"cp -r {self.out}/{id}/{id}/2.ass/{id}_scaffolds.fasta {self.p_fa}/{id}.fa")
            else:
                cmd.append(f"cp -r {self.out}/{id}/{id}/2.ass/{id}_consensus.fa {self.p_fa}/{id}.fa")
            common.cmd2shell(cmd,p_sh)
        return p_sh


    # 进化树
    def tree(self,d_tree):
        p_sh = f"{self.out}/0.shell/tree.sh"
        p_tree = f"{self.out}/Tree"
        p_yaml = f"{p_tree}/tree.yaml"
        os.makedirs(p_tree,exist_ok=True)

        fa_dic={}
        for file in os.listdir(self.p_fa):
            id=file.replace('.fa','')
            fa_dic[id] = os.path.join(self.p_fa,file)

        for dir in os.listdir(d_tree):
            fa_dic[dir] = os.path.join(d_tree,dir)  

        tree_dic = {}
        tree_dic['samples'] = {}
        tree_dic['samples'] = fa_dic
        tree_dic['library'] = 'Tree'
        tree_dic['result_dir'] = p_tree
        tree_dic['kindom'] = self.kind

        #生成yaml
        with open(p_yaml, "w", encoding="utf-8", newline="") as f:
            f.write(yaml.dump(tree_dic, allow_unicode=True))
        
        cmd = []
        cmd.append(f"{python} {TREE}/main.py -p {self.tree_mode} -i {p_yaml}")   #run进化树脚本
        for id in self.id_list: #拷贝结果到各样本的结果目录
            cmd.append(f"mkdir -p {self.out}/{id}/{id}/Tree_Typing/Tree")
            cmd.append(f"cp -r {p_tree}/Tree/Upload/*  {self.out}/{id}/{id}/Tree_Typing/Tree")
        common.cmd2shell(cmd,p_sh)
        return p_sh



    # 分型
    def typing(self,typing):
        cmd = []
        p_sh=f"{self.out}/0.shell/typing.sh"
        for id in self.id_list:
            p_typing = f"{self.out}/{id}/typing"
            p_res = f"{self.out}/{id}/{id}/Tree_Typing/Typing"

            cmd.append(f"""# typing
mkdir -p {p_typing}
mkdir -p {p_res}
{python} {TYPING}/main.py --samples {id} -d {typing} --fastas {self.p_fa}/{id}.fa -o {p_typing}
cp -r {p_typing}/Upload/{id} {p_res}
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh



    # 报告
    def report(self):
        p_sh=f"{self.out}/0.shell/report.sh"
        cmd = []
        for id in self.id_list:
            cmd.append(f"""#report
cp -r {p_GDHR}/src {self.out}/{id}/{id}/Tree_Typing
{perl} {PATH}/lib/report/tree_typing_report.pl {id} {self.out}/{id}/{id}/Tree_Typing

cd {self.out}/{id}
zip -qr {id}.zip {id}
cd {self.out}/{id}/{id}
zip -qr Tree_Typing.zip Tree_Typing

cp -r {self.out}/{id}/{id} {self.upload}
cp -r {self.out}/{id}/{id}.zip {self.upload}
""")
        common.cmd2shell(cmd,p_sh)
        return p_sh




#### 参数
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-table', required=True,type=click.Path(),help="task_id yaml file")
@click.option('-kind', required=True,type=click.Choice(["bacteria","fungi","virus"]),help="pathogen kind")
@click.option('-out', required=True,type=click.Path(),help="out dir")
@click.option('-upload', required=False,type=click.Path(),default='None',help="upload dir")
@click.option('-ref', required=False,type=click.Path(),default='None',help="ref dir")
@click.option('-tree', required=False,type=click.Path(),help="input tree dir")
@click.option('-typing', required=False,type=click.STRING,help="trace species")


def main(table,kind,out,upload,tree,typing,ref):
    logfile=f"{out}/log"
    os.makedirs(logfile,exist_ok=True)

    project = MultiRun(out,upload,kind)
    path_shell_list = project.single_denovo(table,ref)
    common.mul_pool(path_shell_list, logfile,2)

    sh2 = project.prepare_fa()
    common.prun((sh2, logfile))

    if tree != "blank":
        sh3 = project.tree(tree)
        common.prun((sh3, logfile))

    if typing != "blank":
        sh4 = project.typing(typing)
        common.prun((sh4, logfile))

    if tree != "blank" or  typing != "blank":
        sh5 = project.report()
        common.prun((sh5, logfile))


if __name__ == "__main__":
    main()
