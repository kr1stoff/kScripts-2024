#!/usr/bin/env python
import os
import logging
from multiprocessing import Pool
from subprocess import Popen,PIPE,run 
import itertools
import sys

#写shell
def cmd2shell(cmds, sh):
    """
    命令list写到shell脚本里
        cmds:  命令列表
        shell: 输出脚本
    """
    with open(sh, "w", encoding="utf-8", newline="") as w:
        w.write(f"#!/usr/bin/bash \n")
        for cmd in cmds:
            w.write(cmd + "\n")
    os.system(f"chmod 755 {sh}")




def set_log(filename):

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(filename)# 创建一个handler，用于写入日志，输出到文件
    ch = logging.StreamHandler(sys.stdout)    # 创建一个handler，用于输出到屏幕
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M')
 
    # 绑定formatter 到handler上
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
 
    # 绑定handler到logger对象上
    logger.addHandler(fh)
    logger.addHandler(ch)


#返回时间戳
def timestamp(name, statu):
    if statu == "start":
        cmd = f"echo -e [`date +%y-%m-%d.%H:%M:%S`] {name} Start"
    elif statu == "end":
        cmd = f"echo -e [`date +%y-%m-%d.%H:%M:%S`] {name} End\n"
    return cmd


#创建软连接
def creat_link(p_file, p_link):
    if os.path.isfile(p_file) and not os.path.exists(p_link):
        os.symlink(p_file, p_link)
    else:
        logging.debug(f"不存在文件<{p_file}> 或 链接已存在<{p_link}>")



#多进程池
def prun(sh_tuple):
    shfile = sh_tuple[0]
    logdir = sh_tuple[1]

    shname=os.path.basename(shfile)
    p=Popen(shfile,shell=True,encoding='utf-8',stdout=PIPE,stderr=PIPE)
    stdout,stderr = p.communicate()

    open(f"{logdir}/{shname}.o",'w',encoding='utf-8').write(stdout)
    open(f"{logdir}/{shname}.e",'w',encoding='utf-8').write(stderr)



def mul_pool(plist, logdir, num=None):
    """plist: shell文件名"""

    pool= Pool(num) if num else Pool(len(plist))
    argv_list=[(tmp,logdir) for tmp in plist]
    pool.map(prun, argv_list)




# 获取系统的cpu / 线程数
def get_cfg():
    # CPU 核数
    cmd = "cat /proc/cpuinfo| grep 'cpu cores'| uniq |cut -d' ' -f3"
    cpu = int(os.popen(cmd).read())

    # 线程数
    cmd = "cat /proc/cpuinfo| grep 'processor'| wc -l"
    threads = int(os.popen(cmd).read())
    return cpu,threads