
from encodings import utf_8_sig
import sys
import json
from jsonpath import jsonpath 
import os
import pandas as pd

rgi_json=os.path.abspath(sys.argv[1])
outfile=os.path.abspath(sys.argv[2])
merge_file=os.path.abspath(sys.argv[3])

OUT=open(outfile,'w',encoding='utf_8_sig')
OUT.write("比对基因ID\t耐药基因ID\t耐药基因名\t耐药基因种类\t测序深度\t覆盖度%\n")
dir=os.path.dirname(outfile)

if len(sys.argv) != 4:
    print('Usage: python %s <rgi_json> <outfile> <merge_file> '%sys.argv[0])
    sys.exit(-1)


#文件大小
if os.path.getsize(rgi_json) <= int(10) or not os.path.exists(rgi_json):
    
    print("f{rgi_json} is empty")
    exit(-1)


result_dic={}
with open(rgi_json, "r", encoding='utf-8') as f,open(merge_file,'r')as m:
    #json格式字符串转dic
    rgi=json.load(f)

    dic={}
    for line in m:
        line=line.strip().split('\t')
        dic[line[0]]=line[1:]


    for gene_id in rgi.keys():
        #简化id
        id=gene_id.split(' ',1)[0]
        print(id)

        depth=format(float(dic[id][0]),'.1f')
        cov=format(float(dic[id][1]),'.1%')
        #print(cov,type(cov))
    
        #多重字典内容读取
        data=jsonpath(rgi[gene_id], '$..ARO_category')[0]
        
        #提取内容
        for key in data.keys():
            aro_num  = data[key]["category_aro_accession"]
            aro_name = data[key]["category_aro_name"]
            aro_class= data[key]["category_aro_class_name"]
            #aro_des  = data[key]["category_aro_description"]

            #写入内容
            OUT.write(f"{id}\t{aro_num}\t{aro_name}\t{aro_class}\t{depth}\t{cov}\n")


OUT.close()

#生成xlsx

data=pd.read_csv(outfile,sep="\t")
pd.DataFrame(data).to_excel(f"{dir}/detail_drug_resistance.xlsx",index=False,encoding='utf_8_sig')

