
import sys
import json
from jsonpath import jsonpath 
import os

rgi_json=sys.argv[1]
outfile=sys.argv[2]
merge_file=sys.argv[3]

OUT=open(outfile,'w')


if len(sys.argv) != 4:
    print('Usage: python %s <rgi_json> <outfile> <merge_file> '%sys.argv[0])
    sys.exit(-1)



if not os.path.exists(rgi_json):
    os.system(r"touch {}".format(outfile)) 
    exit(-1)


result_dic={}
with open(rgi_json, "r", encoding='utf-8') as f,open(merge_file,'r')as m:
    #json格式字符串转dic
    rgi=json.load(f)
    OUT.write("Gene_id\tAro_accession\tAro_name\tClass_name\tDepth\tCoverage(%)\n")

    dep={}
    coverage={}
    for line in m:
        line=line.strip().split('\t')
        dep[line[0]]=line[1]
        coverage[line[0]]=line[2]


    for gene_id in rgi.keys():
        #简化id
        id=gene_id.split(' ',1)[0]
        #print(id)
        #print (dep[id])

        depth=format(float(dep[id]),'.1f')
        #depth=dep[id]
        cov=format(float(coverage[id]),'.1%')
        #cov=coverage[id]
        #print(cov,type(cov))
    
        #多重字典内容读取
        data=jsonpath(rgi[gene_id], '$..ARO_category')[0]
        
        #提取内容
        for key in data.keys():
            aro_num  = data[key]["category_aro_accession"]
            aro_name = data[key]["category_aro_name"]
            aro_class= data[key]["category_aro_class_name"]
            aro_des  = data[key]["category_aro_description"]

            #写入内容
            OUT.write(f"{id}\t{aro_num}\t{aro_name}\t{aro_class}\t{depth}\t{cov}\n")

OUT.close()

