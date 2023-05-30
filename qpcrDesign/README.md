# qpcrDesign
基于`primer3`的qPCR设计脚本, 替代Windows版`PrimerExpress3.0.1`实现自动化引物&探针设计.

## 设计思路
- qpcr  
  1. 输入FASTA都格式化成大写字母  
  2. 内引物设计, 使用`primer3`
  3. 内引物解析, 自编脚本

- qpcr-batch  
  1. 从BED获取FASTA, 并拆分.
  2. 使用`qpcr`批量跑, 初次设计使用70-110扩增产物长度阈值.
  3. 使用`qpcr`批量跑, 针对70-110阈值不成功的区域再次设计,使用50-150扩增产物长度阈值.

## TODO

## 更新
- [x] 230519 批量脚本, Name列和BED关联
- [x] 230529 `primer3_core`部分更新, 使用Primer3Plus固定参数提高设计敏感性, 仅`PRODUCT_SIZE_RANGE`可变
