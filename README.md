# AsSb-gene项目

## 项目简介

AsSb-cso项目旨在识别和分析As和Sb相关的同源物。通过使用特定的数据库和工具，我们能够有效地进行序列比对和同源物识别。

根据现有文献，aioA是砷专性的，因此本项目目前没有包括。

## 功能描述

- 使用hmmsearch工具进行序列比对，识别AsSb-cso相关的同源物。
- 提取并处理蛋白质序列以生成同源物列表。
- 生成输出文件，包含识别到的AsSb-cso同源物信息。

## 使用说明

1. 确保安装了必要的工具，如hmmsearch和seqtk。
2. 运行脚本以开始同源物识别过程。
3. 检查输出文件以获取识别结果。

## 目录结构

- `dbbase/`: 构架马尔科夫链的原始蛋白质序列。
- `hmmdb/`: 用于比对的砷锑代谢数据库，根据pfam和uniprot等构建。
- `outputs/`: 存储生成的同源物列表和其他结果文件。

## Ref

项目参考了以下资源和文献：

汞甲基化：http://rcees.cas.cn/jz/202212/t20221227_6825428.html

MeHg：https://github.com/elizabethmcd/MEHG

MeHg：https://github.com/ericcapo/marky-coco

数据库构建：http://hmmer.org

蛋白质序列获取：https://www.uniprot.org/uniprotkb?query=arsA

多序列比对：https://mafft.cbrc.jp/alignment/software/macosx.html
