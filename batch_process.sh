#!/bin/bash
# 批量处理

# 检查AS.txt文件是否存在
if [ ! -f "AS.txt" ]; then
    echo "AS.txt文件不存在"
    exit 1
fi

# 创建或清空done文件
> done

# 逐行读取AS.txt文件
while read -r sample; do
    # 去除可能的空白字符
    sample=$(echo $sample | xargs)
    
    # 检查样品名是否为空
    if [ -z "$sample" ]; then
        continue
    fi
    
    # 执行第一个命令
    ./cp.sh $sample
    
    # 执行第二个命令
    setsid ./MeHg_pe_all.sh $sample
    
    # 记录完成
    echo "$sample" >> done
    
done < "AS.txt"

echo "所有样品处理完成"
