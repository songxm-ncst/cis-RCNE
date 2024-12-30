import sys

# 定义两个文件名
ce_tab_file = sys.argv[1]  # 'CE.tab' phastCons
filter_phylop_file = sys.argv[2]  # 'Filter_phyloP_element-scores.txt' phyloP
bed_file = sys.argv[3] # Ath_mulitway_most-cons.bed
# 初始化两个集合来存储两个文件中的元素
ce_elements = set()
filter_phylop_elements = set()

# 读取CE.tab文件并提取第一列元素
with open(ce_tab_file, 'r') as file:
    for line in file:
        columns = line.strip().split('\t')
        ce_elements.add(columns[0])

# 读取Filter_phyloP_element-scores.txt文件并提取第四列元素
with open(filter_phylop_file, 'r') as file:
    for line in file:
        columns = line.strip().split('\t')
        filter_phylop_elements.add(columns[3])

# 找出两个集合的交集
intersection = ce_elements.intersection(filter_phylop_elements)

# 打印交集结果
#print("交集结果：")
#for element in intersection:
#    print(element)

with open(bed_file, 'r') as file:
    for line in file:
        columns = line.strip().split('\t')
        if columns[3] in intersection:
            print(line.strip())  # 打印符合条件的行

