# partition a string
recursion
index of the next character to be processed
the output string
[partition a string](https://www.geeksforgeeks.org/print-ways-break-string-bracket-form/)

[x << y](https://wiki.python.org/moin/BitwiseOperators)

按位运算符

xrange() generator
# On the Design of Oligos for Gene Synthesis

子片段

[GeneArt](https://www.thermofisher.com/cn/zh/home/life-science/cloning/geneart-type-ii-assembly-kits.html?ef_id=Cj0KCQiA2uH-BRCCARIsAEeef3nQo7Ts8HyCmHL_4dk9iRTHvU5vsekEDzNTqLYiyLJngcpewAEUzKcaAvixEALw_wcB:G:s&s_kwcid=AL!3652!3!476186507209!p!!g!!bsai&cid=bid_mol_clo_r01_co_cp1358_pjt0000_bid00000_0se_gaw_nt_pur_con&gclid=Cj0KCQiA2uH-BRCCARIsAEeef3nQo7Ts8HyCmHL_4dk9iRTHvU5vsekEDzNTqLYiyLJngcpewAEUzKcaAvixEALw_wcB)

    Design:
    - DNA短片段 + 限制酶识别序列 + 反向互补序列（使得2-10个碱基在3’端悬挂）

    拆分序列:
    1. 给定长度范围，把目标序列拆分为相互重叠的DNA短片段
    Input: (sequence, min, max)
    Output: (oligoGroups)
    - main body overlap with each other (len_overlap=2)
    determine restriction enzyme:
    - restriction site/avoid site not in the main bodies (res="BtsI"/res="BsrDI")
    构造寡核苷酸序列:
    - division point, 5 oligos, (1,4)(2,3)(3,2)(4,1), positive/negative
    - BtsI/BsrDI site
    - reverse complement sequence length (rc_len=14)
    - unique overhang of each oligo in a group (2)
    其他:
    - deltaG, secondary structure 7-12 1.5kj


x = [[0, 1], [7, 8], [2, 3], [4, 6]]
y = {}
y[0] = [0,1]
y[7] = [7,8]
y[2] = [2,3]
y[4] = [4,6]
y
{0: [0, 1], 2: [2, 3], 4: [4, 6], 7: [7, 8]}
y.keys()
[0, 2, 4, 7]
sorted(y.keys())