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

    1. 给定长度范围，把目标序列拆分为DNA短片段
    2. 确定酶切位点：延长 - 判断 - 选择
DNA短片段相互重叠，酶切位点确定 BtsI/BsrDI
    3. 
    4. 根据division point，5 oligos, (1,4)(2,3)(3,2)(4,1)，确定正负链
    5. 添加酶切位点 
    6. 添加互补序列 (rc_len=14)
    7. 粘端唯一 unique overhang of each oligo in a group (2)

    其他:
    - deltaG, secondary structure 7-12 1.5kj


