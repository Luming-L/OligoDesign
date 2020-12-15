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

```python
def findCombinations(string, index, out):
    if index == len(string):
        print(out)
 
    for i in range(index, len(string), 1):
 
        # append substring formed by str[index, i] to output string
        findCombinations(string, i + 1, out + "(" + string[index:i + 1] + ")")
 
# Driver Code
if __name__ == "__main__":
 
    # input string
    string = "abcd"
    findCombinations(string, 0, "")
 
# This code is contributed by
# sanjeev2552
```