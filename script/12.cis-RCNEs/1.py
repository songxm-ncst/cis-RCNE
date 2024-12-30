file_in = './alternate_cisRCNEs.maf'
file_out1 = open(r'./cisRCNEs.maf', 'w')  #每一部分大于17
file_out2 = open(r'./out2.maf', 'w') #每一不发小于17
with open(file_in,'r') as file_in:
    a = 0
    l = []
    for line in file_in:
        if line.strip() == '':
            if a < 30:
                file_out2.write(''.join(l) + '\n')
            else:
                file_out1.write(''.join(l)+ '\n')
            l = []
            a = 0
        else:
            if line[0] != 'e':
                l.append(line)
            a += 1



