with open('Data/string_identifiers.txt', 'w') as outfile:
    for line in open('Data/9606.protein.links.v11.0.txt', 'r'):
        row = line.strip().split(' ')
        row[0] = row[0].replace('9606.', '')
        row[1] = row[1].replace('9606.', '')
        outfile.write(row[0])
        outfile.write('\n')
        outfile.write(row[1])
        outfile.write('\n')


plt.hist(sizes,bins=50)
plt.yscale('log')
plt.xlabel('Size')
plt.ylabel('Freq')
plt.show()
import pysam
