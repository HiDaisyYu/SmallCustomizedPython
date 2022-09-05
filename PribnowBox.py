## Pribnow box in E.Coli genome

## Define a function to read whome genome sequence of E.Coli
def readFASTA(f1):
    # Construct an empty dictionary to store genome sequence
    dict_genome = {}
    
    # Read file line by line
    for line in f1:
        # Save the bacteria strain info(starts with '>') as the dictionary key
        if line.startswith('>'):
            name = line.rstrip('\n')
            dict_genome[name] = ''
        # Save the sequence(not starts with '>') as value to its strain info key
        else:
            dict_genome[name] = dict_genome[name] + line.rstrip('\n')

    f1.close()
    return dict_genome


## Define a function to read position information
def readPOS(f2):
    # Construct an empty dictionary to store position info
    dict_position = {}
    # Construct an empty list to store values of dictionary
    list_position = []
    
    # Read file line by line
    for line in f2:
        a = line.split()
        # Extract gene ID column
        gene_id = a[0]
        # Extract position column
        gene_pos = a[3]
        
        # To check if each line is complement
        if gene_pos.startswith('c'):
            # Split each line 
            gene_pos_1 = gene_pos.split('(',1)
            # I notice there are joined position genes during analyzing data
            # eg. gene ID b0502
            # so I intend to split these genes info specifically
            if gene_pos_1[1].startswith('join'):
                gene_pos_11 = gene_pos_1[1].split('(', 1)
                gene_pos_12 = gene_pos_11[1].split('..', 1)
                gene_pos_13 = gene_pos.rsplit(')', 2)
                gene_pos_14 = gene_pos_13[0].rsplit('..',1)
                strand = "Reverse"
                reverse_1 = int(gene_pos_14[1])
                reverse_2 = int(gene_pos_12[0])
                # Add info into list
                list_position = [strand, reverse_1, reverse_2]
                # Add each gene's info as key and value to dictionary
                dict_position[gene_id] = list_position

            # non-join position for complement gene
            else:
                gene_pos_2 = gene_pos_1[1].split('..',1)
                gene_pos_3 = gene_pos_2[1].rsplit(')')
                strand = "Reverse"
                # Convert the position info into numeric data
                reverse_1 = int(gene_pos_3[0])
                reverse_2 = int(gene_pos_2[0])
                # Add info into list
                list_position = [strand, reverse_1, reverse_2]
                # Add each gene's info as key and value to dictionary
                dict_position[gene_id] = list_position
                
        # forward strand        
        else:
            # I notice there are joined positions for forward strands during analyzing data
            # eg. gene ID b4659
            # So I intend to split these gene info specifically
            if gene_pos.startswith('join'):
                # Split joined position info
                gene_pos_4 = gene_pos.split('(', 1)
                gene_pos_5 = gene_pos_4[1].split('..', 1)
                gene_pos_6 = gene_pos_4[1].rsplit(')', 1)
                gene_pos_7 = gene_pos_6[0].rsplit('..', 1)
                strand = "Forward"
                # Convert the position info into numeric data
                forward_1 = int(gene_pos_5[0])
                forward_2 = int(gene_pos_7[1])
                list_position = [strand, forward_1, forward_2]
                dict_position[gene_id] = list_position

            # non-join position for forward strand  
            else:
                gene_pos_8 = gene_pos.split('..', 1)
                strand = "Forward"
                # Convert the position info into numeric data
                forward_1 = int(gene_pos_8[0])
                forward_2 = int(gene_pos_8[1])
                # Add info into list
                list_position = [strand, forward_1, forward_2]
                # Add each gene's info as key and value to dictionary
                dict_position[gene_id] = list_position

    f2.close()
    return dict_position


## Define the extract upstream function
def extract_upstream(str_genome, position_region, nucleotide_num):

    # check if this gene position region is reverse
    if position_region[0] > position_region[1]:
        start_pos = int(position_region[0])
        end_pos = int(position_region[0] + nucleotide_num)
    # for forward strand
    else:
        end_pos = int(position_region[0] -1)
        start_pos = int(end_pos - nucleotide_num)
    return str_genome[start_pos : end_pos]

    
## Define the main() function 
def main():
    # Read genome file
    f1 = open('e.coli.genome', 'r')
    dict_genome = readFASTA(f1)

    # Read gene position file
    f2 = open('e.coli.pos', 'r')
    dict_position = readPOS(f2)

    # Get ths whome genome sequence and make it to a string
    str_genome = str(list(dict_genome.values())[0])

    # Store all the gene ID into a list
    gene_ID_list = list(dict_position.keys())

    # Store all the 'reverse' and 'forward' into a list
    position_region = list(dict_position.values())
    count_gene = len(dict_position.keys())
    strand_list = []
    for i in range(0, count_gene):
        # delete the integers to get all the 'reverse' and 'forward'
        strand_list.append(position_region[i][0])

    # Obtain position region and make it to a list
    # which includes lists made of two integers
    for i in range(0, count_gene):
        # delete the 'reverse' and 'forward' element
        del position_region[i][0]

    # Extract all the upstream sequence with extract_upstream() function
    # and store all these upstream sequence into a list
    upstream_list = []
    nucleotide_num = 20
    for i in range(0, count_gene):
        up = extract_upstream(str_genome, position_region[i], nucleotide_num)
        upstream_list.append(up)

    # Print the whome genome information as title 
    # Change the dictionary key to string type
    title_pre = str(list(dict_genome.keys())[0])
    # Delete the '>'
    tittle = title_pre.replace('>', '')
    print(tittle)
    print('GeneID' + '\t' + 'Start' + '\t' + 'End' + '\t' +
          'Strand' + '\t' + '20bp_promoter_seq_5to3')
    

        
    # Check if there is TATAAT in each upstream sequence
    for i in range(0, count_gene):
        # Divide the printing into 'reverse' and 'forward'
        if strand_list[i] == 'Reverse':
            if 'ATTATA' in upstream_list[i]:
                # print result and highlight TATAAT with <>
                print(str(gene_ID_list[i]) + '\t' + str(position_region[i][0]) + '\t' +
                      str(position_region[i][1]) + '\t' + str(strand_list[i]) + '\t' +
                      '*' + str(upstream_list[i]).replace('ATTATA', str('<ATTATA>')))
        if strand_list[i] == 'Forward':
            if 'TATAAT' in upstream_list[i]:
                print(str(gene_ID_list[i]) + '\t' + str(position_region[i][0]) + '\t' +
                      str(position_region[i][1]) + '\t' + str(strand_list[i]) + '\t' +
                      str(upstream_list[i]).replace('TATAAT', str('<TATAAT>')) + '*')

main()    




# I have only obtained 12 sequences which include the TATAAT within 20 bp promoter regions.
# I think there are two steps which influence the result.

## (1) The genome assembly result.
# As these 20bp promoter regions were defined based on the gene positions provided by KEGG POS file,
# I was thinking if these gene positions/annotations are really correct.

# Solution for problem(1).
# In order to check my idea, I went to UCSC genome browser and NCBI geneBank to check different assembly results.
# I found the gene assembly results for this MG1655 E.Coli are roughly the same, except for minor differences.
# 1st difference: The total gene number, there are 4494 genes in our POS file and I notice there are different
# results include 4403 and 4466 genes respectively.
# 2nd difference: The gene position. Although they are almost the same from different assembly data.
# I notice there are genes called 'repeat region' in other assembly data, which do not exist in POS file.

# Taken together, the gene number and position differences are really small.
# I think these minor differences may change the result a bit,
# but can not improve a lot for detecting pribnow box.

## (2) The searching method
# The Pribnow box is usually represented by the consensus sequence TATAAT,
# but not strictly TATAAT/ATTATA in each promoter region.
# There is usually one or more nucleotide difference in each gene's pribnow box.
# Searching exactly the 'TATAAT'/'ATTATA' will definitely miss quite a lot pribnow boxes.

# Solution for problem(2).
# For 20bp promoter regions extracted for each gene, I plan to obtain the position-specific scoring matrix (PSSM)
# of every 6 nucleotides in each promoter region sequence.
# And calculated the Relative Scores of TATAAT(forward strand), ATTATA(reverse strand) in each 6-nucleotide unit.
# For those 6-nucleotide units, whose Relative Scores >= 80%, can be defined as a Pribnow Box.


# I give an example to make it clear:
# e.g.
# 20bp_promoter_seq_5to3(forward): ACTGCATGTACTCTAGCTAC*
# Divide it into 15 6-nucleotide units:
# ACTGCA
# CTGCAT
# TGCATG
# GCATGT
# CATGTA
# ATGTAC
# TGTACT
# GTACTC
# TACTCT
# ACTCTA
# CTCTAG
# TCTAGC
# CTAGCT
# TAGCTA
# AGCTAC

# For the first one, ACTGCA,
# obtain the Position Frequency Matrix (PFM) first,
# A 100001
# C 010010
# G 000100
# T 001000

# Convert it to position-specific scoring matrix (PSSM),
# A:	 0.916	-0.693	-0.693	-0.693	-0.693	0.916	
# C:	-0.693	0.916	-0.693	-0.693	0.916	-0.693	
# G:	-0.693	-0.693	-0.693	0.916	-0.693	-0.693	
# T:	-0.693	-0.693	0.916	-0.693	-0.693	-0.693

# then get the Relative Scores of TATAAT(as this is forward strand),
# Relative Scores = 34%,
# 34% < 80%, so this 6-nucleotide unit ACTGCA, is NOT Pribnow box.

# Repeat calculation for these 15 6-nucleotide units;
# If there are more than one 6-nucleotide unit with Relative Scores>=80%,
# I plan to choose the highest score unit as the Pribnow bow in this gene upstream region.

# Repeat above steps for 4494 20bp_promoter_seq_5to3 sequences,
# to get as many Pribnow box.
    
    
