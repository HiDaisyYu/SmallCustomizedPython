filename = 'sample.vcf'

## A class representing simple SNPs
class SNP:
    def __init__(self, chrname, pos, snpid, refallele, altallele):
        assert refallele != altallele, "Error: ref == alt at pos " + str(pos)
        self.chrname = chrname
        self.pos = pos
        self.snpid = snpid
        self.refallele = refallele
        self.altallele = altallele
 
    # Returns True if refallele/altallele is A/G, G/A, C/T, or T/C
    def is_transition(self):
        if self.refallele == "G" or self.refallele == "A":
            if self.altallele == "G" or self.altallele == "A":
                return True

        if self.refallele == "C" or self.refallele == "T":
            if self.altallele == "C" or self.altallele == "T":
                return True
        return False

## A class representing a chromosome, which has a collection of SNPs
class Chromosome:
    def __init__(self, chrname):
        self.chrname = chrname
        self.locations_to_snps = dict()

    # Returns the chromosome name
    def get_name(self):
        return self.name

    # Given all necessary information to add a new SNP, create
    # a new SNP object and add it to the SNPs dictionary. If
    # a SNP already exists at that location, or
    # the given chrname doesn't match self.chrname, an error is reported.
    def add_snp(self, chrname, pos, snpid, refallele, altallele):
        # If there is already an entry for that SNP, throw an error
        open_location = (pos not in self.locations_to_snps)
        assert open_location, "Duplicate SNP: " + self.chrname + ":" + str(pos)

        # If the chrname doesn't match self.chrname, throw and error
        assert chrname == self.chrname, "Chr name mismatch!"

        # Otherwise, create the SNP object and add it to the dictionary
        newsnp = SNP(chrname, pos, snpid, refallele, altallele)
        self.locations_to_snps[pos] = newsnp

    # Returns the number of snps between l and m
    def snps_region(self, l, m):
        count = 0
        for location in self.locations_to_snps.keys():
            if location >= l and location <= m:
                count = count + 1

        return count

    # Returns the number of transition snps between l and m
    def transitions_region(self, l, m):
        trans_count = 0
        for location in self.locations_to_snps.keys():
            snp = self.locations_to_snps[location]
            if location >= l and location <= m:
                if snp.is_transition():
                    trans_count = trans_count + 1

        return trans_count

    # Returns the position of the last SNP known
    def get_last_snp_position(self):
        locations = list(self.locations_to_snps.keys())
        locations.sort()
        return locations[len(locations) - 1]

    # Given a region size, looks at non-overlapping windows
    # of that size and returns a list of four elements:
    # [start of region, end of region, transition number, snp number]
    def region_answer(self, region_size):
        region_start = 1
        last_snp_position = self.get_last_snp_position()
        # construct a dictionary for adding loop function output
        answer_dict = dict()

        # read each 100,000bp region using while loop
        while region_start < last_snp_position:
            # if region end position is smaller than the last SNP position,
            # continue to read next 100,000bp region,
            # until the region end extends the last SNP postion(see elif code below)
            if region_start + region_size - 1 < last_snp_position:
                region_end = region_start + region_size - 1
                number_of_snp = self.snps_region(region_start, region_end)
                number_of_transition = self.transitions_region(region_start, region_end)
                answer_dict[region_start] = (region_start, region_end, number_of_transition, number_of_snp)
                region_start = region_start + region_size
                
            # if region end positon extends the last SNP position,
            # define the last SNP position as the region end,
            # meanwhile output the result into the dictionary named answer_dict
            elif region_start + region_size - 1 >= last_snp_position:
                region_end = last_snp_position
                number_of_snp = self.snps_region(region_start, region_end)
                number_of_transition = self.transitions_region(region_start, region_end)
                answer_dict[region_start] = (region_start, region_end, number_of_transition, number_of_snp)
                return answer_dict

## Create chrnames_to_chrs dictionary, parse the input file
chrnames_to_chrs = dict()
fhandle = open(filename, "r")

for line in fhandle:
    # Don't attempt to parse header lines
    if line[0] != '#':
        line_list = line.strip().split()

        chrname = line_list[0]
        pos = int(line_list[1])
        snpid = line_list[2]
        refallele = line_list[3]
        altallele = line_list[4]

        # Put the data in the dictionary named chr_obj
        # If this chrname is already in chr_obj, then add snp information
        if chrname in chrnames_to_chrs:
            chr_obj = chrnames_to_chrs[chrname]
            chr_obj.add_snp(chrname, pos, snpid, refallele, altallele)
        # If this chrname is not in chr_obj, then get the charname from Chromosome,
        # before adding SNP information
        else:
            chr_obj = Chromosome(chrname)
            chr_obj.add_snp(chrname, pos, snpid, refallele, altallele)
            chrnames_to_chrs[chrname] = chr_obj

fhandle.close()

## Print the results
# First print the headline
print("chromosome" + "\t" + "region" + "\t" + "percent_transitions" + "\t" + "num_snps")

# For chrnames_to_chrs dictionary which storing the file input information,
for chrname in chrnames_to_chrs:
    chr_obj = chrnames_to_chrs[chrname]
    # perform the self.region_answer to read each 100,000bp region for each chrname
    answer_dict = chr_obj.region_answer(100000)
    # Only retrive the value from answer_dict,
    # and assign these values into answer_list
    for value in answer_dict:
        answer_list = answer_dict[value]
        
        region_start = str(answer_list[0])
        region_end = str(answer_list[1])
        number_of_transition = answer_list[2]
        number_of_snp = answer_list[3]

        if number_of_snp == 0:
            percent_transition = 0  
        elif number_of_snp != 0:
            percent_transition = number_of_transition/number_of_snp
        # Standard the format of output regarding the example
        percent = str('%.4f' % percent_transition)
            
        print(chrname + "\t" + region_start + ".." + region_end + "\t" + percent + "\t" + str(number_of_snp))


    
