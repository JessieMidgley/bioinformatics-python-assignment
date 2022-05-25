class Read_FastA:
    def Read_FastA_Names_And_Sequences(self, filepath):
        print("Reading fastA sequences...")
        self.sequence_names = []
        self.sequences = []
        fasta_file = open(filepath, 'r')
        self.counter = 0
        for line in fasta_file:
            if (line[0] == '>'):
                self.counter += 1
                self.sequence_names.append(line[1:].strip())
                self.sequences.append(str())
            else:
                self.sequences[self.counter - 1] += line.strip()
        fasta_file.close()
        return (self.sequence_names, self.sequences)


class Read_GFF:
    def Get_Gene_Positions(self, list_of_gff_features, filepath, feature):
        gff_file = open(filepath, 'r')
        list_of_gff_features = []
        for gff_line in gff_file:
            if gff_line[0] != '#':
                list_of_column_items = gff_line.split(sep='\t')
                if list_of_column_items[2] == 'gene':
                    lst = (list_of_column_items[0], list_of_column_items[3],
                           list_of_column_items[4], list_of_column_items[6])
                    list_of_gff_features.append(lst)
        return list_of_gff_features


class My_Codons:
    def Make_List_Of_Codons(self):
        self.nucleotides = ['G', 'A', 'T', 'C']
        self.codons = []
        self.tempcodon = ''
        for a in self.nucleotides:
            for b in self.nucleotides:
              for c in self.nucleotides:
                self.tempcodon = a + b + c
                self.codons.append(self.tempcodon)
        return (self.codons)

    def Count_Codons(self, new_sliced_sequence, codons, number_of_occurrences):
        codon_dict = {}
        for x in range(0,64):
            temp_dict = {codons[x]:number_of_occurrences[x]}

            codon_dict.update(temp_dict)
        for gene in new_sliced_sequence:
            for start in range(0, len(gene) + 1, 3):
                test_codon = gene[start:start+3]

                if len(test_codon) == 3:
                    codon_dict[test_codon] = codon_dict.get(test_codon, 0) + 1

        values = codon_dict.values()
        number_of_occurrences = list(values)

        return (number_of_occurrences)


# MAIN CODE
# ==================================================================
path_of_gff_file = 'saccharomyces_cerevisiae_2022.gff'  # change the path string if yours is different
path_of_fasta_file = 'saccharomyces_cerevisiae_2022.fna'  # change the path string if yours is different

# Get the positions and offset of all codings sequences (CDS) in the yeast genome
GFF_file_object = Read_GFF()
list_of_gff_features = []
total_sequence_length = 0
list_of_gff_features = GFF_file_object.Get_Gene_Positions(list_of_gff_features, path_of_gff_file, 'CDS')

# make a list of all 64 possible codons
codon_object = My_Codons()
codons = codon_object.Make_List_Of_Codons()

# Read the chromosome sequences
FASTA_file_object = Read_FastA()
sequence_name, sequences = FASTA_file_object.Read_FastA_Names_And_Sequences(path_of_fasta_file)

# Loop over list_of_gff_features, using one entry at a time
number_of_occurrences = [0] * 64
print('Counting codons...')
new_sliced_sequence = []
for gff_line in list_of_gff_features:
# get chromosome and slice the gene sequence of the chromosome with the calculated index
    chromosome_sequence = sequences[sequence_name.index(gff_line[0])]
    strand = chromosome_sequence[int(gff_line[1]) - 1:int(gff_line[2])]
# reverse transcribe negative strand
    base_pairs = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    positive_strand = ''
    if gff_line[3] == '-':
        reverse_strand = strand[::-1]
        for base in reverse_strand:
            positive_strand += base_pairs.get(base)
        new_sliced_sequence.append(positive_strand)
    else:
        positive_strand = strand
        new_sliced_sequence.append(positive_strand)

number_of_occurrences = codon_object.Count_Codons(new_sliced_sequence, codons, number_of_occurrences)

# Print out the total codons
total_codons = sum(number_of_occurrences)
print("Total codons:", total_codons)

format_str = "{0:<5}    {1:<6}    {2}"
for i in range(0, 64):
    if (i == 0):
        print(format_str.format('Codon', 'Number', 'Frequency (/1000)'))
    print(format_str.format(codons[i], number_of_occurrences[i], 1000 * number_of_occurrences[i] / total_codons))

