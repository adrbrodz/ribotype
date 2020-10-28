import argparse
import os
import random
import string
import subprocess
import sys
import tempfile

from Bio import SeqIO

class Riboswitch:
    def __init__(self, index=0, chromosome='',
                 ribo_class='', strand='',
                 aptamer_start=0, aptamer_stop=0,
                 rbs_position=0, rbs_pattern='',
                 start_codon='', start_codon_position=0,
                 riboswitch_sequence='', rho_independent=False,
                 rho_start='---', rho_stop='---',
                 rho_score='---', rho_seq='---',
                 aptamer_sequence='', constrain='',
                 constrain_full='', aptamer_secondary_structure='',
                 riboswitch_secondary_no_constrain='',
                 riboswitch_secondary_structure='',
                 translation_inhibition=False, translation_promotion=False,
                 shift=0, inhibition_score=0,
                 rho_sequence_position='---', rbs_sequence_position=0,
                 rbs_region_constrain='', rbs_region_no_constrain='',
                 riboswitch_mechanism='unknown', riboswitch_length=0
                 ):
        self.index = index
        self.chromosome = chromosome
        self.ribo_class = ribo_class
        self.strand = strand
        self.aptamer_start = aptamer_start
        self.aptamer_stop = aptamer_stop
        self.rbs_position = rbs_position
        self.rbs_pattern = rbs_pattern
        self.start_codon = start_codon
        self.start_codon_position = start_codon_position
        self.riboswitch_sequence = riboswitch_sequence
        self.rho_independent = rho_independent
        self.rho_start = rho_start
        self.rho_stop = rho_stop
        self.rho_score = rho_score
        self.rho_seq = rho_seq
        self.aptamer_sequence = aptamer_sequence
        self.constrain = constrain
        self.constrain_full = constrain_full
        self.aptamer_secondary_structure = aptamer_secondary_structure
        self.riboswitch_secondary_no_constrain = \
             riboswitch_secondary_no_constrain
        self.riboswitch_secondary_structure = riboswitch_secondary_structure
        self.translation_inhibition = translation_inhibition
        self.translation_promotion = translation_promotion
        self.shift = shift
        self.inhibition_score = inhibition_score
        self.rho_sequence_position = rho_sequence_position
        self.rbs_sequence_position = rbs_sequence_position
        self.rbs_region_constrain = rbs_region_constrain
        self.rbs_region_no_constrain = rbs_region_no_constrain
        self.riboswitch_mechanism = riboswitch_mechanism
        self.riboswitch_length = riboswitch_length

    def cmsearch_tabfile_parser(self, cmsearch_output):
        """Parser for cmsearch output with the tabfile format ('--tblout'). 

        Args:
            cmsearch_output: The output file.

        Returns:
            riboswitches: List of Riboswitch class objects.
        """
        riboswitches = []
        self.index = 0

        with open(cmsearch_output) as cmsearch_output:
            
            for line in cmsearch_output:
                if not line.startswith('#'):
                    line = line.split()
                    self.chromosome = line[0]
                    self.ribo_class = line[2]
                    self.aptamer_start = int(line[7])
                    self.aptamer_stop = int(line[8])
                    self.strand = line[9]
                    self.index += 1
                    
                    riboswitches.append(
                        Riboswitch(
                                   self.index,
                                   self.chromosome,
                                   self.ribo_class,
                                   self.strand,
                                   self.aptamer_start,
                                   self.aptamer_stop
                    ))
    
        return riboswitches
    
    def write(self, output_file):
        riboswitch_data = (
                           '{index}' 
                      '\t' '{family}'
                      '\t' '{chromosome}'
                      '\t' '{riboswitch_start}'
                      '\t' '{riboswitch_stop}'
                      '\t' '{strand}' 
                      '\t' '{riboswitch_length}'
                      '\t' '{riboswitch_sequence}'
                      '\t' '{riboswitch_secondary}'
                      '\t' '{riboswitch_secondary_constrain}'
                      '\t' '{aptamer_start}'
                      '\t' '{aptamer_stop}'
                      '\t' '{aptamer_sequence}'
                      '\t' '{aptamer_secondary_sequence}'
                      '\t' '{constrain}'
                      '\t' '{start_codon}'
                      '\t' '{start_codon_position}'
                      '\t' '{rho_start}'
                      '\t' '{rho_stop}'
                      '\t' '{rho_sequence}'
                      '\t' '{rbs_pattern}'
                      '\t' '{rbs_region}'
                      '\t' '{rbs_region_constrain}'
                      '\t' '{rho_score}'
                      '\t' '{mechanism}'
                      '\n'
        .format(
            index=self.index,
            chromosome=self.chromosome,
            family=self.ribo_class,
            strand=self.strand,
            aptamer_start=self.aptamer_start,
            aptamer_stop=self.aptamer_stop,
            aptamer_sequence=self.aptamer_sequence,
            aptamer_secondary_sequence=self.aptamer_secondary_structure,
            constrain=self.constrain,
            rbs_position=self.rbs_position,
            rbs_sequence_position=self.rbs_sequence_position,
            rbs_pattern=self.rbs_pattern,
            start_codon=self.start_codon,
            start_codon_position=self.start_codon_position,
            rho_start=self.rho_start,
            rho_stop=self.rho_stop,
            rho_sequence=self.rho_seq,
            rho_sequence_position=self.rho_sequence_position,
            riboswitch_length=self.riboswitch_length,
            riboswitch_sequence=self.riboswitch_sequence,
            riboswitch_secondary=self.riboswitch_secondary_no_constrain,
            riboswitch_secondary_constrain=self.riboswitch_secondary_structure,
            rbs_region=self.rbs_region_no_constrain,
            rbs_region_constrain=self.rbs_region_constrain,
            rho_score=self.rho_score,
            translation_inhibition=self.translation_inhibition,
            translation_promotion=self.translation_promotion,
            rho_independent=self.rho_independent,
            mechanism=self.riboswitch_mechanism,
            riboswitch_start=self.aptamer_start,
            riboswitch_stop=self.start_codon_position
        ))
        output_file.write(riboswitch_data)


class RhoTerm:
    def __init__(self, index=0,
                 chromosome='', strand='',
                 rho_start=0, rho_stop=0,
                 rho_score=0, rho_seq=''
                 ):
        self.index = index
        self.chromosome = chromosome
        self.strand = strand
        self.rho_start = rho_start
        self.rho_stop = rho_stop
        self.rho_score = rho_score
        self.rho_seq = rho_seq
    
    def transterm_parser(self, transterm_output):
        """Parser for TranstermHP output file.
        Terminator sequences are converted to RNA. Terminators on the negative
        strand are coverted to reverse complementary sequences.

        Args:
            transterm_output: The output file.

        Returns:
            rho-terminators: List of RhoTerm class objects.
        """
        rho_terminators = []
        index = -1

        with open(transterm_output) as l:
            rhoterm_output = l.readlines()

        for record in SeqIO.parse(args.input, "fasta"):
            chrom_start_line = 0
            chrom_id = str(record.id)
            chrom_seq = str(record.seq)
            self.chromosome = chrom_id

            chrom_start_line = find_chromosome(
                chrom_start_line,
                rhoterm_output, chrom_id
                )
            rho_terminators, index = RhoTerm.extract_rho_terms(
                RhoTerm, index,
                rhoterm_output, chrom_start_line,
                rho_terminators, chrom_seq
                )

        return rho_terminators

    def extract_rho_terms(self, index, output_file, chrom_start_line,
                          rho_terminators, chrom_seq):

        for line in output_file[chrom_start_line:]:

            if "SEQUENCE" in line:
                break

            if 'TERM' in line:
                line = line.split()

                self.index = index
                self.rho_start = int(line[2])
                self.rho_stop = int(line[4])
                self.strand = line[5]
                self.rho_score = line[7]

                if self.strand == '+':
                    self.rho_seq = chrom_seq[self.rho_start-1:self.rho_stop]
                    self.rho_seq = self.rho_seq.replace('T', 'U')

                elif self.strand == '-':
                    rho_seq = chrom_seq[self.rho_stop:self.rho_start]
                    self.rho_seq = reverse_complement(rho_seq)
                    self.rho_seq = self.rho_seq.replace('T', 'U')

                index += 1

                rho_terminators.append(
                    RhoTerm(
                            self.index,
                            self.chromosome,
                            self.strand,
                            self.rho_start,
                            self.rho_stop,
                            self.rho_score,
                            self.rho_seq
                ))

        return rho_terminators, index


class RBS:
    def __init__(self, chromosome='', strand='',
                 rbs_position=0, rbs_pattern='',
                 start_codon='', start_codon_position=0,
                 shift=0
                 ):
        self.strand = strand
        self.chromosome = chromosome
        self.rbs_position = rbs_position
        self.rbs_pattern = rbs_pattern
        self.start_codon = start_codon
        self.start_codon_position = start_codon_position
        self.shift = shift

    def rbs_finder_parser(self, rbs_finder_output):
        """Parser for rbs_finder.pl output file.

        Args:
            rbs_finder_output: The output file.
            chromosome: The chromosome id of input sequence.

        Returns:
            RBSlist: List of RBS class objects.

        """    
        RBSlist = []

        with open(rbs_finder_output) as l:
            rbs_output = l.readlines()

        for record in SeqIO.parse(args.input, "fasta"):
            chrom_start_line = 0
            chrom_id = str(record.id)
            self.chromosome = chrom_id

            chrom_start_line = find_chromosome(chrom_start_line,
                                               rbs_output, chrom_id)
                
            for line in rbs_output[chrom_start_line:]:

                if ">" in line and not self.chromosome in line:
                    break

                if 'orf' in line and '---' not in line:
                    line = line.split()

                    orf_start = int(line[1])
                    orf_stop = int(line[2])

                    self.strand = '+' if orf_start < orf_stop else '-'
                    self.rbs_pattern = line[3]
                    self.rbs_position = int(line[4])
                    self.start_codon = line[5].replace('T', 'U')
                    self.shift = int(line[6])
                    self.start_codon_position = int(line[8])

                    RBSlist.append(
                            RBS(
                                self.chromosome,
                                self.strand,
                                self.rbs_position,
                                self.rbs_pattern,
                                self.start_codon,
                                self.start_codon_position,
                                 self.shift
                    ))

        return RBSlist

class Aptamer:
    def __init__(self, chromosome='', aptamer_start=0,
                 aptamer_sequence='', constrain=''
                 ):
        self.aptamer_start = aptamer_start
        self.chromosome = chromosome
        self.aptamer_sequence = aptamer_sequence
        self.constrain = constrain

    def cmsearch_parser(self, cmsearch_output):
        """Parser for cmsearch result.

        Args:
            cmsearch_output: The output file.

        Returns:
            aptamers: List of Aptamer class objects.
        
        """
        aptamers = []
        x = 0

        with open(cmsearch_output) as l:
            record = l.readlines()
        
        for line in record:
            x += 1
            if "[]" in line:
                self.chromosome = record[x+5].split()[0]
                self.aptamer_start = int(line.split()[9])

                aptamer_seq1 = record[x+5]
                aptamer_seq2 = record[x+12]
                aptamer_seq3 = record[x+19]

                constrain = Aptamer.extract_constrain(Aptamer, x, record)
                
                # If non-conserved nucleotides are found use remove_gap.
                if "*[" in aptamer_seq1:
                    aptamer_seq1, constrain = Aptamer.remove_gap(
                    Aptamer, aptamer_seq1, constrain
                    )
                if "*[" in aptamer_seq2:
                    aptamer_seq2, constrain = Aptamer.remove_gap(
                    Aptamer, aptamer_seq2, constrain
                    )
                if "*[" in aptamer_seq3:
                    aptamer_seq3, constrain = Aptamer.remove_gap(
                    Aptamer, aptamer_seq3, constrain
                    )
                aptamer_sequence = (
                  Aptamer.seq_parse(Aptamer, aptamer_seq1)
                + Aptamer.seq_parse(Aptamer, aptamer_seq2)
                + Aptamer.seq_parse(Aptamer, aptamer_seq3)
                )
                self.constrain = constrain

                # If gaps are found in sequence the constrain is
                # modified with the correct_constrain function.
                if "-" in aptamer_sequence:
                    self.constrain = Aptamer.correct_constrain(
                    Aptamer, aptamer_sequence, constrain
                    )

                # Removing gaps from the sequence.
                self.aptamer_sequence = (aptamer_sequence
                                        .upper()
                                        .replace("-", "")
                                        )
                constrain = Aptamer.constrain_remove_non_canonical(Aptamer, self.aptamer_sequence, self.constrain)
                self.constrain = constrain
                aptamers.append(Aptamer(self.chromosome, self.aptamer_start, self.aptamer_sequence, self.constrain))

        return aptamers

    def extract_constrain(self, record_line, output):
        x = record_line
        record = output
        long_constrain = True

        constrain1 = Aptamer.parse_constrain(Aptamer, record[x+2])
        string_start = record[x+2].find('CS')-len(constrain1)-1
        constrain1_NC = record[x+1][string_start : string_start
                                                 + len(constrain1)]

        if 'HMM' not in record[x+9] and 'rank' not in record[x+9]:
            constrain2 = Aptamer.parse_constrain(Aptamer, record[x+9])
            constrain2_NC = record[x+8][string_start : string_start
                                                     + len(constrain2)]
        else:
            constrain2 = ''
            constrain3 = ''
            constrain2_NC = ''
            constrain3_NC = ''
            long_constrain = False

        if 'HMM' not in record[x+16] and 'rank' not in record[x+16] \
                 and long_constrain == True:
            constrain3 = Aptamer.parse_constrain(Aptamer, record[x+16])
            constrain3_NC = record[x+15][string_start : string_start
                                                      + len(constrain3)]
        else:
            constrain3 = ''
            constrain3_NC = ''

        constrain = (constrain1 + constrain2 + constrain3)
        constrain = constrain.replace("…", "...")
        constrain_NC = constrain1_NC + constrain2_NC + constrain3_NC
        constrain = Aptamer.remove_nc(Aptamer, constrain, constrain_NC)
    
        return constrain

    def parse_constrain(self, constrain):
        constrain = (
                     constrain.replace(" ","").replace("CS","")
                              .replace("_", ".").replace("\n","")
                              .replace(":", ".").replace("(","<")
                              .replace(")",">").replace(",", ".")
                              .replace("-", ".").replace("[", "<")
                              .replace("]",">").replace("{", "<")
                              .replace("}",">")
        )
        return constrain

    def remove_nc(self, constrain, constrain_NC):
        c = 0
        for symbol in constrain_NC:
            c += 1
            if symbol == "v":
                constrain = constrain[:c-1]+'.'+constrain[c:]

        return constrain

    def remove_gap(self, seq, constrain):
        """Marks non-conserved nucleotides from cmsearch output as "X".
        
        Non-conserved region can differ in length and is marked as "*[ n]*"
        where n is the number of non-canonnical nucleotides. During folding
        "X" is interpreted as non-cannonical nucleotide and coverted to "_". 
        Non-conserved sites in the corresponding constrain are marked as "~"
        This modifications allow construction of input that is properly read
        by RNAfold.

        Args:
            seq: The string with aptamer sequence.
            constrain: The string with corresponding constrain.

        Returns:
            seq: String with modified sequence.
            constrain: String with modified correspodning constrain.

        """
        # Extracting inforation about amount of non-conserved nucleotides
        n = int(
                seq[seq.find("*[")+2 : seq.find("*[")+4]
                .replace("]", "").replace(" ", "")
        )
        # Indicating non-conserved nucleotides as "X"
        seq = (
                seq.replace("*[ ", "")
                   .replace("*[", "")
                   .replace("]*", "")
                   .replace(str(n), n*"X")
        )
        # Modyfing correspodning constrain
        constrain = (
                     constrain[: constrain.find("~")]
                               + n*"."
                               + constrain[constrain.find("~"):]
                               .replace("~", "")
        )
        
        return seq, constrain

    def seq_parse(self, seq):
        """Extracts sequence from string if not empty."""
        if not seq == '\n' and not 'Query' in seq and not 'local' in seq:
            seq = seq.split()
            seq = seq[2]
        else:
            seq = ''
        return seq

    def correct_constrain(self, aptamer_sequence, constrain):
        """Corrects gaps in corresponding constrain."""
        z = 0
        for symbol in aptamer_sequence:
            z += 1
            if symbol == "-":
                constrain = constrain[:z-1]+constrain[z:]
                z+=-1

        return constrain
        
    def constrain_remove_non_canonical(self, sequence, constrain):
        pairs = {}
        pstack = []
        canonical_pairs = ['A:U','U:A','U:G','G:U','C:G','G:C']

        for i, c in enumerate(constrain):
            if c == '<':
                pstack.append(i)
            elif c == '>':
                if len(pstack) == 0:
                    raise IndexError("No matching closing parens at: "
                                     + str(i))
                pairs[pstack.pop()] = i

        if len(pstack) > 0:
            raise IndexError("No matching opening parens at: "
                             + str(pstack.pop()))

        for pair in pairs:
            if (sequence[pair] + ':' + sequence[pairs[pair]]) \
                    not in canonical_pairs:
                constrain = constrain[:pair] + "." + constrain[pair+1:]
                constrain = constrain[:pairs[pair]] \
                            + "." + constrain[pairs[pair]+1:]

        return constrain

def find_chromosome(chrom_start_line, output_file, chrom_id):
    chrom_start_line = chrom_start_line
    for line in output_file:
        chrom_start_line += 1
        if chrom_id in line:
            return chrom_start_line

def rnafold(handle, sequence, constrain = None, index = None):
    """Folds the input sequence with optional constrain using RNAfold.

    Args:
        handle: Identifyer for writing input and output files.
        sequence: The string with aptamer sequence.
        constrain: The string with corresponding constrain (optional).

    Returns:
        secondary_structure: String with the generated secondary structure.

    """
    with open(handle, 'w') as input_file:
        if constrain == None:
            input_file.write(sequence)
        else:
            constrain = constrain+"."*(len(sequence)-len(constrain))
            input_file.write(sequence+"\n"+constrain)
    if constrain == None:
        subprocess.run(['RNAfold', '--noPS', "-o"+handle+'_output', handle])
    else:
        if '<>' not in constrain and '<.>' not in constrain and '<..>' not in constrain:
            subprocess.run(['RNAfold', '-C','--enforceConstraint', '--noPS', "-o"+handle+'_output', handle])
            if index != None:
                print("Aptamer (" + str(index) + ") folded using --enforceConstraint.")
        else:
            subprocess.run(['RNAfold', '-C', '--noPS', "-o"+handle+'_output', handle])
            if index != None:
                print("Aptamer (" + str(index) + ") folded using -C.")
    with open(handle+'_output') as output_file:
        output_file.readline()
        secondary_structure = output_file.readline().split()[0]
    return secondary_structure

def complement(dna):
    complement_dna = []
    for nucleotide in dna:
        if nucleotide == 'A':
            complement_dna.append('T')
        elif nucleotide == 'T':
            complement_dna.append('A')
        elif nucleotide == 'C':
            complement_dna.append('G')
        elif nucleotide == 'G':
            complement_dna.append('C')
    complement_dna = "".join(complement_dna)
    return complement_dna

def reverse_complement(dna):
    reverse_complement_dna = complement(dna)
    reverse_complement_dna = list(reverse_complement_dna)
    reverse_complement_dna.reverse()
    reverse_complement_dna = "".join(reverse_complement_dna)
    return reverse_complement_dna

def generate_prefix():
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join((random.choice(lettersAndDigits) for i in range(8)))


### ArgumentParser for command line options
parser = argparse.ArgumentParser()
parser.add_argument("input", help="input file with sequence in FASTA format")
parser.add_argument("output", help="output file name")
parser.add_argument("-t", "--temp", action="store_true", help="keep temporary files", default=False)
parser.add_argument("-l", "--maxlength", help="set the maximum length of riboswitch sequence (default is 1000)", default =1000)
parser.add_argument("-r", "--range", help="modify the 'range' option for RBSfinder (default is 25)", default='25')
parser.add_argument("-rs", "--rhoscore", help="set the minimum threshold for rho-independent terminators confidence value (default is 80)", default=80)
parser.add_argument("-ti", "--inhibition", help="set the threshold for translation inhibition score (default is 60)", default=20)
parser.add_argument("-tp", "--promotion", help="set the threshold for translation promotion score (default is 60)", default=20)
args = parser.parse_args()

print('\n\nRiboType\n\nRiboType searches for riboswitches in the input '
      'FASTA file and tries to identify their molecular mechanims of function '
      'based on presence of rho-indepdendant terminators and changes in the '
      'secondary structure of RBS regions\n\n')

prefix = generate_prefix()+'_'
print('Prefix generated for this run: ' + prefix + '\n')
print('Searching for the input FASTA file...')
try:
    with open(args.input) as input_fasta:
        header = input_fasta.readline()
        if ">" in header:
            pass
        else:
            print("ERROR: Header line should start with '>'.\n")
            sys.exit()
except FileNotFoundError:
    print('ERROR: FASTA file not found.\n')
    sys.exit()
print('*SUCCESSFUL*\n')

print("Searching for aptamers with Infernal 'cmsearch'...")
### Localizing riboswitch aptamers on both strands in the input sequence
if not args.temp:
    with tempfile.NamedTemporaryFile() as rb, tempfile.NamedTemporaryFile() as rb_detail:
        try:
            subprocess.check_output(['cmsearch','-g', '--cut_ga', '--rfam', '--nohmmonly', '--notrunc', '--tblout', rb.name, 'rfam.cm', args.input])
            subprocess.check_output(['cmsearch','-g', '--cut_ga', '--rfam', '--nohmmonly', '--notrunc', '-o', rb_detail.name, 'rfam.cm', args.input])
        except subprocess.CalledProcessError:
            print('ERROR: cmsearch subprocess failed.\n')
            sys.exit()
        riboswitches = Riboswitch.cmsearch_tabfile_parser(Riboswitch, rb.name)
        aptamers = Aptamer.cmsearch_parser(Aptamer, rb_detail.name)
        print('*SUCCESSFUL*\n')
else:
    try:
        subprocess.check_output(['cmsearch', '-g', '--rfam', '--cut_ga', '--nohmmonly', '--notrunc', '--tblout', prefix+args.output+'_riboswitches', 'rfam.cm', args.input])
        subprocess.check_output(['cmsearch', '-g', '--rfam', '--cut_ga', '--nohmmonly', '--notrunc', '-o', prefix+args.output+'_riboswitches_detail', 'rfam.cm', args.input])
    except subprocess.CalledProcessError:
        print('ERROR: cmsearch subprocess failed.\n')
# Identifying the aptamer positions in genome with cmsearch.
    riboswitches = Riboswitch.cmsearch_tabfile_parser(Riboswitch, prefix+args.output+'_riboswitches')
# Assigning the constrains from cmsearch to aptamer sequences.
    aptamers = Aptamer.cmsearch_parser(Aptamer, prefix+args.output+'_riboswitches_detail')
    print('*SUCCESSFUL*\n')

print('Searching for putative genes with Glimmer3...')
### Localizing RBS pattern and position in the input sequence
# glimmer3 subprocess for generating gene predictions required by RBSfinder
devnull = open(os.devnull, 'w')
try:
    subprocess.run(['long-orfs', args.input, prefix+'long-orfs'], stdout=devnull, stderr=devnull)
except subprocess.CalledProcessError:
    print("Glimmer3 'long-orfs' failed.")
    sys.exit()
try:
    subprocess.run(['extract', args.input, prefix+'long-orfs'], stdout=open(prefix+'train', 'w'), stderr=devnull)
except subprocess.CalledProcessError:
    print("Glimmer3 'extract' failed.")
    sys.exit()
try:
    subprocess.run(['build-icm', prefix+'icm'], stdin=open(prefix+'train'))
except subprocess.CalledProcessError:
    print("Glimmer3 'build-icm' failed.")
    sys.exit()
try:
    subprocess.run(['glimmer3', args.input, prefix+'icm', prefix], stdout=devnull, stderr=devnull)
except subprocess.CalledProcessError:
    print("Glimmer3 'glimmer3' failed.")
    sys.exit()
print("*SUCCESSFUL*\n")

print('Searching for RBS regions with rbs_finder.pl...')
try:
    subprocess.check_output(['perl', 'rbs_finder.pl', args.input, prefix+'.predict', prefix+args.output+'_rbs', args.range])
except subprocess.CalledProcessError:
    print("Perl scirpt 'rbs_finder.pl' failed.")
    sys.exit()

RBSlist = RBS.rbs_finder_parser(RBS, prefix+args.output+'_rbs')
if not args.temp:
    os.remove(prefix+args.output+'_rbs')
os.remove(prefix+'long-orfs')
os.remove(prefix+'icm')
os.remove(prefix+'train')
os.remove(prefix+'.predict')
os.remove(prefix+'.detail')
print("*SUCCESSFUL*\n")

print('Searching for rho-independent terminators with TranstermHP...')
### Identification of rho-independent terminators
# Creating fake.coords file for TranstermHP
with open(prefix+'fake.coords', 'w') as fake:
    for record in SeqIO.parse(args.input, "fasta"):
        fake.write("\tfakegene1\t1 2\t{}\n\tfakegene2\t{} {}\t{}\n".format(str(record.id), str(len(str(record.seq))),str(len(str(record.seq))-1), str(record.id)))
# TranstermHP subprocces for localisation of rho-depentent terminators in the input sequence
if not args.temp:
    with tempfile.NamedTemporaryFile() as transterm_output:
        try:
            subprocess.run(['transterm','-c 0', '-p', 'expterm.dat', args.input, prefix+'fake.coords'], stdout=open(transterm_output.name, 'w'), stderr=devnull)
        except subprocess.CalledProcessError:
            print("TranstermHP 'transterm' failed.")
            sys.exit()
        rho_terminators = RhoTerm.transterm_parser(RhoTerm, transterm_output.name)
        print("*SUCCESSFUL*\n")
else:
    try:
        subprocess.run(['transterm','-c '+str(args.rhoscore), '-p', 'expterm.dat', args.input, prefix+'fake.coords'], stdout=open(prefix+args.output+'_transterm', 'w'), stderr=devnull)
    except subprocess.CalledProcessError:
        print("TranstermHP 'transterm' failed.")
        sys.exit()
    rho_terminators = RhoTerm.transterm_parser(RhoTerm, prefix+args.output+'_transterm')
    print("*SUCCESSFUL*\n")
os.remove(prefix+'fake.coords')

# Extracting genomic sequences from the input file
sequences = {}
for record in SeqIO.parse(args.input, "fasta"):
    sequences[str(record.id)] = str(record.seq)

for r in riboswitches:
    sequence = sequences[r.chromosome]
    ### Assigning RBS and start codons to identyfied aptamer sequences.
    for u in RBSlist:
        if r.chromosome == u.chromosome:
            # For positive strand breaking the loop after assigning first RBS position following start of aptamer sequence.
            if r.strand == u.strand == '+' and u.rbs_position > r.aptamer_stop:
                r.rbs_position = u.rbs_position
                r.rbs_sequence_position = r.rbs_position - r.aptamer_start
                r.rbs_pattern = u.rbs_pattern
                r.start_codon = u.start_codon
                r.start_codon_position = u.start_codon_position
                r.shift = u.shift
                break
            # For negative strand iterating over RBS list as long as the RBS position precedes start of aptamer sequence.
            elif r.strand == u.strand == '-' and u.rbs_position < r.aptamer_stop:
                r.rbs_position = u.rbs_position
                r.rbs_sequence_position = r.aptamer_start - r.rbs_position
                r.rbs_pattern = u.rbs_pattern
                r.start_codon = u.start_codon
                r.start_codon_position = u.start_codon_position
                r.shift = u.shift
    
    ### Constructing a full riboswitch sequence based on aptamer position and assigned start codons.
    if r.strand == '+':
        r.riboswitch_sequence = sequence[r.aptamer_start:r.start_codon_position+2+r.shift+20].replace('T', 'U')
    elif r.strand == '-':
        r.riboswitch_sequence = reverse_complement(sequence[r.start_codon_position-3-r.shift-20:r.aptamer_start]).replace('T', 'U')
    r.riboswitch_length = len(r.riboswitch_sequence)

    if r.riboswitch_length < args.maxlength:

        ### Iterating over the list of rho-independent terminators and trying assign the one with highest confidence value to constructed riboswitches.
        rho_terminators_in_sequence = {}
        for t in rho_terminators:
            if r.chromosome == t.chromosome:
                if r.strand == '+':
                    if t.rho_seq in r.riboswitch_sequence:
                        #if t.rho_start in range(r.aptamer_stop, r.start_codon_position):
                        if t.rho_start in range(r.aptamer_start, r.start_codon_position):
                            rho_terminators_in_sequence[t.index] = int(t.rho_score)
                            max_score = max(rho_terminators_in_sequence.values())
                            best_rho_index = (list(rho_terminators_in_sequence.keys())[list(rho_terminators_in_sequence.values()).index(max_score)])
                        try:
                            if t.index == best_rho_index:
                                r.rho_start = t.rho_start
                                r.rho_stop = t.rho_stop
                                r.rho_score = int(t.rho_score)
                                r.rho_seq = t.rho_seq
                        except:
                            pass
                elif r.strand == '-':
                    if t.rho_seq in r.riboswitch_sequence:
                        #if t.rho_start in range(r.start_codon_position, r.aptamer_stop):
                        if t.rho_start in range(r.start_codon_position, r.aptamer_start):
                            rho_terminators_in_sequence[t.index] = int(t.rho_score)
                            max_score = max(rho_terminators_in_sequence.values())
                            best_rho_index = (list(rho_terminators_in_sequence.keys())[list(rho_terminators_in_sequence.values()).index(max_score)])
                        try:
                            if t.index == best_rho_index:
                                r.rho_start = t.rho_start
                                r.rho_stop = t.rho_stop
                                r.rho_score = int(t.rho_score)
                                r.rho_seq = t.rho_seq
                        except:
                            pass
            if not r.rho_start == '---':
                r.rho_sequence_position = r.riboswitch_sequence.find(r.rho_seq)

        ### Correcting the non-conserved nucleotides from cmsearch result in aptamer sequence to represent genomic sequence.
        for a in aptamers:
            if r.chromosome == a.chromosome:
                if r.aptamer_start == a.aptamer_start:
                    r.aptamer_sequence = a.aptamer_sequence
                    r.constrain = a.constrain
                # If non-conserved nucleotides found in aptamer sequence use the corresponding genomic sequence to fill the gaps.
                if "X" in r.aptamer_sequence:
                    aptamer_seq1 = r.aptamer_sequence[:r.aptamer_sequence.find("X")]
                    aptamer_seq1 = aptamer_seq1 + r.riboswitch_sequence[r.riboswitch_sequence.find(aptamer_seq1)+len(aptamer_seq1):r.riboswitch_sequence.find(aptamer_seq1)+len(aptamer_seq1)+r.aptamer_sequence.count("X")]
                    aptamer_seq2 = r.aptamer_sequence[r.aptamer_sequence.find("X")+r.aptamer_sequence.count("X"):]
                    r.aptamer_sequence = aptamer_seq1 + aptamer_seq2
                r.riboswitch_sequence = r.aptamer_sequence + r.riboswitch_sequence[r.riboswitch_sequence.find(r.aptamer_sequence)+len(r.aptamer_sequence):] 
            
        ### Gnerating aptamer and riboswitch secondary sequences with RNAfold
        r.aptamer_secondary_structure = rnafold(prefix+'aptamer_fold', r.aptamer_sequence, r.constrain, r.index)
        r.riboswitch_secondary_structure = rnafold(prefix+'riboswitch_fold', r.riboswitch_sequence, r.constrain)
        r.riboswitch_secondary_no_constrain = rnafold(prefix+'riboswitch_fold_no_constrain', r.riboswitch_sequence)
        os.remove(prefix+'aptamer_fold')
        os.remove(prefix+'aptamer_fold_output')
        os.remove(prefix+'riboswitch_fold') 
        os.remove(prefix+'riboswitch_fold_output')
        os.remove(prefix+'riboswitch_fold_no_constrain')
        os.remove(prefix+'riboswitch_fold_no_constrain_output')

        ### Extracting secondary structure of the RBS regions for comparison
        r.rbs_region_no_constrain = r.riboswitch_secondary_no_constrain[r.rbs_sequence_position:r.rbs_sequence_position+5]
        r.rbs_region_constrain = r.riboswitch_secondary_structure[r.rbs_sequence_position:r.rbs_sequence_position+5]
        
        ### Identyfing the riboswitch mechanism based on changes in RBS regions after applying constrain
        # If there are no changes there is no translational inhibition mechanism
        if r.rbs_region_no_constrain == r.rbs_region_constrain:
            r.inhibition_score = 0
        else:
            RBS_bind_nc = False
            RBS_bind_c = False

            if "..." in r.rbs_region_no_constrain or r.rbs_region_no_constrain == "))).." or r.rbs_region_no_constrain == "..(((":
                RBS_bind_nc = True


            if "..." in r.rbs_region_constrain or r.rbs_region_constrain == "))).." or r.rbs_region_constrain == "..(((":
                RBS_bind_c = True
            
            if RBS_bind_nc == True and RBS_bind_c == False:
                r.translation_inhibition = True

            if RBS_bind_nc == False and RBS_bind_c == True:
                r.translation_promotion = True
    
        if r.rho_score == '---':
            if r.translation_inhibition == True:
                r.riboswitch_mechanism = 'translation inhibition'
            if r.translation_promotion == True:
                r.riboswitch_mechanism = 'translation promotion'
        else:
            if r.rho_score >= 95:
                r.rho_independent = True
                r.riboswitch_mechanism = 'rho-independent termination'
            else:
                if r.translation_inhibition == True:
                    r.riboswitch_mechanism = 'translation inhibition'
                if r.translation_promotion == True:
                    r.riboswitch_mechanism = 'translation promotion'
            if not r.translation_inhibition and not r.translation_promotion and r.rho_score >= args.rhoscore:
                r.rho_independent = True
                r.riboswitch_mechanism = 'rho-independent termination'            

with open(args.output+'.csv', 'w') as ribotype_output:
    ribotype_output.write("{index}\t{family}\t{chromosome}\t{riboswitch_start}\t{riboswitch_stop}\t{strand}\t{riboswitch_length}\t{riboswitch_sequence}\t{riboswitch_secondary}\t{riboswitch_secondary_constrain}\t{aptamer_start}\t{aptamer_stop}\t{aptamer_sequence}\t{aptamer_secondary_sequence}\t{constrain}\t{start_codon}\t{start_codon_position}\t{rho_start}\t{rho_stop}\t{rho_sequence}\t{rbs_pattern}\t{rbs_region}\t{rbs_region_constrain}\t{rho_score}\t{mechanism}\n"
            .format(index='Index',chromosome='Chromosome',family='Family',strand='Strand',aptamer_start='Aptamer Start',aptamer_stop='Aptamer Stop',aptamer_sequence='Aptamer Sequence',
                    aptamer_secondary_sequence='Aptamer Secondary Structure',constrain='Constrain',rbs_position='RBS Position in Genome',rbs_sequence_position='RBS Position in Sequence',
                    rbs_pattern='RBS Pattern',start_codon='Start Codon',start_codon_position='Start Codon Position',rho_start='Rho Start',rho_stop='Rho Stop',rho_sequence='Rho Sequence',
                    rho_sequence_position='Rho Position in Sequence', riboswitch_sequence='Riboswitch Sequence',riboswitch_secondary='Riboswitch Secondary Structure',
                    riboswitch_secondary_constrain='Riboswitch Secondary Structure Constrained',rbs_region='RBS Region',rbs_region_constrain='RBS Region Constrained',
                    rho_score='Rho Score', translation_inhibition='Translation Inhibition', translation_promotion='Translation Promotion',
                    rho_independent='Rho-independent Termination', mechanism='Riboswitch Mechanism',riboswitch_length='Riboswitch Length',riboswitch_start='Riboswitch Start',
                    riboswitch_stop='Riboswitch Stop'))
    for r in riboswitches:
        if r.riboswitch_length < args.maxlength:
            r.write(ribotype_output)

print("\n*DONE*\n\nCheck " + str(args.output) +  ".csv for results.\n")