import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# primer pattern: ATTACTAGGTAGAGGAGGAGNGG

bases = ["A", "T", "G", "C"]

def generate_binding_site(primer_sequence, max_mismatches):
    global bases
    num_mismatches = random.randint(0, max_mismatches)
    positions = random.sample(range(0, len(primer_sequence)), k=num_mismatches)
    primer_sequence = list(primer_sequence)

    for pos in positions:
        new_base = random.choice(bases)
        primer_sequence[pos] = new_base
    
    target_sequence = f"{''.join(primer_sequence)}{random.choice(bases)}GG"

    return target_sequence

def generate_read(target_length, length_variability, binding_site_probability):
    global primer_sequence, max_mismatches, bases

    trim_length = random.randint(0, int(length_variability * target_length))
    actual_length = target_length - trim_length
    add_binding_site = random.random() < binding_site_probability
    
    if add_binding_site:
        binding_site = generate_binding_site(primer_sequence, max_mismatches)
        if random.choice((0, 1)):
            binding_site = str(Seq(binding_site).reverse_complement())
            print(binding_site)
        else:
            binding_site = str(Seq(binding_site).complement())
            print(binding_site)
        start_length = random.randint(0, actual_length - len(binding_site))
        end_length = actual_length - start_length - len(binding_site)
        start_seq = "".join([random.choice(bases) for i in range(start_length)])
        end_seq = "".join([random.choice(bases) for i in range(end_length)])
        seq = f"{start_seq}{binding_site}{end_seq}"
        return seq
    else:
        seq = "".join([random.choice(bases) for i in range(actual_length)])
        return seq

def generate_read_set(num_reads, prefix):
    global target_length, length_variability, binding_site_probability

    seqs = []

    for i in range(num_reads):
        read_id = f"{prefix}_{i+1}"
        seq = Seq(generate_read(target_length, length_variability, binding_site_probability))
        rec = SeqRecord(seq, id=read_id, name="", description="")
        seqs.append(rec)

    return seqs

primer_sequence = "ATTACTAGGTAGAGGAGGAG"
max_mismatches = 5
target_length = 500000
length_variability = 0.1
binding_site_probability = 0.75
num_reads = 20
prefix = "DUMMY"

dummy_reads = generate_read_set(num_reads, prefix)
SeqIO.write(dummy_reads, "dummy_genome.fasta", "fasta")