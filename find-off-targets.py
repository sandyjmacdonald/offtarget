from search_off_target import search_genome

primers = ["ATTACTAGGTAGAGGAGGAG"]

input_file = "dummy_genome.fasta"

for p in primers:
    output_file = f"{p}_offtarget.bed"
    search_genome(input_file, p, output_file, max_mismatches=5)