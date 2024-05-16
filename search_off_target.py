import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Look for primer sequence AS IS (reverse strand) plus NGG
# and complement of primer sequence (forward strand) plus NCC


def hamming_distance(s1, s2):
    """
    Returns the Hamming distance between two strings of equal length.

        Parameters:
            s1 (str): First string to check
            s2 (str): Second string to check

        Returns:
            h_dist (int): Hamming distance between s1 and s2
    """
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")

    h_dist = sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

    return h_dist


def search_genome(
    input_file, primer_sequence, output_file, max_mismatches=5, append=False, expand=expand
):
    """
    Takes an input fasta reference sequence file (e.g. a genome), a primer
    sequence, and optional maximum number of mismatches, and returns a bed
    file with off-target mutation sites where the primer could bind.

        Parameters:
            input_file (str): Filename of the input fasta file
            guide primer_sequence (str): A guide primer sequence (nucleotide) to search for
            output_file (str): Filename of the output bed file

        Returns:
            A bed file with off-target mutation sites in the input reference
    """
    if primer_sequence[-3:] == "NGG":
        primer_sequence = primer_sequence[:-3]

    search_length = len(primer_sequence) + 3

    if append:
        mode = "a"
    else:
        mode = "w"

    with open(output_file, mode) as out:
        if not append:
            out.write("track itemRgb=On\n")
        for seq in SeqIO.parse(input_file, "fasta"):
            for i in range(0, len(seq) - search_length):
                start = i
                end = i + search_length
                this_seq = str(seq[start:end].seq)
                match = False

                for direction in ("fwd", "rev"):
                    if direction == "rev":
                        search_seq = primer_sequence
                        mismatches = hamming_distance(this_seq[:-3], search_seq)
                        if mismatches <= max_mismatches and this_seq[-2:] == "GG":
                            match = True
                            strand = 2
                    elif direction == "fwd":
                        search_seq = str(Seq(primer_sequence).reverse_complement())
                        mismatches = hamming_distance(this_seq[3:], search_seq)
                        if mismatches <= max_mismatches and this_seq[:2] == "CC":
                            match = True
                            strand = 1

                if match:
                    name = "."
                    if max_mismatches > 0:
                        score = int(1000 - ((mismatches / max_mismatches) * 1000))
                    else:
                        score = 1000

                    green = int((score / 1000) * 255)
                    red = 255 - green
                    blue = 0
                    rgb = f"{red},{green},{blue}"

                    if strand == 1:
                        strand_str = "+"
                    else:
                        strand_str = "-"
                    out.write(
                        f"{seq.id}\t{start-expand}\t{end+expand}\t{name}\t{score}\t{strand_str}\t{start-expand}\t{end+expand}\t{rgb}\n"
                    )


if __name__ == "__main__":
    # Set up command line arguments:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", "--input", type=str, help="input fasta filename", required=True
    )
    parser.add_argument(
        "-o", "--output", type=str, help="output bed file name", required=True
    )
    parser.add_argument(
        "-p",
        "--primer",
        type=str,
        help="guide primer sequence (*excluding* NGG)",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--mismatches",
        type=int,
        help="maximum number of base mismatches in primer sequences",
        default=5,
    )
    parser.add_argument(
        "-e",
        "--expand",
        type=int,
        help="number of bases to expand out from match",
        default=0,
    )
    parser.add_argument("-a", "--append", action="store_true")

    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    primer_sequence = args.primer
    max_mismatches = args.mismatches
    expand = args.expand
    append = args.append

    search_genome(
        input_file,
        primer_sequence,
        output_file,
        max_mismatches=max_mismatches,
        expand=expand,
        append=append,
    )
