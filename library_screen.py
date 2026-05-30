
#---------------------------------------------------------------------------------------------------------------------------------------
# IMPORT PACKAGES
#---------------------------------------------------------------------------------------------------------------------------------------

# Import Packages 
from Bio import SeqIO
from Bio.Seq import Seq

#---------------------------------------------------------------------------------------------------------------------------------------
# CREATE SEQUENCE FILE
#---------------------------------------------------------------------------------------------------------------------------------------

# Create FASTA sequence file
# The sequence and its header
header = ">Cetuximab_scFv"
CDR3_H = "QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSINKDNSKSQVFFKMNSLQSNDTAIYYCARALTYYDYEFAYWGQGTLVTVSA"
CDR3_L = "DILLTQSPVILSVSPGERVSFSCRASQSIGTNIHWYQQRTNGSPRLLIKYASESISGIPSRFSGSGSGTDFTLSINSVESEDIADYYCQQNNNWPTTFGAGTKLELK"
Linker = "GGGGSGGGGSGGGGS"

dna_seq = CDR3_H + Linker + CDR3_L


# Create and save the FASTA file
with open("scFv_Cetuximab.fasta", "w") as file:
    file.write(f"{header}\n{dna_seq}\n")

# Mask CDR3 Regions for Mutations
fasta_file = "scFv_Cetuximab.fasta"

motif_CDR3_L_target = "NNNW"

motif_CDR3_H_target = "TYYDY"

motif_CDR2_H = "WSGGN"

def get_CDR_positions(fasta_file, motif):
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        start = seq.find(motif)

        if start == -1:
            raise ValueError("Motif not found")
        
        end = start + len(motif)
        
        print("0-based Python indices", (start, end))
        print("1-based Python indices", (start + 1, end + 1))

        return start, end

CDR3_L_pos = get_CDR_positions(fasta_file, motif_CDR3_L_target)
CDR3_H_pos = get_CDR_positions(fasta_file, motif_CDR3_H_target)
CDR2_H_pos = get_CDR_positions(fasta_file, motif_CDR2_H)

CDR_ranges = {
    "VH_CDR3": CDR3_H_pos,
    "VH_CDR2": CDR2_H_pos,
    "VL": CDR3_L_pos
}

#---------------------------------------------------------------------------------------------------------------------------------------
# MASK FOR MUTATIONS 
#---------------------------------------------------------------------------------------------------------------------------------------

# Define Mask Function
def mask_sequence(seq, start, end, mask_char="X"):
    seq_list = list(seq)
    seq_list[start:end] = [mask_char] * (end - start)
    return "".join(seq_list)



def mask_fasta(input_fasta, output_fasta, cdr_ranges):
    masked_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq)

        vh_cdr3_start, vh_cdr3_end = CDR_ranges["VH_CDR3"]
        vh_cdr2_start, vh_cdr2_end = CDR_ranges["VH_CDR2"]
        vl_start, vl_end = CDR_ranges["VL"]

        seq = mask_sequence(seq, vh_cdr3_start, vh_cdr3_end)
        seq = mask_sequence(seq, vh_cdr2_start, vh_cdr2_end)
        seq = mask_sequence(seq, vl_start, vl_end)

        record.seq = Seq(seq)   

        masked_records.append(record)

    SeqIO.write(masked_records, output_fasta, "fasta")

mask_fasta(
    "scFv_Cetuximab.fasta",
    "scFv_Cetuximab_masked.fasta",
    CDR_ranges
)
