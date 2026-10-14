import pyhmmer

alphabet = pyhmmer.easel.Alphabet.amino()

protein_file = "/Users/karthikmandigiri/Geneweaver_v3/Structural_annotation/test_data.faa"
with pyhmmer.easel.SequenceFile(protein_file, digital=True, alphabet=alphabet) as sf:
    proteins = sf.read_block()

Pfam_db = "/Users/karthikmandigiri/Geneweaver_v3/Structural_annotation/databases/Pfam-A.hmm"
with pyhmmer.plan7.HMMFile(Pfam_db) as hf:
    hmms = list(hf)

with open("pfam_hits.tsv", "w") as out:
    out.write("protein\tpfam_family\tpfam_accession\tevalue\tscore\n")
    for hits in pyhmmer.hmmsearch(hmms, proteins, cpus=0, bit_cutoffs="gathering"):
        for hit in hits.included:
            out.write(f"{hit.name}\t{hits.query.name}\t{hits.query.accession}\t{hit.evalue:.2e}\t{hit.score:.1f}\n")