from Bio import SeqIO
import sys

annotation_file = sys.argv[1]
annotation_prefix = ".".join(annotation_file.split(".")[:-1])
rnt_file = annotation_prefix + ".rnt"
ptt_file = annotation_prefix + ".ptt"
fasta_file = annotation_prefix + ".fna"
r = SeqIO.parse(annotation_file, "gb")
for record in r:
   fasta_file = open(fasta_file, "a")
   SeqIO.write(record, fasta_file, "fasta")
   fasta_file.close()
   record.features = [f for f in record.features if f.type == "rRNA" or f.type == "tRNA"]
   fout = open(rnt_file, "a")
   fout.write("{0} - 0..{1}\n".format(record.description, len(record)))
   fout.write("{0} RNAs\n".format(len(record.features)))
   fout.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
   strand = {1:'+', -1:'-'}
   for f in record.features:
        fout.write("{0}\n".format("\t".join([str(f.location.start).replace("<", "")+".."+str(f.location.end).replace(">", ""),strand[f.strand],str(abs(f.location.start-f.location.end)),'-',f.qualifiers["locus_tag"][0],f.qualifiers["locus_tag"][0],f.qualifiers["locus_tag"][0],"-",f.qualifiers["product"][0]])))
   fout.close()

r = SeqIO.parse(annotation_file, "gb")
for record in r:
   record.features = [f for f in record.features if f.type == "CDS"]
   fout = open(ptt_file, "a")
   fout.write("{0} - 0..{1}\n".format(record.description, len(record)))
   fout.write("{0} proteins\n".format(len(record.features)))
   fout.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
   for f in record.features:
       fout.write("{0}\n".format("\t".join([str(f.location.start).replace("<", "")+".."+str(f.location.end).replace(">", ""),strand[f.strand],str(abs(f.location.start-f.location.end)),'-',f.qualifiers["locus_tag"][0],f.qualifiers["locus_tag"][0],f.qualifiers["locus_tag"][0],"-",f.qualifiers["product"][0]])))
   fout.close()
