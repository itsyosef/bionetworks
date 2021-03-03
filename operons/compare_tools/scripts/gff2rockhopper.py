# from Bio import SeqIO
import gffutils
import sys
import os

annotation_file = sys.argv[1]
annotation_prefix = ".".join(annotation_file.split(".")[:-1])
rnt_file = annotation_prefix + ".rnt"
ptt_file = annotation_prefix + ".ptt"

db = gffutils.create_db("562.47719.gff", "temp_db")
db = gffutils.FeatureDB("temp_db")
r = """list of chromosomes"""
for record in r:
   fout = open(rnt_file, "a")
   fout.write("{0} - 0..{1}\n".format(record.description, len(record)))
   fout.write("{0} RNAs\n".format(len(record.features)))
   fout.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
   strand = {1:'+', -1:'-'}
    #ToDo: Edit the following to use correct syntax in getting the info from a gffutils.Feature and correxponding .attributes info (has same information as below)
   for f in db.features_of_type(["rRNA", "tRNA"]):
        fout.write("{0}\n".format("\t".join([str(f.location.start).replace("<", "")+".."+str(f.location.end).replace(">", ""),strand[f.strand],str(abs(f.location.start-f.location.end)),'-',f.qualifiers["locus_tag"][0],f.qualifiers["locus_tag"][0],f.qualifiers["locus_tag"][0],"-",f.qualifiers["product"][0]])))
   fout.close()

for record in r:
   fout = open(ptt_file, "a")
   fout.write("{0} - 0..{1}\n".format(record.description, len(record)))
   fout.write("{0} proteins\n".format(len(record.features)))
   fout.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
   for f in db.features_of_type(["CDS"]):
       fout.write("{0}\n".format("\t".join([str(f.location.start).replace("<", "")+".."+str(f.location.end).replace(">", ""),strand[f.strand],str(abs(f.location.start-f.location.end)),'-',f.qualifiers["locus_tag"][0],f.qualifiers["locus_tag"][0],f.qualifiers["locus_tag"][0],"-",f.qualifiers["product"][0]])))
   fout.close()
