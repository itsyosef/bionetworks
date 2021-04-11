##seqkit sample --two-pass -n 3000 msa_0310.fasta | seqkit subseq -r 22200:22250 | seqkit grep -s -i -v -p N | seqkit sample -n 1000 >./small_msa_pos_22k.fasta
##./harvesttools -m ./small_msa_pos_22k.fasta -B ./small_msa_pos_22k.backbone -F ./small_msa_pos_22k.ref.fasta
##python gengraphTool.py --out_file_name small_cv19_pos_22k.graph --out_file_path ~/projects/git_repos/cid_work/viral_graphs/march_11_2021/msa_0310/ --input ~/projects/git_repos/cid_work/viral_graphs/march_11_2021/msa_0310/small_cv19_pos_22k.msa.fasta make_graph_from_fasta

```
grep -v "^>" ./small_cv19_pos_22k.msa.fasta | sort | uniq
ATAGCGCAACATTACCTAAAGGCATAATG--ATGAATGTCGCAAAATATAC
ATAGTGCAACATTACCTAAAGGCATAATG--ATGAATGTCGCAAAATACAC
ATAGTGCAACATTACCTAAAGGCATAATG--ATGAATGTCGCAAAATATAC
ATAGTGCAATATTACCTAAAGGCATAATG--ATGAATGTCGCAAAATATAC
```
