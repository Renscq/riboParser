grep -vE 'Annotated|Internal|Truncated|Extended|CDSFrameOverlap|Known' ribotish/gmx4_all.txt | sed 's/gene-//g' > test.txt
