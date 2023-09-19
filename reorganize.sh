for f in $(ls olddata/Stephen/aBSREL_Input/*.fasta); do
  echo "--- $f ---"
  bn=$(basename $f);
  echo "basename - $bn"
  phylo=${bn%%_*}
  echo "phylo - $phylo"
  gene_with_extension=${bn#*_}
  gene=${gene_with_extension%.fasta}
  echo "gene - $gene"
  mkdir -p data/$gene/$phylo
  cp $f data/$gene/codons.fasta
  newick=${f%.fasta}.nwk
  echo "newick - $newick"
  cp $newick data/$gene/$phylo/tree.nwk
  bn_noext=${bn%.fasta}
  oldabsrel=olddata/absrel/$bn_noext.ABSREL.json
  newabsrel=data/$gene/$phylo/absrel.json
  cp $oldabsrel $newabsrel
done
