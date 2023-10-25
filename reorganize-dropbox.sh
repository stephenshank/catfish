echo Rearranging data from $OLD_BASE to $NEW_BASE...

reference_directory=$OLD_BASE/Neoarius_graeffei_reference_sequences
tree_filepath=$OLD_BASE/trees.txt
traits=("B2P" "M2F" "M2F_NOEURY" "S2E")
methods=("aBSREL" "BUSTED" "BUSTED-PH" "RELAX")
for reference_gene_fasta_filepath in $(ls $reference_directory/C*); do
  reference_gene_fasta_basename=$(basename $reference_gene_fasta_filepath)
  reference_gene_id=${reference_gene_fasta_basename%.fasta}
  gene_directory=$NEW_BASE/$reference_gene_id
  mkdir -p $gene_directory

  old_reference_filepath=$reference_directory/$reference_gene_fasta_basename
  new_reference_filepath=$gene_directory/reference.fasta
  cp -n $old_reference_filepath $new_reference_filepath

  unaligned_fasta_directory=$OLD_BASE/FASTA_FILES
  old_unaligned_fasta_filepath=$unaligned_fasta_directory/$reference_gene_fasta_basename
  new_unaligned_fasta_filepath=$gene_directory/unaligned.fasta
  cp -n $old_unaligned_fasta_filepath $new_unaligned_fasta_filepath

  for tree in $(cat tables/trees.txt); do
    root_to_tip_directory=$gene_directory/root_to_tip/$tree
    mkdir -p $root_to_tip_directory
    old_root_to_tip_newick=$OLD_BASE/aBSREL_root-to-tip/${tree}_$reference_gene_id.nwk
    new_root_to_tip_newick=$root_to_tip_directory/tree.nwk
    cp -n $old_root_to_tip_newick $new_root_to_tip_newick
  done

  for trait in "${traits[@]}"; do
    for method in "${methods[@]}"; do
      old_newick=$OLD_BASE/$trait/$method/allraxml_$reference_gene_id.nwk
      new_dir=$gene_directory/$method/$trait
      new_newick=$new_dir/tree.nwk
      mkdir -p $new_dir
      cp $old_newick $new_newick
    done
  done

done

echo "...done!"
