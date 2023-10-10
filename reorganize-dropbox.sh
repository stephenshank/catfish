echo Rearranging data from $OLD_BASE to $NEW_BASE...

reference_directory=$OLD_BASE/Neoarius_graeffei_reference_sequences
tree_filepath=$OLD_BASE/trees.txt
for reference_gene_fasta_filepath in $(ls $reference_directory/C*); do
  reference_gene_fasta_basename=$(basename $reference_gene_fasta_filepath)
  reference_gene=${reference_gene_fasta_basename%.fasta}
  gene_directory=$NEW_BASE/$reference_gene
  mkdir -p $gene_directory

  old_reference_filepath=$reference_directory/$reference_gene_fasta_basename
  new_reference_filepath=$gene_directory/reference.fasta
  cp -n $old_reference_filepath $new_reference_filepath

  unaligned_fasta_directory=$OLD_BASE/FASTA_FILES
  old_unaligned_fasta_filepath=$unaligned_fasta_directory/$reference_gene_fasta_basename
  new_unaligned_fasta_filepath=$gene_directory/unaligned.fasta
  cp -n $old_unaligned_fasta_filepath $new_unaligned_fasta_filepath
done

echo "...done!"
