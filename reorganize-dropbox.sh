echo Rearranging data from $OLD_BASE to $NEW_BASE...

tree_filepath=data/trees.txt
traits=("B2P" "M2F" "M2F_NOEURY" "S2E")
methods=("aBSREL" "BUSTED-E" "BUSTED-PH" "RELAX")
for reference_gene_id in $(cat tables/genes.txt); do
  echo "  for gene $reference_gene_id..."
  gene_directory=$NEW_BASE/$reference_gene_id
  mkdir -p $gene_directory

  for tree in $(cat tables/trees.txt); do
    root_to_tip_directory=$gene_directory/root_to_tip/$tree
    mkdir -p $root_to_tip_directory
    old_root_to_tip_json=$OLD_BASE/aBSREL-ROOT-TO-TIP/${tree}_${reference_gene_id}_aBSREL_R2T.json
    new_root_to_tip_json=$root_to_tip_directory/absrel.json
    cp -n $old_root_to_tip_json $new_root_to_tip_json
  done

#  for trait in "${traits[@]}"; do
#    for method in "${methods[@]}"; do
#      old_json=$OLD_BASE/$trait/$method/${reference_gene_id}_${trait//_/}_${method//-/}.json
#      new_dir=$gene_directory/$method/$trait
#      new_json=$new_dir/hyphy.json
#      mkdir -p $new_dir
#      cp $old_json $new_json
#    done
#  done

done

echo "...done!"
