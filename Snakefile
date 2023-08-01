from catfish import mean_pss


rule extract_fasta_names:
  input:
    'data/fasta'
  output:
    'data/fasta_names.txt'
  shell:
    'ls {input} | sed "s/\.fasta//g" > {output}'

rule absrel:
  input:
    alignment='data/Stephen/aBSREL_Input/{name}.fasta',
    tree='data/Stephen/aBSREL_Input/{name}.nwk'
  output:
    'data/absrel/{name}.ABSREL.json' 
  shell:
    'mpirun -np 16 HYPHYMPI absrel --alignment {input.alignment} --tree {input.tree} --branches Foreground --output {output}'

rule extract:
  input:
    rules.absrel.output[0]
  output:
    'data/absrel/{name}.extract.json' 
  run:
    mean_pss(input[0], output[0])
