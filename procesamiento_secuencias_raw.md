# Procesamiento de secuencias 16S rRNA

## Ensamblado

```bash
#!/usr/bin/env bash
for i in $(cat muestras.txt);                                      
do                                                                 
  archivo="L$i"                                                    
  # Descomprimir y sustituir                                       
  zcat Rawdata/*${archivo}*_1.fq.gz |\                             
    sed '/^@/ s/$/:23:C5A1KACXX:1:1101:1:2 1:N:0:ATCACG/' >\       
    ${archivo}_forward.fq                                          
                                                                   
  zcat Rawdata/*${archivo}*_2.fq.gz |\                             
    sed '/^@/ s/$/:23:C5A1KACXX:1:1101:1:2 2:N:0:ATCACG/' >\       
    ${archivo}_reverse.fq                                          
                                                                   
  # Ensamblar                                                     
  ./pandaseq -f ${archivo}_forward.fq -r ${archivo}_reverse.fq -B \
    -t .9 -w ensamblados/ensamblado_${archivo}.fa -l 400 -L 500 \  
    -G ensamblados/log_${archivo}.bz2
                                                                   
  # Longitudes secuencias ensambladas                              
  infoseq ensamblados/ensamblado_${archivo}.fa | \                 
    awk -v nom="$archivo" 'NR > 1 {print $3,$6,nom}' >> \          
    ensamblados/metadatos_raw_final.tsv                            
                                                                   
  rm ${archivo}_forward.fq ${archivo}_reverse.fq                   
done                                                               
                                                                   
```

## Agrupar en OTUs

Cambiar nombres:

```bash

for i in $(cat muestras.txt); do 
  nombre=L$i; perl header.fasta.numbers.pl $nombre ensamblado_$nombre.fa
done

```

Agrupar secuencias al 97% de identidad

```bash

qsub -N cdhit_leo -b y -j y -cwd -V \
  "cd-hit-est -c 0.97 -T 20 -M 0 -i leo_tomate_raw_completo.fa -o tomate_raw"

```

Generar tabla de OTUS

```bash

perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g; s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' \
  tomate_raw.clstr > raw.otu

```

Otener secuencias representativas:

```bash

pick_rep_set.py -i raw.otu -f leo_tomate_raw_completo.fa -o rep_set_raw.fa

```

## Filtrado de OTUs

### Singletons

Buscar singletons en los clusters

```bash

awk 'NF < 3' raw.otu | cut -f1 > ids_singletons.txt
wc -l ids_singletons.txt
cut -f1 raw.otu > ids_raw.txt
```

Quitar ids de singletons. En total son 498, 648 OTUs, 
338, 446 son singletons.

```bash
awk 'NR==FNR { values[$0]; next } !($0 in values)' ids_singletons.txt \
ids_raw.txt  > ids_singletons_removed.txt
wc -l ids.txt
```
Quedan 160, 202 OTUs

Generar fasta sin singleton

```bash

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' \
ids_singletons_removed.txt rep_set_raw.fa > singletons_removed.fa

grep ">" -c singletons_removed.fa

```


### Contaminantes

Primer búsqueda de contaminantes contra Greengenes

```bash

parallel_assign_taxonomy_blast.py -i singletons_removed.fa -o not_16S_screen \
-r /qiime/gg_otus-13_8-release/rep_set/70_otus.fasta \
-t /qiime/gg_otus-13_8-release/taxonomy/70_otu_taxonomy.txt \
-X contaminantes_1 -O 20&

```

758 OTUs no dieron hit contra Green Genes

```bash
grep "No blast hit" not_16S_screen/singletons_removed_tax_assignments.txt |\
  cut -f1 > ids_not_16S_screen_1.txt

wc -l ids_not_16S_screen_1.txt
```

Verificar OTUs sin hit contra SILVA

```bash

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' \
ids_not_16S_screen_1.txt singletons_removed.fa > not_16S_screen_1.fa

grep ">" -c not_16S_screen_1.fa

parallel_assign_taxonomy_blast.py -i not_16S_screen_1.fa -o not_16S_screen_2 \
-r /databases/silva/138_2/silva.fna \
-t /home/ltrueba/16S_tomate/silva_tax_leo.txt -O 20 -X contaminantes_2&

grep "No blast hit" -c not_16S_screen_2/not_16S_screen_1_tax_assignments.txt

grep "No blast hit" not_16S_screen_2/not_16S_screen_1_tax_assignments.txt |\
  cut -f1 > ids_not_16S_screen_2.txt

wc -l ids_not_16S_screen_2.txt # quedan 42823 OTUs
```
Solamente 1 OTU no dio contra SILVA

Crear fasta con secuencias sin contaminantes y sin singletons, 
son 160, 201 OTUs

```bash

awk 'NR==FNR { values[$0]; next } !($0 in values)' ids_not_16S_screen_2.txt \
ids_singletons_removed.txt > ids_contaminants_singletons_removed.txt

wc -l ids_contaminants_singletons_removed.txt

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' \
  ids_contaminants_singletons_removed.txt singletons_removed.fa \
  > cont_sing_removed.fna

grep ">" -c cont_sing_removed.fna

```

### Quimeras

Alinear para identificar quimeras

```bash

/home/qiime/bin/parallel_align_seqs_pynast.py -i cont_sing_removed.fna \
  -o chimera_align -O 20 -X chim_align&

/home/qiime/bin/parallel_identify_chimeric_seqs.py -m blast_fragments \
  -i cont_sing_removed.fna -a chimera_align/cont_sing_removed_aligned.fasta \
  -o chimera_list.txt -O 20 -X chimerablast \
  --id_to_taxonomy_fp /qiime/gg_otus-13_8-release/taxonomy/97_otu_taxonomy.txt \
  -r /qiime/gg_otus-13_8-release/rep_set/97_otus.fasta&

awk '{print $1}' chimera_list.txt > ids_chimeras.txt

wc -l ids_chimeras.txt

```

Quitar quimeras, son 117, 568. Quedan 42 633 OTUs

```bash

awk 'NR==FNR { values[$0]; next } !($0 in values)' ids_chimeras.txt \
  ids_contaminants_singletons_removed.txt > ids_cont_sing_chim_removed.txt

wc -l ids_cont_sing_chim_removed.txt

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' \
  ids_cont_sing_chim_removed.txt cont_sing_removed.fna \
  > cont_sing_chim_removed.fna

grep ">" -c cont_sing_chim_removed.fna

```

## Asignación taxonómica

Anotar contra silva

```bash

parallel_assign_taxonomy_blast.py -i cont_sing_chim_removed.fna -o taxonomy \
  -r /databases/silva/138_2/silva.fna \
  -t /home/ltrueba/16S_tomate/silva_tax_leo.txt -O 20 -X taxleo&

```

Limpiar asignación taxonómica, quitando secuencias sin hit, 
secuencias de cloroplasto y mitocondria. En total son 1474 que remover.

```bash

grep -i -c "Chloroplast" taxonomy/cont_sing_chim_removed_tax_assignments.txt
grep -i -c "Mitochondria" taxonomy/cont_sing_chim_removed_tax_assignments.txt
grep "No blast hit" -c taxonomy/cont_sing_chim_removed_tax_assignments.txt


grep "No blast hit" taxonomy/cont_sing_chim_removed_tax_assignments.txt |\
  cut -f1 > ids_remove.txt
grep -i "Chloroplast" taxonomy/cont_sing_chim_removed_tax_assignments.txt |\
  cut -f 1 >> ids_remove.txt
grep -i "Mitochondria" taxonomy/cont_sing_chim_removed_tax_assignments.txt |\
  cut -f 1 >> ids_remove.txt

wc -l ids_remove.txt

```

Tabla de taxonomía limpia filtrada:

```bash

grep -e "No blast hit" -e "Chloroplast" -e "Mitochondria" -v\
    taxonomy/cont_sing_chim_removed_tax_assignments.txt > taxonomy/tax_limpia.txt

wc -l taxonomy/tax_limpia.txt

```

Generar nueva tabla deOTUs sin singletons, contaminantes, cloroplasto y mitocondria:
Quedan 41, 159 OTUs.

```bash

awk 'NR==FNR { values[$0]; next } !($0 in values)' \
  ids_remove.txt ids_cont_sing_chim_removed.txt > ids_final.txt

awk 'NR==FNR { values[$0]; next } $1 in values { print }'\
  ids_final.txt raw.otu > raw_filtrado.otu

wc -l raw_filtrado.otu

```

Formatear tabla de otus
```bash

grep ">" rep_set_raw.fa | tr -d ">" | sed 's/_/\t/g' |\
 cut -d " " -f 2 | cut -f1 | sort | uniq > samples.txt

python3 mod_leo.py samples.txt raw_filtrado.otu

```


Fasta con datos filtrados

```bash

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' \
  ids_final.txt cont_sing_chim_removed.fna > final.fna

grep ">" -c final.fna

```

## Árbol filogenético OTUs filtrados

Preparar para el alineamiento

```bash

ssu-prep final_para_arbol.fna arbol_raw 20 prefix_align.txt -f

```

Alinear 
```bash

bash arbol_raw.ssu-align.sh

```

Limpiar alineamiento

```bash

ssu-mask --stk2afa arbol_raw

perl -i.bak -pe 's/\./-/g;s/U/T/gi' arbol_raw/arbol_raw.bacteria.afa

sed -i 's/>OTT/>OTU/' arbol_raw/arbol_raw.bacteria.afa

```

Calcular árbol

```bash

fasttreeMP -nt -fastest arbol_raw/arbol_raw.bacteria.afa > tomate_raw.nwk

```
