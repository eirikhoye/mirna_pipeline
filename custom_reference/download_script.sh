while read p
do
  echo "Download Mature and Star ref from MirGeneDB.org for $p"
  wget -qO- https://mirgenedb.org/fasta/$p?mat=1 > $p.mature_star_ref.fas
  wget -qO- https://mirgenedb.org/fasta/$p?star=1 >> $p.mature_star_ref.fas
done < species_id.txt
