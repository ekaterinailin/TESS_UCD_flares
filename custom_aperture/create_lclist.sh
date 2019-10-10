END=194
for i in $(seq 1 $END);
    do
        ls $i*tess*.fits -t | head -1 >> lclist.txt
    done

while read path
    do
        echo python $path >> ../cluster_scripts/tessjobs.sh
    done < lclist.txt

rm lclist.txt