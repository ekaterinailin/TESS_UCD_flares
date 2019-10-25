END=194
for i in $(seq 0 $END);
    do
        ls "$i"ca_tess*.fits -t | head -n 1 >> lclist.txt
    done

while read path
    do
        echo python /home/eilin/TESSUCDs/src/single_tess_injrec.py $path >> ../cluster_scripts/tessjobs.sh
    done < lclist.txt

rm lclist.txt
