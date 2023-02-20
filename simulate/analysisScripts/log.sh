i=0;for seed in 333 666 999; do for file in data/*.xml;do if [[ $i -lt 1000 ]] && [[ ! -f $seed/$(basename $file | sed "s/.xml//").log ]]; then sbatch -c 1 -p htc -t 0-2 -o ologs/flipflopBeast_$(basename ${file})_${seed}-%A_%a.out ./flipflop.sh $seed $file;i=$(( $i+1 ));fi;done;done

i=0;for seed in 333 666 999; do for file in data/*.xml;do if [[ $i -lt 1000 ]] && [[ ! -f $seed/$(basename $file | sed "s/.xml//").log ]]; then sbatch -c 1 -t 1-0 -o ologs/flipflopBeast_$(basename ${file})_${seed}-%A_%a.out ./flipflop.sh $seed $file;i=$(( $i+1 ));fi;done;done

for n in 333 666 999; do for i in ${n}/*.log; do nlines=$(wc -l $i | sed "s/ .*$//");echo $nlines;done | sort | uniq -c > ${n}.lines;done

runthething() { Rscript ../scripts/flipFlopRWTY.R 0.1 1 analysis/ $1 666/$(basename $1) 999/$(basename $1); }; export -f runthething; parallel -j 8 runthething ::: 333/*.trees
