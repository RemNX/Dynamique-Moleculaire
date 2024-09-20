N=1000

for opti in -O0 -O1 -O2 -O3 -ffast-math
do
    total=0
    for (( i=1; i<=$N; i++))
    do
        start=$(date +%s%N)
        cc $opti Particules_mieux_ppar.c -lm && ./a.out
        end=$(date +%s%N)

        total=$(( total + end - start ))
    done
    total_s=$(echo "scale=6; $total / $N / 1000000000" | bc -l)

    echo $opti : $total_s
done
