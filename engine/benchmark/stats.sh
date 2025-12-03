nums=$(grep "load" log.txt | grep -o '\[[0-9]\+ ms\]' | sed -E 's/\[([0-9]+) ms\]/\1/')

echo "+++ Load +++"
echo "Mode:   $(echo "$nums" | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')"
echo "Min:    $(echo "$nums" | sort -n | head -n1)"
echo "Max:    $(echo "$nums" | sort -n | tail -n1)"
echo "Median: $(echo "$nums" | sort -n | awk '{a[NR]=$1} END {if (NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2}}')"
echo "Mean:   $(echo "$nums" | awk '{sum+=$1} END {if(NR>0) print sum/NR}')"


nums=$(grep "filter" log.txt | grep -o '\[[0-9]\+ ms\]' | sed -E 's/\[([0-9]+) ms\]/\1/')

echo "+++ Filter +++"
echo "Mode:   $(echo "$nums" | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')"
echo "Min:    $(echo "$nums" | sort -n | head -n1)"
echo "Max:    $(echo "$nums" | sort -n | tail -n1)"
echo "Median: $(echo "$nums" | sort -n | awk '{a[NR]=$1} END {if (NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2}}')"
echo "Mean:   $(echo "$nums" | awk '{sum+=$1} END {if(NR>0) print sum/NR}')"


nums=$(grep "mined" log.txt | grep -o '\[[0-9]\+ ms\]' | sed -E 's/\[([0-9]+) ms\]/\1/')

echo "+++ Mine +++"
echo "Mode:   $(echo "$nums" | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')"
echo "Min:    $(echo "$nums" | sort -n | head -n1)"
echo "Max:    $(echo "$nums" | sort -n | tail -n1)"
echo "Median: $(echo "$nums" | sort -n | awk '{a[NR]=$1} END {if (NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2}}')"
echo "Mean:   $(echo "$nums" | awk '{sum+=$1} END {if(NR>0) print sum/NR}')"
echo "Stddev: $(echo "$nums" | awk '{x[NR]=$1; s+=$1} END {if (NR>1){m=s/NR; for(i=1;i<=NR;i++) ss+=($1=x[i]-m)^2; print sqrt(ss/(NR-1))}}')"


nums=$(grep "memcpy" log.txt | grep -o '\[[0-9]\+ ms\]' | sed -E 's/\[([0-9]+) ms\]/\1/')

echo "+++ Transfer +++"
echo "Mode:   $(echo "$nums" | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')"
echo "Min:    $(echo "$nums" | sort -n | head -n1)"
echo "Max:    $(echo "$nums" | sort -n | tail -n1)"
echo "Median: $(echo "$nums" | sort -n | awk '{a[NR]=$1} END {if (NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2}}')"
echo "Mean:   $(echo "$nums" | awk '{sum+=$1} END {if(NR>0) print sum/NR}')"


nums=$(grep "wrote" log.txt | grep -o '\[[0-9]\+ ms\]' | sed -E 's/\[([0-9]+) ms\]/\1/')

echo "+++ Store +++"
echo "Mode:   $(echo "$nums" | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')"
echo "Min:    $(echo "$nums" | sort -n | head -n1)"
echo "Max:    $(echo "$nums" | sort -n | tail -n1)"
echo "Median: $(echo "$nums" | sort -n | awk '{a[NR]=$1} END {if (NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2}}')"
echo "Mean:   $(echo "$nums" | awk '{sum+=$1} END {if(NR>0) print sum/NR}')"


















