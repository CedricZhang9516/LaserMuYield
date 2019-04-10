
i=-15

while [ "$i" != "20" ]
do
    p=$(( 1000000 * $i ))

	/Users/zhangce/WorkArea/LaserMuYield/LaserSim Detail_output_$p.root $p
    
    i=$(($i+5))
done
