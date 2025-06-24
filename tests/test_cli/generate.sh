for i in $( ls `pwd`/../test_topologies | grep top |grep -v lambda)
do
	filename=$(basename $i .top) 
       	python3 ../../src/repex_topology_parser.py --verbose -p `pwd`/../test_topologies/$i -O topol_rest2_$filename -P ./ -m rest2 -H 0 -n 20 $(if [[ "$i" == *sys1* ]]; then echo "-T 300 -M 450 "; else echo "-T 300 -M 500"; fi) 
done

for i in $( ls `pwd`/../test_topologies | grep top |grep -v lambda)
do
	filename=$(basename $i .top) ; python3 ../../src/repex_topology_parser.py --verbose -p `pwd`/../test_topologies/$i -O topol_ssrest3_$filename -P ./ -m ssrest3 -H 0 -n 20 $(if [[ "$i" == *sys1* ]]; then echo "-T 300 -M 450 "; else echo "-T 300 -M 500"; fi) -k 1.1 --kappa-atomtypes OW_tip4pd 
done
