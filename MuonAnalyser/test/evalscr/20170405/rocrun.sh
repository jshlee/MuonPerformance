namesuf="moreeta"
echo -n "" > isovaluecutlist_$namesuf.txt

for id in `echo "Tight Loose LooseMod TightModNoIP"`
do for pu in `echo "PU0 PU140 PU200"`
do for type in `echo "1 2"`
    do python evalscr/20170405/roccurve_moreeta.py $pu $type $id >> isovaluecutlist_$namesuf.txt & 
done
done
done

while true; do if [ "`ps -ef| grep 8144| grep pytho[n]| wc -l`" = "0" ]; then break; fi; done
