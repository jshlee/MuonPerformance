echo -n "" > isovaluecutlist_ME0.txt

for id in `echo "Loose"`
do for pu in `echo "PU200"`
do for type in `echo "1 2"`
    do python evalscr/20170417/roccurve_me0.py $pu $type $id >> isovaluecutlist_ME0.txt & 
done
done
done

while true; do if [ "`ps -ef| grep 8144| grep pytho[n]| wc -l`" = "0" ]; then break; fi; done
