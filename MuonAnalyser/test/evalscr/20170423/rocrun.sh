date=20170423

python evalscr/$date/generate_roccuves.py

for name in `ls jsonconf/$date/roccurves_*/*/*.json`; do
  while true; do 
    if [ "`ps -ef| grep 8144| grep pytho[n]| wc -l`" -lt "20" ]; then 
      break
    fi
    
    sleep 3
  done
  
  python roccurve_withjson.py $name & 
done

while true; do if [ "`ps -ef| grep 8144| grep pytho[n]| wc -l`" = "0" ]; then break; fi; done


