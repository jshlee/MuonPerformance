python evalscr/20170405/generaterates_newenv.py

for name in `ls jsonconf/20170405/effbkgplots_*/*/*/*/*.json`; do 
    while true; do 
        if [ "`ps -ef| grep 8144| grep pytho[n]| wc -l`" -lt "20" ]; then
            break
        fi
        
        sleep 3
    done
    
    python ../python/universaldrawer.py $name & 
done

while true; do if [ "`ps -ef| grep 8144| grep pytho[n]| wc -l`" = "0" ]; then break; fi; done


