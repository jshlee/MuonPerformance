mkdir -p jsonconf/20170805/roc_various_pt

python evalscr/20170805/generate_template_various_pt.py && python bendingJSON.py evalscr/20170805/template_roc_various_pt.json

if [ "$?" = "0" ]; then
    for name in `ls jsonconf/20170805/roc_various_pt/*_00_24_* jsonconf/20170805/roc_various_pt/*_24_28_*` ; do
      while true; do 
        if [ "`ps -ef| grep $PSID | grep pytho[n]| wc -l`" -lt "20" ]; then 
          break
        fi
        
        sleep 3
      done
      
      python roccurve_withjson.py $name & 
    done

    while true; do if [ "`ps -ef| grep $PSID | grep pytho[n]| wc -l`" = "0" ]; then break; fi; done
fi


