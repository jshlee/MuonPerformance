date=20170805

for name in `ls jsonconf/20170805/roc_various_pt/*` ; do
  while true; do 
    if [ "`ps -ef| grep $PSID | grep pytho[n]| wc -l`" -lt "20" ]; then 
      break
    fi
    
    sleep 3
  done
  
  python roccurve_withjson.py $name & 
done

while true; do if [ "`ps -ef| grep $PSID | grep pytho[n]| wc -l`" = "0" ]; then break; fi; done


