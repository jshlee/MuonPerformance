for name in `cat listroclist_run2roc.txt` ; do python evalscr/20171026/roccurvenew_comp_p8.py $name & done
while true ; do if [ "`ps -ef | grep $UID | grep pytho[n] | grep roccurve | wc -l`" = "0" ] ; then break ; fi ; sleep 3 ; done
