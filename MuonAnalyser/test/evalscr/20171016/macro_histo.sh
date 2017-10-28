for name in `cat listroclist.txt` ; do python evalscr/20171016/makerochistos.py $name Sig & python evalscr/20171016/makerochistos.py $name Bkg & done
while true ; do if [ "`ps -ef | grep $UID | grep pytho[n] | grep makeroc | wc -l`" = "0" ] ; then break ; fi ; sleep 3 ; done
