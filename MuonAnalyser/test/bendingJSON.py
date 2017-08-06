import sys, json, copy


strBendHead    = "###bend"
strBendMul     = "###multibend"
strBendMulList = "###mullist"
strBendMulVal  = "###mulval"
strBendSave    = "###savedic"
strBendJSONInf = "###jsoninfo"


def myCreateDir(strPath):
  strDirCurr = ""
  listPartPath = strPath.split("/")
  
  for strPart in listPartPath[ 0:-1 ]: 
    strDirCurr = strDirCurr + strPart + "/"
    if not os.path.isdir(strDirCurr): os.mkdir(strDirCurr)


def myBendJSON(dicIn, dicStack, dicRead, dicSaveInfo, listOutput):
  nIsSplited = 0
  
  for itemIn in dicIn.keys():
    valItem = dicIn[ itemIn ]
    
    if type(dicIn[ itemIn ]) == dict and strBendHead in dicIn[ itemIn ]: 
      strID = dicIn[ itemIn ][ strBendHead ]
      
      if strID in dicRead:
        valItem = dicIn[ itemIn ][ dicRead[ strID ] ]
      else: 
        for strCase in dicIn[ itemIn ].keys():
          if strCase == strBendHead: continue
          
          dicInSub = copy.deepcopy(dicIn)
          dicStackSub = copy.deepcopy(dicStack)
          dicReadSub = copy.deepcopy(dicRead)
          
          dicInSub.pop(itemIn)
          dicStackSub[ itemIn ] = dicIn[ itemIn ][ strCase ]
          dicReadSub[ strID ] = strCase
          
          myBendJSON(dicInSub, dicStackSub, dicReadSub, dicSaveInfo, listOutput)
        
        nIsSplited = 1
        break
    elif type(dicIn[ itemIn ]) == dict and strBendMul in dicIn[ itemIn ]: 
      nIsAllNeed = 1
      
      for itemNeed in dicIn[ itemIn ][ strBendMul ]: 
        if itemNeed not in dicRead: 
          nIsAllNeed = 0
          break
      
      if nIsAllNeed == 0: 
        continue
      
      for dicMulVal in dicIn[ itemIn ][ strBendMulList ]: 
        nIsLookFor = 1
        
        for itemNeed in dicIn[ itemIn ][ strBendMul ]: 
          if dicMulVal[ itemNeed ] != dicRead[ itemNeed ]: 
            nIsLookFor = 0
            break
        
        if nIsLookFor == 1: 
          valItem = dicMulVal[ strBendMulVal ]
          break
        
    #elif type(dicIn[ itemIn ]) == dict
    
    dicStack[ itemIn ] = valItem
    dicIn.pop(itemIn)
  
  if nIsSplited == 0: 
    strItemSave = ""
    strOutputFile = ""
    
    for itemIn in dicStack.keys(): 
      if type(dicStack[ itemIn ]) == dict and strBendSave in dicStack[ itemIn ]: 
        strItemSave = itemIn
        #break
      if itemIn == strBendJSONInf and "jsonname" in dicStack[ strBendJSONInf ]: 
        strOutputFile = dicStack[ strBendJSONInf ][ "jsonname" ]%dicRead
    
    if strItemSave != "": 
      dicStack[ strItemSave ] = dicStack[ strItemSave ][ strBendSave ]%dicRead
    
    if strOutputFile == "":
      strOutputFile = dicSaveInfo[ "outputfile" ]%dicRead
    
    if strBendJSONInf in dicStack: dicStack.pop(strBendJSONInf)
    
    
    json.dump(dicStack, open(strOutputFile, "w"), indent=2)
    listOutput.append(strOutputFile)


dicBendJSON = json.load(open(sys.argv[ 1 ], "r"))
dicOutputInfo = {}

if len(sys.argv) >= 3: 
  dicOutputInfo = json.loads(sys.argv[ 2 ])

listResJSON = []
myBendJSON(dicBendJSON, {}, {}, dicOutputInfo, listResJSON)

if strBendJSONInf in dicBendJSON and "cmd" in dicBendJSON[ strBendJSONInf ]: 
  strMulCmd = ""

  for strJSON in listResJSON: 
    #strMulCmd += "python ../python/universaldrawer.py " + strJSON + " & "
    strMulCmd += dicBendJSON[ strBendJSONInf ][ "cmd" ]%strJSON + " & \n"

  strMulCmd += "while true ; do if [ \"`ps -ef | grep $PSID | grep pytho[n] | wc -l`\" = \"0\" ] ; then break ; fi ; sleep 0.5 ; done"
  print strMulCmd


