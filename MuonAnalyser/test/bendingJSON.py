import sys, json, copy


strBendHead    = "###bend"
strBendMul     = "###multibend"
strBendMulList = "###mullist"
strBendMulVal  = "###mulval"
strBendSave    = "###savedic"


def myBendJSON(dicIn, dicStack, dicRead, dicSaveInfo):
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
          
          myBendJSON(dicInSub, dicStackSub, dicReadSub, dicSaveInfo)
        
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
    
    for itemIn in dicStack.keys(): 
      if type(dicStack[ itemIn ]) == dict and strBendSave in dicStack[ itemIn ]: 
        strItemSave = itemIn
        break
    
    if strItemSave != "": 
      dicStack[ strItemSave ] = dicStack[ strItemSave ][ strBendSave ]%dicRead
    
    json.dump(dicStack, open(dicSaveInfo[ "outputfile" ]%dicRead, "w"), indent=2)


myBendJSON(json.load(open(sys.argv[ 1 ], "r")), {}, {}, json.loads(sys.argv[ 2 ]))


