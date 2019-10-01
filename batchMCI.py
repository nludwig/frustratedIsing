##existing bugs/issues
#if launch an OS with long name, say
#job-kT10-11- ...
#and subsequently an OS with shorter
#name that would grep long OS name, say
#job-kT10
#will get an incorrect dependency of
#short name on long name. Since I don't
#plan to use this script in that manner,
#I won't fix it for now

#functions and imports
import sys
import subprocess as sp
import shlex
import os
from itertools import product
from rngSeedTable import RngSeedTable

def getArgs(nargs,errStr):
  err=sys.stderr #alias
  print('nargs:{}\tcmdline args:{}\targs:'.format(nargs,len(sys.argv)),file=err)
  for i in range(0,len(sys.argv)):
    print('{}'.format(sys.argv[i]),end=' ',file=err)
  print('',file=err)
  if len(sys.argv) != nargs:
    print('usages: python3 {}'.format(errStr), file=err)
    exit(1)
  arg=iter(sys.argv[1:])
  return iter(sys.argv[1:])

def grep(string,lines):
  lstring=len(string)
  for line in lines:
    same=0
    for lett in line:
      if lett==string[same]: same+=1
      else:                  same=0
      if same==lstring:
        yield line
        break

def getNandNsu(f):
  if type(f)==str: f=open(f,'r')
  lines=[line.strip() for line in f]
  f.close()
  N=-1; Nsu=-1
  for line in lines:
    a=line.split()
    if a[0]=='N':     N=a[1]
    elif a[0]=='Nsu': Nsu=a[1]
    if N!=-1 and Nsu!=-1:      break
  return N,Nsu

def totalWorkers(var,varkeys):
  nWorkers=1
  for key in varkeys:
    nWorkers*=len(var[key])
  return nWorkers

def totalOverseers(var,varkeys,maxCpuPerNode):
  nWorkers=totalWorkers(var,varkeys)
  nOverseers=nWorkers//maxCpuPerNode
  if nWorkers%maxCpuPerNode!=0:
    nOverseers+=1
  return nOverseers

def uniqueValuesByIndexAmongstTupleList(tl):
  l=len(tl)
  ll=len(tl[0]) #assume all same length
  uniques=[[] for i in range(0,ll)]
  for i in range(0,l):
    assert len(tl[i])==ll
    for j in range(0,ll):
      new=True
      for k in range(0,len(uniques[j])):
        if tl[i][j]==uniques[j][k]:
          new=False
          break
      if new is True:
        uniques[j].append(tl[i][j])
  return uniques

def nUniqueValuesByIndexAmongstTupleList(tl):
  uniques=uniqueValuesByIndexAmongstTupleList(tl)
  lUniques=[len(i) for i in uniques]
  return lUniques

#assumes sorted
def unique(l):
  ll=len(l)
  luni=[-1]*ll
  ln=[0]*ll
  i=0;j=1;k=0
  while j<ll:
    if l[i]==l[j]:
      j+=1
    else:
      luni[k]=l[i]
      ln[k]=j-i
      i=j
      j+=1
      k+=1
  if l[i]==l[-1]:
    luni[k]=l[i]
    ln[k]=j-i
    k+=1
  return luni[:k],ln[:k]

#format of lines returned from '/home/nickludwig/bashScripts/sqn.sh':
#[jobN,serverPfix,userNm,jobNm,userNm,status,runT,nodes,nodeNm]
def findRecentDep(cmd,grepstrings,col=0,disqual=None,ret='N'):
  if type(grepstrings)!=list:
    try: grepstrings=[grepstrings]
    except: exit(1)
  cmd=shlex.split(cmd)
  try: out=sp.check_output(cmd).decode('UTF-8').split(sep='\n')
  except sp.CalledProcessError:
    print('findRecentDep: error calling bash fn {}'.format(cmd))
    exit(1)
    #return -1,''
  for gs in grepstrings: out=grep(gs,out) #grep each gs in sequence
  out=list(out)
  lout=len(out)
  if lout==0:
    print('findRecentDep: grep returns nothing from lines from fn {}'.format(cmd))
    return -1,''   #no such jobs running
  out=[i.split() for i in out]
  i=0
  while i<lout:
    if out[i][col][0]!=grepstrings[0][0]:
      del out[i]
      lout-=1
    else: i+=1

  if lout==0:
    print('findRecentDep: grep2 returns nothing from lines from fn {}'.format(cmd))
    return -1,''   #no such jobs running
      
  steps=[0]*lout
  for j in range(0,lout):
    i=-1
    if disqual is not None:
      foundDisqual=False
      for k in range(0,len(out[j][col])):
        if out[j][col][k]==disqual:
          foundDisqual=True
          break
      if foundDisqual:
        steps[j]=-1
        continue    #don't consider entries w/ disqual
    while out[j][col][i]!='-': i-=1
    k=0
    while(True):
      if k!=0:
        try: steps[j]=int(out[j][col][i+1:k])
        except ValueError: k-=1
        else: break
      else:
        try: steps[j]=int(out[j][col][i+1:])
        except ValueError: k-=1
        else: break
  stepm=max(steps)
  if len(steps)>1:
    allSame=True
    for i in range(0,len(steps)):
      if stepm!=steps[i]:
        allSame=False
        break
    if allSame is True:
      return -1,''
  i=0
  while(steps[i]!=stepm): i+=1
  if type(ret)==int:
    return out[i][ret],out[i][col]    #jobN,jobNm
  elif ret.lower()=='n':
    return stepm,out[i][col]

class JobLauncher:
  def __init__(self,N,Nsu,name,var,varkeys,paraTemplate,umbrK,rotation,
                    lStr,sigma,ewaldIn,enFreq,dumpFreq,dataFreq,
                    suComFreq,umbrFreq,eSteps,pSteps,neutralOverall):
    self.N=N
    self.Nsu=Nsu
    self.name=name
    self.var=var
    self.varkeys=varkeys
    self.paraTemplate=paraTemplate
    self.umbrK=umbrK
    self.rotation=rotation
    self.lStr=lStr
    self.sigma=sigma
    self.ewaldIn=ewaldIn
    self.enFreq=enFreq
    self.dumpFreq=dumpFreq
    self.dataFreq=dataFreq
    self.suComFreq=suComFreq
    self.umbrFreq=umbrFreq
    self.eSteps=eSteps
    self.pSteps=pSteps
    self.neutralOverall=neutralOverall

  #N,Nsu,var,varkeys
  def makeWorkersOutPfx(self):
    nWorkers=totalWorkers(self.var,self.varkeys)
    outprefix=['' for i in range(0,nWorkers)]
    string='N{}-Nsu{}'
    for key in self.varkeys:
      string+='-'+key+'{}'
            
    itervar=[iter(self.var[key]) for key in self.varkeys]
    #for i,tup in enumerate(product(*itervar)):
    #  outprefix[i]=string.format(self.N,self.Nsu,*tup)

    tuples=list(product(*itervar))
    tuni,tn=unique(tuples)
    luni=len(tuni)
    print('tuples:{}\ntuni:{}\ntn:{}'.format(tuples,tuni,tn))
    ct=0
    for i in range(0,luni):
      for j in range(0,tn[i]):
        outprefix[ct]=string.format(self.N,self.Nsu,*(tuni[i]))+'-v{}'.format(j)
        ct+=1
    return outprefix

    #itervari=[iter(self.var[key]) for key in self.varkeys]
    #itervarj=[iter(self.var[key]) for key in self.varkeys]
    #dupes=[] #multiple trjs w/ same paras?
    #for i,tupi in enumerate(product(*itervari)):
    #  for j,tupj in enumerate(product(*itervarj)):
    #    if i>=j:
    #      continue
    #    if tupi==tupj:
    #      
    #      dupes.append([j,i])
    #ldupes=len(dupes)
    #if ldupes==0:
    #  itervar=[iter(self.var[key]) for key in self.varkeys]
    #  for i,tup in enumerate(product(*itervar)):
    #    outprefix[i]=string.format(self.N,self.Nsu,*tup)
    #else:
    #  dupeNum=[[d[0],0] for d in dupes]
    #  for i in range(1,ldupes):
    #    for j in range(0,i):
    #      if dupes[i][1]==dupes[j][1]:

  #N,Nsu,varkeys
  def makeOverseersOutPfx(self,maxCpuPerNode):
    nWorkers=totalWorkers(self.var,self.varkeys)
    osWorkers,_=self.assignWorkersToOverseers(maxCpuPerNode)
    nOverseers=len(osWorkers)
    outprefix=['' for i in range(0,nOverseers)]
    for i in range(0,nOverseers):
      string='OS-N{}-Nsu{}-'
      nWorkersi=len(osWorkers[i])
      unique=uniqueValuesByIndexAmongstTupleList(osWorkers[i])
      print(unique)
      for item in unique:
        item.sort()
      for j,key in enumerate(self.varkeys):
        string+=key
        for k in range(0,len(unique[j])):
          string+='{}-'.format(unique[j][k])
      string=string[:-1]
      outprefix[i]=string.format(self.N,self.Nsu)
    return outprefix

  def assignWorkersToOverseers(self,maxCpuPerNode,wpfix=None):
    nWorkers=totalWorkers(self.var,self.varkeys)
    nOverseers=totalOverseers(self.var,self.varkeys,maxCpuPerNode)
    itervar=[iter(self.var[key]) for key in self.varkeys]
    tlist=[]
    osWorkers=[]
    wpfixByOS=[]
    for i,tup in enumerate(product(*itervar)):
      tlist.append(tup)
      if (i+1)%maxCpuPerNode==0:
        osWorkers.append(tlist)
        tlist=[]
        if wpfix is not None:
          ind=(i+1)//maxCpuPerNode-1
          wpfixByOS.append(wpfix[ind:ind+maxCpuPerNode])
          print('i:{}\tind:{}\tind+maxCpu:{}\twpfixByOS:{}'.format(i,ind,ind+maxCpuPerNode,wpfixByOS))
    if nWorkers%maxCpuPerNode!=0: #add stragglers to another, non-full node
      osWorkers.append(tlist)
      tlist=[]
      if wpfix is not None:
        ind=(i+1)//maxCpuPerNode-1
        wpfixByOS.append(wpfix[ind+maxCpuPerNode:])
        print('post loop i:{}\tind:{}\tind+maxCpu:{}\twpfixByOS:{}'.format(i,ind,ind+maxCpuPerNode,wpfixByOS))
    if len(osWorkers)!=nOverseers:
      print('assignWorkersToOverseers: losWorkers!=nOverseers: {}!={}'.format(len(osWorkers),nOverseers),file=sys.stderr)
    if wpfix is not None:
      return osWorkers,wpfixByOS
    else:
      return osWorkers,None

  #var,varkeys,paraTemplate,name,umbrK,umbrR
  def makeParas(self,umbrR):
    #make .para file from template
    #read in
    if type(self.paraTemplate)==str:
      self.paraTemplate=open(self.paraTemplate,'r')
    lines=[line.strip() for line in self.paraTemplate]
    self.paraTemplate.close()
    nlines=len(lines)
 
    nkeys=len(self.varkeys)
    lineN=[0]*nkeys
  
    for i in range(0,nlines):
      for j in range(0,nkeys):
        try:
          if lines[i].split()[0]==self.varkeys[j]:
            lineN[j]=i
            break
        except IndexError:
          continue

    iter_opfix=iter(self.makeWorkersOutPfx())
  
    itervar=[0]*nkeys
    for i,key in enumerate(self.varkeys):
      itervar[i]=iter(self.var[key])
  
    replLine=[0]*nkeys
    for i,tup in enumerate(product(*itervar)):
      for j in range(0,nkeys):
        if self.varkeys[j]!='L':
          replLine[j]='{} {}'.format(lines[lineN[j]].split()[0],tup[j])
        else:
          d=-1
          for line in lines:
            if line.split()[0]=='d':
              d=int(lines.split()[1])
          if d!=-1:
            replLine[j]='{}'.format(lines[lineN[j]].split()[0])
            for i in range(0,d): replLine[j]+=' {}'.format(tup[j])
          else:
            replLine[j]='{} {}'.format(lines[lineN[j]].split()[0],tup[j])
      o='{}.{}'.format(next(iter_opfix),'para')
      if self.name is not None:  o='{}-{}'.format(self.name,o)
      if umbrR is not None: o='{}-{}'.format(umbrR,o)
  
      linesout=['' for i in range(0,nlines)]
      for i in range(0,nlines):
        l=lines[i]
        for j in range(0,nkeys):
          if i==lineN[j]:
            l=replLine[j]
            break
        linesout[i]=l
      if umbrR is not None:
        #replace K_VAL
        K_VALUElineNs=[]
        for i,line in enumerate(linesout):
          if len(line.split(sep='K_VALUE'))>1:
            K_VALUElineNs.append(i)
        K_VALUElines=['' for i in range(0,len(K_VALUElineNs))]
        for i,j in enumerate(K_VALUElineNs):
          tmp=linesout[j].split(sep='K_VALUE')
          if len(tmp)==2:
            K_VALUElines[i]='{}{}{}'.format(tmp[0],self.umbrK,tmp[1])
          else:
            print('got more than one K_VALUE in line: {}'.format(linesout[j]),file=sys.stderr)
            exit(1)
        #replace R0_VAL
        for i in range(0,len(K_VALUElines)):
          tmp=K_VALUElines[i].split(sep='R0_VALUE')
          if len(tmp)==2:
            K_VALUElines[i]='{}{}'.format(tmp[0],umbrR)
          else:
            print('got more than one R0_VALUE in line: {}'.format(K_VALUElines[i]),file=sys.stderr)
            exit(1)
        for i,j in enumerate(K_VALUElineNs):
          linesout[j]=K_VALUElines[i]
      with open(o,'w') as f:
        for l in linesout:
          print('{}'.format(l),file=f)

  #lines,jobN,pfixi,rnList,execInParentDir,rotation,data,N,lStr,sigma,ewaldIn,sweepsMult,kTeff,enFreq,dumpFreq,dataFreq,suComFreq,umbrFreq,eSteps,pSteps,neutralOverall
  def makeSbatch(self,lines,jobN,ospfixi,rnList,execInParentDir,data,wpfixi,sweepsMult,kTeff,umbrR=None):
    ewaldIn='../'+self.ewaldIn if execInParentDir is True else self.ewaldIn
    #prep
    rnStr=''
    for rn in rnList:
      rnStr+='{},'.format(rn)
    rnStr=rnStr.rstrip(',')
  
    oldpfix=ospfixi+'-'+str(jobN-1)
    nextpfix=ospfixi+'-'+str(jobN)
    wpfixi=[self.name+'-'+w for w in wpfixi]
    oldwpfix=[w+'-'+str(jobN-1) for w in wpfixi]
    nextwpfix=[w+'-'+str(jobN) for w in wpfixi]
  
    #modify
    NAMElineNs=[]
    for i,line in enumerate(lines):
      if len(line.split(sep='NAME'))>1:
        NAMElineNs.append(i)
    NAMElines=['' for i in range(0,len(NAMElineNs))]
    for i,j in enumerate(NAMElineNs):
      tmp=lines[j].split(sep='NAME')
      if len(tmp)==2:
        NAMElines[i]='{}{}{}'.format(tmp[0],nextpfix,tmp[1])
      else:
        print('got more than one NAME in line: {}'.format(lines[j]),file=sys.stderr)
        exit(1)
  
    #prep workers in to para, data, trj
    paraVal=''
    for i in range(0,len(wpfixi)):
      if jobN==0:
        paraVal+=wpfixi[i]+'.para,' #hack
      else:
        paraVal+=oldwpfix[i]+'.para,'
    paraVal=paraVal.rstrip(',')
    if umbrR is not None:
      paraVal='{}-{}'.format(umbrR,paraVal)
    if data is not None:
      dataVal=data
    else:
      dataVal=''
      for w in oldwpfix:
          dataVal+=w+'.data,'
      dataVal=dataVal.rstrip(',')
    lammpsVal=''
    for i,w in enumerate(nextwpfix):
      lammpsVal+=w+'.fitrj,'
    lammpsVal=lammpsVal.rstrip(',')
  
    o='{}.{}'.format(nextpfix,'sbatch')
    #output
    nlines=len(lines)
    with open(o,'w') as f:
      for i in range(0,nlines):
        n=0
        for j in range(0,len(NAMElineNs)):
          if i!=NAMElineNs[j]:
            n+=1
          else:
            l=NAMElines[j]
            del NAMElines[j]
            del NAMElineNs[j]
            n=-1
            break
        if n==len(NAMElineNs):
          l=lines[i]
        print('{}'.format(l),file=f)
  
      cmd='./MCIsingFIv28_icc_sa 1 0 {} {} {} {} {:.4f} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(self.neutralOverall,self.rotation,self.N,self.lStr,self.sigma,self.ewaldIn,paraVal,dataVal,sweepsMult,kTeff,self.enFreq,self.dumpFreq,self.dataFreq,self.suComFreq,self.umbrFreq,self.eSteps,self.pSteps,rnStr,lammpsVal)
      if execInParentDir is True:
        cmd='.'+cmd
      print('{}'.format(cmd),file=f)
    return o

  def submitBatchLoop(self,nWorkers,umbrR,jobNStart,jobIterations,jobMax,dataIn,lines,seeds,execInParentDir,sweepsMultByOS,kTeffByOS,maxCpuPerNode):
    wpfix=self.makeWorkersOutPfx()
    ospfix=self.makeOverseersOutPfx(maxCpuPerNode)
    osWorkers,wpfixByOS=self.assignWorkersToOverseers(maxCpuPerNode,wpfix=wpfix)
    #print('wpfix: {}\nosWorkers: {}\nwpfixByOS: {}'.format(wpfix,osWorkers,wpfixByOS))
    nOverseers=totalOverseers(self.var,self.varkeys,maxCpuPerNode)
    print('nOverseers: {}\tospfix: {}'.format(nOverseers,ospfix),file=sys.stderr)
    for i in range(0,nOverseers):
      ospfixi=ospfix[i]
      print(ospfixi)
      if self.name is not None:
        ospfixi='{}-{}'.format(self.name,ospfixi)
      if umbrR is not None:
        ospfixi='{}-{}'.format(umbrR,ospfixi)

      N,Nm=findRecentDep('/home/nickludwig/bashScripts/sqn.sh',[ospfixi],col=3,ret=0)
      if N==-1:
        for j in range(0,len(osWorkers[i])):
          Nj=[0]*len(osWorkers[i])
          jobNStart,Nm=findRecentDep('ls',[osWorkers[i][j],'.data'],disqual='_')
          Nm=Nm[:-5]
        for j in range(0,len(osWorkers[i])):
          for k in range(j+1,len(osWorkers[i])):
            if Nj[j]!=Nj[k]:
              print('submitBatchLoop: not all .data files have same number',file=sys.stderr)
              exit(1)
        jobNStart=Nj[0]
      else:
        ind=-1
        while Nm[ind]!='-':
          ind-=1
        jobNStart=int(Nm[ind+1:])
      if jobNStart!=0:
        jobNStart+=1
      if N!=-1:
        slurmJNum=N
      else:
        slurmJNum=None

      for jobN in range(jobNStart,jobNStart+jobIterations):
        print('jobN:{}\tjobNStart:{}\tjobMax:{}\tdataIn:{}\tNm:{}'.format(jobN,jobNStart,jobMax,dataIn,Nm))
        if jobN>=jobMax:
          print('jobN>=jobMax; skipping')
          continue
  
        makeSbatchArgList=lambda d: [lines,jobN,ospfixi,seeds[i][jobN-jobNStart],\
                execInParentDir,d,wpfixByOS[i],sweepsMultByOS[i],kTeffByOS[i],umbrR]

        print('in loop, handing makeSbatch wpfixByOSi:{}'.format(wpfixByOS[i]))
        if jobN!=jobNStart:
          o=self.makeSbatch(*makeSbatchArgList(None))
        else:
          if dataIn is not None:
            o=self.makeSbatch(*makeSbatchArgList(dataIn))
          #elif Nm!='': #if want this for multithread, need automated scheme to pull .data's
          #  o=makeSbatch(*makeSbatchArgList(Nm+'.data'))
          else:
            o=self.makeSbatch(*makeSbatchArgList('none'))
  
        if slurmJNum is not None:
          sbatchcall='sbatch --dependency=afterok:{} {}'.format(slurmJNum,o)
        else:
          sbatchcall='sbatch {}'.format(o)
  
        print('submitting: {}'.format(sbatchcall))
        sbatchcall=shlex.split(sbatchcall)
        #sp.call(sbatchcall)
        slurmJNum=sp.check_output(sbatchcall).split()
        slurmJNum=int(slurmJNum[-1])    #slurm output: "blah blah <job #>"

def main():
  err=sys.stderr
  nargs=26
  arg=getArgs(nargs,'batchMCI.py jobNameOrNone dataInOrNone l0,l1,l2 sigma ewaldInOrNone sweepsMultOS0,1,..(or0.0) kTeffOS0,1,..(or0.0) neutralOverall?(0/1) rotation template.para template.sbatch jobNStart jobIterations jobMax maxCpuPerNode var0Nm,var00,var01,..:var1Nm,var10,var11,..:.. umbrK(orNone) umbrR0,1,...(orNone) enFreq dumpFreq dataFreq suComFreq umbrFreq equiStep prodSteps')
  
  name=next(arg)
  if name.lower()=='none': name=None
  dataIn=next(arg)
  if dataIn.lower()=='none': dataIn=None
  #N=int(next(arg))
  lStr=next(arg)
  sigma=float(next(arg))
  ewaldIn=next(arg)
  if ewaldIn.lower()=='none':
    ewaldIn='none'
  sweepsMult=[float(i) for i in next(arg).split(sep=',')]
  kTeff=[float(i) for i in next(arg).split(sep=',')]
  neutralOverall=next(arg)
  rotation=next(arg)
  paraIn=next(arg)
  sbatchIn=next(arg)
  jobNStart=int(next(arg))
  jobIterations=int(next(arg))
  jobMax=int(next(arg))
  maxCpuPerNode=int(next(arg))
  var=next(arg).split(sep=':')
  umbrK=next(arg)
  umbrR=next(arg)
  enFreq=next(arg)
  dumpFreq=next(arg)
  dataFreq=next(arg)
  suComFreq=next(arg)
  umbrFreq=next(arg)
  eSteps=next(arg)
  pSteps=next(arg)

  var=[item.split(sep=',') for item in var]
  varkeys=[item[0] for item in var] #keep order of input intact
  nkeys=len(varkeys)
  var={ item[0] : item[1:] for item in var }
  for key in varkeys:
    print('{}:{}'.format(key,var[key]))

  if umbrK.lower()=='none': umbrK=None
  else:                     umbrK=float(umbrK)
  if umbrR.lower()=='none': umbrR=None
  else:                     umbrR=[float(i) for i in umbrR.split(sep=',')]

  N,Nsu=getNandNsu(paraIn)

  with open(sbatchIn,'r') as f:
    lines=[line.strip() for line in f]
  nlines=len(lines)

  if umbrR is not None:
    paraIn='../'+paraIn
    ewaldIn='../'+ewaldIn
  jl=JobLauncher(N,Nsu,name,var,varkeys,paraIn,umbrK,
                 rotation,lStr,sigma,ewaldIn,enFreq,dumpFreq,
                 dataFreq,suComFreq,umbrFreq,eSteps,pSteps,neutralOverall)

  rngTable=RngSeedTable()
  numRNs=2    #hardcode: need 2/run for pcg32: seed & seq.

  nWorkers=totalWorkers(var,varkeys)
  nOverseers=totalOverseers(var,varkeys,maxCpuPerNode)

  if sweepsMult[0]==0.0 and len(sweepsMult)==1:
    sweepsMult=[0.0 for i in range(0,nOverseers)]
  if kTeff[0]==0.0 and len(kTeff)==1:
    kTeff=[0.0 for i in range(0,nOverseers)]

  if umbrR is not None:
    paraIn='../'+paraIn
    seeds=[[[rngTable.getSeed(numRNs) for i in range(0,jobIterations)] for i in range(0,nWorkers)] for i in range(0,len(umbrR))]
    execInParentDir=True
    for i,umbrRi in enumerate(umbrR):
      try:
        os.mkdir('{}'.format(umbrRi))
      except FileExistsError:
        pass
      os.chdir('{}'.format(umbrRi))
      jl.makeParas(umbrRi)
      jl.submitBatchLoop(nWorkers,umbrRi,jobNStart,jobIterations,jobMax,dataIn,lines,seeds[i],execInParentDir,sweepsMult,kTeff,maxCpuPerNode)
      os.chdir('..')
  else:
    seeds=[[rngTable.getSeed(numRNs) for i in range(0,jobIterations)] for i in range(0,nWorkers)]
    execInParentDir=False
    jl.makeParas(umbrR)
    jl.submitBatchLoop(nWorkers,umbrR,jobNStart,jobIterations,jobMax,dataIn,lines,seeds,execInParentDir,sweepsMult,kTeff,maxCpuPerNode)


if __name__ == "__main__":
  main()
