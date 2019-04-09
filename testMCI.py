import sys
import subprocess as sp
import shlex
from getArgs import getArgs

def grep(string,lines):
  lstring=len(string)
  for line in lines:
    on=False
    same=0
    for lett in line:
      if lett==string[same]:
        same+=1
        on=True
      else:
        same=0
        on=False
      if same==lstring:
        break
    if on is True: yield line

def buildOptForJobSet(jobName,lStr,sigma,ewaldIn,sweepsMultOSStr,kTeffOSStr,neutOverall,rot,tempPara,tempSbatch,jobIter,
                      jobMax,maxCpuPerNode,paraList,umbrK,umbrR,enFreq,dumpFreq,suComFreq,umbrFreq,equiSteps,prodSteps):
  opt={}
  opt['jobName']=jobName
  opt['lStr']=lStr
  opt['sigma']=sigma
  opt['ewaldIn']=ewaldIn
  opt['sweepsMultOSStr']=sweepsMultOSStr
  opt['kTeffOSStr']=kTeffOSStr
  opt['neutOverall']=neutOverall
  opt['rot']=rot
  opt['tempPara']=tempPara
  opt['tempSbatch']=tempSbatch
  opt['jobIter']=jobIter
  opt['jobMax']=jobMax
  opt['maxCpuPerNode']=maxCpuPerNode
  opt['paraList']=paraList
  opt['umbrK']=umbrK
  opt['umbrR']=umbrR
  opt['enFreq']=enFreq
  opt['dumpFreq']=dumpFreq
  opt['suComFreq']=suComFreq
  opt['umbrFreq']=umbrFreq
  opt['equiSteps']=equiSteps
  opt['prodSteps']=prodSteps
  return opt

#batchMCI template:
#python3 batchMCI.py jobNameOrNone dataInOrNone neutralOverall?(0/1) template.para template.sbatch jobNStart jobIterations "var0Nm,var00,var01,..;var1Nm,var10,var11,..;.. enFreq dumpFreq dataFreq umbrFreq equiStep prodSteps'
def launchJobSet(opt):
  err=sys.stderr
  jobName=opt['jobName']
  lStr=opt['lStr']
  sigma=opt['sigma']
  ewaldIn=opt['ewaldIn']
  sweepsMultOSStr=opt['sweepsMultOSStr']
  kTeffOSStr=opt['kTeffOSStr']
  neutOverall=opt['neutOverall']
  rot=opt['rot']
  tempPara=opt['tempPara']
  tempSbatch=opt['tempSbatch']
  jobIter=opt['jobIter']
  jobMax=opt['jobMax']
  maxCpuPerNode=opt['maxCpuPerNode']
  paraList=opt['paraList']
  umbrK=opt['umbrK']
  umbrR=opt['umbrR']
  enFreq=opt['enFreq']
  dumpFreq=opt['dumpFreq']
  suComFreq=opt['suComFreq']
  umbrFreq=opt['umbrFreq']
  equiSteps=opt['equiSteps']
  prodSteps=opt['prodSteps']

#./MCIsing verbose?(0/1) spoof?(0/in.lammpstrj) neutralOverallOnSoluteInsertion?(0/1) rotation?(0/1) N L0,L1,L2 sigma ewaldIn(or'none') parametersIn0,1,.. dataIn0,1,..(or'none') sweepsMult0,1,..(or0.0) +kTeff0,1,..(or0.0) enFreq dumpFreq dataFreq suComFreq umbrFreq equiSteps prodSteps rngSeed,rngSeq confDump0,1,..

#usages: python3 batchMCI.py jobNameOrNone dataInOrNone l0,l1,l2 sigma ewaldInOrNone sweepsMultOS0,1,..(or0.0) kTeffOS0,1,..(or0.0) neutralOverall?(0/1) rotation template.para template.sbatch jobNStart jobIterations jobMax maxCpuPerNode var0Nm,var00,var01,..:var1Nm,var10,var11,..:.. umbrK(orNone) umbrR0,1,...(orNone) enFreq dumpFreq dataFreq suComFreq umbrFreq equiStep prodSteps

  cmd='python3 batchMCI.py {} none {} {} {} {} {} {} {} {} {} 0 {} {} {} {} {} {} {} {} 0 {} {} {} {}' \
      .format(jobName,lStr,sigma,ewaldIn,sweepsMultOSStr,kTeffOSStr,neutOverall,rot,tempPara,tempSbatch, \
      jobIter,jobMax,maxCpuPerNode,paraList,umbrK,umbrR,enFreq,dumpFreq,suComFreq,umbrFreq,equiSteps,prodSteps)
  #cmd='python3 batchMCI.py {} none {} {} {} 0 {} {} {} {} 0 {} {} {}' \
  #      .format(jobName,neutOverall,tempPara,tempSbatch,jobIter,paraList,enFreq,dumpFreq,umbrFreq,equiSteps,prodSteps)
  cmdcall=shlex.split(cmd)
  try:
    ret=sp.run(cmdcall,check=True)
  except sp.CalledProcessError:
    print('cmdcall failed\n{}'.format(cmdcall),file=err)
    #print('launchJobSet: got CalledProcessError',file=err)
    #print('returncode: {}'.format(ret.returncode),file=err)
    #print('cmd: {}'.format(ret.cmd),file=err)
    #print('output: {}'.format(ret.output),file=err)
    #print('stdout: {}'.format(ret.stdout),file=err)
    #print('stderr: {}'.format(ret.stderr),file=err)
    exit(1)

def main():
  err=sys.stderr
  nargs=8
  arg=getArgs(nargs,'testMCI.py tempParaIPfx(orNone) tempParaFIPrfx(orNone) tempSbatch(orNone) l0,l1,l2 ewaldInF(orNone) umbrKIn(orNone) umbrRIn(orNone)')
  print('{}'.format(sys.argv), file=err)

  tempParaI=next(arg)
  tempParaFI=next(arg)
  tempSbatch=next(arg)
  l=[int(i) for i in next(arg).split(sep=',')]
  ewaldInF=next(arg)
  umbrKIn=next(arg)
  umbrRIn=next(arg)

  lStr=','.join(str(i) for i in l)
  N=1
  for item in l: N*=item

  #defaults
  if tempParaI[-4:]=='para' or tempParaFI[-4:]=='para':
    print('tempPara should be PFXs. attempting to proceed by chopping .para',file=err)
    tempParaI=tempParaI[:-5]
    tempParaFI=tempParaFI[:-5]
  if tempParaI.lower()=='none':
    tempParaI='template-i'
  if tempParaFI.lower()=='none':
    tempParaFI='template-fi'
  if tempSbatch.lower()=='none':
    tempSbatch='template.sbatch'

  #cases:
  ##Pure Ising
  ###0 su, swap v clst
  ###1 su
  ####no move s v c
  ####yes move s v c
  ###2 su
  ####no move s v c
  ####yes move s v c
  ####umbr s v c
  ###F vs T trending to cmp. to theory (??)

  ##Frustrated Ising
  ###"

  isingTypes=['I','FI']
  suTypes=['0','1','2','u']
  suMvTypes=['0','1']
  mvTypes=['swap','clst']

  import datetime
  now=datetime.datetime.now()
  yr=str(now.year)[2:]
  mo=str(now.month)
  da=str(now.day)
  if len(mo)==1: mo='0'+mo
  if len(da)==1: da='0'+da
  nowPfix=''.join([yr,mo,da])

  for iType in isingTypes:
    for mvType in mvTypes:
      for suType in suTypes:
        for suMvType in suMvTypes:
          if suType=='0' and suMvType!='0':
            continue
          #for init'l testing
          if iType!='FI': continue
          if suType=='2' or suType=='u': continue

          sigma='1.0'

          if iType=='I':
            ewaldIn='none'
            sweepsMultOSStr='0.0'
            kTeffOSStr='0.0'
            tempPara=tempParaI
          else:
            ewaldIn=ewaldInF
            sweepsMultOSStr='0.0'
            kTeffOSStr='0.0'
            tempPara=tempParaFI

          neutOverall=1
          rot=0

          umbrK='none'
          umbrR='none'
          if suType=='0':
            tempPara+='-0'
          elif suType=='1':
            tempPara+='-1'
          elif suType=='2':
            tempPara+='-2'
          elif suType=='u':
            tempPara+='-u'
            umbrK=umbrKIn
            umbrR=umbrRIn
          if suType!='0' and suType!='u':
            if suMvType=='0':
              tempPara+='-0'
            elif suMvType=='1':
              tempPara+='-1'
            elif suMvType=='2':
              tempPara+='-2'
            elif suMvType=='3':
              tempPara+='-3'
          tempPara+='.para'

          jobIter=2
          jobMax=jobIter
          maxCpuPerNode=1
          paraList='kT,10.0'
          #paraList='kT,10.0,11.0,12.0,13.0,14.0'
          desiredNOut=1000
          baseNumEquiSteps=500
          equiSteps=baseNumEquiSteps if mvType=='swap' else 0
          prodSteps=N*baseNumEquiSteps/100 if mvType=='clst' else 0
          enFreq=int(prodSteps/desiredNOut) if mvType=='clst' else int(equiSteps*N/desiredNOut)
          dumpFreq=enFreq*100
          suComFreq=0 if suTypes!='u' else enFreq
          umbrFreq=0 if suTypes=='0' else enFreq

          jobName='-'.join([nowPfix,iType,mvType,suType,suMvType])

          opt=buildOptForJobSet(jobName,lStr,sigma,ewaldIn,sweepsMultOSStr,kTeffOSStr,neutOverall,rot,tempPara,tempSbatch,jobIter,
                                jobMax,maxCpuPerNode,paraList,umbrK,umbrR,enFreq,dumpFreq,suComFreq,umbrFreq,equiSteps,prodSteps)
          #opt=buildOptForJobSet(jobName,neutOverall,tempPara,tempSbatch,jobIter,paraList,enFreq,dumpFreq,umbrFreq,equiSteps,prodSteps)
          print('\n\n\n\nattempting launch jobset iType {}, mvType {}, su{}, suMv{}\n\n\n'.format(iType,mvType,suType,suMvType),file=err)
          launchJobSet(opt)
          print('successful launch',file=err)

if __name__ == "__main__":
  main()
