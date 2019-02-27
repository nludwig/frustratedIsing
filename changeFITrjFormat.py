import sys

def getArgs(nargs,errStr):
  err=sys.stderr #alias
  print('nargs:{}\tcmdline args:{}\nargs:{}'.format(nargs,len(sys.argv),sys.argv),file=err)
  if len(sys.argv) != nargs:
    print('usages: python3 {}'.format(errStr), file=err)
    exit(1)
  arg=iter(sys.argv[1:])
  return iter(sys.argv[1:])

def py2round(f):
  if round(f+1)-round(f) != 1:
    return f+abs(f)/f/2
  return round(f)

def prepForReadShared(inp):
  #get N
  for i in range(0,3):
    inp.readline()
  N=int(inp.readline())
  #determine number of iterations
  inp.seek(0,0)
  numLines=0
  for line in inp:
    numLines+=1
  header=9
  itera=py2round(numLines/(N+header))
  #find N+, N-, N0
  inp.seek(0,0)
  hi=0
  for line in inp:
    hi+=1
    if hi==header: break
  return N,itera,header

def siteFIToSiteL(typeFI,su,ctt):
  if su==-1:
    ctt[typeFI]+=1
    return ctt[typeFI]-1
  else:
    ctt[su+3]+=1
    return ctt[su+3]-1

def prepForReadLammps(inp):
  N,itera,header=prepForReadShared(inp)
  Np=0;Nm=0;Nz=0;suList=[]
  ni=0
  for line in inp:
    l=line.split()
    su=int(l[2])
    if su==-1:
      t=int(l[1])
      if t==0:   Nm+=1
      elif t==1: Np+=1
      elif t==2: Nz+=1
    else:
      suList.append(su)
    ni+=1
    if ni==N: break
  try:               Nsu=max(suList)+1
  except ValueError: Nsu=0
  Nt=[Nm,Np,Nz,Nsu]
  if Nsu>0:
    Nsites=[0 for i in range(0,Nsu)]
    inp.seek(0,0)
    hi=0
    for line in inp:
      hi+=1
      if hi==header: break
    ni=0
    for line in inp:
      l=line.split()
      su=int(l[2])
      if su!=-1:
        Nsites[su]+=1
      ni+=1
      if ni==N: break
    Nt+=Nsites
  inp.seek(0,0)
  return N,itera,header,Nt

def buildInitCtTForLammps(Nt):
  ctt=[0 for i in range(0,Nt[3]+3)]
  Nsites=Nt[4:]
  ctt[0]=0
  ctt[1]=Nt[0]
  ctt[2]=Nt[0]+Nt[1]
  s=Nt[0]+Nt[1]+Nt[2]
  for i in range(0,Nt[3]):
    ctt[3+i]=s
    s+=Nsites[i]
  return ctt

##fi format:
#site# typeFI su# x y z
#with typeFI:
#number type
#0      -
#1      +
#2      0
#3      hydrophobic
##lammps format:
#site# typeL x y z
#with typeL:
#number type
#0      -
#1      +
#2      0
#3+     su#-3
def fiLineToLammpsLine(line,ctt):
  siteFI,typeFI,su,x,y,z=[int(item) for item in line.split()]
  siteL=siteFIToSiteL(typeFI,su,ctt)
  typeL=typeFIToTypeLForLammps(typeFI,su)
  lout=[siteL,typeL,x,y,z]
  return '{} {} {} {} {}'.format(*lout)

def typeFIToTypeLForLammps(typeFI,su):
  if su==-1 and typeFI==3:
    print('typeFIToTypeL got su==-1 and typeFI==3. exiting',file=sys.stderr)
    exit(1)
  if su==-1: return typeFI
  else:      return su+3

def prepForReadHydroVMD(inp):
  N,itera,header=prepForReadShared(inp)
  Np=0;Nm=0;Nz=0;Nh=0
  ni=0
  for line in inp:
    l=line.split()
    t=int(l[1])
    if t==0:   Nm+=1
    elif t==1: Np+=1
    elif t==2: Nz+=1
    elif t==3: Nh+=1
    ni+=1
    if ni==N: break
  Nt=[Nm,Np,Nz,Nh]
  inp.seek(0,0)
  return N,itera,header,Nt

def buildInitCtTForHydroVMD(Nt):
  ctt=[0 for i in range(0,4)]
  ctt[0]=0
  ctt[1]=ctt[0]+Nt[0]
  ctt[2]=ctt[1]+Nt[1]
  ctt[3]=ctt[2]+Nt[2]
  return ctt

##fi format:
#site# typeFI su# x y z
#with typeFI:
#number type
#0      -
#1      +
#2      0
#3      hydrophobic
##lammps format:
#site# typeL x y z
#with typeL:
#number type
#0      -
#1      +
#2      0
#3+     su#-3
def fiLineToHydroVMDLine(line,ctt):
  siteFI,typeFI,su,x,y,z=[int(item) for item in line.split()]
  siteL=siteFIToSiteL(typeFI,-1,ctt)
  typeL=typeFI
  lout=[siteL,typeL,x,y,z]
  return '{} {} {} {} {}'.format(*lout)

def main():
  err=sys.stderr #alias
  nargs=4
  arg=getArgs(nargs,'python3 changeFITrjFormat.py fileIn.fitrj toLammps(0/1) toHydroVMD(0/1)')
  fNm=next(arg)
  toLammps=True if next(arg).lower()=='1' else False
  toHydroVMD=True if next(arg).lower()=='1' else False
  if (toLammps is True and toHydroVMD is True) or (toLammps is False and toHydroVMD is False):
    print('set ONE AND ONLY ONE of toLammps or toHydroVMD to 1',file=err)
    exit(1)
  lammpsEndOfHeader='ITEM: ATOMS id type x y z'
  with open(fNm,'r') as f:
    if toLammps is True:
      N,itera,header,Nt=prepForReadLammps(f)
      cttBase=buildInitCtTForLammps(Nt)
      for i in range(0,itera):
        hi=0
        for line in f:
          if hi<header-1: print('{}'.format(line),end='')
          else:           print('{}'.format(lammpsEndOfHeader))
          hi+=1
          if hi==header: break
        ctt=list(cttBase) #local copy
        ni=0
        for line in f:
          lline=fiLineToLammpsLine(line,ctt)
          print('{} '.format(lline))
          ni+=1
          if ni==N: break
    elif toHydroVMD is True:
      N,itera,header,Nt=prepForReadHydroVMD(f)
      cttBase=buildInitCtTForHydroVMD(Nt)
      for i in range(0,itera):
        hi=0
        for line in f:
          if hi<header-1: print('{}'.format(line),end='')
          else:           print('{}'.format(lammpsEndOfHeader))
          hi+=1
          if hi==header: break
        ctt=list(cttBase) #local copy
        ni=0
        for line in f:
          lline=fiLineToHydroVMDLine(line,ctt)
          print('{} '.format(lline))
          ni+=1
          if ni==N: break
      

if __name__ == "__main__":
  main()
