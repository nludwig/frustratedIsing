import sys
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

#example/template main
def main():
  err=sys.stderr #alias
  nargs=3      #hardcode
  arg=getArgs(nargs,'NAME.py inF out')
  inF=next(arg)
  out=next(arg)

  print('not implemented!',file=err)
  exit(1)

  with open(inF,'r') as f:
    inp=[line.strip() for line in f]
  for i,line in enumerate(inp):
    with open(line.strip(),'r') as f:
      pass
      #get filenames from file
  pass
  #do other stuff

if __name__ == "__main__":
  main()
