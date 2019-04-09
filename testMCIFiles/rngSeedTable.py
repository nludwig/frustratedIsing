#rngSeedTable.py
#written by Nicholas B. Ludwig
#started on 16/10/03 (yy/mm/dd)
#first working version on 16/10/12
#current version 16/10/18

import sys,fcntl

class RngSeedTable:
  """
  Creates and maintains tables of random numbers to be used for RNG seeds.

  Master file contains list of all tables of RNG seeds, with most recent
  at end. Tables are formatted as lines of numbers separated by spaces
  with a first line which defines the current position in the file to
  read from.

  Master file has lines of the form:
  tableName.txt SEED0 SEED1 ...
  where SEEDi's are the seeds used to create that table. 
  Tables have first lines of the form:
  curLine numLines curWord numWords masterFile

  Call getSeed(n) in order to get n new seeds from table; table header
  will be updated so next call will get new seeds. If reach last line of
  table file, new table file will be created using RNG and seeds based on
  last line of current table file. Master file will then be updated with
  the new table file name.

  Currently assumes requested number of seeds will be small (ie. 1 or 2)
  and so with reasonable table file line lengths (default 255), there is
  no risk of reading from more than 2 lines. In particular, ...
  getSeeds() and updateTableFirstLine() and haveEnoughSeedsLeft()
  make use of this assumption. Generalization should be possible without
  too much hair pulling.

  Note that since table files have a header for the first line, these
  files should have numLines+1 lines in total. Thus we should have a one
  line grace which may help to save from coding errors. I don't expect
  this to be such a waste of space that it is worth the arcane coding
  that would be necessary to make use of the extra line.

  firstRunArgs are chosen arbitrarily (one table small [<1MB], but expect
  it should be sufficient for a single project [~10^5 seeds]); may need
  tuning.

  Note that I get a bit paranoid with fcntrl flocks. I plan to use this
  class to pull RNGs with potentially simulataneous calls, and it is my
  hope that these flocks allow such simultaneity safely. It may be that
  this is overkill, but since speed is not effected, will keep it for now.
  """
  def __init__(self,overwrite=False,masterFile=None,prng='pcg',firstRunArgs=[1,255,0,255]):
    """
    Sets up seed table if none exists; else loads existing table.
    """
    if prng=='pcg':     self.prng=prng
    else:
      print('RngSeedTable: PRNG {} not supported; try "pcg". Exiting.'.format(prng),file=sys.stderr)
    while True:
      if masterFile is not None:
        self.masterF=masterFile
        mfRead=open(self.masterF,'r')
        fcntl.flock(mfRead, fcntl.LOCK_EX)
        for line in mfRead:
          pass
        fcntl.flock(mfRead, fcntl.LOCK_UN)
        mfRead.close()
        self.tableF=line.split()[0]
        #line0 format:
        #curLine numLines curWord numWords masterFile
        try:
          tfRead=open(self.tableF,'r')
        except FileNotFoundError:
          self.curLine=firstRunArgs[0]
          self.numLines=firstRunArgs[1]
          self.curWord=firstRunArgs[2]
          self.numWords=firstRunArgs[3]
          self.makeNewTable()
          break
        fcntl.flock(tfRead, fcntl.LOCK_EX)
        line=tfRead.readline().split()
        fcntl.flock(tfRead, fcntl.LOCK_UN)
        tfRead.close()
        self.curLine=int(line[0])
        self.numLines=int(line[1])
        self.curWord=int(line[2])
        self.numWords=int(line[3])
        if self.curLine==self.numLines:
          self.appendNewTableToMaster()
          self.makeNewTable()
        else:   break
      else:
        masterFile=self.makeMaster(overwrite=overwrite)

  #primary support function: return requested
  #number of random numbers from table & update
  #so next request will be handled smoothly.
  #possible updates include:
  ##change header of table file
  def getSeed(self,nNums):
    """
    Gets nNums seeds from table and returns them.

    If not enough seeds in current table, creates new table
    and updates master list.
    """
    if self.haveEnoughSeedsLeft(nNums):
      tfRead=open(self.tableF,'r')
      fcntl.flock(tfRead, fcntl.LOCK_EX)
      i=0
      for line in tfRead:  #read header + curLine lines
        i+=1
        if i==self.curLine+1:   break
      line=line.split()
      if self.curWord+nNums<=self.numWords:
        ret=[int(i) for i in line[self.curWord:(self.curWord+nNums)]]
      else:
        #assumes max one line overlap => small number seeds;
        #could generalize with some sort of loop
        ret=[int(i) for i in line[self.curWord:self.numWords]]
        line=tfRead.readline().split()
        for i in range(0,self.curWord+nNums-self.numWords):
          ret.append(int(line[i]))
      fcntl.flock(tfRead, fcntl.LOCK_UN)
      tfRead.close()
      self.updateTableFirstLine(nNums)
      return ret
    else:
      self.appendNewTableToMaster()
      self.makeNewTable()
      return self.getSeed(nNums)

  def haveEnoughSeedsLeft(self,nNums):
    """
    Checks if have nNums seeds left in table. Returns boolean based on findings.

    Makes certain that plenty of seeds are left in table for the creation of a
    new table.
    """
    if self.curLine==self.numLines:
      return False
    elif self.curLine==self.numLines-1  and  self.curWord+nNums-self.numWords >= 0:
      return False
    else: 
      return True

  def updateTableFirstLine(self,nNums):
    """
    Updates the first line of the current table to reflect that nNums
    seeds have been read from it.

    Creates new table if necessary.
    """
    tfRead=open(self.tableF,'r')
    fcntl.flock(tfRead, fcntl.LOCK_EX)
    lines=[line.strip() for line in tfRead]
    fcntl.flock(tfRead, fcntl.LOCK_UN)
    tfRead.close()
    #line0 format:
    #curLine numLines curWord numWords masterFile
    if self.curWord+nNums<self.numWords:
      self.curWord+=nNums
      lines[0]='{} {} {} {} {}'.format(self.curLine,self.numLines,self.curWord,self.numWords,self.masterF)
    else:
      self.curLine+=1
      self.curWord+=(nNums-self.numWords)
      lines[0]='{} {} {} {} {}'.format(self.curLine,self.numLines,self.curWord,self.numWords,self.masterF)
    tfWrite=open(self.tableF,'w')
    fcntl.flock(tfWrite, fcntl.LOCK_EX)
    for line in lines:
      print('{}'.format(line),file=tfWrite)
    fcntl.flock(tfWrite, fcntl.LOCK_UN)
    tfWrite.close()
    if not self.haveEnoughSeedsLeft(nNums):
      self.appendNewTableToMaster()
      self.makeNewTable()

  #scroll through defined file to line,word and
  #return nNums
  def fetchSeed(self,fName,line,word,nNums):
    """
    Outputs nNums seeds from fName starting from line, word. Does not
    update table first line.

    Use case is to find previously used seeds based on arg.s or to
    get new seeds without updating first line of table. If you want
    seeds for your RNG, use getSeed instead!
    """
    tfRead=open(fName,'r')
    fcntl.flock(tfRead, fcntl.LOCK_EX)
    i=0
    for l in tfRead:
      i+=1
      if i==line: break
    l=[int(i) for i in l.split()[word:word+nNums]]
    fcntl.flock(tfRead, fcntl.LOCK_UN)
    tfRead.close()
    return l

  #wraps fetchSeed to print out requested seeds to fOut
  def printSeed(self,fName,line,word,nNums,fOut=sys.stdout):
    """
    Wraps fetchSeed to print nNums seeds from fName starting at line,word to fOut.

    Use case is to find previously used seeds based on arg.s or to
    get new seeds without updating first line of table. If you want
    seeds for your RNG, use getSeed instead!
    """
    seeds=self.fetchSeed(fName,line,word,nNums)
    for seed in seeds:  print('{}'.format(seed),end=' ',file=fOut)
    print('',file=fOut)   #\n

  #use parameters on last line of current table to make a new table
  #get table prefix & suffix from master; printout first line, then
  #generate numbers from chosen RNG, default pcg
  def makeNewTable(self):
    """
    Creates new table using self.prng seeded from previous table.

    If this is used to make the first table, read seeds from master.
    These seeds are currently (161018) pulled from /dev/urandom .
    """
    mfRead=open(self.masterF,'r')
    fcntl.flock(mfRead, fcntl.LOCK_EX)
    l0=mfRead.readline().split()
    for i in range(1,len(l0)):
      l0[i]=int(l0[i])
    prefix=l0[0].split(sep='0')[0]      #assumes first entry->'prefix'+'0'+'.txt'

    mfRead.seek(0)      #in case only one line in mfRead
    for line in mfRead:
      pass
    fcntl.flock(mfRead, fcntl.LOCK_UN)
    mfRead.close()
    #lines of form:
    #"tableName#.txt            seed0 seed1 ... seed(n-1)"
    suffix=line.split()[0]
    suffix=suffix.split(sep='.')[0]
    suffix=suffix.split(sep=prefix)[-1]

    if self.prng=='pcg':
      try:      import pcg32_cy as pcg32
      except ImportError:
                import pcg32
      n=2
      rng=pcg32.rng()

    if suffix!='0':
      tfRead=open(self.tableF,'r')
      fcntl.flock(tfRead, fcntl.LOCK_EX)
      for line in tfRead:
        pass
      line=[int(i) for i in line.split()[:n]]
      if self.prng=='pcg':
        pcg32.srandom_r(rng,*line)
      self.appendTableSeedsToMaster(line)
      fcntl.flock(tfRead, fcntl.LOCK_UN)
      tfRead.close()
    else:       #ie. first table
      if self.prng=='pcg':
        pcg32.srandom_r(rng,*l0[1:])

    self.tableF='{}{}.txt'.format(prefix,suffix)
    tfWrite=open(self.tableF,'w')
    fcntl.flock(tfWrite, fcntl.LOCK_EX)
    self.curLine=1
    self.curWord=0
    #line0 format:
    #curLine numLines curWord numWords masterFile
    print('{} {} {} {} {}'.format(self.curLine,self.numLines,self.curWord,self.numWords,self.masterF),file=tfWrite)
    for i in range(0,self.numLines):
      for j in range(0,self.numWords):
        print('{}'.format(pcg32.random_r(rng)),end=' ',file=tfWrite)
      print('',file=tfWrite) #\n
    fcntl.flock(tfWrite, fcntl.LOCK_UN)
    tfWrite.close()
      
  #append new table to end of master file
  def appendNewTableToMaster(self):
    """
    Updates master file to reflect new table.

    Works in concert with appendTableSeedsToMaster in a hacky way.
    This function called first to append file name; appendTableSeeds...
    called second which edits the appended line to also include
    the seeds used to generate the new table.
    """
    mfRead=open(self.masterF,'r')
    fcntl.flock(mfRead, fcntl.LOCK_EX)
    lines=[]
    for line in mfRead:
      lines.append(line.strip())
    fcntl.flock(mfRead, fcntl.LOCK_UN)
    mfRead.close()
    lprefix=len(lines[0].split()[0].split(sep='0')[0])
    tmp=lines[-1].strip().split()[0].split(sep='.')      #split from .txt
    num=int(tmp[0][lprefix:])+1        #assumes no '.' in fname
    lines.append('{}{}.{}'.format(tmp[0][:lprefix],num,tmp[1]))
    mfWrite=open(self.masterF,'w')
    fcntl.flock(mfWrite, fcntl.LOCK_EX)
    for line in lines:
      print('{}'.format(line),file=mfWrite)
    fcntl.flock(mfWrite, fcntl.LOCK_UN)
    mfWrite.close()

  #append new table to end of master file
  def appendTableSeedsToMaster(self,nums):
    """
    Updates master file to reflect new table.

    Works in concert with appendNewTableToMaster in a hacky way.
    That function called first to append file name; this one
    called second to edit the appended line to also include
    the seeds used to generate the new table.
    """
    print('appendTableSeedsToMaster: input {}'.format(nums))
    lines=[]
    mfRead=open(self.masterF,'r')
    fcntl.flock(mfRead, fcntl.LOCK_EX)
    for line in mfRead:
      lines.append(line.strip())
    fcntl.flock(mfRead, fcntl.LOCK_UN)
    mfRead.close()
    s=''
    for num in nums: s+=str(num)+' '
    lines[-1]+='\t\t{}'.format(s)
    mfWrite=open(self.masterF,'w')
    fcntl.flock(mfWrite, fcntl.LOCK_EX)
    for line in lines:
      print('{}'.format(line),file=mfWrite)
    fcntl.flock(mfWrite, fcntl.LOCK_UN)
    mfWrite.close()

  #make new master file. Return name of master file.
  def makeMaster(self,overwrite=False,fname='seedMaster.txt'):
    """
    Creates a master file to manage seed tables.

    Lines of the masterfile will be of the form:
    tableName.txt       SEED0 SEED1 ...
    where the SEEDi's are the seeds used to generate
    that table.
    Note that seeds for the first table are created in this 
    function using /dev/urandom . Further tables are seeded
    from previous tables.
    """
    if overwrite is False:
      try:
        mfWrite=open(fname,'x')
      except FileExistsError:
        return fname      #if default filename file exists, don't overwrite
    else:
      mfWrite=open(fname,'w')
    fcntl.flock(mfWrite, fcntl.LOCK_EX)
    #roll & record initial seed(s) for prng
    if self.prng=='pcg':
      nNum=2
      bytesPerNum=4
    from os import urandom
    bytes=urandom(nNum*bytesPerNum)
    nums=[0]*nNum
    if bytesPerNum==4:
      for i in range(0,nNum):
        nums[i]=int.from_bytes(bytes[i*bytesPerNum:(i+1)*bytesPerNum],byteorder=sys.byteorder)

    prefix=fname.split(sep='Master') #assumes fname->'prefix'+'Master'+'.txt'
    if len(prefix)==2:
      print('{}{}{}{}'.format(prefix[0],'Table',0,prefix[1]),end='\t\t',file=mfWrite)
    for num in nums:
      print('{}'.format(num),end=' ',file=mfWrite)
    print('',file=mfWrite)   #\n
    fcntl.flock(mfWrite, fcntl.LOCK_UN)
    mfWrite.close()
    return fname

if __name__ == "__main__":
  nargs=4      #hardcode
  print('nargs:{}\tcmdline args:{}'.format(nargs,len(sys.argv)),file=sys.stderr)
  if len(sys.argv) != nargs:
    print('usage: python3 rngSeedTable.py testRun? "none"Or"overwrite" nSeeds', file=sys.stderr)
    sys.exit()
  print('{}'.format(sys.argv), file=sys.stderr)

  arg=iter(sys.argv)
  next(arg)     #skip BLAH.py
  testRun=int(next(arg))
  master=next(arg)
  nSeeds=int(next(arg))
  overwrite=False
  if master.lower()=='none':
    master=None
  elif master.lower()=='overwrite' or master.lower()=='ow':
    master=None
    overwrite=True
  if testRun==1:
    from time import perf_counter
    tstart=perf_counter()
    table=RngSeedTable(overwrite=overwrite,masterFile=master,firstRunArgs=[1,10,0,10])
    for i in range(0,1000):
      print('{}'.format(table.getSeed(nSeeds)))
    tstop=perf_counter()
    print('runtime: {:.5f} s total'.format(tstop-tstart))
  else:
    table=RngSeedTable(overwrite=overwrite,masterFile=master)
    print('{}'.format(table.getSeed(nSeeds)))
