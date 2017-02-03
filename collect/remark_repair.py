import os
import re

pattern = re.compile(r'{(?P<key>[^\{\}]*)\:(?P<value>[^\{\}]*)}')

def parse_remarks(filename,filedir=None):
    if filedir is None:
        fileloc=filename
    else:
        fileloc = os.path.join(filedir,filename)

    remarks=[]
    with open(fileloc,'r') as o:
        filelines= o.readlines()
        for eachline in filelines:
            if 'Remark:' in eachline:
                #print(eachline)
                rdict ={}
                for m in pattern.finditer(eachline):
                    result =m.groupdict()
                    rdict[result['key']]=result['value']

                remarks.append(rdict)
    return remarks

if __name__=="__main__":
    a=parse_remarks('sample.mol')
    print(a)
