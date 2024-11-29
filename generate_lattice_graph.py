# The main entry point for the script
import argparse
import sys
from math import floor
#
#   Arguments and parsing
#
parser = argparse.ArgumentParser(description='Generate a latice network')
parser.add_argument('-row', help='number of rows', type=int,default=10)
parser.add_argument('-col', help='number of cols', type=int,default=10)
parser.add_argument('-srlg', help='the number of height each srlg', type=int,default=2)


args = parser.parse_args()

row=args.row
col=args.col
srlg=args.srlg
original_stdout = sys.stdout # Save a reference to the original standard output

filename='net/latice_row'+str(row)+'_col'+str(col)+'_srlg'+str(srlg)+'.grf'
print(filename)
with open(filename, 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    print(row*col+2)
    id=0
    for c in range(col):
        for r in range(row):
            print(c+1,r,id)
            id+=1
    sid=id
    tid=id+1
    print(0,int(floor(r/2)),sid)
    print(c+2,int(floor(r/2)),tid)
    print(sid,tid)
    print(' ')
    id=0
    print((row-1)*col+row*(col-1)+2*row)
    for c in range(col-1):
        for r in range(row):
            print(r+c*row,r+(c+1)*row,id)
            id+=1
    for c in range(col):
        for r in range(row-1):
            print(r+c*row,r+1+c*row,id)
            id+=1
    for r in range(row):
        print(sid,r,id)
        id+=1
        print(r+(col-1)*row,tid,id)
        id+=1
    print(' ')
    srlgnum=(row-1)*(col-1)+1
    print(srlgnum)
    for s in range(srlgnum):
        strsrlg=""
        for ss in range(srlg):
            strsrlg+=str(s+ss)+' '
        print(srlg,strsrlg)
sys.stdout = original_stdout
