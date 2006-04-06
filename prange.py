#!/usr/bin/env python

def prange(N):
    if N>0:
	yield 0,0
    i = 1
    
    for exp in range(30,0,-1):
	base = 1<<exp
	for pi in range(base >> 1,N,base):
	    yield i,pi
	    i += 1

def piter(list):
    for i,p in prange(len(list)):
	yield list[p]

if __name__ == "__main__":
    for n,pn in prange(27):
	print n,pn
