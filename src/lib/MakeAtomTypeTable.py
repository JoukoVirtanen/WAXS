import re, os, sys

def readlines(file_path):
	file = open(file_path, 'r')
	lines = file.readlines()
	file.close()
	return lines

def main():
	print "hello"
	os.system("awk \'{if ($1==\"const\") print }\' /home2/jouko/project/HeaderFiles/AtomCode.h> /home/jouko/temp.txt")
	lines=readlines('/home2/jouko/temp.txt')
	NumLines=len(lines)
	print NumLines
	for i in range(NumLines):
		code=lines[i].split('=')
		residue=(code[0])[10:13]
		atom=(code[0])[14:len(code[0])]
                atom=atom.replace('_','\'')
		print residue+'\t'+atom 

if (__name__== '__main__'):
	main()
