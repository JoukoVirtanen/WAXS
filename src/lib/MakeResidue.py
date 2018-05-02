import re, os, sys

def readlines(file_path):
	file = open(file_path, 'r')
	lines = file.readlines()
	file.close()
	return lines

def main():
	print "hello"
	os.system("awk \'{if ($1==\"const\") print }\' /home2/jouko/HeaderFiles/AtomCode.h> /home/jouko/temp.txt")
	lines=readlines('/home2/jouko/temp.txt')
	NumLines=len(lines)
	print NumLines
	for i in range(NumLines):
		code=lines[i].split('=')
		code1=(code[0])[10:len(code[0])]
		print 'ResidueName[' + code1 + ']=\"' + (lines[i])[10:13] + '\";'

	for i in range(NumLines):
		code=lines[i].split('=')
		code1=(code[0])[10:len(code[0])]
		code2=(code[0])[14:len(code[0])]
		print 'AtomName[' + code1 + ']=\"' + code2 + '\";'
if (__name__== '__main__'):
	main()
