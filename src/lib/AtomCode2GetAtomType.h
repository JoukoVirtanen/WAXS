

def main():
        file='/home2/jouko/project/HeaderFiles/AtomCode.h'
        lines=IOUtils.readlines(file)
        for i in range(len(lines)):
                if lines[i][0:4]=='const':
                        res=lines[i].split(' ')[2]
                        res=res[i].split('_')[0]
                        underscore=lines[i].find('_')
                        equals=lines[i].find('=')
                        
