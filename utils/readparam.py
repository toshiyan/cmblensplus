
class config:

    def __init__ (self,filename=''):

        self.filename = filename
        self.pdict = {}

    def pdict(self):
    
        f = open(self.filename)

        for line in f:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            s = line.split('#')
            line = s[0]
            s = line.split('=')
            if len(s) != 2:
                print("Error parsing line:")
                print(line)
                continue
            key = s[0].strip()
            val = s[1].strip()
            self.pdict[key] = val

        f.close()


    def getstr(self,p,val=None):
        
        return self.pdict[p]


    def getint(self,p,val=None):

        return int(self.pdict[p])



