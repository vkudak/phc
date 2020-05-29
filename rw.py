class group(object):
    """Group class"""
    def __init__(self, Date=[], Time=[], dt=0, B=[], V=[], c=0):
        self.Date = []
        self.Time = []
        self.dt = 0.0
        self.B = []
        self.V = []
        self.c = 0


class star(object):
    '''Star object'''
    def __init__(self, B=[], RMSb=0, V=[], RMSv=0, ImpB=0, ImpV=0,
                 mB=0, mV=0, bmv=0, c=0, UT=0, RA=0, DEC=0, Mz=0):
        self.B = []
        self.RMSb = 0
        self.V = []
        self.RMSv = 0
        self.ImpB = 0
        self.ImpV = 0
        self.mB = 0
        self.mV = 0
        self.bmv = 0
        self.c = 0
        self.UT = 0
        self.RA = 0
        self.DEC = 0
        self.Mz = 0


def read(path):
    '''
    Read Mastugin photometry file and return GroupList(grList)
    GroupCount is len(grList)-1 becouse groups starts from 1.
    Each group is a group objets'''
    f = open(path, "r")
    f.readline()
    ### Date
    date = f.readline()
    date = date.split(':', 1)
    date = date[1]
    #print date
    ### TIME
    Time = f.readline()
    Time = Time.split("Begin time :")
    Time = Time[1]
    Time = Time.split(' ')
    #print Time
    ###Groups count
    #GrN = len(Time) - 1
    #print GrN
    #  Dt
    dt = f.readline()
    dt = dt.split(' ')
    #print dt
    # f.readline()
    culm_time = f.readline()
    culm_time = culm_time.replace(" ", "") #.strip()
    f.readline()
    f.readline()  # pass
    f.readline()
    f.readline()
    ###GrNumber
    f.readline()
    N = 1
    grList = []
    grList.append(0)  # fictive!!! just for start from #1
    grList.append(group())
    #f.readline
    for line in f:
        #line = lines
        #print line
        if line.split()[0] == 'GROUP' or line[0] in ['    ', '#'] or len(line) < 12:
            pass
        elif int(line.split(',')[1]) == 999999:
            grList[N].Time = Time[N - 1]
            grList[N].Date = date
            grList[N].dt = float(dt[N - 1]) / 1000
            grList[N].c = len(grList[N].B)
            #print N,grList[N].c
            N = N + 1
            grList.append(group())
        else:
            l = line.split('.')
            l = l[1:]
            l = l[0].split(',')
            grList[N].B.append(int(l[0]))
            grList[N].V.append(int(l[1]))
    f.close()
    grList.pop()
    return grList, culm_time


def test():
    test_file = "test.txt"
    L = read(test_file)
    print 'Total Group count=', len(L) - 1
    print 'Group #3 time=', L[3].Time
    print 'Group #3 dt=', L[3].dt
    print 'Group #3 count=', L[3].c
    print 'Group #3 V array=', L[3].V
