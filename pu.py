'''Additional routins
for photometry standartisation
'''

import numpy as np
from scipy import interpolate


def interp(A, an, bn):
    tck = interpolate.splrep(range(an), A, s=0)
    #d = 0.0
    d = float(an) / float(bn)
    #print an
    #print bn
    #print 'd=',d
    xnew = []
    i = 0
    while i < an:
        xnew.append(i)
        i = i + d
    A2 = interpolate.splev(xnew, tck, der=0)
    return A2


def chunks(l, n):
    return [l[i:i + n] for i in range(0, len(l), n)]


def RMS_mean(t_list, gate):
    l = chunks(t_list, 10)
    #print l
    for i in range(len(l)):
        rms = np.sqrt(np.mean(l[i] ** 2))
        m = np.mean(l[i])
        #print 'RMS=',rms, 'Mean=',m
        for j in range(len(l[i])):
            if l[i][j] > 3 * rms:
                l[i][j] = m
    return t_list


def RMS_del(A, value):
    '''Delete elements of array A until A.RMS>value'''
    A = np.array(A)
    while A.std(axis=0) > value:
        #rms = A.std(axis=0)
        mean = A.mean(axis=0)
        d = []  # X-mean
        maxx = 0
        for i in range(len(A)):
            d.append(abs(A[i] - mean))
            if d[i] > maxx:
                maxx = d[i]
                imax = i
        A = np.delete(A, imax)
    return A


def lsqFit(y, x):
    '''
    y=ax+c
    Return a, c, residual
    '''
    x = np.array(x)
    y = np.array(y)
    A = np.vstack([x, np.ones(len(x))]).T
    res = []
    # linearly generated sequence
    if y != []:
        wb = np.linalg.lstsq(A, y)  # obtaining the parameters
        a, c = wb[0]
        #residual = wb[1][0]
    #else:
    #    residual = 999
    for i in range(0, len(x)):
        res.append(abs(y[i] - (a * x[i] + c)))
    res = np.array(res)
    res_max = np.max(res)
    res_min = np.min(res)
    if res_max > abs(res_min):
        ind = np.argmax(res)
    else:
        ind = np.argmin(res)
        res_max = abs(min)
    return a, c, res_max, ind


def test_rms_mean():
    test_list = np.arange(22)
    print 'Start=', test_list
    test_list = RMS_mean(test_list, 10)
    print 'Stop=', test_list
    print 'Yeah!'
