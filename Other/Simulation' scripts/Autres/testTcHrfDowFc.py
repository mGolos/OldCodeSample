# #!/usr/bin/env python
# ''' ???
# @author: gmoaltos@hotmail.com
# '''

'##################################################  LIBRARIES  #################################################'
from pylab import *
from numpy import load as np_load
from time import time

def compareMatrices(lMat,lTitl, aspect=None, switch=0):
    nb = len(lMat)
    nbb= nb + nb%2
    nc = int(ceil(sqrt(nbb)))
    nl = int(floor(nbb / nc))
    figure('Comparaison entre Matrices')
    for i in xrange(nb):
        if switch: subplot(nc,nl,i+1)
        else:      subplot(nl,nc,i+1)
        if aspect:
            imshow(lMat[i], interpolation = 'nearest', aspect='auto')
        else:
            imshow(lMat[i], interpolation = 'nearest')
        title(lTitl[i])
        colorbar()
    show()

def HRF(ms = 1.):
    ''' The Heamodynamic response function
    '''
    T     = 10
    tau_s = 0.8
    tau_f = 0.4
    scale = 1000. / ms
    v_time= linspace(0., T, scale * T)
    sqrt_tfts = sqrt(1./ tau_f - 1./ (4.* tau_s ** 2))
    exp_ts    = exp(-0.5 * (v_time / tau_s))
    h         = exp_ts * sin(v_time * sqrt_tfts) / sqrt_tfts
    return h


TC = zeros((20*1000,4))
#TC = np_load('../FilesIN/Time_Courses/test_TC_998.npy')
BT = zeros(TC.shape)
tmax,N = TC.shape
for n in range(N):
    TC[500:1000+int(tmax*n/float(N)), n] = 1
Hrv = HRF(1)

t0 = time()
for k in range(N):
    t1 = time()
    BT[:,k] += convolve(TC[:,k], Hrv)[:tmax]
    print k,time()-t1
print 'tot :',time()-t0


from mpl_toolkits.mplot3d import axes3d, Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
X, Y = meshgrid(linspace(0,tmax/1000.,tmax/100),arange(N))
fig = figure()
fig.canvas.set_window_title('fgh')
ax = p3.Axes3D(fig)
ax.scatter3D(X,Y,TC[::100].T*BT.max(), c='b')
ax.scatter3D(X,Y,BT[::100].T, c='r')
ax.set_xlabel(r'$Time\ (s)$', fontsize=15)
ax.set_zlabel(r'$BOLD\ Signal$', fontsize=15)
#compareMatrices([ TC.T, BT.T], \
                #['TC','Btot'], aspect='auto')
#subplot(221); plot(TC[:,0])
#subplot(221); plot(TC[:,1])
#subplot(222); plot(BT[:,0])
#subplot(222); plot(BT[:,1])
xlim(0,20)
ylim(0,3)
show()

#BS = zeros(TC.shape)
#s1 = int(tmax/3.)
#s2 = int(2*tmax/3.)
#TC_= [TC[:s1], TC[s1:s2], TC[s2:]]
#div = linspace(0,tmax,3,endpoint=False).astype(int32)
#t0 = time()
#for k in range(N):
    #for d in range(len(div)):
        #BS[:,k] += convolve(TC_[d][:,k], append(zeros(div[d]), HRF()))[:tmax]
        #print k,d,time()-t0
#print 'sum :',time()-t0
#compareMatrices([ TC.T, BS.T, BT.T], \
                #['TC','Bsum','Btot'], aspect='auto')



from pylab import *
from scipy.stats import pearsonr, spearmanr, ss
from numpy import load as np_load, save as np_save
from time import time

def loadTimeCourse(fileAdress):
    ofi = open(fileAdress, 'r')
    TA = array([line.split() for line in ofi]).astype(float)
    ofi.close()
    return TA.T

def compareMatrices(lMat,lTitl, aspect=None):
    nb = len(lMat)
    nbb= nb + nb%2
    nc = int(ceil(sqrt(nbb)))
    nl = int(floor(nbb / nc))
    figure('Comparaison entre Matrices')
    for i in xrange(nb):
        subplot(nl,nc,i+1)
        if aspect:
            imshow(lMat[i], interpolation = 'nearest', aspect='auto')
        else:
            imshow(lMat[i], interpolation = 'nearest')
        title(lTitl[i])
        colorbar()
    show()

def compareFunctionnalConnectivity(FC1, FC2):
    if FC1.shape != FC2.shape:
        print 'Functional Connectomes shape are not the same'
        return 0
    FCC = zeros(FC1.shape)
    for i in range(FCC.shape[0]):
        FCC[i,i:] = FC1[i,i:]
        FCC[i:,i] = FC2[i:,i]
    return FCC

def test1(TC):  t0 = time();    corr = cov(TC.T);           print 'Covariance   ', time() - t0;  return corr
def test2(TC):  t0 = time();    corr = corrcoef(TC.T);      print 'Kaushik Ghose', time() - t0;  return corr
def test3(TC):  t0 = time();    corr = spearmanr(TC)[0];    print 'Spearman     ', time() - t0;  return corr
def test4(TC):
    t0  = time()
    ms  = TC.mean(axis=0)[(slice(None,None,None),None)]
    TCm = TC.T - ms
    TCss= sqrt(ss(TCm, axis=1))
    N   = TC.shape[1]
    corr= zeros((N, N))
    for i in xrange(N):
        corr[i,i:]  = dot(TCm[i:], TCm[i].T)
        corr[i,i:] /= TCss[i:] * TCss[i]
        corr[i:,i]  = corr[i,i:]
    print 'Pearson      ', time() - t0
    return corr
def test5(TC):
    ''' Attention for TC which is modified by the function
    '''
    t0  = time()
    ms  = TC.mean(axis=0)[(slice(None,None,None),None)]
    TC -= ms.T
    TCss= sqrt(ss(TC, axis=0))
    N   = TC.shape[1]
    corr= zeros((N, N))
    for i in xrange(N):
        corr[i,i:]  = dot(TC[:,i:].T, TC[:,i])
        corr[i,i:] /= TCss[i:] * TCss[i]
        corr[i:,i]  = corr[i,i:]
    print 'Pearson      ', time() - t0
    return corr

#TCa  = loadTimeCourse('../FilesIN/Time_Courses/TC_998_1a.dat') # 66
#TCa += loadTimeCourse('../FilesIN/Time_Courses/TC_998_2a.dat') # 66
#TCa += loadTimeCourse('../FilesIN/Time_Courses/TC_998_3a.dat') # 66
#TCa += loadTimeCourse('../FilesIN/Time_Courses/TC_998_4a.dat') # 66
#TCa += loadTimeCourse('../FilesIN/Time_Courses/TC_998_5a.dat') # 66
#TCb  = loadTimeCourse('../FilesIN/Time_Courses/TC_998_1b.dat') # 66
#TCb += loadTimeCourse('../FilesIN/Time_Courses/TC_998_2b.dat') # 66
#TCb += loadTimeCourse('../FilesIN/Time_Courses/TC_998_3b.dat') # 66
#TCb += loadTimeCourse('../FilesIN/Time_Courses/TC_998_4b.dat') # 66
#TCb += loadTimeCourse('../FilesIN/Time_Courses/TC_998_5b.dat') # 66
#corr1 = (test1(TCa / 5.) + test1(TCb / 5.)) / 2.
#corr2 = (test2(TCa / 5.) + test2(TCb / 5.)) / 2.
#corr3 = (test3(TCa / 5.) + test3(TCb / 5.)) / 2.
#corr4 = (test4(TCa / 5.) + test4(TCb / 5.)) / 2.
#corr5 = (test5(TCa / 5.) + test5(TCb / 5.)) / 2.

#TC = np_load('../FilesIN/Time_Courses/test_TC_66.npy')
TC = np_load('../FilesIN/Time_Courses/test_TC_998.npy')[:2000]
#corr1 = test1(TC)
#corr2 = test2(TC)
#corr3 = test3(TC)
#corr4 = test4(TC)
corr5 = test5(TC)

compareMatrices([corr1,corr2,corr3,corr4],['cov','corrcoef','spearmanr','pearson'])







from pylab import *
from time import time
import scipy.sparse as sparse
def matRef(TC,I):
    N = len(I)
    J = dot(ones((N,1)), [arange(N)]).astype(int)
    M = []
    for i in range(N):
        M.append([])
        for j in range(N):
            M[i].append(TC[I[i,j], J[i,j]])
    return M
def thresholding(x,theta):
    return 0.5 * (1. + tanh(30. * (x - theta)))
def twoDimOne(mat,N):
    vec=zeros(len(mat)*len(mat))
    for i in range(N):
        for j in range(N):
            vec[i*N+j] = N * mat[i,j] + j
    return vec.astype(int)
def boucle1(A0, tmax):
    global Ii, Ij, tau, W, theta
    Tx = zeros((Ii.max()+1, N))
    Tx[0] = A0.copy()
    t0 = time()
    for t in range(tmax):
        print t,Tx.shape, Ii.shape, Ij.shape
        TE = Tx[Ii, Ij]
        A = thresholding(TE)
        Tx = append([Tx[0] + dt / taux * (- Tx[0] + Psi* (W * A).sum(axis=1) + 0.3 * rand(N))], Tx[:-1], axis=0)
        theta += dt / tauT * ( - theta + Omega * (A.sum() / N + 0.6 * rand(N)) )
        #for n in range(4):
            #Tx[:,n] = append([TE + tau * (- TE + 1.* W.dot(thresholding(TE)) + 0.01)], Tx[:-1])
    print 'boucle 1 takes %.2fs' %(time()-t0)
    return Tx
def boucle2(A0, tmax):
    global Ii, Ij, tau, W, theta, N, Psi, Omega
    II = twoDimOne(Ii, N)
    TA = zeros((Ii.max()+1, N))
    TA[0] = A0.copy()
    for t in range(1,len(TA)):
        TA[t] = TA[0].copy()
    t1 = time()
    x  = rand(N)
    for t in range(tmax):
        print t,TA.shape, Ii.shape, Ij.shape
        t0=time()
        #TAE = TA.flatten()[II].reshape(N,N)
        TAE = take(TA, II).reshape(N,N)
        print (time()-t0)*1000
        t0=time()
        x     += dt / taux * (- x     + Psi   * (W.multiply(TAE).sum(1).A1 + 0.3 * rand(N)))
        print (time()-t0)*1000
        t0=time()
        theta += dt / tauT * (- theta + Omega * (TA[0].sum() / N + 0.8 * rand(N)) )
        print (time()-t0)*1000
        t0=time()
        #TA[1:] = TA[:-1]
        TA = roll(TA, 1, axis=0)
        print (time()-t0)*1000
        t0=time()
        TA[0]  = thresholding(x,theta)
        print (time()-t0)*1000
    print 'boucle 3 takes %.2fs' %(time()-t1)
    return TA
def boucle3(A0, tmax):
    global Ii, tau, W, theta, N, Psi, Omega, V
    TA = zeros((Ii.max()+1, N, 1))
    TA[0] = A0.copy()[:, newaxis]
    for t in range(1,len(TA)):
        TA[t] = TA[0].copy()
    t1 = time()
    x  = rand(N)
    TAE= matRef(TA,Ii)
    TAE0 = array(TAE)
    for t in range(tmax):
        print 'TAE', TAE[3][:3]
        print t,TA.shape, Ii.shape
        t0=time()
        #x     += dt / taux * (- x     + Psi   * (W.multiply(array(TAE)[:,:,0]).sum(1).A1                 + 0.3 * rand(N)))  # 180
        #x     += dt / taux * (- x     + Psi   * ((V * TAE).sum(1).T[0]                 + 0.3 * rand(N)))  # 180
        x     += dt / taux * (- x     + Psi   * ((V * TAE).sum(1).flatten()                 + 0.3 * rand(N)))  # 180
        #x     += dt / taux * (- x     + Psi   * (W.multiply(array(TAE)[:,:,0]).sum(1).A1                 + 0.3 * rand(N)))  # 180
        #x     += dt / taux * (- x     + Psi   * (array(W.multiply(array(TAE)[:,:,0]).sum(1)).T[0]       + 0.3 * rand(N)))  # 180
        #x     += dt / taux * (- x     + Psi   * ((W.toarray()[:,:,newaxis] * TAE).reshape(N,N).sum(1)   + 0.3 * rand(N)))  # 180
        #x     += dt / taux * (- x     + Psi   * ((array(W)[:,:,newaxis] * TAE).reshape(N,N).sum(1)      + 0.3 * rand(N)))  # 180
        #x     += dt / taux * (- x     + Psi   * (array(W.multiply(array(TAE).reshape(N,N)).sum(1)).T[0] + 0.3 * rand(N)))  # 180
        print (time()-t0)*1000
        t0=time()
        theta += dt / tauT * (- theta + Omega * (TA[0].sum() / N + 0.8 * rand(N)) )
        print (time()-t0)*1000
        t0=time()
        TA[1:] = TA[:-1]
        print (time()-t0)*1000
        t0=time()
        TA[0]  = thresholding(x,theta)[:,newaxis]
        print (time()-t0)*1000
    TAEF=array(TAE)
    print 'boucle 3 takes %.2fs' %(time()-t1), ' TAE identique: ', (TAE0==TAEF).all()
    return TA

N    = 998
v    = 5      # m/s
dt   = 0.01   # ms
taux = 1      # ms
tauT = 6      # ms
Psi  = 30
Omega= 30
theta= 0.5
dens = 0.75
L    = 1000 * rand(N,N)
Ii   = floor(L / v / dt / 10.).astype(int)
Ij   = dot(ones((N,1)), [arange(N)]).astype(int)
A0   = rand(N)
V    = rand(N, N)
for t in range(N):
    V[t,permutation(N)[range(int(dens * N))]] = zeros(int(dens * N))
    V[:,t]= V[t,:].copy()
    V[t,t]= 0
V   /= norm(V)
W    = sparse.csr_matrix(V)
V    = V[:,:,newaxis]
tmax = 100
#TA1  = boucle1(A0, tmax)
TA2  = boucle2(A0, tmax)
#TA3  = boucle3(A0, tmax)
#imshow(TA3.T,aspect='auto'); show()
#print 'same ->', (Tx1==Tx2).all()




from pylab import *
def twoDimOne(mat):
    N = len(mat)
    vec=zeros(N * N)
    for i in range(N):
        for j in range(N):
            vec[i*N+j] = N * mat[i,j] + j
    return vec.astype(int)

#def listRef():#TC, Ind):
    ##N  = len(Ind)
    #N = 3
    #ll = []
    #for i in range(N):
        #ll.append([])
        #for j in range(N):
            #ll[i].append(1)
    #return ll

    #TC[Ind[i,j]

N = 3
T = 5
a = arange(N*T).reshape(T,N)
b = (T * rand(N,N)).astype(int)
c = dot(ones((N,1)), [arange(N)]).astype(int)
d = a[b,c]
e = a.base[twoDimOne(b)]
f = a[:,:,newaxis].copy()
g = f[b,c]


def matRef(TC,I):
    N = len(I)
    J = dot(ones((N,1)), [arange(N)]).astype(int)
    M = []
    for i in range(N):
        M.append([])
        for j in range(N):
            M[i].append(TC[I[i,j], J[i,j]])
    return M

