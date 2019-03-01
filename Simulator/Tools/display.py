#!/usr/bin/env python
#-*- coding:Utf-8 -*-


from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.pyplot import cm
from matplotlib.collections import LineCollection
from pylab import array, log, where, isnan, get_cmap, figure, zeros, rand, axes, close, vstack
from pylab import arange, axis, sqrt, ceil, fill, hstack, fill_diagonal, fft, sin, cos
from pylab import linspace, ones, disp, meshgrid, exp, c_, r_, savefig, deg2rad, newaxis
from Tools.matrices import sortBy
from Tools.lists import sortByIrregularDimension
from Tools.functions import gaussian
from time import time
import copy


def printOptions(linewidth=142, **kwa):
    from numpy import set_printoptions
    set_printoptions(linewidth=linewidth)


def inLign(data, t=None, labels=None, crowd=1., xlbl=None, ttl=None,
           ax=None, fig=None, fs=None, colors='k', offsets=None, **kwa):
    '''Put the different subplots data[0,:] to data[n,:] in different rows.
    '''
    numRows, numSamples = data.shape[0], data.shape[-1]
    
    if t == None:
        t = arange(numSamples)
    
    if labels == None:
        labels = arange(numRows)
    
    if ax == None:
        if fig == None:
            fig = figure(figsize=fs)
        ax = fig.add_subplot(111)
        ax.set_title(ttl)
        
    dmin, dmax = data.min(), data.max()
    dr = (dmax - dmin) * crowd
    y0, y1 = dmin, (numRows - 1) * dr + dmax
    ax.set_ylim(y0, y1)
    ax.set_xlim(t.min(), t.max())
    #ax.set_xticks(arange(10))
    ticklocs, segs = [], []
    
    if len(data.shape) == 2:
        for i in range(numRows):
            segs.append(hstack((t[:, newaxis], data[i, :, newaxis])))
            ticklocs.append(i * dr)
    elif len(data.shape) == 3:
        for j in range(data.shape[1]):
            for i in range(numRows):
                segs.append(hstack((t[:, newaxis], data[i, j, :, newaxis])))
                ticklocs.append(i * dr)

    if offsets == None:
        offsets = zeros((numRows*data.shape[1], 2), dtype=float)
        offsets[:, 1] = ticklocs

    lines = LineCollection(segs, offsets=offsets, transOffset=None, **kwa)
    lines.set_color(colors)
    ax.add_collection(lines)
    ax.set_yticks(ticklocs)
    ax.set_yticklabels(labels)
    if xlbl != None:
        ax.set_xlabel(xlbl)


def powerSpectra(TC, sr=0.5, nf=2, sw=False, avg=False, 
                 fig=None, ax=None, fs=(5,5), xlim=None, mean=False, half=True, **kwa):
    """Time has to be the first dimension.
    sr: sampling rate (Hz)
    """
    tx = TC.shape[0]
    t  = 2.* arange(tx)
    tmp = TC - avg * TC.mean(0)
    FT = []
    for i in range(TC.shape[1]):
        FT.append( abs(fft(tmp[:,i], n=nf*tx)[:tx] / tx)**2 )
    FT = array(FT)
    fq = linspace(0, sr, tx, endpoint=False) /nf
    
    if sw or fig != None:
        mx = len(fq) / 2.
        if fig == None:
            fig = figure(figsize=fs)
            ax = fig.add_subplot(111)
        if xlim == None:
            ax.set_xlim(0, fq[:mx].max())
        else:
            ax.set_xlim(*xlim)
        if mean:
            ax.plot(fq[:mx], FT[:,:mx].mean(0), 'k', **kwa)
        else:
            ax.plot(fq[:mx], FT[:,:mx].T, 'k', **kwa)
        ax.set_xlabel('Freq (Hz)')
        ax.set_ylabel('Amplitude')
        if sw:
            close(fig)
            return fig
    else:
        if half:
            h = len(fq) / 2
            return fq[:h], FT[:, :h]
        else:
            return fq, FT
    

def noteToSlides(adress, ie=False):
    from os import system, path
    system("ipython nbconvert reveal %s --to slides" %adress)
    if ie:
        system("start %s.slides.html" %path.realpath(adress.rstrip('ipynb')))
    

def tic(ret=False):
    global time_t0
    time_t0 = time()
    
    if ret == True:
        return time_t0


def tac(ret=False, t0=None, st=""):
    global time_t0
    
    if t0 != None: tUsed = t0
    else:          tUsed = time_t0
        
    if ret == True:
        return time() - tUsed
    
    elif ret == 'h':
        Te = time() - tUsed
        print 'Time elapsed: %ih%.2imin%.2is' %((Te/3600)%24, (Te/60)%60, Te%60) + st
        
    else:
        print "Time elapsed: %.2fs" %(time() - tUsed) + st


def estimAndPercent(i, Tmax, avg=100, t0=None):
    global time_t0
    
    if t0 != None: tUsed = t0
    else:          tUsed = time_t0
    
    if i == 0:
        return tic()
        
    t1 = tac(ret=True, t0=tUsed)
        
    if not((i+1)%(Tmax/10)):
        print '(%.2i%%), time elapsed: %.2fs' %(100*(i+1)/Tmax, t1)
       
    elif i == avg:
        Te = Tmax * t1 / avg
        disp('Time estimated: %ih%.2imin%.2is' %((Te/3600)%24, (Te/60)%60, Te%60))


def focus(tc, vmin, vmax):
    '''Time has to be the last dimension.
    '''
    nruter  = array(tc)
    nruter -= nruter.mean(-1)
    nruter /= max(nruter.max(-1), -nruter.min(-1))
    nruter *= (vmax-vmin)/2.
    nruter += (vmax-vmin)/2.
    return nruter


def mapMatrices(lMat, lTitl=None, lX=None, lY=None, fig=None, fs=None, labels=None, lxy=None,
                multp=False, disp='imshow', grid='subplot', x=None, ncl=None, shAll=True,
                label_mode='L', lign=False, transp=0, axV=True, axL=True, xpad=0.1, 
                cbar=False, cbarBin=False, fmt='%.2f', lvl=420, lxystep=[2,2], rotations=[0,0],
                cmapmax=[], tight=False, fontsize=20, **kwa):
    '''Show specific displays (imshow, plot) for list/nparrays as a grid
    depending on the size of the input lMat.
    TODO Trier tout ça
    '''
    from lists import dimensions
    
    # Options
    if type(lxystep) != list:
        lxystep = [lxystep, lxystep]
    
    # Create a figure if no inputed
    if fig == None:
        fig = figure(facecolor='white', figsize=fs)
      
    
    # Define the case depending on dimensions and the type of display.
    type1 = ['imshow','stackedHistograms','spy','contourf']
    type2 = ['plot']
    dims  = dimensions(lMat)
    ndim  = len(dims)
    
    
    if ncl != None:
        ncl = ncl[::-1]
        case = 4
    
    elif (ndim == 4 and disp in type1)\
    or (ndim == 3 and disp in type2 and not multp)\
    or (ndim == 4 and disp in type2 and multp):
        case = 1
        
    elif (ndim == 3 and disp in type1)\
    or   (ndim == 2 and disp in type2 and not multp)\
    or   (ndim == 3 and disp in type2 and multp):
        if lign:
            case = 2
        else:
            case = 3
    else:
        print "dimensions problem"
        return 0
              
    # Initialize number of raws/colums and reshape lMat and lTitl.
    # Raws and columns defined by the first dimensions.
    if case == 1:
        nc = len(lMat[0])
        nl = len(lMat)
        nb = nc * nl
        sh = dims[2:]
        sh.insert(0, nb)

        if transp:
            mats = lMat.reshape(*sh)
            xs = lX.reshape(*sh)
            if lTitl != None:
                tits = lTitl.reshape(*sh)
            else:
                tits = None
        else:
            mats = lMat.swapaxes(0,1).reshape(*sh)
            xs = lX.swapaxes(0,1).reshape(*sh)
            if lTitl != None:
                tits = lTitl.swapaxes(0,1).reshape(*sh)
            else:
                tits = None
                
        if lign:
            nl = nb
            nc = 1
        
    # Only one raw
    elif case == 2:
        nl = len(lMat)
        nb = nl
        nc = 1
        mats = lMat
        tits = lTitl
        xs = lX
        
    # Raws and columns generated by square root formula.
    elif case == 3:
        from numpy import ceil, sqrt
        nb  = len(lMat)
        nl  = int(ceil(sqrt(nb)))
        nc  = int(ceil(nb / float(nl)))
        mats = lMat
        tits = lTitl
        xs = lX
        
    # Number of raw and column defined in ncl
    elif case == 4:
        nc, nl = ncl
        nb = len(lMat)
        mats = lMat
        tits = lTitl
        xs = lX
        
    # Change raws for columns if transpose if True
    if transp:
        tmp = nc
        nc = nl
        nl = tmp
        
        order = hstack(arange(nc*nl).reshape(nl, nc).T).tolist()
        for i in range(nb, nc*nl):
            order.remove(i)
    else:
        order = arange(nb)
        
    # Function definition
    if grid.lower() == 'subplot':
        from matplotlib.axes import Axes
        fType = Axes
    elif grid.lower() == 'axesgrid':
        from mpl_toolkits.axes_grid1 import AxesGrid
        fType = AxesGrid._defaultLocatableAxesClass
        
    if disp == 'imshow':
        func = fType.imshow.im_func
    elif disp == 'contourf':
        func = fType.contourf.im_func
    elif disp == 'stackedHistograms':
        func = stackedHistograms
    elif disp == 'spy':
        func = Axes.spy.im_func
    elif disp == 'plot':
        if ndim != 2:
            #TODO to generalize
            mats = mats.swapaxes(-1,-2)
            if lX != None:
                xs = xs.swapaxes(-1,-2)
        func = fType.plot.im_func
        

    # Display depending on the type of grid
    if grid.lower() == 'subplot':            
        for i in range(nb):
            ax = fig.add_subplot(nc,nl,i+1)
            #im = ax.func(mats[order[i]], **kwa)
            #im = func(ax, mats[order[i]], **kwa)
            
            if cmapmax != []:
                mymap = betweenMap(cmapmax, mats[order[i]].min(), mats[order[i]].max())
                kwa.update({'cmap': mymap})
            
            if disp == 'stackedHistograms':
                im = func(mats[order[i]], ax=ax, **kwa)
            elif disp == 'contourf':
                if lX != None and lY != None:
                    im = func(ax, lX[order[i]], lY[order[i]], mats[order[i]],
                              levels=linspace(mats[order[i]].min(),mats[order[i]].max(),lvl), **kwa)
                else:
                    im = func(ax, mats[order[i]], levels=linspace(mats[order[i]].min(),mats[order[i]].max(),lvl), **kwa)
            elif x != None:
                try:
                    x[0][0]
                    im = func(ax, x[order[i]], mats[order[i]], **kwa)
                except:
                    im = func(ax, x, mats[order[i]], **kwa)
            elif lX != None:
                im = func(ax, xs[order[i]], mats[order[i]], **kwa)
            else:
                im = func(ax, mats[order[i]], **kwa)
                
            if labels != None:
                if i%nl == 0:
                    if type(labels[1]) != list: 
                        ax.set_ylabel(labels[1], fontsize=fontsize)
                    else:
                        ax.set_ylabel(labels[1][i/nl], fontsize=fontsize)
                if i >= nl * (nc-1):
                    if type(labels[0]) != list:
                        ax.set_xlabel(labels[0], fontsize=fontsize)
                    else: 
                        ax.set_xlabel(labels[0][i%nl], fontsize=fontsize)
                    
            if lxy != None:
                ax.set_xticks([])
                ax.set_yticks([])
                if i%nl == 0:
                    ax.set_yticks(list(range(len(lxy[1])))[::lxystep[1]])
                    ax.set_yticklabels(lxy[1][::lxystep[1]], rotation=rotations[1])
                if i >= nl * (nc-1):
                    ax.set_xticks(list(range(len(lxy[0])))[::lxystep[0]])
                    ax.set_xticklabels(lxy[0][::lxystep[0]], rotation=rotations[0])
                    
            
            if tits != None:
                ax.set_title(tits[order[i]])
            
            if not axV:
                #ax.set_visible(False)
                #ax.set_axis_off()
                ax.set_xticks([])
                ax.set_yticks([])
                
            if cbar:
                cb = fig.colorbar(im, ax=ax, format=fmt)
                if cbarBin:
                    cb.set_ticks([mats[order[i]].min(), mats[order[i]].max()])
                
            if not axL:
                axis('off')
    
    
    elif grid.lower() == 'axesgrid':
        if disp in ['imshow','contourf']:
            kwG  = {'cbar_location':'right',
                    'cbar_mode':'single'}
        
        elif disp in ['stackedHistograms','plot']:
            kwG  = {'aspect':False}
            
        ax = AxesGrid(fig, 111,
                      nrows_ncols = (nc,nl),
                      share_all = shAll,  # Share the same colorbar
                      label_mode = label_mode,
                      axes_pad = xpad,
                      **kwG)

        for i in range(nb):
            
            if cmapmax != []:
                mymap = betweenMap(cmapmax, mats[order[i]].min(), mats[order[i]].max())
                kwa.update({'cmap': mymap})

            #im = ax[i].func(mats[order[i]], **kwa)
            if disp == 'stackedHistograms':
                im = func(mats[order[i]], ax=ax[i], **kwa)
            elif x != None and disp == 'plot':
                try:
                    x[0][0]
                    im = func(ax[i], x[i], mats[order[i]], **kwa)
                except:
                    im = func(ax[i], x, mats[order[i]], **kwa)
            elif lX != None:
                im = func(ax[i], xs[order[i]], mats[order[i]], **kwa)
            else:
                im = func(ax[i], mats[order[i]], **kwa)
            
            if disp == 'imshow' and x!= None:
                im.set_extent([x.min(),x.max(),x.min(),x.max()])
            
            if labels != None:
                if i%nl == 0:
                    ax[i].set_ylabel(labels[1])
                if i >= nl * (nc-1):
                    ax[i].set_xlabel(labels[0])
                    
            if lxy != None:
                ax[i].set_xticks([])
                ax[i].set_yticks([])
                if i%nl == 0:
                    ax[i].set_ylabel(r'$\sigma_x$', fontsize=fontsize)
                    ax[i].set_yticks(list(range(len(lxy[1])))[::lxystep[1]])
                    ax[i].set_yticklabels(lxy[1][::lxystep[1]], rotation=rotations[1])
                if i >= nl * (nc-1):
                    ax[i].set_xlabel(r'$P$', fontsize=fontsize)
                    ax[i].set_xticks(list(range(len(lxy[0])))[::lxystep[0]])
                    ax[i].set_xticklabels(lxy[0][::lxystep[0]], rotation=rotations[0])
                    
            if tits != None:
                ax[i].set_title(tits[order[i]])
                
            if axV == False:
                #ax[i].set_visible(False)
                #ax[i].set_axis_off()
                ax[i].set_xticks([])
                ax[i].set_yticks([])
                
            if not axL:
                axis('off')
                
    if disp in ['imshow','contourf'] and grid.lower() != 'subplot':
        ax.cbar_axes[0].colorbar(im)

    fig.show(warn='off')
    if tight:
        fig.tight_layout()



def htmlAnimation(fig, func, frames=100, interval=20, blit=True):
    from IPython.display import HTML
    from matplotlib.animation import FuncAnimation
    from tempfile import NamedTemporaryFile
    
    VIDEO_TAG = """<video controls>
    <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
    Your browser does not support the video tag.
    </video>"""

    def anim_to_html(anim):
        if not hasattr(anim, '_encoded_video'):
            with NamedTemporaryFile(suffix='.mp4') as f:
                anim.save(f.name, fps=20, extra_args=['-vcodec', 'libx264'])
                video = open(f.name, "rb").read()
            anim._encoded_video = video.encode("base64")
        
        return VIDEO_TAG.format(anim._encoded_video)

    def display_animation(anim):
        close(anim._fig)
        return HTML(anim_to_html(anim))


    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = FuncAnimation(fig=fig,
                         func=func,
                         frames=frames,
                         interval=interval,
                         blit=blit)

    # call our new function to display the animation
    vid = display_animation(anim)
    
    return anim, vid


def stackedHistograms(lY, lX=None, X=None, color=None, rnorm=False, ax=None, order=None):
    '''Staked histograms.'''
    m = len(lY)
    
    if order == None:
        order = list(range(m))
    if lX == None or type(lX[0][0]) != int:
        from functions import indicesFrom2DList
        lX = indicesFrom2DList(lY)
    
    if X == None:
        allX = []
        for lXi in lX:
            allX.extend(lXi)
        X = list(set(allX))
    
    if ax == None:
        fig = figure()
        ax  = fig.add_subplot(111)
        
    if color == None:
        color = rand(m, 3)
        
    cusu0 = zeros(len(X))
    cusu1 = cusu0.copy()
    for i in order:
        if rnorm:
            cusu1[lX[i]] += (lY / lY.sum(axis=0))[i]
        else:
            cusu1[lX[i]] += lY[i]
            
        ax.fill_between(X, cusu0, cusu1, color=color[i], alpha=.7)
        cusu0 = cusu1.copy()
        
    ax.set_xlim((min(X), max(X)))
    ax.set_ylim(bottom=0)
    if rnorm:
        ax.set_ylim((0, 1))


def cmap_powerlaw_adjust(cmap, a):
    if a < 0.:
        return cmap
    cdict = copy.copy(cmap._segmentdata)
    fn = lambda x : (x[0]**a, x[1], x[2])
    for key in ('red','green','blue'):
        cdict[key] = map(fn, cdict[key])
        cdict[key].sort()
        assert (cdict[key][0]<0 or cdict[key][-1]>1), \
            "Resulting indices extend out of the [0, 1] segment."
    return LinearSegmentedColormap('colormap', cdict, 1024)


def cmap_cutter(cmap, part):  
    cdict = copy.copy(cmap._segmentdata)
    n = len(cdict['red']) / 2 + 1
    ndict = {}
    for key in ('red','green','blue'):
        ndict[key] = zeros((n,3))#[[0.]*3]*n
        for i in range(n):
            ndict[key][i][0]  = i / float(n-1)
            if part in ['up','U',0]:
                ndict[key][i][1:] = cdict[key][i+n-1][1:]
            elif part in ['down','D',1]:
                ndict[key][i][1:] = cdict[key][i][1:]
    return LinearSegmentedColormap('colormap', ndict, 1024)


def cmap_center_adjust(cmap, vmin, vmax, center=0, relative=0):    
    center_ratio = (center - vmin) / (vmax - vmin)
    if center_ratio <= 0.1:#\
    #and center >= vmax:
        cmap.set_under('w', 1.0)
        if relative:
            return cmap_powerlaw_adjust( cmap_cutter(cmap, part='up'), vmax)
        else:
            return cmap_cutter(cmap, part='up')
    elif 0.9 <= center_ratio:#\
    #and center <= vmin:
        cmap.set_over('w', 1.0)
        if relative:
            return cmap_powerlaw_adjust( cmap_cutter(cmap, part='down'), -1./ vmin)
        else:
            return cmap_cutter(cmap, part='down')
    else:
        a = log(abs(center_ratio)) / log(0.5)
        return cmap_powerlaw_adjust(cmap, a)
 

def initParcellation(adress='Regions/Images', org='Deco'):
    from pylab import imread
    lbls, down, avg = parcelOrganization(org)
    lbls = lbls[:33]
    
    int_0 = imread(adress + '/Int.png')      # Contours of internal part
    ext_0 = imread(adress + '/Ext.png')      # Contours of external part
    int_N = zeros((33,513,837))              # Internal part
    ext_N = zeros((33,513,837))              # External part
    for i in range(33):
        try:    int_N[i] = imread(adress + '/Int%s.png' %lbls[i][1:])
        except: pass
        try:    ext_N[i] = imread(adress + '/Ext%s.png' %lbls[i][1:])
        except: pass
    return int_0, ext_0, int_N, ext_N


def hagmannParcellGrid(fig=None, axes_pad=0., cbar_location='right', **kwa):
    '''Return an AxesGrid made for an hagmannParcellation.'''
    from mpl_toolkits.axes_grid1 import AxesGrid
    if fig == None:
        fig = figure(facecolor='white')
    
    return AxesGrid(fig, 111, nrows_ncols = (2,2),
                    share_all = True,
                    label_mode = 'L',
                    cbar_mode = 'single',
                    axes_pad = axes_pad,
                    cbar_location = cbar_location,
                    **kwa)


def hagmannParcellation(vec, color='RdBu_r', Ims=None, 
                        maxs=None, grid=None, fig=None, org='Deco', 
                        disp=True, saveN=None, dpi=300):
    '''Shows the Hagmann parcellation depending on the 66 or 998 weights of "vec"".'''
        
    # Shape test    
    if len(vec) not in [66,998]:
        print 'Vector shape not in [66, 998] as Hagmann parcellation'
        return None
        
    # Create a grid if no inputed
    if grid == None:
        # Create a figure if no inputed
        if fig == None:
            fig = figure(facecolor='white')
            
        from mpl_toolkits.axes_grid1 import AxesGrid
        fig.clear()
        grid = AxesGrid(fig, 111,
                        nrows_ncols = (2,2),
                        axes_pad = 0.0,
                        share_all = True,
                        label_mode = 'L',
                        cbar_location = 'right',
                        cbar_mode = 'single')
    
    elif fig == None:
        fig = grid._divider._fig
    
        
    # Downsampling if 998 nodes
    if len(vec) == 998:
        vec = downsample(vec, org)
        
    # Define maximi for imshow
    #if maxs[0] < (vec.max() - vec.min()) / 2.
    if maxs != None:
        vmin = maxs[0]
        vmax = maxs[1]
        bord = vmax
    else:
        vmin = vec.min()
        vmax = vec.max()
        if vmin == vmax:
            vmin = vmin * (vmax < 0)
            vmax = vmax * (vmax > 0) + 1.* (vmax==0)
        cRatio= - vmin / (vmax - vmin)
        bord  = where(cRatio < 0.9, vmax, vmin)
        
    # Load the image parts if not inputed
    if Ims == None:
        int_0, ext_0, int_N, ext_N = initParcellation(org=org)
    else:
        int_0, ext_0, int_N, ext_N = Ims
        
    # Generate the image from weightened parts
    int_R = bord * int_0 + sum([int_N[i]    * (vec[i])    for i in range(33)])
    ext_R = bord * ext_0 + sum([ext_N[i]    * (vec[i])    for i in range(33)])
    if org in ['D', 'Deco']:
        int_L = bord * int_0 + sum([int_N[-1-i] * (vec[33+i]) for i in range(33)])
        ext_L = bord * ext_0 + sum([ext_N[-1-i] * (vec[33+i]) for i in range(33)])
    elif org in ['Hagmann', 'H', 'Raw', 'R']:
        int_L = bord * int_0 + sum([int_N[i] * (vec[33+i]) for i in range(33)])
        ext_L = bord * ext_0 + sum([ext_N[i] * (vec[33+i]) for i in range(33)])
    

    # Display
    color = get_cmap(color)
    kwa = {'cmap':cmap_center_adjust(color, vmin, vmax, center=0, relative=0),
           'norm':Normalize(vmin=vmin, vmax=vmax, clip = False),
           }
    #color.set_under('w', 0.0)
    #kwa = {'cmap':color,
           #'norm':Normalize(vmin=vmin, vmax=vmax, clip = False),
           ##'alpha':0.5
           #}
    grid[0].clear()
    grid[1].clear()
    grid[2].clear()
    grid[3].clear()
    im1 = grid[0].imshow(int_L[:,::-1], **kwa)
    im2 = grid[1].imshow(int_R,         **kwa)
    im3 = grid[2].imshow(ext_L[:,::-1], **kwa)
    im4 = grid[3].imshow(ext_R,         **kwa)
    grid.cbar_axes[0].colorbar(im1)
    grid.axes_llc.set_xticks([])
    grid.axes_llc.set_yticks([])
    
    if saveN:
        fig.tight_layout()
        fig.savefig(saveN, dpi=dpi)
    if disp:
        fig.show()
    else:
        return fig


def downsample(vec, org='Deco', allAxis=False, nan0=True, diag=None, sum=False):
    '''Return the downsampled numpy array(s) of 66 nodes from the 998 nodes numpy array(s)
    with the 'org' organization.'''
    
    if allAxis:
        vec = downsample(vec.T, org=org, allAxis=False, nan0=nan0, sum=sum).T
    
    if type(org) != str:
        avg = array([(org == i).sum() for i in range(max(org)+1)])
    else:
        org, avg = parcelOrganization(org)[1:]
    
    L1, L2 = len(vec), len(avg)
    
    if sum:
        avg = 1.
    
    if vec.ndim == 1:
        vec66  = zeros(L2)
        for n in range(L1):
            vec66[org[n]] += vec[n]
        nruter = vec66 / avg
    
    elif vec.ndim == 2:
        if vec.shape[0] == L1:
            N, M = vec.shape
            vec66  = zeros((L2, M))
            for m in range(M):
                for n in range(N):
                    vec66[org[n], m] += vec[n, m]
                vec66[:, m] /= avg
            nruter = vec66
        
        elif vec.shape[1] == L1:
            M, N = vec.shape
            vec66  = zeros((M, L2))
            for m in range(M):
                for n in range(N):
                    vec66[m, org[n]] += vec[m, n]
                vec66[m] /= avg
            nruter = vec66
        
    else:
        print 'Dimension higher than 2.'
        return 0
    
    if allAxis and diag != None:
        fill_diagonal(nruter, diag * ones(L2))
    
    if nan0:
        return where(isnan(nruter), 0, nruter)
    else:
        return nruter


def TwoTriangle(mat1, mat2, posit=False, seuil=1):
    ''' Show two functional connectivity matrices using only one square matrix
    mat1 is upper right and mat2 bottom left
    '''
    if mat1.shape != mat2.shape:
        print 'Functional Connectomes shape are not the same'
        return 0
    mat12 = zeros(mat1.shape)
    for i in range(mat12.shape[0]):
        mat12[i,i:]   = seuil * mat2[i,i:]
        mat12[i+1:,i] = mat1[i+1:,i]
    if posit:
        mat12 = where(isfinite(mat12), mat12, 0.)
        return where(mat12 < 0, 0., mat12)
    else:     
        return mat12


def parcelAxis(axes, org='Deco', size=18, angle='vertical', l=1):
    axes.set_xticks(arange(66)*l)
    axes.set_yticks(arange(66)*l)
    lbls = parcelOrganization(org='Deco')[0]
    axes.set_xticklabels(lbls, fontsize=size, rotation=angle)
    axes.set_yticklabels(lbls, fontsize=size)
    return axes
    

def parcelOrganization(org='Deco'):
    if org in ['Hagmann', 'H']:
        return hagmannOrganization()
    elif org in ['Deco', 'D']:
        return decoOrganization()
    elif org in ['Raw', 'R']:
        return rawOrganization()
    

def hagmannOrganization():
    lbls = array(['rBTST', 'rCAC', 'rCMF', 'rCUN', 'rENT', 'rFP', 'rFUS', 'rIP', \
        'rIT', 'rISTC', 'rLOCC', 'rLOF', 'rLING', 'rMOF', 'rMT', 'rPARC', \
        'rPARH', 'rPOPE', 'rPORB', 'rPTRI', 'rPCAL', 'rPSTC', 'rPC', \
        'rPREC', 'rPCUN', 'rRAC', 'rRMF', 'rSF', 'rSP', 'rST', 'rSMAR', \
        'rTP', 'rTT', 'lBTST', 'lCAC', 'lCMF', 'lCUN', 'lENT', 'lFP', \
        'lFUS', 'lIP', 'lIT', 'lISTC', 'lLOCC', 'lLOF', 'lLING', 'lMOF', \
        'lMT', 'lPARC', 'lPARH', 'lPOPE', 'lPORB', 'lPTRI', 'lPCAL', \
        'lPSTC', 'lPC', 'lPREC', 'lPCUN', 'lRAC', 'lRMF', 'lSF', 'lSP', \
        'lST', 'lSMAR', 'lTP', 'lTT'], dtype='|S5')
    
    down = array([ \
        0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,\
        4,  4,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,\
        7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,\
        8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11,\
        11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,\
        13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15,\
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 19, 19,\
        19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,\
        21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,\
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,\
        24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,\
        26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27,\
        27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,\
        28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29,\
        29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 32, 32, 32, 33, 33, 33, 33, 33, 34, 34, 34, 34, 35,\
        35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36, 36, 36, 37, 37, 37, 38, 38, 39, 39, 39, 39, 39, 39, 39, 39, 39,\
        39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,\
        40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 42, 42, 42, 43, 43, 43, 43, 43,\
        43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44,\
        44, 44, 44, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47,\
        47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 49, 50,\
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 54, 54,\
        54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, 55,\
        55, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,\
        56, 56, 56, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 59, 59, 59, 59,\
        59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60,\
        60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 61, 61, 61,\
        61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62,\
        62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63,\
        63, 63, 63, 63, 64, 64, 64, 64, 65, 65, 65, 65])

    avg = array([ \
        7,  4, 13, 10,  2,  2, 22, 28, 19,  8, 19, 19, 17, 12, 20, 12,  6,\
        10,  6,  8, 10, 31,  7, 36, 23,  4, 22, 46, 27, 28, 16,  3,  3,  5,\
        4, 13,  8,  3,  2, 22, 25, 17,  8, 22, 20, 16, 12, 19, 11,  6, 11,\
        6,  7,  9, 30,  7, 36, 23,  4, 19, 50, 27, 29, 19,  4,  4])
    
    return lbls, down, avg


def decoOrganization():
    lbls = array(['rENT', 'rPARH', 'rTP', 'rFP', 'rFUS', 'rTT', 'rLOCC', 'rSP', 'rIT', \
        'rIP', 'rSMAR', 'rBTST', 'rMT', 'rST', 'rPSTC', 'rPREC', 'rCMF', \
        'rPOPE', 'rPTRI', 'rRMF', 'rPORB', 'rLOF', 'rCAC', 'rRAC', 'rSF', \
        'rMOF', 'rLING', 'rPCAL', 'rCUN', 'rPARC', 'rISTC', 'rPCUN', 'rPC', \
        'lPC', 'lPCUN', 'lISTC', 'lPARC', 'lCUN', 'lPCAL', 'lLING', 'lMOF', \
        'lSF', 'lRAC', 'lCAC', 'lLOF', 'lPORB', 'lRMF', 'lPTRI', 'lPOPE', \
        'lCMF', 'lPREC', 'lPSTC', 'lST', 'lMT', 'lBTST', 'lSMAR', 'lIP', \
        'lIT', 'lSP', 'lLOCC', 'lTT', 'lFUS', 'lFP', 'lTP', 'lPARH', 'lENT'],  \
        dtype='|S5')

    down = array([ \
        0,  0,  1,  1,  1,  1,  1,  1,  2,  2,  2,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,\
        4,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,\
        7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,\
        8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10,\
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,\
        12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14,\
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15,\
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16,\
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19,\
        19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,\
        21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,\
        24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25,\
        25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28,\
        28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31,\
        31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 34, 34, 34,\
        34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36,\
        36, 36, 36, 36, 36, 37, 37, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39,\
        39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,\
        41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 42, 42,\
        42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 45, 45, 45, 45, 45, 45, 46, 46,\
        46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48,\
        48, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,\
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51,\
        51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52,\
        52, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 55, 55, 55,\
        55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,\
        56, 56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58,\
        58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59,\
        59, 59, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 62, 62, 63,\
        63, 63, 63, 64, 64, 64, 64, 64, 64, 65, 65, 65])

    avg = array([ \
        2,  6,  3,  2, 22,  3, 19, 27, 19, 28, 16,  7, 20, 28, 31, 36, 13,\
        10,  8, 22,  6, 19,  4,  4, 46, 12, 17, 10, 10, 12,  8, 23,  7,  7,\
        23,  8, 11,  8,  9, 16, 12, 50,  4,  4, 20,  6, 19,  7, 11, 13, 36,\
        30, 29, 19,  5, 19, 25, 17, 27, 22,  4, 22,  2,  4,  6,  3])

    return lbls, down, avg


def rawOrganization():
    lbls = array(['rLOF', 'rPORB', 'rFP', 'rMOF', 'rPTRI', 'rPOPE', 'rRMF', 'rSF', \
        'rCMF', 'rPREC', 'rPARC', 'rRAC', 'rCAC', 'rPC', 'rISTC', 'rPSTC', \
        'rSMAR', 'rSP', 'rIP', 'rPCUN', 'rCUN', 'rPCAL', 'rLOCC', 'rLING', \
        'rFUS', 'rPARH', 'rENT', 'rTP', 'rIT', 'rMT', 'rBTST', 'rST', 'rTT', \
        'lLOF', 'lPORB', 'lFP', 'lMOF', 'lPTRI', 'lPOPE', 'lRMF', 'lSF', \
        'lCMF', 'lPREC', 'lPARC', 'lRAC', 'lCAC', 'lPC', 'lISTC', 'lPSTC', \
        'lSMAR', 'lSP', 'lIP', 'lPCUN', 'lCUN', 'lPCAL', 'lLOCC', 'lLING', \
        'lFUS', 'lPARH', 'lENT', 'lTP', 'lIT', 'lMT', 'lBTST', 'lST', 'lTT'],  \
        dtype='|S5')

    down = array([ \
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  2,  2,  3,  3,  3,  3,  3,  3,  3,\
        3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,\
        6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,\
        7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,\
        8,  8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,\
        9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14,\
        14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,\
        15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,\
        17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,\
        18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20,\
        20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 23,\
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,\
        24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 26, 26, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,\
        29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, 31,\
        31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33,\
        33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 34, 34, 34, 34, 34, 34, 35, 35, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 37, 37, 37, 37,\
        37, 37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 40,\
        40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,\
        40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 42,\
        42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 43, 43, 43, 43,\
        43, 43, 43, 43, 43, 43, 43, 44, 44, 44, 44, 45, 45, 45, 45, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48,\
        48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 49, 49, 49,\
        49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,\
        50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52,\
        52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54,\
        54, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,\
        56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 59,\
        59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62,\
        62, 62, 62, 62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,\
        64, 64, 64, 64, 64, 64, 64, 64, 65, 65, 65, 65])

    avg = array([ \
        19,  6,  2, 12,  8, 10, 22, 46, 13, 36, 12,  4,  4,  7,  8, 31, 16,\
        27, 28, 23, 10, 10, 19, 17, 22,  6,  2,  3, 19, 20,  7, 28,  3, 20,\
        6,  2, 12,  7, 11, 19, 50, 13, 36, 11,  4,  4,  7,  8, 30, 19, 27,\
        25, 23,  8,  9, 22, 16, 22,  6,  3,  4, 17, 19,  5, 29,  4])

    return lbls, down, avg


def plot4D(x,y,z, cmap='hsv', kwaFig={}, kwaSca={}, kwaLin={'lw':.5, 'ls':'-','c':'k','alpha':0.4},
           fig=None, line=False, lim=None, rev=None, cbar=False, view=None, ret=False):
    '''
    -lim: (xmin, xmax, ymin, ymax, zmin, zmax)
    '''
    
    from mpl_toolkits.mplot3d import Axes3D
    if fig == None:
        fig = figure(**kwaFig)
    ax = fig.add_subplot(111, projection='3d', aspect='auto')
    
    if line:
        ax.plot(x,y,z, **kwaLin)
    ax.scatter(x,y,z, c=cm.get_cmap(cmap)(linspace(0,255,len(x)).astype(int)), **kwaSca)
    
    if lim != None:
        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[2], lim[3])
        ax.set_zlim(lim[4], lim[5])
        
    if rev != None:
        if rev == 'rand':
            rev = (-1)**(rand(3)<0.5)
        ax.set_xlim(ax.get_xlim()[::rev[0]])
        ax.set_ylim(ax.get_ylim()[::rev[1]])
        ax.set_zlim(ax.get_zlim()[::rev[2]])
       
    if cbar:
        #ax = fig.add_axes()
        #ax.set_axis_off()
        sm = cm.ScalarMappable(cmap=cmap)
        sm._A = [0,len(x)]
        cbar = fig.colorbar(sm, shrink=0.5)
        #cbar.ax.set_yticklabels([])
    
    if view  != None:
        ax.view_init(*view)
        
    if ret:
        return ax, fig, cbar
    else:
        fig.show()
    


def _blob(x,y,area,colour):
    """
    Draws a square-shaped blob with the given area (< 1) at
    the given coordinates.
    """
    hs = sqrt(area) / 2
    xcorners = array([x - hs, x + hs, x + hs, x - hs])
    ycorners = array([y - hs, y - hs, y + hs, y + hs])
    fill(xcorners, ycorners, colour, edgecolor=colour)


def hinton(W, maxWeight=None):
    """
    Draws a Hinton diagram for visualizing a weight matrix.
    Temporarily disables matplotlib interactive mode if it is on,
    otherwise this takes forever.
    """
    from pylab import isinteractive, ioff, clf, title, axis, ion
    reenable = False
    if isinteractive():
        ioff()
    clf()
    height, width = W.shape
    if not maxWeight:
        maxWeight = 2**ceil(log(abs(W).max())/log(2))

    fill(array([0,width,width,0]),array([0,0,height,height]),'gray')
    title(r'$Hinton\ diagram$')
    axis('off')
    axis('equal')
    for x in range(width):
        for y in range(height):
            _x = x+1
            _y = y+1
            w = W[y,x]
            if w > 0:
                _blob(_x - 0.5, height - _y + 0.5, min(1,w/maxWeight),'white')
            elif w < 0:
                _blob(_x - 0.5, height - _y + 0.5, min(1,-w/maxWeight),'black')
    if reenable:
        ion()
        

def parcelFromCenters(x, centers, s=2, intext=False, lign=True,
                      tle=None, tlx=1.05, tly=0.55,
                      vmin=None, vmax=None,  cbar=True, cmap='Spectral_r',
                      fig=None, fs=None, fsave=None, dpi=300, fcol=None,
                      **kwa):
    '''Visualization of an array using dots of different size for their proximity with the viewer and colors for their values.
    x : spacial pattern
    centers: [x,y,z] coordinates
    s : dots size
    intext : internal visualization 
    lign : visualization over one lign (false: two ligns two columns)
    tle : title, with tlx and tly to change position if any problems
    fig : existing figure
    fs : figure size
    fsave : adress to save the figure
    fcol : background color.
    '''
    cmap = get_cmap(cmap)
    if fig == None:
        if fs == None:
            if lign: fs = (7,1)
            else:    fs = (5,3.1)
        fig = figure(figsize=fs)
        
    if vmin == None or vmax == None:
        vmin, vmax = x.min(), x.max()
        X = (x - vmin) / (vmax - vmin) * 256
    else:
        X = x.copy()
        if vmin != None:
            X = where(X > vmin, X, vmin)
        if vmax != None:
            X = where(X < vmax, X, vmax)
        X = (x - vmin) / (vmax - vmin) * 256
       
    if lign:
        nl, nc = 1, 4
        a1, a2, a3, a4 = 1, 2, 3, 4
    else:
        nl, nc = 2, 2
        a1, a2, a3, a4 = 1, 3, 4, 2
    
    x = centers[:,1] - (centers[:,1].max()+centers[:,1].min())/2.
    y = centers[:,2] - (centers[:,2].max()+centers[:,2].min())/2.
    z = centers[:,0] - (centers[:,0].max()+centers[:,0].min())/2.
    xm, xM = x.min(), x.max()
    ym, yM = y.min(), y.max()
    zm, zM = z.min(), z.max()
    
    def form(ax,inv=False):
        ax.axis('equal')
        ax.set_xlim((xm)*1.1, (xM)*1.1)
        ax.set_ylim((ym)*1.1, (yM)*1.1)
        ax.set_axis_off()
        return ax
    
    si = sortBy(z, inverse=True)[0]
    ax = fig.add_subplot(nl,nc,a1)
    ax.scatter(x[si], y[si], c=cmap(X[si].astype(int)), s=s * (-z[si]), **kwa)
    ax = form(ax)
    ax.invert_xaxis()
    if tle != None: 
        ax.text(tlx * max(ax.get_xlim()), tly * max(ax.get_ylim()), tle, fontsize=18)
    
    si = sortBy(z, inverse=False)[0]
    ax = fig.add_subplot(nl,nc,a4)
    ax.scatter(x[si], y[si], c=cmap(X[si].astype(int)), s=s * (z[si]), **kwa)
    ax = form(ax)
    
    if intext:
        si = sortBy(z, inverse=False)[0][:len(z)/2]
        ax = fig.add_subplot(nl,nc,a2)
        ax.scatter(x[si], y[si], c=cmap(X[si].astype(int)), s=s * (z[si]-zm), **kwa)
        ax = form(ax)
        
        si = sortBy(z, inverse=True)[0][:len(z)/2]
        ax = fig.add_subplot(nl,nc,a3)
        ax.scatter(x[si], y[si], c=cmap(X[si].astype(int)), s=s * (zM-z[si]), **kwa)
        ax = form(ax)
        ax.invert_xaxis()
    
    else:
        si = sortBy(x, inverse=True)[0]
        ax = fig.add_subplot(nl,nc,a2)
        ax.scatter(z[si], y[si], c=cmap(X[si].astype(int)), s=s * (xM - x[si]) / 3., **kwa)
        ax = form(ax)
        ax.set_xlim((zm)*1.2, (zM)*1.2)
        ax.set_ylim((ym)*1.2, (yM)*1.2)
        
        si = sortBy(y, inverse=False)[0]
        ax = fig.add_subplot(nl,nc,a3)
        ax.scatter(z[si], x[si], c=cmap(X[si].astype(int)), s=s * (y[si] + abs(ym)) / 2., **kwa)
        ax = form(ax)
        ax.set_xlim((zm)*1.1, (zM)*1.1)
        ax.set_ylim((xm)*1.1, (xM)*1.1)
     
    if cbar:
        ax = fig.add_axes([0.87, 0.05, 0.1, 0.87])
        #ax = fig.add_axes([0.9, 0., 0.1, 1.])
        ax.set_axis_off()
        sm = cm.ScalarMappable(cmap=cmap)
        sm._A = [vmin, vmax]
        cbar = fig.colorbar(sm)
        #cbar = fig.colorbar(sm, aspect=42)
        cbar.set_ticks(sm._A)
        #cbar.ax.set_yticklabels(sm._A, fontsize=8);
        #cbar.ax.set_yticklabels(['0']+9*['']+[vmax])
        fig.subplots_adjust(left=0, bottom=0, right=0.95, top=1, hspace=0., wspace=0.)
    else:
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0., wspace=0.)

    if fcol:
        fig.patch.set_facecolor(fcol)
    
    if fsave != None:
        fig.savefig(fsave, dpi=dpi, bbox_extra_artists=(cbar,))
        #fig.savefig(fsave, dpi=300)
    fig.show()


def linesFromCenters(x, centers, norm=True, s=2, vmin=None, vmax=None, fsave=None, tle=None, lines=None,
                     fig=None, cmap='Spectral_r', intext=False, lign=True, cbar=True, fcol=None, 
                     alphal=1, disp=True, **kwa):
    '''x [space]
    centers [x,y,z]
    '''
    cmap = get_cmap(cmap)
    if fig == None:
        if lign:
            fig = figure(figsize=(7,1))
        else:
            fig = figure(figsize=(5,3.1))
        
    if vmin == None or vmax == None:
        vmin, vmax = x.min(), x.max()
        X = (x - vmin) / (vmax - vmin) * 256
        #if norm:
            #X = x / x.max() * 256
            #vmax = '%.2f'%x.max()
        #else:
            #X = x * 256
            #vmax = '1'
    else:
        X = x.copy()
        if vmin != None:
            X = where(X > vmin, X, vmin)
        if vmax != None:
            X = where(X < vmax, X, vmax)
        X = (x - vmin) / (vmax - vmin) * 256
        #X = where(x < vmax, x, vmax) / vmax * 256
        #vmax = '%.2f'%vmax
       
    if lign:
        nl, nc = 1, 4
        a1, a2, a3, a4 = 1, 2, 3, 4
    else:
        nl, nc = 2, 2
        a1, a2, a3, a4 = 1, 3, 4, 2
        
    #f = lambda x,i,si: x[:,i][si].max() /2. - x[:,i][si]
    
    x = centers[:,1] - (centers[:,1].max()+centers[:,1].min())/2.
    y = centers[:,2] - (centers[:,2].max()+centers[:,2].min())/2.
    z = centers[:,0] - (centers[:,0].max()+centers[:,0].min())/2.
    xm, xM = x.min(), x.max()
    ym, yM = y.min(), y.max()
    zm, zM = z.min(), z.max()
    
    def form(ax,inv=False):
        ax.axis('equal')
        ax.set_xlim((xm)*1.1, (xM)*1.1)
        ax.set_ylim((ym)*1.1, (yM)*1.1)
        ax.set_axis_off()
        return ax    
    
    si = sortBy(z, inverse=True)[0]
    ax = fig.add_subplot(nl,nc,a1)
    if lines != None:
        ax.plot(x[lines], y[lines], alpha=alphal, c='k')
    ax.scatter(x[si], y[si], c=cmap(X[si].astype(int)), s=s * (-z[si]), **kwa)
    ax = form(ax)
    ax.invert_xaxis()
    if tle != None: 
        ax.text(1.05* max(ax.get_xlim()), .55* max(ax.get_ylim()), tle, fontsize=18)
    
    si = sortBy(z, inverse=False)[0]
    ax = fig.add_subplot(nl,nc,a4)
    if lines != None:
        ax.plot(x[lines], y[lines], alpha=alphal, c='k')
    ax.scatter(x[si], y[si], c=cmap(X[si].astype(int)), s=s * (z[si]), **kwa)
    ax = form(ax)
    
    if intext:
        si = sortBy(z, inverse=False)[0][:len(z)/2]
        ax = fig.add_subplot(nl,nc,a2)
        if lines != None:
            ax.plot(x[lines], y[lines], alpha=alphal, c='k')
        ax.scatter(x[si], y[si], c=cmap(X[si].astype(int)), s=s * (z[si]-zm), **kwa)
        ax = form(ax)
        
        si = sortBy(z, inverse=True)[0][:len(z)/2]
        ax = fig.add_subplot(nl,nc,a3)
        if lines != None:
            ax.plot(x[lines], y[lines], alpha=alphal, c='k')
        ax.scatter(x[si], y[si], c=cmap(X[si].astype(int)), s=s * (zM-z[si]), **kwa)
        ax = form(ax)
        ax.invert_xaxis()
    
    else:
        si = sortBy(x, inverse=True)[0]
        ax = fig.add_subplot(nl,nc,a2)
        if lines != None:
            ax.plot(z[lines], y[lines], alpha=alphal, c='k')
        ax.scatter(z[si], y[si], c=cmap(X[si].astype(int)), s=s * (xM - x[si]) / 3., **kwa)
        ax = form(ax)
        ax.set_xlim((zm)*1.2, (zM)*1.2)
        ax.set_ylim((ym)*1.2, (yM)*1.2)
        
        si = sortBy(y, inverse=False)[0]
        ax = fig.add_subplot(nl,nc,a3)
        if lines != None:
            ax.plot(z[lines], x[lines], alpha=alphal, c='k')
        ax.scatter(z[si], x[si], c=cmap(X[si].astype(int)), s=s * (y[si] + abs(ym)) / 2., **kwa)
        ax = form(ax)
        ax.set_xlim((zm)*1.1, (zM)*1.1)
        ax.set_ylim((xm)*1.1, (xM)*1.1)
     
    if cbar:
        ax = fig.add_axes([0.87, 0.05, 0.1, 0.87])
        ax.set_axis_off()
        sm = cm.ScalarMappable(cmap=cmap)
        sm._A = [vmin, vmax]
        cbar = fig.colorbar(sm)
        cbar.set_ticks(sm._A)
        #cbar.ax.set_yticklabels(sm._A, fontsize=8);
        #cbar.ax.set_yticklabels(['0']+9*['']+[vmax])
        fig.subplots_adjust(left=0, bottom=0, right=0.95, top=1, hspace=0., wspace=0.)
    else:
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0., wspace=0.)

    if fcol:
        fig.patch.set_facecolor(fcol)
    
    if fsave != None:
        fig.savefig(fsave, dpi=300)#bbox_extra_artists=(cbar,))
    if disp:
        fig.show()
    else:
        close(fig)
        return fig
        
        
def plot3DfromCenters(x, centers, fig=None, cmap='Spectral_r', s=50):
    from mpl_toolkits.mplot3d import Axes3D
    cmap = get_cmap(cmap)
    if fig == None:
        fig = figure(figsize=(10,10))
        
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.axis('equal')
    ax.scatter(centers[:,0], centers[:,1], centers[:,2], 
               c=cmap((x / x.max() * 256).astype(int)), s=s)
    ax.set_ylim(centers[:,1].min()*1.1, centers[:,1].max()*1.1) 
    ax.set_axis_off()
    fig.tight_layout()
    fig.show()
    
    
def histoParcel(lY, lX, X, iA, AP, AT, 
                m=-1, order=None, maxs=None, gauss=None, tend=0, color=None, figsize=(5,5)):
    
    # Parameters
    L = len(lY[:m])
    if color == None:
        color = rand(L, 3)
    if order not in [None, False, 0]:
        if order in [True, 1]:
            order = sortByIrregularDimension(lY, inverse=1)[0][:m]
        LX = [lX[o] for o in order]
        LY = [lY[o] for o in order]
        IA = [iA[o] for o in order]
    else:
        order = None
        LX = lX[:m]
        LY = lY[:m]
        IA = iA[:m]
    
    # Stacked Histograms
    fig = figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    ax.set_xlabel('P', fontsize=20)
    ax.set_ylabel('Predominance', fontsize=20)
    ax.set_ylim((0, 1))
    fig.tight_layout()
    ax.tick_params(labelsize=15)
    ax.tick_params(labelsize=15)
    
    stackedHistograms(lY=LY, lX=LX, X=X, ax=ax, color=color)
        
    # Generate Global Patterns
    pattIs = zeros((L, len(AP[0][0])))
    for gloP in range(L):
        for pos in range(len(IA[gloP])):
            for patt in range(len(IA[gloP][pos])):
                
                pattIs[gloP] += AP [LX[gloP][pos]] [IA[gloP][pos][patt]]
                
                if tend:
                    pattIs[gloP] *= AT [LX[gloP][pos]] [IA[gloP][pos][patt]]
                    
            if gauss != None:
                pattIs[gloP] *= gaussian(arange(200), *gauss)[LX[gloP][pos]]
                    
    # Show Patterns
    Ims = initParcellation()
    for i in range(len(pattIs)):
        f = figure(figsize=figsize)
        hagmannParcellation(pattIs[i]/pattIs[i].max(), Ims=Ims, maxs=maxs, fig=f)
        #ax = f.add_subplot(11,9,50)
        ax=f.add_axes([0.445,0.455,0.1,0.1], polar=1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.patch.set_facecolor(color[i])
    
    
def imshow3D(dim3, X=None, Y=None, Z=None, zdir='z',
             fig=None, fs=(5,5), lbls=None, elev=None, azim=None, cbar=False, cmap='ocean', **kwa):
    '''Return a 3d imshow figure of dim3 array where its first dimension is the set of images.
    '''
    
    from mpl_toolkits.mplot3d import Axes3D
    import mpl_toolkits.mplot3d.art3d as art3d
    sh = dim3.shape

    if fig == None:
        fig = figure(figsize=fs)
        ax = fig.add_subplot(111, projection='3d')

    if X == None or Y == None:
        X, Y = meshgrid(linspace(0,1,sh[2]), linspace(0,1,sh[1]))
    elif X.ndim == 1 or Y.ndim==1:
        X, Y = meshgrid(X, Y)
        
    if Z == None:
        Z = linspace(0,1,sh[0])
       
    if type(cmap) == str:
        cmap = cm.get_cmap(cmap)
        
    for i in range(sh[0]):
        #p = ax.pcolor(X, Y, dim3[i].T, **kwa)
        #art3d.poly_collection_2d_to_3d(p, Z[i], zdir)
        p = ax.plot_surface(X, Y, Z[i]*ones((sh[2],sh[1])), rstride=1, cstride=1, facecolors=dim3[i], 
                            shade=False, vmin=dim3[i].min(), vmax=dim3[i].max())

    lims = array([[X.min(), X.max()], [Y.min(), Y.max()], [Z.min(), Z.max()]])
    if zdir == 'x':
        z = [2,1,0]
    elif zdir == 'y':
        z = [0,1,2]
    elif zdir == 'z':
        z = [0,2,1]
        
    ax.set_xlim(lims[z[0]])
    ax.set_zlim(lims[z[1]])
    ax.set_ylim(lims[z[2]])
    ax.view_init(elev,azim)
    
    if lbls != None:
        ax.set_xlabel(lbls[z[0]])
        ax.set_ylabel(lbls[z[2]])
        ax.set_zlabel(lbls[z[1]])
    
    if cbar:
        fig.colorbar(p)

    return fig


def nbChiHist(x):
    '''Return the number of bins recquired for a chi distribution.
    From Otnes, Enockson, Digital time series analysis, Wiley, 1972.
    '''
    return exp(0.626 + 0.4*log(len(x) - 1))


def tableShow(dataframe, cmap='Spectral', ind=None, uns=None, spe={}, nspe={},
              sort=False, ival=None, asc=None, fmt='%.2f', nan=None, maxs=None, justDF=False):
    from IPython.display import HTML
    
    if ind == None:
        return 'Need to specify index.'
    
    df = dataframe.copy()
    
    if len(spe) or len(nspe):        
        cols = df.columns.tolist()
        endcol = list(range(len(cols)))
        
        for k in spe.keys():
            for val in spe[k]:
                df = df[getattr(df, k).values == val]
            #if len(spe[k]) == 1:
            #endcol.remove(cols.index(k))
        for k in nspe.keys():
            for val in nspe[k]:
                df = df[getattr(df, k).values != val]
        if ival != None:
            endcol.remove(ival)
            endcol.append(ival)
        df = df[df.columns[endcol]]
    
    if uns != None and sort:
        if asc == None:
            asc = [1] * len(uns)
        df = df.sort(uns, inplace=0)#ascending=asc)
        
    cma = cm.get_cmap(cmap)
    val = array(df.val).astype(float)
    if nan != None:
        val = where(isnan(val), 0., val)
        
    if maxs == None:
        vco = 255 * cma( (val - val.min()) / (val.max() - val.min()) )
    else:
        vco = 255 * cma( (val - maxs[0]) / (maxs[1] - maxs[0]) )
        
    df.val = [("<div style='background-color:#%02x%02x%02x;' >"+fmt)%(c[0],c[1],c[2], v) for c,v in zip(vco,val)]
    
    df = df.set_index(ind)
    if uns != None:
        df = df.unstack(uns)
        df.columns = df.columns.droplevel(0)
    if nan != None:
        df = df.fillna(nan)
    df = df.sort(inplace=0)#ascending=asc)
    
    if justDF:
        return df
    else:
        return HTML(df.to_html(escape=False))


def circular(links, node_names=None, node_angles=None, node_colors=None, down=None, cmap='bone_r',
             facecolor=(0,0,0,0.), dpi=300, save=None, nwidth=2, con=None, avg=None, fig=None, sym=-1,
             node_edgecolor=(0,0,0,0.), fcol='white', falpha=0., cbar=1, node_linewidth=1, linewidth=1,
             hemi=True, rotation=-90, direction=-1, disp=True, vmin=None, vmax=None, **kwa):
    '''Return a circular connectivity plot using mne library.
    direction=-1 or None
    '''
    from mne.viz import circular_layout, plot_connectivity_circle
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning) 
    
    if node_names == None:
        node_names = arange(links.max()+1)

    if down != None:
        
        # If no renormalisation is asked
        if avg == None:
            avg = ones(86)
            
        indices = c_[down[links[0]], down[links[1]]].tolist()
        tmp, CON = [],[]
        
        # If connectivity strengh is given or not
        if con == None:
            con = ones((len(down),len(down)))
            
        for l in range(len(indices)):
            # indices i et j de 0 a 85 des ROIs sous forme ij = i,j for l links
            ij = indices[l]
            IJ = links[:,l]
            # density of links across the ROI in and out
            #we = 1./ (avg[ij[0]] * avg[ij[1]])
                #CON.append(1./ avg[ij].max())
            we = con[IJ[0],IJ[1]] / (avg[ij].prod())
            
            if ij not in tmp:
                tmp.append(ij)
                CON.append(we)
            else:
                CON[tmp.index(ij)] += we
        indices = array(tmp).T
        CON = array(CON)
        #node_names = [lbl[sorted(lbl.keys())[k]] for k in range(len(lbl))]
        
        # Left hemisphere from top left to bottom left
        # then right hemisphere from top right to bottom right
        #node_angles = c_[linspace(100,260,43), linspace(80,-80,43)].flatten()
        if node_angles == None:
            if hemi:
                node_angles = r_[linspace(100,260,43), 
                                linspace(-80,80,43)[::sym]].flatten()
            else:
                node_angles = (linspace(0,360,len(node_names)) -rotation)[::direction]
        
    else: 
        indices = c_[links, [0,0]]
        if con == None:
            con = r_[ones(indices.shape[1]-1), 0]
        CON = con
        #if lbl:
            #node_names = [lbl[sorted(lbl.keys())[k]] for k in range(len(lbl))]
            
        # Left hemisphere from top left to bottom left
        # then right hemisphere from top right to bottom right
        #node_angles = c_[linspace(100,260,len(node_names)/2.), linspace(80,-80,len(node_names)/2.)].flatten()
        if node_angles == None:
            if hemi:
                node_angles = r_[linspace(100,260,len(node_names)/2.), 
                                linspace(-80,80,len(node_names)/2.)[::sym]].flatten()
            else:
                node_angles = (linspace(0,360,len(node_names)) -rotation)[::direction]

    fig,_ = plot_connectivity_circle(CON, node_names, indices=indices, node_angles=node_angles,
                                     node_width=nwidth, colorbar=cbar, facecolor=facecolor,
                                     node_edgecolor=node_edgecolor, node_colors=node_colors,
                                     node_linewidth=node_linewidth, vmin=vmin, vmax=vmax, 
                                     linewidth=linewidth,
                                     textcolor='black', colormap=cmap, fig=fig, show=disp, **kwa)
    
    fig.patch.set_facecolor(fcol)
    fig.patch.set_alpha(falpha)
    if save:
        savefig(save, dpi=dpi)
    if disp:
        fig.show()
    else:
        close(fig)
        return fig
  
  
def square(links, node_names=None, down=None, cmap='bone_r', cbar=1, avg=None, 
             facecolor='white', dpi=300, save=None, nwidth=2.5, **kwa):
    '''Return a square connectivity plot.
    '''

    if down != None:
        
        # If no renormalisation is asked
        if avg == None:
            avg = ones(86)
            
        indices = c_[down[links[0]], down[links[1]]].tolist()
        tmp = []
        con = zeros((86,86))
        
        for l in range(len(indices)):
            # indices i et j de 0 a 85 des ROIs sous forme ij = i,j for l links
            ij = indices[l]
            # density of links across the ROI in and out
            #we = 1./ (avg[ij[0]] * avg[ij[1]])
                #con.append(1./ avg[ij].max())
            we = 1./ (avg[ij].prod())
            con[ij[0],ij[1]] += we
            con[ij[1],ij[0]] += we
        
    fig = figure(figsize=(7,6))
    ax = fig.add_subplot(111)
    ax.imshow(con, cmap='bone_r', interpolation='nearest')
        #con, node_names, indices=indices, node_angles=node_angles, node_width=nwidth,
                             #colorbar=cbar, facecolor=facecolor, textcolor='black', colormap=cmap, **kwa)
    
    if save:
        fig.savefig(save, dpi=dpi)
        
    
def linksToNet(links, N, alpha=0.01, ax=None, lbl=None, name=''):
    import networkx as ne
    G = ne.Graph()
    pos = {}
        #node_angles = c_[linspace(100,260,43), linspace(80,-80,43)].flatten()
    angles = deg2rad(r_[linspace(-20,-160,N/2.), linspace(20,160,N/2.)].flatten())
    #angles = deg2rad(c_[linspace(-20,-160,N/2.), linspace(20,160,N/2.)].flatten())
    for n in range(N):
        pos[n] = (sin(angles[n]), cos(angles[n]))
        G.add_node(n)

    for l in links.T:
        if l[0] != l[1]:
            G.add_edge(l[0], l[1], attr_dict={'alpha':alpha})
       
    if ax == None:
        fig = figure(figsize=(7,6))
        ax = fig.add_subplot(111)
        
    #ne.draw_networkx_edges(G, pos=pos, with_labels=False, alpha=alpha, ax=ax,
                           #edge_color=1.+rand(G.number_of_edges())*0.00001, 
                           #edge_cmap=speCmap('hot',r=(0.0,0.9), su={'color':'w', 'alpha':alpha}, so={'color':'w', 'alpha':alpha}))
    ne.draw_networkx_nodes(G, pos=pos, node_color='w', node_size=0.01, ax=ax)
    ne.draw_networkx_edges(G, pos=pos, alpha=alpha, ax=ax)
    
    ax.set_title(name + " (%i)"%len(links.T))
    ax.axis('off')
    ax.set_aspect('equal')
    cut = 1.05
    xmax = cut * max(xx for xx,yy in pos.values())
    ymax = cut * max(yy for xx,yy in pos.values())
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-ymax, ymax)    


def transpCmap(cmap, alpha=0.5):
    '''Return a cmap from a string or a actual cmap using a scalar or list of alphas.
    Liste de la forme : np.linspace(0, 1, cmap.N)
    '''
    from matplotlib.colors import ListedColormap
    if type(cmap) == str:
        cmap = cm.get_cmap(cmap)
    my_cmap = cmap(arange(cmap.N))
    my_cmap[:,-1] = alpha #np.linspace(0, 1, cmap.N)
    my_cmap = ListedColormap(my_cmap)
    return my_cmap

    
def speCmap(cmap, r=(0,1), su=None, so=None):
    '''Return a part of a colorbar using 'r' as a range of it in percent.
    su: set a color for values under vmin using a string 'k' or a dict {'color':'k', 'alpha':1.0}
    so: same for values over vmax.
    '''
    if type(cmap) == str:
        cmap = cm.get_cmap(cmap)
        
    cmaplist = [cmap(i) for i in range(cmap.N)]
    for i in range(int(r[1]*cmap.N), cmap.N)[::-1]:
        cmaplist.pop(i)
    for i in range(int(r[0]*cmap.N))[::-1]:
        cmaplist.pop(i)
        
    cmap = cmap.from_list('Custom cmap', cmaplist)
    
    if su != None:
        try:     cmap.set_under(**su)
        except : cmap.set_under(su)
    if so != None:
        try:     cmap.set_under(**so)
        except : cmap.set_under(so)
        
    return cmap


def betweenMap(cmaps, vmin ,vmax, err=1e-2):
    '''Return a cmap1 surrounded by two cmap2 using vmin and vmax for the proportions.
    vmin and vmax has to be in [0;1].
    err is for not having the maximum value of the color of the cmap2.
    The cmap is only limited to a single color for now.
    cmaps takes:
    - a string value cmap1=cmap2=cmap3(-"_r")
    - a list of two string cmap1, cmap2=cmap3
    - a list of three strings cmap1, cmap2, cmap3
    '''
    try:
        cmap1, cmap2, cmap3 = cmaps
    except:
        try:
            cmap1, cmap2 = cmaps
            cmap3 = cmap2
        except:
            cmap1 = cmaps
            cmap2 = cmap1
            if cmap1[-1] == 'r':
                cmap3 = cmap1[:-2]
            else:
                cmap3 = cmap1 + '_r'
    colors1 = cm.get_cmap(cmap2)(linspace(0., 0., int(256*vmin)))
    colors2 = cm.get_cmap(cmap1)(linspace(0., 1., ceil(256*(vmax+err-vmin)).astype(int)))
    colors3 = cm.get_cmap(cmap3)(linspace(0., 0., int(256*(1-vmax))))
    colors = vstack((colors1, colors2, colors3))
    return LinearSegmentedColormap.from_list('my_colormap', colors)


def figToHtml(fig):
    from IPython.core.pylabtools import print_figure 
    from base64 import b64encode
    return "<img src='data:image/png;base64,%s'>" % b64encode(print_figure(fig)).decode("utf-8")
   
   
def new_section(title):
    from IPython.display import HTML
    style = "text-align:center;background:#66aa33;padding:30px;color:#ffffff;font-size:3em;"
    return HTML('<div style="{}">{}</div>'.format(style, title))


def htmlTable(tab, ind=True, col=None, lign=None, fmt='%.2f', angle='', 
              html=True, transf=False, width=100, cmap=None, maxs=None):
    if transf:
        for i in range(len(tab)):
            for j in range(len(tab[0])):
                tab[i][j] = figToHtml(tab[i][j])
    
    if ind or col or lign:
        if col == None:
            col = arange(len(tab[0])).astype(str)
        if lign== None:
            lign = arange(len(tab)).astype(str)
    
    
    if cmap != None:
        cma = cm.get_cmap(cmap)
        vco = tab.copy()
        if maxs == None:
            vco = 255 * cma( (tab - tab.min()) / (tab.max() - tab.min()) )
        else:
            vco = 255 * cma( (tab - maxs[0]) / (maxs[1] - maxs[0]) )

        
    table = '<table style="width:%i%%">'%width
    
    if ind:
        table += '<tr>'
        table += '<td>'+ angle +'</td>'
        for txt in col:
            table += '<td>'+ txt +'</td>'
        table += '</tr>'
    
    for i in range(len(tab)):
        table += '<tr>'
        if ind:
            table += '<td>'+ lign[i] +'</td>'
        
        for j in range(len(tab[0])):
            if cmap != None:
                try:
                    table += '<td>'+ "<div style='background-color:#%02x%02x%02x;' >" %(vco[i][j][0],vco[i][j][1],vco[i][j][2]) + fmt % tab[i][j] +'</td>'
                except:
                    import pdb; pdb.set_trace()
            else:
                table += '<td>'+ fmt % tab[i][j] +'</td>'
        table += '</tr>'
    table += '</tr></table>'
    
    if html:
        from IPython.display import HTML
        return HTML(table)
    else:
        return table
    
    
def plots(data, time=None, fig=None, ax=None, crowdf=0.7, disp=False, fs=None,
          ticklocs=None, xlim=None, xticks=None, xlabel=None, **kwa):
    '''Plots the time series "data" on different rows.
    The dimensions of "data"  have to be [Time x Nodes].
    '''
    from pylab import subplots, hstack, c_, newaxis
    from matplotlib.collections import LineCollection
    if ax == None:
        if fig == None: fig, ax = subplots(1, figsize=fs)
        else:           ax = fig.add_subplot(111)
            
    numRows = len(data[0])
    dmin, dmax = data.min(), data.max()
    dr = (dmax - dmin) * crowdf
    ym = (numRows - 1) * dr + dmax
    if time == None:
        time = arange(len(data))
            
    if ticklocs == None:
        ticklocs = []
        for i in range(numRows):
            ticklocs.append(i * dr)
    
    segs = []
    for i in range(numRows):
        segs.append(hstack((time[:, newaxis], data[:, i, newaxis])))
        
    lines = LineCollection(segs, offsets=c_[zeros(numRows), ticklocs], transOffset=None, **kwa)
    
    ax.add_collection(lines)
    ax.set_yticks(ticklocs)
    ax.set_yticklabels(list(range(numRows)))
    ax.set_ylim(dmin, ym)
    if xlim != None: ax.set_xlim(xlim)
    else:            ax.set_xlim(time.min(), time.max())
    if xticks != None: ax.set_xticks(xticks)
    if xlabel != None: ax.set_xlabel(xlabel)
    fig.tight_layout()
    if disp:
        fig.show()
    else:
        close(fig)
        return fig
        

def waves(data, time=None, fig=None, ax=None, crowdf=4, cutoff=0,
          colors=["#D3D820", "#C9CC54", "#D7DA66", "#FDFE42"], 
          disp=False, fs=None, bgcolor='black', axes=False,
          ticklocs=None, xlim=None, xticks=None, xlabel=None, **kwa):
    '''Plots the time series "data" on different rows as waves with cutoff.
    The dimensions of "data" have to be [Time x Nodes].
    '''
    from itertools import cycle
    from pylab import figure, fill_between
    if ax == None:
        if fig == None: 
            fig = figure(figsize=fs)
        ax = fig.add_subplot(111, axisbg=bgcolor)

    ax.xaxis.set_visible(axes)
    ax.yaxis.set_visible(axes)
    face_colors = cycle(colors)

    k = data.shape[-1]
    if time == None:
        time = linspace(0, 1, data.shape[0])
    TC = data.T.copy()

    displace = TC.max() / crowdf
    TC[TC<cutoff] = None     # Add a cutoff

    for n,y in enumerate(TC):
        # Vertically displace each plot
        y0 = ones(y.shape) * n * displace
        y1 = y + n*displace

        ax.fill_between(time, y0, y1, lw=1, 
                     facecolor=face_colors.next(), zorder=len(TC)-n)
    if disp:
        fig.show()
    else:
        close(fig)
        return fig
    
    
def boundaryCmap(cols=['beige','seagreen','darkorange','firebrick'], bounds=[0,0.01,0.05,0.1,1], alphas=[1,1,1,1], dic=True):
    '''Return a cmap and norm.
    '''
    from matplotlib import colors
    crgbas = []
    for ic, c in enumerate(cols):
        chex = colors.cnames[c]
        crgbas.append( list(colors.hex2color(chex)) + [alphas[ic]] )
    cmap = colors.ListedColormap(crgbas)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    if dic:
        return {'cmap':cmap, 'norm':norm}
    else:
        return cmap, norm
    
    
cmapStr = [['BrBG','bwr','coolwarm','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','seismic'],
           ['afmhot','autumn','bone','cool','copper','gist_heat','gray','hot','pink','spring','summer','winter'],
           ['Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples',
            'RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd'],
           ['Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3'],
           ['gist_earth','terrain','ocean','gist_stern','brg','CMRmap','cubehelix','gnuplot','gnuplot2','gist_ncar',
            'nipy_spectral','jet','rainbow','gist_rainbow','hsv','flag','prism']]
