from os import path
from Tools.ext import humanSize, array2data, data2array
from numpy import array
  
    
class database(dict):
    '''Create a database using the base of dictionaries.
    The goal is to have adress and description of all the static objets used in projects
    and be able to find and display them easily.
    '''
    defaultDir = '/home/golos/Main/Tools/DB.npy'
    
    def __init__(_, new=False, size=False, dir=defaultDir):
        _.dir = dir
        if not new:
            from Tools.ext import data2array
            _.update(data2array(_.dir, dic=True))
            #_.delete(adress=['CFSS','Norm'])
            if size: 
                print _.size()
            

    def add(_, name, adress, project='None', description='None', type='None', **kwa):
        '''Add a generic entry for an object.
        '''
        _[name] = {'name': name,
                   'project': project,
                   'adress': adress, 
                   'description': description,
                   'type': type,
                   'ftype': adress.split('.')[-1]}
        _[name].update(kwa)
        
        try:    _[name].update({'size': humanSize(path.getsize(adress))})
        except: pass
        try:    _[name].update({'dir': '/'.join(_[name]['adress'].split('/')[:-1])})   
        except: pass


    def find(_, subk=None, **kwa):
        '''Find the different entries using a dictionary.
        The value entered can be part of the real description.
        subk hav to be a list.
        Optimisable : premieres boucles sur kwa en decrementant nruter.
        '''
        nruter = {}
        for k in _.keys():
            try:
                cop = True
                for attr in kwa.keys():
                    
                    if type(kwa[attr]) == list:
                        for aa in kwa[attr]:
                            
                            if aa[0] == '!':
                                if aa[1:].lower() in _[k][attr].lower():
                                    cop = False
                            elif aa.lower() not in _[k][attr].lower():
                                cop = False
                        
                    elif type(kwa[attr]) == str:
                        if kwa[attr][0] == '!':
                            if kwa[attr][1:].lower() in _[k][attr].lower():
                                cop = False
                        elif kwa[attr].lower() not in _[k][attr].lower():
                            cop = False

                if cop:
                    if type(subk) == list:
                        nruter[k] = {sk:_[k][sk] for sk in subk if sk in _[k]}
                    elif type(subk) == str:
                        nruter[k] = _[k][subk]
                    else:
                        nruter[k] = dict(_[k])
            except:
                pass
        return nruter


    def save(_, dir_=None):
        '''Save the database in the directory .dir by default or specified during the initialization.
        If dir_ is specified, it change the directory of the database while saving it.
        '''
        try:
            size = humanSize(path.getsize(_.dir))
            if dir_: 
                _.dir = dir_
            array2data(_, _.dir)
            print 'Saved here : %s (%s, %i items)' %(_.dir, size, _.len())
            
        except:
            print 'Problem with the adresses.'
            return None
        
        
    def delete(_, **kwa):
        '''Remove the objects using "find" function.
        '''
        items = _.find(**kwa)
        for k in items.keys():
            del _[k]
            
            
    def elementsDiff(_, *elements, **kwa):
        '''Return the different kinds of elements as lists. 
        For exemple: 
        db.elementsDiff("project","type") returns the two lists of the different projects and the different file types.
        '''
        L = len(elements)
        items = _.find(**kwa)
        nruter = [[] for i in range(L)]
        for l in range(L):
            for k in items.keys():
                try:
                    tmp = items[k][elements[l]]

                    if tmp not in nruter[l]:
                        nruter[l].append(tmp)
                except:
                    pass
                    
        if L == 1 : return nruter[0]
        else:       return nruter
        
        
    def display(_, titles=None, dispOpt={}, transp=False, sGM=False, **kwa):
        '''Display from finding function using display tools.
        Profites of the loading function to add the files shape in the database.
        '''
        from Tools.display import mapMatrices
        items = _.find(**kwa)
        mats = []
        if titles == None:
            titles = []
            upti = True
        else:
            upti = False
            
        for k in sorted(items.keys()):
            tmp = data2array(items[k]['adress'])
            
            try:    _[k]['shape']
            except: _[k]['shape'] = tmp.shape
            if sGM:
                tmp -= tmp.mean()
                
            if transp: mats.append(tmp.T)
            else:      mats.append(tmp)
            
            if upti:
                titles.append( items[k]['name'] )  
        if upti:
            titles = dimTitles(titles)
            
        mapMatrices(mats, lTitl= titles, **dispOpt)
        
        
    def size(_, items=None, **kwa):
        '''Print the total size of the files found by options.
        items can be a dictionnary obtained by .find function.
        '''
        if items == None:
            items = _.find(**kwa)
        sumSize = 0
        for k in sorted(items.keys()):
            sumSize += path.getsize(items[k]['adress'])
        return humanSize(sumSize) + ' for %i items'%len(items.keys())
        
        
    def len(_):
        return len(_.keys())
    
    
    def merge(_, dbase):
        '''Merge two databases.
        Not implemented yet.
        '''
        pass
    
        
    def load(_, transp=False, slices=slice(None), mmap_mode=None, sprint=True, **kwa):
        '''Load files to a dictionnary.
        mmap_mode='r' for memory access.
        '''
        items = _.find(subk=['adress'], **kwa)
        mats, titles = [], []
        for k in sorted(items.keys()):
            if slices != None:
                tmp = array(data2array(items[k]['adress'], mmap_mode='r')[slices])
            else:
                tmp = data2array(items[k]['adress'], mmap_mode=mmap_mode)
            
            try:    _[k]['shape']
            except: _[k]['shape'] = tmp.shape
                
            if transp: mats.append(tmp.T)
            else:      mats.append(tmp)
                
            titles.append( k )  
            
        nruter = {}
        for k in range(len(mats)):
            nruter[titles[k]] = mats[k]
            
        if sprint: print _.size(items)
        return nruter
    
    
    def commonAdress(_, rmBeg=None, **kwa):
        '''Find the common adress, spliting it where the parameters enters.
        rmBeg(have to be negative) remove the end of the begining part.
        '''    
        items = _.find(subk='adress', **kwa)
        begin = [dimTitles(items.values(), beg=True) [:rmBeg]]
        follow = dimTitles(items.values())[0]
        interm = follow.split('_')[1::2]
        for i in range(len(interm)):
            interm[i] = '_'+ interm[i] +'_'
        end = ['.' + items.values()[0].split('.')[-1]]

        return begin + interm + end
        
    
    def test(_):
        '''Test the function by printing an example and a research on it.
        '''
        db = database(new=True)
        db.add('1', project='a')
        db.add('2', adress='test', project='a')
        db.add('3', project='b')
        print 'Objet'
        print db
        print 'Research (.find(project="a", adress="test")'
        print db.find(project='a', adress='test')
        

def dimTitles(liste, beg=False):
    '''Find the identical begining string of each objects of the list.
    Remove it or save the other part (if end is True) in a new list.
    "liste" needs to have at least two elements.
    '''
    i = 0
    while liste[0][:i+1] in liste[1][:i+1]:
        i += 1
        
    for k in range(1, len(liste)):
        while liste[k][:i] not in liste[0][:i]:
            i -= 1
        
    if beg:
        return liste[0][:i]
    else:
        if liste[0].split('.')[-1] in ['txt','npy','dat']:
            return ['.'.join(liste[k][i:].split('.')[:-1]) for k in range(len(liste))]
        else:
            return [liste[k][i:] for k in range(len(liste))]
    
    
def commonString(liste=[''], fileType=True):
    '''Find the common string.
    '''    
    begin = [dimTitles(liste, beg=True)]
    follow = dimTitles(liste)[0]
    interm = follow.split('_')[1::2]
    for i in range(len(interm)):
        interm[i] = '_'+ interm[i] +'_'
    
    if fileType:
        end = ['.' + liste[0].split('.')[-1]]
    else:
        end = ['']
    
    return begin + interm + end
