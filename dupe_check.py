# Checking for duplicates
'''This script searches two datasets of two point features for replication of
data between the two.  The GNDB set is the control set, while the External set
will update source history on duplicated names, and conflate unique names
to the features to which they belong, if there are any.
'''

# GNDB set from GNS or Geomedia Output? NOT GNS! Too many extraneous fields
# ...or Oracle PL/SQL?

__metaclass__ = type # still have to deal with old vs. new style classes in 2.6

from math import radians, cos, sin, asin, sqrt
import sys, os
# tested to see if load time for modules was any less if I didn't import all of arcpy
from arcpy import SearchCursor, Point,InsertCursor,UpdateCursor, env,\
     Exists, FeatureClassToGeodatabase_conversion,GetCount_management,\
     ListFields, PointGeometry, ListFields

#for exception handling
import traceback, datetime
import time

class aFactory:
    def newOne(self, name):
        '''
        initializes the second object for use in the comparison methods
        '''
        return featureClass(name)
        

class featureClass:
    """

    This is a dataset with the GNDB output schema. This can be initialized for any case in which
    such a dataset is needed.  NOTE: two instances of this class are needed for most of the
    methods to function properly.

    Inputs:

        --- name - to initialize, indicate the full path of the file to use

    Attributes:

        --- name - the file name (no path)

        --- path - the directory path to the file

        --- workspace - the location in the directory tree idenifying the directory scope of operations

        --- gdb - the ESRI geodatabase where outputs are to be stored

        --- namePairs - the pairs of names (and other info) found as a result of the proximity search
                    method

        --- insGNDB - a list of GNDB features to be edited, and therefore added to the set of external
                    names.  Initialized as an empty list but populated with the result of self.demote_gndb()

    Methods:

        --- getGNDBRecords - a generator function that iterates through the GNDB dataset and holds pertinent info from
                    one record at a time

                    input: name - the file name (with path) of a feature class

                    output: rec - a list of one record's relevant features:
                                        [UFI, FULL_NAME, LON, LAT, DSG, UNI, NT]

        --- handle_dsg - handler for designation code issues

                    inputs:
                        x = the designation as recorded in the external set
                        y = the designation as recorded in the GNDB set

                    outputs:
                        None, to signal if the designations are different, indicating that
                        the two features are unrelated and should not be picked up in the
                        list as potential duplicates

        --- proximitySearch - iterates through the matrix of both feature classes and performs the
                        Haversine formula to determine the distance, around a curved surface, between
                        two points given as geographic (LAT/LON) coordinates. Uses the above methods
                        in its definition

                    inputs:
                        name = the name of the GNDB feature class (the external class's info is already
                            indicated as attributes of the featureClass class)

                    outputs:
                        no return value, but builds the namePairs attribute (a list)

                    
            <still building>            
    """
    def __init__(self, name):
                
        self.name = os.path.basename(name)
        self.path = os.path.dirname(name)
        if 'gdb' or 'mdb' in name:
            
            self.workspace = os.path.dirname(os.path.dirname(name))
            self.gdb = os.path.dirname(name).split('/')[-1][:-4]
        else:
            self.workspace = os.path.dirname(name)
            self.gdb = os.path.dirname(name)
              
    def getGNDBRecords(self, name):
        ''' iterates through the GNDB dataset and holds pertinent info from one record
            at a time.
        '''
        try:      
            d0 = time.clock()
            sc = SearchCursor(name)
            row = sc.next()
            
            while row:
                rec = []
                rec.append(row.getValue('UFI'))
                rec.append(row.getValue('FULL_NAME'))
                rec.append(row.getValue('LON'))
                rec.append(row.getValue('LAT'))
                rec.append(row.getValue('DSG'))
                rec.append(row.getValue('UNI'))
                rec.append(row.getValue('NT'))

                yield rec
                row = sc.next()

        except Exception, e:
            print "Error."
            e = sys.exc_info()
            log_file = open("H:/Programming/PyErrorLog.txt", 'a')
            log_file.write(str(datetime.datetime.now().time()) + "\n")
            log_file.write(str(e[0]) + "\n")
            log_file.write(str(e[1]) + "\n")
            traceback.print_tb(e[2],None,log_file)
            log_file.write('\n')
            log_file.close()
           

    def handle_dsg(self,x,y):
        '''
        handler for designation code issues

        inputs:
            x = the designation as recorded in the external set
            y = the designation as recorded in the GNDB set

        outputs:
            None, to signal if the designations are different, indicating that
            the two features are unrelated and should not be picked up in the
            list as potential duplicates

        '''
        ### handle PPLAs
        if x == 'PPLA' and y[:3] == 'PPL':
            y = x

        ### handle MT vs HLL
        mt_dsgs = ['MT','HLL','PK','MTS']
        if x in mt_dsgs and y in mt_dsgs:
            y = x

        ### handle LKI vs LK
        elif x[:2] ==  'LK':
            y = x

        ### handle PPL vs FRM
        elif x == 'PPL' and y == 'FRM':
            y = x
            
        if x != y:

            return None
        else:
            return x,y

    # name2 = the name of the GNDB set    
    def proximitySearch(self, name2):
        try:
            proxFile = 'H:/proxFile.txt'
            pF = open(proxFile,'w')
            
            count= 1
            d0 = time.clock()
            self.name2 = os.path.basename(name2)
            self.path2 = os.path.dirname(name2)
            if self.path != self.path2:
                if Exists(self.path+'/'+self.name2) == False:
                    FeatureClassToGeodatabase_conversion(self.path2+'/'+\
                                                           self.name2,self.path)
                    
            factory = aFactory()
            self.obj2 = factory.newOne(self.name2)

            c1 = GetCount_management(self.path+'/'+self.name).getOutput(0)
            
            ### search cursor for external dataset
            
            sc = SearchCursor(self.path+'/'+self.name)
            rowZ = sc.next()
            self.namePairs = []
            while rowZ:
                # gather relevant info for each value in the external dataset
                if count == c1:     # once the cursor is at the end of the set...
                    break           # ...exit the loop
                ufi = rowZ.getValue('UFI')
                nym = rowZ.getValue('FULL_NAME')
                lat = rowZ.getValue('LAT')
                lon = rowZ.getValue('LON')
                dsg = rowZ.getValue('DSG')
                if rowZ.getValue('GENERIC'):
                    generic = True
                else:
                    generic = False
                nt = rowZ.getValue('NT')

                rec = self.getGNDBRecords(name2)
                r = rec.next()      # GNDB record: [UFI, Name, Lon, Lat, Dsg, UNI, NT]

                if count % 10 == 0:
                    print count

                while r:                    
                    try:
                        # don't consider proximate features if they're not the same type

                        ########################
                        # handle dsg codes separately
    
                        d = self.handle_dsg(dsg, r[4])
                        if d == None:
                            r = rec.next()
                            continue
                                                                       
                        #
                        ########################
                        
                        #following lines converts to radians
                        lonDB,latDB,lonEx,latEx = map(radians,[r[2],r[3],lon,lat])

                        # get X,Y pairs distance, in km
                        distLon = (lonDB - lonEx)
                        distLat = (latDB - latEx)
                        a = sin(distLat/2)**2+cos(latEx)*cos(latDB)*sin(distLon/2)**2
                        c = 2 * asin(sqrt(a))
                        km = 6367*c
                    
                        threshold = 5
                        pF.write('%s (%s): %s (%s) == %0.5f\n' % (r[1].strip().encode('utf8'),r[4].encode('utf8')\
                                                                  ,nym.strip().encode('utf8'),dsg.encode('utf8'),km))
                        if km < threshold:
                            # store GNDB UFI, GNDB Name, External UFI...,
                            # External Name, GNDB UNI, generic, NT
                            pairInfo = [r[0],r[1].strip(),ufi,nym.strip(),r[5],generic, nt]
                            self.namePairs.append(pairInfo)
                        r = rec.next()
                            
                    except StopIteration: 
                        rowZ = sc.next()
                        count += 1
                        break

            del sc
     
            print "Proximity search complete."
            pF.close()
            print len(self.namePairs),'pairs found within %d km.\n' % threshold
            d1 = time.clock()
            delta = d1-d0
            print "%s took %f s to complete." % (sys._getframe().f_code.co_name,\
                                               delta)
            return self.namePairs # returns list of (GNDBUFI, GNDBName, exUFI, exName, GNDB_UNI) tuples  
            
        except Exception, e:
            print "Error."
            e = sys.exc_info()
            log_file = open("H:/Programming/PyErrorLog.txt", 'a')
            log_file.write(str(datetime.datetime.now().time()) + "\n")
            log_file.write(str(e[0]) + "\n")
            log_file.write(str(e[1]) + "\n")
            traceback.print_tb(e[2],None,log_file)
            log_file.write('\n')
            log_file.close()

    def genericFilter(self,**kwargs):
        try:
            x = kwargs['name']
            x = x.split(',')

            if len(x[0])== 0:
                pass
            for p in [self.path]:
                uc = UpdateCursor(p)
    
                for row in uc:
                    fn = row.getValue('FULL_NAME')
                    for j in x:
                        if j in fn:
                            fn = fn.replace(j,'')
                        row.setValue('FULL_NAME', fn)
                    uc.updateRow(row)
            del uc

        except:
            print "Error in %s." % sys._getframe().f_code.co_name
            e = sys.exc_info()
            log_file = open("H:/Programming/PyErrorLog.txt", 'a')
            log_file.write(str(datetime.datetime.now().time()) + "\n")
            log_file.write(str(e[0]) + "\n")
            log_file.write(str(e[1]) + "\n")
            traceback.print_tb(e[2],None,log_file)
            log_file.write('\n')
            log_file.close()
            
    def LevDist(self, a):
        '''
        Take the UFI/Name pairs and determine the Levenshtein distance (LD)
        between them.

        input: "a" is the list of lists that each contain (gndbUFI, gndbName,
          exUFI, exName,generic flag[T/F]) where:
          
          exUFI = the ufi of the feature from the extrnal dataset
          exName = that feature's N name
          GNDBUFI = the ufi of the feature from the dataset exported from GNDB
          GNDBName = that feature's N Name
          generic = whether there is a generic in the name or not
          nt = name type

        output:
                output : tuple of 6 items: dup/var flag, gndb_ufi,gndb_name,exUFI,exName
                       gndb_uni
        '''

        try:
            LevFile = os.path.dirname(self.path)+'/LevDist.txt'
            f = open(LevFile, 'w')

            last_ufi = -1
            count = 0
            counter = 1
            for item in a:
                # check to see if the pair is from the same UFI as the last iteration; if not, (re)set
                # dupe_flag to 0. dupe_flag needed to make sure that duplicate find is not overwritten
                # by variant find
                if item[0] != last_ufi:
                    dupe_flag = 0
                # x is gndbName
                x = item[1].strip()
                # if compound, may contain generics & could mean no match
                #  - but excluding generics means no match if generics exist in
                # both names, so check both versions - with and without generic
                xList = []
                xList.append(x)         # append the full name
                if item[5] == True:               
                    x = x.split(' ')[0]     
                    xList.append(x)     # append just the pre-space, truncated name
                    
                # y is externalName
                y = item[3].strip()

                for x in xList:
                    
                    n, m = len(x), len(y)
                    if n > m:
                        # Make sure n <= m, to use O(min(n,m)) space
                        x,y = y,x
                        n,m = m,n
        
                    current = range(n+1)

                    for i in range(1,m+1):
                        previous, current = current, [i]+[0]*n
                        for j in range(1,n+1):
                            add, delete = previous[j]+1, current[j-1]+1
                            change = previous[j-1]
                            if x[j-1] != y[i-1]:
                                change = change + 1
                            current[j] = min(add, delete, change)
  
                # the new name is an exact match, and should be considered a dupe                        
                    if current[n] == 0:
                        # output dupeFlag,gUFI,gName,xUFI,xName,gUNI
                        it = ('DUPLICATE',item[0],x, item[2], y,item[4], item[6]) # flag the ufi
                        counter += 1
                        f.write("Feature pair %d: %d ==> %s, %s -- Distance %d.\n" \
                                % (counter, item[2],x.encode('utf8'), y.encode('utf8'), current[n]))
                        dupe_flag = 1
                        last_ufi = item[0]
                        yield it

                    elif current[n] > 3: # these are truly new names?
                        # don't write a new record to the dataset if a duplicate name for this ufi has
                        # already been found
                        if dupe_flag == 1:
                            continue
                        else:
                            counter += 1
                            f.write("Feature pair %d: %d ==> %s, %s -- Distance %d.\n" \
                                % (counter, item[2], x.encode('utf8'), y.encode('utf8'), current[n]))
                                    
                    else: # exName is close to GNDBName, and is considered a variant
                        # output varFlag,gUFI,gName,xUFI,xName,gUNI

                        # don't write a new record to the dataset if a duplicate name for this ufi has
                        # already been found
                        if dupe_flag == 1:
                            continue
                        it = ('VARIANT',item[0],x, item[2], y, item[4], item[6]) #flag exName as variant
                        counter += 1
                        # log the Lev Distance for each pair.
                        f.write("Feature pair %d: %d ==> %s, %s -- Distance %d.\n" \
                            % (counter, item[2], x.encode('utf8'), y.encode('utf8'), current[n]))
                        yield it
            
            f.close()    
        except Exception:
            e = sys.exc_info()
            log_file = open("H:/Programming/PyErrorLog.txt", 'a')
            log_file.write(str(datetime.datetime.now().time()) + "\n")
            log_file.write(str(e[0]) + "\n")
            log_file.write(str(e[1]) + "\n")
            traceback.print_tb(e[2],None,log_file)
            log_file.write('\n')
            log_file.close()

    def handle_variants(self,n,row):
        try:
            uni = 2
            
            row.setValue('MF','M')
            row.setValue('GEONAME_NOTE', n[0])
            row.setValue('UFI',n[1])
            row.setValue('UNI',uni)
            row.setValue('SOURCE_MOD_TYPE1','B')
            if row.getValue('NT') == 'NS':
                row.setValue('NT','NS')
            else:
                row.setValue('NT', 'N')
            return row

        except Exception:
            print 'Error handling generics.'
            e = sys.exc_info()
            log_file = open("H:/Programming/PyErrorLog.txt", 'a')
            log_file.write(str(datetime.datetime.now().time()) + "\n")
            log_file.write(str(e[0]) + "\n")
            log_file.write(str(e[1]) + "\n")
            traceback.print_tb(e[2],None,log_file)
            log_file.write('\n')
            log_file.close()            

    def demote_gndb_rec(self, gndb_ufi,gndb_uni):
        '''
        retrieves a gndb record to external dataset in order to demote its NT field

        inputs:
            gndb_ufi: have to retrieve the row before it can be processed or
            appended; a cursor will query and search

        outputs:
            the row object, which will be appended to the external feature class 
        '''
        try:
            sc = SearchCursor(self.path+'/'+self.name2,"\"UFI\" = %d" % gndb_ufi)

            row_vals = {}
            row = sc.next()
            while row:
                # skip variant names
                if row.getValue('NT')=='N':
                    for field in ListFields(self.path+'/'+self.name2)[1:-1]:
                        row_vals[field.name] = row.getValue(field.name)
                        row_vals["NT"]='V'
                        row_vals['USID1']='FILE_MAINTENANCE'
                                    
                row = sc.next()
            
            del sc
            return row_vals
        except Exception:
            print gndb_ufi, type(gndb_ufi)

    def updateFC(self):
        '''
        Updates the feature class based on the distances found in the other
        two functions

        input:
        lstNames - 'b' from the LevDist function
        '''
        try:
            d0 = time.clock()
            z = self.LevDist(self.namePairs)
            
            dupeCount = 0
            modCount = 0
            ctr = 0
            self.insGNDB = []
            while z:
                # n =  varFlag,gUFI,gName,xUFI,xName,gUNI
                n = z.next()
                       
                # cursor must apply to external names feature class; else will throw exception
                # as fields such as "MF" do not exist
                uc = UpdateCursor(self.path+'/'+self.name)
                
                gUFI = n[1]
                gN = n[2]
                xUFI = n[3]
                xN = n[4]
                gUNI = n[5]

                # initialize an empty list to put GNDB row objects for records to be demoted
                # from N to V NT
                
                for row in uc:
                    try:
                        uni = 1
                        # print row.getValue('FULL_NAME'), xN
                        uc_name = row.getValue('FULL_NAME')
                        uc_ufi = row.getValue('UFI')
                        if uc_name == xN and uc_ufi == xUFI:
                            if n[0] == 'DUPLICATE':
                                dupeCount += 1
                                row.setValue('GEONAME_NOTE',n[0])
                                # set MF column to M to capture coordinates
                                # or other feature mods
                                row.setValue('MF','M')

                                ##################
                                #
                                # need to change this with new sources
                                row.setValue('USID1','TU-GAZETTEER-09')
                                #
                                ##################
                                row.setValue('SOURCE_MOD_TYPE1','F')
                                row.setValue('UFI', gUFI)
                                row.setValue('UNI', gUNI)
                                uc.updateRow(row)
                                ctr += 1
                                # if exName is duplicate of GNDB Variant, find and
                                # demote GNDB N name
                                if n[6] == 'V':
                                    rr = self.demote_gndb_rec(n[1],gUNI)

                            elif n[0] == 'VARIANT':
                                # Turkey Gazetteer is considered authoritative
                                # so discrepancies favor the gazetteer
                                # handles the external record
                                vRow = self.handle_variants(n,row)
                                uc.updateRow(vRow)
                                # demote the GNDB NT to 'V'
                                rr = self.demote_gndb_rec(n[1], gUNI) # = gndb row
                                self.insGNDB.append(rr)
                                modCount +=1
                                ctr +=1

                            print n[1]

                    except StopIteration:
                        print 'uc done'
                
                        
        except StopIteration:
            print "ctr = %d" % ctr
            print "Feature class modified:\n\t%d duplicate features found \
                    \n\t%d variant names added" % (dupeCount,modCount)
            d1 = time.clock()
            delta = d1-d0
            print "%s took %f to complete." % (sys._getframe().f_code.co_name,\
                                               delta)

            ins_cur = InsertCursor(self.path+'/'+self.name)
            ins_num = 0
            for rec in self.insGNDB:
                if rec == None:
                    continue
                pt = Point()
                row = ins_cur.newRow()
                for k,d in rec.iteritems():
                    
                    row.setValue(k,d)
                    row.setValue('MF','M')
                    row.setValue('NT','V')
                    if k == 'LON':
                        pt.X = d
                    elif k == 'LAT':
                        pt.Y = d
                row.setNull('NAME_RANK')
                row.setValue('SOURCE_MOD_TYPE1','N')
                pt_geo = PointGeometry(pt)
                row.setValue('Shape',pt_geo)
                ins_cur.insertRow(row)
                ins_num += 1

            print "%d GNDB name records demoted and added to feature class." % \
                  ins_num
            
            del ins_cur

        except Exception, e:
            e = sys.exc_info()
            log_file = open("H:/Programming/PyErrorLog.txt", 'a')
            log_file.write(str(datetime.datetime.now().time()) + "\n")
            log_file.write(str(e[0]) + "\n")
            log_file.write(str(e[1]) + "\n")
            traceback.print_tb(e[2],None,log_file)
            log_file.write('\n')
            log_file.close()

    def fc_cleanup(self):
        import arcpy
        num=1
        sort_fc = '/'.join([x.path,x.name+'_sort'+num])
        while arcpy.Exists(sort_fc) == True:
            num += 1
            sort_fc = sort_fc.replace(sort_fc[-1],num)                     
        arcpy.Sort_management(x.path+'/'+x.name,x.path+'/'+x.name+'_sort',[["UFI","ASCENDING"]])                        
        uc = UpdateCursor(sort_fc)
        ufi = -1
        row = uc.next()
        feat_dict = {}
        lst_unis = []
        for row in uc:
            last_ufi,ufi = ufi,row.getValue('UFI')
            nt = row.getValue('NT')
            if ufi == last_ufi and nt in ('V',u'V'):
                uni = row.getValue('UNI')
                lst_unis.append(row.getValue('UNI'))
            else:
                feat_dict[ufi] = lst_unis
                lst_unis = []
            uc.updateRow(row)
            row = uc.next()
        del uc
        num += 1
        uc2 = UpdateCursor(sort_fc)
        row = uc.next()
        for row in uc2:
            ufi = row.getValue('UFI')
            uni = row.getValue('UNI')
            if ufi in feat_dict:
                lst_unis = feat_dict[ufi]
                if uni in lst_unis:
                    uc2.deleteRow(row)       
            
if __name__ == '__main__':
    x = featureClass('H:/Turkey_Gazetteer.gdb/Lurkey_39N__425E')
    y = x.proximitySearch('H:/Turkey_Gazetteer.gdb/GNDB_39N__42_5E')
    x2 = x.obj2
    u = x.updateFC()
    x.fc_cleanup()
