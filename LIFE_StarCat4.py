#StarCat4 code
import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables
import importlib #reloading external functions after modification
import math


#self created modules
import helperfunctions as hf
importlib.reload(hf)#reload module after changing it

def StarCat_creation(plx_cut,stars=[],details=False,query=True):

    def query(link,query,catalogs=[]):
        """
        Performs a query via TAP on the service given in the link parameter.
        If a list of tables is given in the catalogs parameter,
        those are uploaded to the service.
        :param link: Service access URL.
        :param query: Query to be asked of the external database service in ADQL.
        :param catalogs: List of astropy tables to be uploaded to the service.
        :return cat: Astropy table containing the result of the query.
        """
        #defining the vo service using the given link
        service = vo.dal.TAPService(link)
        #without upload tables
        if catalogs==[]:
            result=service.run_async(query.format(**locals()), maxrec=160000)
        #with upload tables
        else:
            tables={}
            for i in range(len(catalogs)):
                tables.update({f"t{i+1}":catalogs[i]})
            result = service.run_async(query,uploads=tables,timeout=None,
                                       maxrec=160000)
        cat=result.to_table()
        return cat
    
    #ra, dec, plx, distance, name, sptype, coo_gal_l, coo_gal_b, Teff, R, M, sep_phys, binary_flag, mag_i, mag_j
    adql_query="""
    SELECT o.main_id, sb.coo_ra, sb.coo_dec, sb.plx_value, sb.dist_st_value, sb.sptype_string, sb.coo_gal_l, 
    sb.coo_gal_b, sb.teff_st_value, sb.mass_st_value, sb.radius_st_value, sb.binary_flag,
    sb.mag_i_value, sb.mag_j_value,  sb.class_lum, sb.class_temp, o_parent.main_id AS parent_main_id, sb_parent.sep_ang_value
    FROM life_td.star_basic AS sb
    JOIN life_td.object AS o ON sb.object_idref=o.object_id
    LEFT JOIN life_td.h_link AS h ON o.object_id=h.child_object_idref
    LEFT JOIN life_td.object AS o_parent ON h.parent_object_idref=o_parent.object_id
    LEFT JOIN life_td.star_basic AS sb_parent ON o_parent.object_id=sb_parent.object_idref
    WHERE o.type = 'st' AND sb.dist_st_value < """+str(plx_cut) 
    #we are only interested in object type stars, up to a distance cut and well defined luminocity class (to sort out objects not around main sequence)
    service='http://dc.zah.uni-heidelberg.de/tap'
    
    if query==True:
        print('Querying life_td...')
        cat=query(service,adql_query)
        hf.save([cat],['StarCat4_raw'])
    else:
        print('Loading StarCat4_raw...')
        [cat]=hf.load(['StarCat4_raw'])
    
    print('Criteria: object type must be star, distance must be smaller than',plx_cut)
    if stars!=[]:
        stars=np.array(stars)
        stars=hf.object_contained(stars,cat['main_id'],details)
        
    print('Removing objects which no temperature class which in life_td means not in OBAFGKM...')
    cat=cat[np.where(cat['class_temp']!='')]
    cat.remove_rows(cat['mass_st_value'].mask.nonzero()[0])
    
    if len(stars)>0:
        stars=hf.object_contained(stars,cat['main_id'],details)
    
    print('Removing objects where luminocity class not IV,V,VI or NaN (where we will assume V)...')
    #print(cat.group_by('class_lum').groups.keys)
    cat=cat[np.where(cat['class_lum']!='I')]
    cat=cat[np.where(cat['class_lum']!='II')]
    cat=cat[np.where(cat['class_lum']!='III')]
    #print(cat.group_by('class_lum').groups.keys)
    
    #separating sample
    singles=cat[np.where(cat['binary_flag']=='False')]#2254
    multiples=cat[np.where(cat['binary_flag']=='True')]#1177

    print('Removing higher order multiples (objects where parent object is a child object as well)')
    #query h_link
    adql_query2="""
    SELECT o.main_id as child_main_id,o.object_id
    FROM life_td.object AS o
    JOIN life_td.h_link AS h on o.object_id=h.child_object_idref
    """
    if query==True:
        print('Querying life_td for h_link...')
        h_link=query(service,adql_query2)
        hf.save([cat],['h_link'])
    else:
        print('Loading h_link...')
        [h_link]=hf.load(['best_h_link'])
    #remove those where parent_main_id in h_link as child_main_id
    higher_order_multiples=np.in1d(multiples['parent_main_id'],h_link['child_main_id'])
    multiples=multiples[np.invert(higher_order_multiples)]#all without higher order multiples
    
    if len(stars)>0:
        stars=hf.object_contained(stars,ap.table.vstack([singles,multiples])['main_id'],details)
    
    print('Removing higher order multiples (objects where multiple parent object given or more than two components)')
    double=[]
    grouped=multiples.group_by('main_id')
    ind=grouped.groups.indices
    for i in range(len(ind)-1):
        if ind[i+1]-ind[i]!=1:
        #print(ind[i],ind[i+1])
            double.append(grouped['main_id'][ind[i]])
    multiples=grouped[np.where(np.invert(np.in1d(grouped['main_id'],double)))]

    if len(stars)>0:
        stars=hf.object_contained(stars,ap.table.vstack([singles,multiples])['main_id'],details)
        
    #critical sep:
    print('Removing those binaries where no separation value given...')
    multiples=multiples[np.where(multiples['sep_ang_value'].mask==False)].copy()#180
    
    if len(stars)>0:
        stars=hf.object_contained(stars,ap.table.vstack([singles,multiples])['main_id'],details)
    
    print('Transforming from angular separation into physical one...') 
    #transform into physical separations the remaining ones

    multiples['sep_phys_value']=multiples['sep_ang_value']
    multiples['sep_phys_value'].unit=ap.units.AU
    for i in range(len(multiples)):
        multiples['sep_phys_value'][i]=np.round(multiples['sep_ang_value'][i]*multiples['dist_st_value'][i],1)
    
    
    print('Removing those binaries where not both components in the list...') 
    print('(meaning either higher order multiples or one of the components does not have the required spectral types)') 
    # meaning having the parameters we need for critical sep computation (main sequence & sep given)

    grouped_multiples=multiples.group_by('parent_main_id')
    ind=grouped_multiples.groups.indices
    #print(ind)
    result=grouped_multiples[:0].copy()

    for i in range(len(ind)-1):
        l=ind[i+1]-ind[i]
        if l==2:
            result.add_row(grouped_multiples[ind[i]])
            result.add_row(grouped_multiples[ind[i]+1])
    #print(result)#130
    
    if len(stars)>0:
        stars=hf.object_contained(stars,ap.table.vstack([singles,result])['main_id'],details)
    
    
    print('Removing those that are too close together to have stable planets < 10AU...')
    def crit_sep(eps,mu,a_bin):
        """
        Computes critical semimajor-axis for planet orbit stability in binary 
        system as described in Holman and Wiegert 1999.
        :param eps: Binary orbit excentricity.
        :param mu: mass fraction with mu=m_s/(m_p+m_s), with m_s the mass 
            of the star considered as perturbing binary companion and m_p the
            mass of the star the planet is orbiting.
        :param a_bin: semimajor-axis of the binary stars.
        :return a_crit_s: Critical separation beyond which a planet on a S-type 
            orbit (circumstellar) is not stable any more.
        :return a_crit_p: Critical separation below which a planet on a P-type 
            orbit (circumbinary) is not stable any more.
        """
        a_crit_s=(0.464-0.38*mu-0.631*eps+0.586*mu*eps+0.15*eps**2\
                  -0.198*mu*eps**2)*a_bin
        a_crit_p=(1.6+5.1*eps-2.22*eps**2+4.12*mu-4.27*eps*mu-5.09*mu**2\
                  +4.61*eps**2*mu**2)*a_bin
        return a_crit_s,a_crit_p

    
    result['a_crit_s']=result['sep_phys_value']

    for i in range(len(result)):
        m_p=result['mass_st_value'][i]
        if i % 2 == 0:
            m_s=result['mass_st_value'][i+1]
        else:
            m_s=result['mass_st_value'][i-1]
        mu=m_s/(m_p+m_s)
        result['a_crit_s'][i]=crit_sep(0,mu,result['sep_phys_value'][i])[0]
        #assumed circular orbit and sep_phys = a_bin

    final=result[:0].copy()
    #wait, didn't I already define this? -> was before removing some
    ind=result.group_by('parent_main_id').groups.indices
    a_max=10.
    
    for i in range(len(ind)-1):
        if a_max < min(result['a_crit_s'][ind[i]],result['a_crit_s'][ind[i]+1]):
            final.add_row(result[ind[i]])
            final.add_row(result[ind[i]+1])
    
    if len(stars)>0:
        stars=hf.object_contained(stars,ap.table.vstack([singles,final])['main_id'],details)
    
    #print(final)
    
    #does not look right. a_crit is 8 but stayed in sample??

    #what about the inner limit of 0.5??
    
    StarCat4=ap.table.vstack([singles,final])
    
    print('Adding architecture parameter')
    #flag any object whose declination is contained within the region between 
    #-(23.4+45)*sin(RA) and +(23.4+45)*sin(RA) with the object's RA in degrees.

    
    def ecliptic(ang,ra,dec):
        """
        This function computes if the given position is within a certain angle
        from the ecliptic.
        :param ang: Angle in degrees.
        :param ra: Array of right ascention in degrees.
        :param dec: Array of declination in degrees.
        :return flag: Array of flags.
        """
        ecliptic=(23.4)*np.sin(2*np.pi*ra/360)
        flag=['True' if dec[j] > -ang+ecliptic[j] and dec[j] < ang+ecliptic[j] else 'False' \
              for j in range(len(ra))]
        return flag

    StarCat4['ecliptic_pm45deg']=ecliptic(45,StarCat4['coo_ra'],StarCat4['coo_dec'])
    hf.save([StarCat4],['StarCat4'])
    
    if len(stars)>0:
        print('All criteria were met by:',len(stars))
        if details:
            print(stars)
    
    return StarCat4
    