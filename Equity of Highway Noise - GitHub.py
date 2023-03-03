root = r""

#root should contain:
# State_median_housing values.csv    [must have 2 columns labeled 'State' and 'Median Housing Value', download here: https://factfinder.census.gov/faces/nav/jsf/pages/searchresults.xhtml?refresh=t#none]
# ACS 2019    [folder with unzipped .gdb for each state, download here: https://www2.census.gov/geo/tiger/TIGER_DP/]
# state shp   [folder with unzipped .shp and associated files for each state. This program uses a combination of preprocessed HPMS 2016 and 2017, download here: https://www.fhwa.dot.gov/policyinformation/hpms/shapefiles_2017.cfm]

### program will output 3 csv files into the root folder


import geopandas as gpd
import time
import pandas
from math import log10
import gc
import os
import geopy.distance
import warnings
import math
import utm
from pyproj import CRS
warnings.filterwarnings('ignore')  #this is set to 'ignore' due to a section with frequent pandas merging that draws a warning


house_value = pandas.read_csv(root+"\\State_median_housing values.csv")
#change in hv per dBA from indirect valuation, Deluchi and Hsu (1998) and median housing value see methodology section in paper
hv = {}
for i in range(len(house_value)):
    #annual change per dBA by state = median housing value * annualitzation rate * near road adjustment factor * % reduction per dBA
    hv[house_value.iloc[i]['State']] = house_value.iloc[i]['Median Housing Value']*.0650*.95*.0085

files = os.listdir(root+r"\ACS 2019")


#This list is used for reference between state FIPS, ABB, name, and UTM code
STATES_DATA = [['01', 'AL', 'Alabama', '16'],
               ['02', 'AK', 'Alaska', '5'],
               ['04', 'AZ', 'Arizona', '12'],
               ['05', 'AR', 'Arkansas', '15'],
               ['06', 'CA', 'California', '11'],
               ['08', 'CO', 'Colorado', '13'],
               ['09', 'CT', 'Connecticut', '18'],
               ['10', 'DE', 'Delaware', '18'],
               ['11', 'DC', 'District of Columbia', '18'],
               ['12', 'FL', 'Florida', '17'],
               ['13', 'GA', 'Georgia', '17'],
               ['15', 'HI', 'Hawaii', '4'],
               ['16', 'ID', 'Idaho', '11'],
               ['17', 'IL', 'Illinois', '16'],
               ['18', 'IN', 'Indiana', '16'],
               ['19', 'IA', 'Iowa', '15'],
               ['20', 'KS', 'Kansas', '14'],
               ['21', 'KY', 'Kentucky', '16'],
               ['22', 'LA', 'Louisiana', '15'],
               ['23', 'ME', 'Maine', '19'],
               ['24', 'MD', 'Maryland', '18'],
               ['25', 'MA', 'Massachusetts', '19'],
               ['26', 'MI', 'Michigan', '16'],
               ['27', 'MN', 'Minnesota', '15'],
               ['28', 'MS', 'Mississippi', '16'],
               ['29', 'MO', 'Missouri', '15'],
               ['30', 'MT', 'Montana', '12'],
               ['31', 'NE', 'Nebraska', '14'],
               ['32', 'NV', 'Nevada', '11'],
               ['33', 'NH', 'New Hampshire', '19'],
               ['34', 'NJ', 'New Jersey', '18'],
               ['35', 'NM', 'New Mexico', '13'],
               ['36', 'NY', 'New York', '18'],
               ['37', 'NC', 'North Carolina', '17'],
               ['38', 'ND', 'North Dakota', '14'],
               ['39', 'OH', 'Ohio', '17'],
               ['40', 'OK', 'Oklahoma', '14'],
               ['41', 'OR', 'Oregon', '10'],
               ['42', 'PA', 'Pennsylvania', '18'],
               #['72', 'PR', 'Puerto Rico', '19'],
               ['44', 'RI', 'Rhode Island', '19'],
               ['45', 'SC', 'South Carolina', '17'],
               ['46', 'SD', 'South Dakota', '14'],
               ['47', 'TN', 'Tennessee', '16'],
               ['48', 'TX', 'Texas', '14'],
               ['49', 'UT', 'Utah', '12'],
               ['50', 'VT', 'Vermont', '18'],
               ['51', 'VA', 'Virginia', '17'],
               ['53', 'WA', 'Washington', '10'],
               ['54', 'WV', 'West Virginia', '17'],
               ['55', 'WI', 'Wisconsin', '16'],
               ['56', 'WY', 'Wyoming', '13']]
# s = pandas.DataFrame(STATES_DATA,columns=['FIPS','ST'])
# s.to_csv(r"")

#vehicle class distribution from FHWA, see Tehcnical Report p. 106 Federal Highway Administration (FHWA), Vehicle Classification and Weight Distribution System, 2017, unpublished.
#this is used to distribute the traffic counts from passenger vehicle, single unit truck, and combination truck into 20 different vehicle classes based on weight and axle configuration
veh_class = {'single': [0.798344485,0.172881753,0.028773762],
            'combination': [0.012392679,0.050868298,0.737375027,0.084669144,0.023194662,0.004540104,0.0000250647895381742,0.000942812,0.008944798,0.03569885,0.01332327,0.000131528,0.000276449,0.006229683,0.021387625]}
#NPCE adjustment calculated as ratio of composite noise level for a given vehicle class to composite noise level for passenger cars. Concept from Hokanson (1981), see Technical Report p. 24-25 for more details
#this is used to convert the traffic count for different vehicle classes for each roadway into a single count for easier noise calculation
speed_npce = {'single':{10:[130.2,230.2,364.3],
                        15:[67.24776174,118.886678,188.1257003],
                        25:[29.24585834,51.70347465,81.81532649],
                        35:[16.89957076,29.87659033,47.27657103],
                        45:[11.21941438,19.83469592,31.38632623],
                       55:[8.08938308,14.30114337,22.63005962],
                       65:[6.161093716,10.89213895,17.23566763],
                       75:[6.161093716,10.89213895,17.23566763],
                       85:[6.161093716,10.89213895,17.23566763]}, 
              'combination':{10:[234.0524646,319.4180744,439.586625,445.1743278,538.6709547,628.7931163,344.6212488,437.8395829,593.9066438,442.8460081,539.8044702,620.0304067,720.0891769,833.025738,234.148796],
                             15:[120.8605222,164.941802,226.9947003,229.8800904,278.1600826,324.697561,177.9562722,226.0925589,306.6828082,228.6777876,278.7454098,320.1726539,371.8412199,430.1596477,120.910266],
                             25:[52.56189381,71.73271574,98.71934288,99.97419076,120.9710206,141.2100362,77.3926714,98.32700418,133.375472,99.45131273,121.2255777,139.2421671,161.712678,187.0751947,52.58352725],
                            35:[30.37262348,41.4503856,57.04447108,57.76957856,69.90249009,81.59750251,44.72096224,56.81776015,77.07033936,57.46743614,70.04958463,80.46037934,93.44484998,108.1004515,30.38512425],
                            45:[20.16400615,27.51839434,37.87111332,38.35250314,46.40739187,54.17156484,29.68968941,37.72060276,51.16603766,38.15191455,46.5050461,53.41664295,62.03687116,71.76654231,20.17230526],
                            55:[14.53858148,19.84121684,27.30569822,27.65278823,33.46049605,39.058593,21.40675645,27.19717762,36.89155825,27.50816054,33.53090636,38.51428184,44.72960876,51.74486237,14.54456527],
                            65:[11.07297826,15.11160927,20.79676089,21.06111407,25.48442197,29.74808456,16.30396673,20.71410866,28.09761208,20.9509617,25.53804838,29.33352241,34.06728408,39.41029164,11.07753568],
                            75:[11.07297826,15.11160927,20.79676089,21.06111407,25.48442197,29.74808456,16.30396673,20.71410866,28.09761208,20.9509617,25.53804838,29.33352241,34.06728408,39.41029164,11.07753568],
                            85:[11.07297826,15.11160927,20.79676089,21.06111407,25.48442197,29.74808456,16.30396673,20.71410866,28.09761208,20.9509617,25.53804838,29.33352241,34.06728408,39.41029164,11.07753568]}}
#Dictionary for referring to ACS variables more easily
ACS = {'white':'B02001e2','black':'B02001e3','native':'B02001e4',
      'asian':'B02001e5','pacific':'B02001e6','other':'B02001e7',
      'poverty':'B17021e2','nonwhite':'nonwhite','nonpoverty':'nonpoverty'}
#      'renter':'B25008e3','nonrenter':'nonrenter'}
#B25008e3	TOTAL POPULATION IN OCCUPIED HOUSING UNITS BY TENURE: Total: Renter occupied: Total population in occupied housing units -- (Estimate)
demos = list(ACS.keys())   

# #FIPS <-> names only needed for Tableau plotting
# county_names = pandas.read_csv(root+"\\FIPS names.csv")
# county_names['county_FIPS'] = county_names.apply(lambda row: str(row['county_FIPS']), axis=1 )
# county_names['county_FIPS'] = county_names.apply(lambda row: '0'+row['county_FIPS'] 
#                                                  if len(row['county_FIPS'])<5 
#                                                  else row['county_FIPS'], axis=1 )


#converts AADT for passenger, combination, and single trucks for a road
#into number of passenger car equivalents (NPCE) based on speed on road
def npce(aadt_c,aadt_s,aadt_p,speed):
    
    #Match the speed in HPMS to the closest speed in the npce dictionary
    key = [i for i in speed_npce['single'].keys()]
    key_dif = [abs(i-speed) for i in key]
    speed = key[key_dif.index(min(key_dif))]
    
    aadt_c_npce = [aadt_c*i*j for i,j in 
                   zip(veh_class['combination'],speed_npce['combination'][speed])]
    aadt_s_npce = [aadt_s*i*j for i,j in 
                   zip(veh_class['single'],speed_npce['single'][speed])]
    aadt_npce = sum([aadt_p]+aadt_c_npce+aadt_s_npce)
    return(aadt_npce)


#calculate the distance each noise_lvl reaches  (from the ECAT excel tool Noise!AM29)
#this equation is derived from Haling and Cohen (1996)
def distance(noise_lvl,aadt_npce,speed):
    if aadt_npce==0:
        return(0)
    else:
        return( 10**( (noise_lvl - 38.1*log10(speed*1.60934) - 10*log10(aadt_npce/24/speed) -4.095)/-15))
                                                                    #speed seemingly should be *1.60934 ? 
                                                                    #in Haling and Cohen (1996) it is specified this way, the second speed should be in mph
#aggregation functions for pandas groupby
#div is used to divide the noise damage cost between multiple CBG intersecting a roadway
#gath is used to store all of the CBG that intersect a roadway for when aggregating results to a higher regional level
def div(x):
    return(sum(x)/len(x))
def gath(x):
    return[i for i in set(x)]

#loads housing data for census blocks
#constructs x_meter buffer around road segments
#finds intersecting census block groups for each road buffer
#calculates number of houses in buffer based on (intersection area / census block group area)
#adds column for number of houses and road_buffer area to road.shp geopanda_dataframe
def est_houses(road,x):
    global house_shp
    gdb = root+"\\ACS 2019\\"+acs
    house = gpd.read_file(gdb, layer="X25_Housing_Characteristics"); house = house.drop('geometry',1)

    house = house.rename(columns = {"GEOID":"GEOID_Data"})
    house_shp = pandas.merge(house,c_shp,how='inner',on=['GEOID_Data'])
    house_shp = gpd.GeoDataFrame(house_shp)
    house_shp['Area'] = house_shp['geometry'].area
    keep = ['GEOID_Data','B25001e1','Area','geometry'] #'Shape_Area',
    house_shp = house_shp[keep]
    
    dist = str(int(x/.3048))  #x is passed in meters, name converted to ft
    road_bufferX = gpd.GeoDataFrame(road.buffer(x, cap_style=2)) 
    road_bufferX = road_bufferX.rename(columns={0:'geometry'})
    road_bufferX['id'] = road['id']
    road_bufferX['rb'+dist+'_area_m'] = road_bufferX['geometry'].area  #area is based on the coordinate units, which is m

    c_ = gpd.overlay(road_bufferX,house_shp, how='intersection') 
    c_['area'] = c_['geometry'].area
    c_['houses'+dist] = c_['area']/c_['Area']*c_['B25001e1']  #area is in sq_m
    #c_['ratio'+dist] = c_['area']/c_['Area']
    
    #c_ has more observations than road as some segments intersect multiple census block groups
    #then groupby 'Route_ID' summing houses, merge with road by 'id'
    c = c_.groupby(c_['id']).aggregate({('houses'+dist):'sum',('rb'+dist+'_area_m'):'first'})#,('ratio'+dist):'first'})
    road = pandas.merge(road,c, how='left',on=['id'])
    return(road)

def latlon(c):    #originally data is in lonlat, some functions require latlon
    return((c[1],c[0]))

#some HPMS segments have different geometry types (line, multipart, polygon)
#coordinate conversion method needs different reference for geometry types
def tryconvert(row):
    try: 
        return(int(convert_wgs_to_utm(latlon(row['geometry'].coords[0]))))  #for line
    except:
        try:
            return(int(convert_wgs_to_utm(latlon(row['geometry'].geoms[0].coords[0])))) #for multipart
        except:
            return(int(convert_wgs_to_utm(latlon(row['geometry'].centroid.coords[0]))))  #for polygon


# c is an index for the State number in an alphabetical list
def equity(c):
    global inte, inte_road, county, census, road, inte_census, inte_county, race_shp, road, road_write  #global for testing purposes
    global acs, c_shp  #these are global because est_houses uses them, could pass them to the function instead
    gc.collect()  #resolves memory issue, failure to overwrite road file for each State
    t0=time.time()
    acs = files[c]      #the ACS State nomenclature (space, no space, or _) is different than HPMS, so this method requires folders for all states be in ACS
    state = house_value['State'][c]   #get the state name as it is refered to in various places
    state_utm_zone = STATES_DATA[c][3]
    crs_utm = CRS.from_string('epsg:326' + state_utm_zone)  #to facilitate more accurate distances each State's coordinates are converted to a utm zone
    
    ###################
    # 1. Preparing Data
    print(state,"reading files", end=', ')
    #read in polygons of census block groups (cbg)
    gdb = root+"\\ACS 2019\\"+acs
    c_shp = gpd.read_file(gdb, layer=acs[:-4])
    c_shp = c_shp.to_crs(crs_utm)
    
    #race and poverty demographics of cbg
    #load and merge with geometry in c_shp
    race = gpd.read_file(gdb, layer="X02_Race"); race = race.drop('geometry',1)
    poverty = gpd.read_file(gdb, layer="X17_Poverty"); poverty = poverty.drop('geometry',1)
    race = pandas.merge(race,poverty,how='inner',on='GEOID')
    race['nonwhite'] = race['B02001e1'] - race[ACS['white']]
    race['nonpoverty'] = race['B02001e1'] - race[ACS['poverty']] #could automate these non_steps depending on how many
    race = race.rename(columns = {"GEOID":"GEOID_Data"})
    race_shp = pandas.merge(race,c_shp,how='inner',on=['GEOID_Data'])
    race_shp['county_FIPS'] = race_shp.apply(lambda row: row['GEOID'][0:5], axis=1)
    race_shp['tract_FIPS'] = race_shp.apply(lambda row: row['GEOID'][2:-1], axis=1)
    
    
#     house = gpd.read_file(gdb, layer="X25_Housing_Characteristics"); house = house.drop('geometry',1)
#     house = house.rename(columns = {"GEOID":"GEOID_Data"})
#     keep = ['GEOID_Data','B25001e1']
#     house = house[keep]
#     race_shp = pandas.merge(race_shp,house,how="inner",on=['GEOID_Data'])
#     state_hh_size = sum(race_shp[])
        
    #road file for the state
    #gdb = root+"\\state shp\\"+state+".shp"
    gdb = root+"\\state shp\\"+state+".shp"
    road = gpd.read_file(gdb)
    road = road[road['geometry'] != None]
    road['id'] = road.index
    renames = {i:i.upper() for i in road.columns}    #HPMS 2016 and 2017 use different naming conventions
    renames['geometry'] = 'geometry'
    
    #some cleanup on column types
    road['ROUTE_ID'] = road.apply(lambda row: str(row['ROUTE_ID']), axis=1)
    change = ['F_SYSTEM','AADT','AADT_SINGL','AADT_COMBI','SPEED_LIMI']
    for i in change:
        road[i] = road[i].apply(int)
    
    road['AADT_PASSENGER'] = road['AADT'] - road['AADT_COMBI'] - road['AADT_SINGL']
    
    #sets speed limit for roads without entry. based on avg speed for these road types
    typ_sp = {0:45,1:65,2:55,3:45,4:45,5:35,6:30,7:25}
    typ = list(set(road['F_SYSTEM']))
    typ.sort()
    road['SPEED_LIMI'] = road['SPEED_LIMI'].fillna(0)
    for i in typ:
        road['SPEED_LIMI'][(road['F_SYSTEM']==i) & 
                           ((road['SPEED_LIMI']==0) | (road['SPEED_LIMI']>90))] = typ_sp[i]
        
    #attribute AADT to vehicle classes
    #as their NPCE varies with speed differently
    road['AADT_NPCE'] = road.apply(lambda row : npce(row['AADT_COMBI'],
                        row['AADT_SINGL'],row['AADT_PASSENGER'],
                        row['SPEED_LIMI']), axis=1)
        
    road = road.to_crs(crs_utm)
    road['length'] = road.length/1609.34  #length is in miles and crs is in meters
    
    ############################
    #2. Estimate Housing Density
    print("housing density", end=', ')
    
    ## Count the number of houses in each buffer strip from the road
    dist = [50,300,1000]  #in ft, UTM is in m
    for d in dist:
        #t0=time.time()
        road = est_houses(road,d*.3048)
        road['houses'+str(d)] = road['houses'+str(d)].fillna(0)
        road['rb'+str(d)+'_area_m'] = road['rb'+str(d)+'_area_m'].fillna(1)
        #print(d,time.time()-t0)
        
    #convert the buffer strip areas from meters to feet
    road['rb1000_area_ft'] = road['rb1000_area_m']/(.3048**2)
    road['rb300_area_ft'] = road['rb300_area_m']/(.3048**2)
    road['rb50_area_ft'] = road['rb50_area_m']/(.3048**2)
    
    #subtract the houses/area from smaller buffer, since the larger contains the smaller
    #want the number of houses/area in the buffer_strip
    road['h1000_all'] = road['houses1000']
    road['houses1000'] = road['houses1000'] - road['houses300']
    road['houses300'] = road['houses300'] - road['houses50']
    road['rb1000_area_ft'] = road['rb1000_area_ft'] - road['rb300_area_ft']
    road['rb300_area_ft'] = road['rb300_area_ft'] - road['rb50_area_ft']
    
    for d in dist:
        road['rb'+str(d)+'_area_ft'] = road['rb'+str(d)+'_area_ft'].replace(0,1)
        #some values have zero area after the subtraction due to buffer overlay being the same area
        #alternatively could replace null or inf density values with 0
        #replacing it with 1 is not the true buffer area, but removes the issue of dividing by 0
        #this should only arise for buffers on edge or that don't intersect any ACS, like roads outside State boundary

    #housing density is houses per square mile in the buffer_strip
    ## we then take the width affected by each noise_range
    # and use it to estimate how many houses are affected at that distance
    road['HD_50'] = road['houses50'] / (road['rb50_area_ft']/5280**2)
    road['HD_300'] = road['houses300'] / (road['rb300_area_ft']/5280**2)
    road['HD_1000'] = road['houses1000'] / (road['rb1000_area_ft']/5280**2)


    ########################
    #3. Estimate Noise Level
    print("noise ", end=', ')
    
    #iterates over noise levels and calculates the max distance those levels affect given the traffic and speed
    for i in range(55,85,5):
        road[str(i)+"_dBA_dist"] = road.apply(lambda row : distance(i,row['AADT_NPCE'],row['SPEED_LIMI']),axis=1)
    
    # calculate buffer_strip for each maximum noise level
    # the buffer_strip is the distance between each dBA_dist, the width of the area affected by that noise range
    for i in range(60,85,5):
        road["buffer_"+str(i-5)+"_"+str(i)] = road.apply(lambda row: row[str(i-5)+'_dBA_dist']
                                                       - row[str(i)+'_dBA_dist'],axis=1)
    road["buffer_>80"] = road.apply(lambda row: row['80_dBA_dist']-30 if row['80_dBA_dist']>30 else 0, axis=1)

    road = road.fillna(0)

    #calculate the number of housing units per road mile at each max noise level
    # noise_distance_ft * both_sides / (ft/mile) * housing_density_sq_mile = houses/mile
    for i in range(60,85,5):
        road['houses_per_mile_'+str(i-5)+"_"+str(i)] = road.apply(lambda row:
    row["buffer_"+str(i-5)+"_"+str(i)]*2/5280*row['HD_50'] if row["buffer_"+str(i-5)+"_"+str(i)]<50 else ( 
    row["buffer_"+str(i-5)+"_"+str(i)]*2/5280*row['HD_300'] if row["buffer_"+str(i-5)+"_"+str(i)]<300 else 
    row["buffer_"+str(i-5)+"_"+str(i)]*2/5280*row['HD_1000'] ), axis=1)

    road['houses_per_mile_>80'] = road.apply(lambda row:
    row["buffer_>80"]*2/5280*row['HD_50'] if row["buffer_"+str(i-5)+"_"+str(i)]<50 else ( 
    row["buffer_>80"]*2/5280*row['HD_300'] if row["buffer_"+str(i-5)+"_"+str(i)]<300 else 
    row["buffer_>80"]*2/5280*row['HD_1000'] ), axis=1)
    
    # calculate the noise damage cost
    # houses/mile * noise_effect * housing_value_change_per_noise * road_miles = housing_value_change
    for i in range(60,85,5):
        road['dmg_'+str(i-5)+"_"+str(i)] = road.apply(lambda row: 
                row['houses_per_mile_'+str(i-5)+"_"+str(i)]*(i-57.5)*hv[state]*row['length'], axis=1)
    road['dmg_>80'] = road.apply(lambda row: row['houses_per_mile_>80']*27.5*hv[state]*row['length'], axis=1)
    road['dmg'] = road.apply(lambda row: row['dmg_55_60']+row['dmg_60_65']+row['dmg_65_70']+
                            row['dmg_70_75']+row['dmg_75_80']+row['dmg_>80'],axis=1)

    #############################
    #4. Attribute to Demographics
    #     State and National
    print("attribute to demos ")
    #attribute dmg to each demographic for intersecting census block groups
    #keep track of sum of dmg to each demographic
    #at end compare % of dmg received to population share

    road = gpd.GeoDataFrame(road)
    
    #This method of intersection doesn't account for census block group intersecting the buffer but not the road
    inte = gpd.sjoin(race_shp,road, op='intersects')
    
    agg_functions ={'B02001e1':'sum','dmg':'first','county_FIPS':gath,'tract_FIPS':gath,'GEOID':gath,
                   'URBAN_CODE':'first'}    
    for i in demos:
        agg_functions[ACS[i]] = 'sum'
    inte_road = inte.groupby(inte['index_right']).aggregate(agg_functions)
    #this method ensures each road is counted once, while summing up the populations of the intersecting CBGs
    #essentially forming 'super CBGs' of all the CBGs intersecting a road in order the calculate the share of each demographic population
    #later, when aggregating over larger areas the dmg for each road segment is divided by the number of intersecting CBGs to ensure all dmg is only counted once

    inte_road = inte_road.fillna(0)  
    
    #split the noise dmg for each road segment into demographics based on share of population
    for i in demos:
        inte_road[i+'_dmg'] = inte_road[ACS[i]] / inte_road['B02001e1'] * inte_road['dmg']
    
    inte_road = inte_road.fillna(0)
    
    ################### For writing road.shp with noise dmg estimates ###################
#     keep = ['YEAR_RECOR','STATE_CODE','ROUTE_ID','BEGIN_POIN','END_POINT','ROUTE_NUMB','ROUTE_NAME',
#        'F_SYSTEM','FACILITY_T','NHS','THROUGH_LA','SPEED_LIMI','AADT','AADT_COMBI','AADT_SINGL','AADT_NPCE',
#        'Shape_Leng','length','geometry','HD_50','HD_300','HD_1000']
#     keep = ['STATE_CODE','ROUTE_ID','ROUTE_NUMB','F_SYSTEM','FACILITY_T','URBAN_CODE',
#             'SPEED_LIMI','AADT','AADT_PASSENGER','AADT_COMBI','AADT_SINGL','AADT_NPCE',
#             'length','geometry','HD_50','HD_300','HD_1000','h1000_all']
#     road_write = road[keep]
#     road_write['index_right'] = road_write.index
#     road_write = pandas.merge(inte_road,road_write,how='inner',on=['index_right'])
#     road_write = gpd.GeoDataFrame(road_write)
#     road_write['county_FIPS'] = road_write.apply(lambda row: str(row['county_FIPS']), axis=1 )
#     road_write['tract_FIPS'] = road_write.apply(lambda row: str(row['tract_FIPS']), axis=1 )
#     road_write['GEOID'] = road_write.apply(lambda row: str(row['GEOID']), axis=1 )
#     rname = {'B02001e1':'total_pop','URBAN_CODE_y':'URBAN_CODE'}
#     for i in demos:
#         rname[ACS[i]] = i+'_pop'
#     road_write = road_write.rename(columns=rname)
#     #Having issues with to_file(*.shp) filling out too many decimal places
#     road_write['dmg_length'] = road_write['dmg']/road_write['length']/10  #dmg per 1/10 mile
#     road_write['dmg_length'] = road_write.apply(lambda row: int(row['dmg_length']), axis=1)
#     #keep = ['ROUTE_NUMB','F_SYSTEM','dmg_length','geometry']
    
#     keep_csv = ['STATE_CODE','county_FIPS','tract_FIPS','GEOID','URBAN_CODE',
#                 'ROUTE_NUMB','F_SYSTEM','FACILITY_T', 
#                 'AADT_PASSENGER','AADT_SINGL','AADT_COMBI','AADT_NPCE','SPEED_LIMI', 
#                 'dmg','dmg_length','length','total_pop','h1000_all',
#                 'HD_50','HD_300','HD_1000']
#     for i in demos:
#         keep_csv.append(i+"_pop")
#         keep_csv.append(i+"_dmg")
#     road_write_csv = pandas.DataFrame(road_write[keep_csv])
#     for i in demos:
#         road_write_csv[i+"_house"] = road_write_csv['h1000_all']*road_write_csv[i+"_pop"]/road_write_csv['total_pop']
#         road_write_csv[i+"_dmg_per_house"] = road_write_csv[i+"_dmg"]/road_write_csv[i+"_house"]
#     root_out_csv =r""
# #     road_write_csv.to_csv(root_out_csv+"\\"+state+"_noise_dmg.csv")  #use separate program to combine States
    
#     keep = ['STATE_CODE','county_FIPS','tract_FIPS','GEOID','URBAN_CODE',
#                 'ROUTE_NUMB','F_SYSTEM','FACILITY_T', 
#                 'AADT_PASSENGER','AADT_SINGL','AADT_COMBI','AADT_NPCE','SPEED_LIMI', 
#                 'dmg','dmg_length','length','total_pop','h1000_all',
#                 'HD_50','HD_300','HD_1000','geometry']
#     for i in demos:
#         keep.append(i+"_pop")
#         keep.append(i+"_dmg")
#     road_write = road_write[keep]
#     road_write = road_write.to_crs("epsg:4326")
#     root_out = r""
#     road_write.to_file(root_out+"\\"+state+"_noise_dmg.shp")
    ###################                                               ###################
    
    
    result['state'].append(state)
    result['state_pop'].append( sum(race_shp['B02001e1']))
    result['state_dmg'].append( sum(inte_road['dmg']) )
    for i in demos:
        result[i+"_pop"].append( sum(race_shp[ACS[i]]) )
        result[i+'_dmg'].append( sum(inte_road[i+'_dmg']))
    result['rural_dmg'].append( sum(inte_road[inte_road['URBAN_CODE']==99999]['dmg']))
    result['urban_dmg'].append( sum(inte_road[inte_road['URBAN_CODE']!=99999]['dmg']))
    
    
    #This section prints out the noise-equity ratios for the state
    print("\nState: ", state, 
          " Time: ",round(time.time() - t0,2), sep='')
    print("-------------------")
    print("White")
    print(round(sum(inte_road['white_dmg']) / sum(inte_road['dmg']),4),'dmg')
    print(round(sum(race_shp['B02001e2']) / sum(race_shp['B02001e1']),4),'pop')

    print("Black")
    print(round(sum(inte_road['black_dmg']) / sum(inte_road['dmg']),4))
    print(round(sum(race_shp['B02001e3']) / sum(race_shp['B02001e1']),4))

    print("Native")
    print(round(sum(inte_road['native_dmg']) / sum(inte_road['dmg']),4))
    print(round(sum(race_shp['B02001e4']) / sum(race_shp['B02001e1']),4))

    print("Asian")
    print(round(sum(inte_road['asian_dmg']) / sum(inte_road['dmg']),4))
    print(round(sum(race_shp['B02001e5']) / sum(race_shp['B02001e1']),4))
    
    print("Pacific")
    print(round(sum(inte_road['pacific_dmg']) / sum(inte_road['dmg']),4))
    print(round(sum(race_shp['B02001e6']) / sum(race_shp['B02001e1']),4))
    
    print("Other")
    print(round(sum(inte_road['other_dmg']) / sum(inte_road['dmg']),4))
    print(round(sum(race_shp['B02001e7']) / sum(race_shp['B02001e1']),4))
    

    #############################
    #5. Attribute to Demographics
    #    County and Census Tract
    
    #below avoids multi-counting for roads that intersect multiple counties/census tracts, 
    # dmg needs to be divided since the summation is across counties/census tracts when reaggregated
    # assumes dmg is split between intersecting counties/census tracts [in revision will split based on proportion of road intersecting]
    inte_road['dmg'] = inte_road['dmg'] / inte_road.apply(lambda row: len(row['GEOID']), axis=1)
    for i in demos:
        inte_road[i+'_dmg'] = inte_road[i+'_dmg'] / inte_road.apply(lambda row: len(row['GEOID']), axis=1)
    inte_road['index_right'] = inte_road.index
    inte_road.index.name = 'id'
    inte_road = inte_road.rename(columns={'dmg':'total_dmg'})   #to distinguish it for remerge
    keep = ['index_right','total_dmg']#'nonwhite_dmg'
    for i in demos:
        keep.append(i+'_dmg')   #don't need 'road_block' population counts, those were used only for proportion of dmg allocation
    inte_road = inte_road[keep]
    ######  ** ##
    keep_inte = ['county_FIPS','tract_FIPS','GEOID','index_right']
    inte = inte[keep_inte]   #to simplify dataframe for later steps
    inte_ = pandas.merge(inte,inte_road, how='inner',on=['index_right']) 
    
    agg_functions = {'total_dmg':'sum'}#, 'nonwhite_dmg':'sum'} #just to count the damage groups
    for i in demos:
        agg_functions[i+'_dmg'] = 'sum'
    
    inte_county = inte_.groupby(inte_['county_FIPS']).aggregate(agg_functions)
    inte_county['county_FIPS'] = inte_county.index
    inte_county.index.names = ['id']
    
    agg_pop = {'B02001e1':'sum'} #just to count the population groups
    for i in demos:
        agg_pop[ACS[i]] = 'sum'
    county_shp = race_shp.groupby(race_shp['county_FIPS']).aggregate(agg_pop)
    
    #inte_county has the correct damage for each, county_shp has the correct population
    inte_county = pandas.merge(inte_county,county_shp, how='inner',on=['county_FIPS'])
    
    rname = {'B02001e1':'total_pop'}
    for i in demos:
        rname[ACS[i]] = i+'_pop'
    inte_county = inte_county.rename(columns=rname)
    inte_county = inte_county.fillna(0)
    
    #calculate the noise-equity ratios (called ndp here for noise% divided by population%)
    for i in demos:
        inte_county[i+'_ndp'] = inte_county[i+'_dmg']/inte_county['total_dmg'] / (inte_county[i+'_pop']/inte_county['total_pop'])
    inte_county = inte_county.fillna(0)
    keep = ['county_FIPS','name']#,'white_ndp','black_ndp','native_ndp','asian_ndp','pacific_ndp','other_ndp']
    for i in demos:
        keep.append(i+'_ndp')  
        keep.append(i+'_dmg')
        keep.append(i+'_pop')
    #keep.append('nonwhite_ndp')
    inte_county = pandas.merge(inte_county,county_names,how='inner',on='county_FIPS')
    
    inte_county_ = inte_county[keep]
#     root_out = r""
#     inte_county_.to_csv(root_out+"\\"+state+"_noise_dmg_equity_county.csv")
    
    #add county names
    
    inte_county2 = pandas.melt(inte_county_,id_vars=['county_FIPS','name'],var_name='metrics',value_name='values')
    inte_county2["Demographic"] = inte_county2.apply(lambda row: row['metrics'][:-4], axis=1)
    inte_county2['metrics'] = inte_county2.apply(lambda row: row['metrics'][-3:], axis=1)
    if c==0:
        county= pandas.DataFrame(columns=list(inte_county2.columns))
    county = pandas.concat([county,inte_county2])
    
    ###################
    ##### Census  #####
    inte_census = inte_road.groupby(inte_['tract_FIPS']).aggregate(agg_functions)
    inte_census['tract_FIPS'] = inte_census.index
    inte_census.index.names = ['id']
    
    census_shp = race_shp.groupby(race_shp['tract_FIPS']).aggregate(agg_pop)
    
    #inte_census has the correct damage for each, census_shp has the correct population
    inte_census = pandas.merge(inte_census,census_shp, how='inner',on=['tract_FIPS'])
    
    inte_census = inte_census.rename(columns=rname)
#     inte_census = inte_census.rename(columns={'B02001e1':'total_pop','B02001e2':'white_pop','B02001e3':'black_pop',
#                       'B02001e4':'native_pop','B02001e5':'asian_pop','B02001e6':'pacific_pop',
#                       'B02001e7':'other_pop', 'dmg':'total_dmg' })
    inte_census = inte_census.fillna(0)
    
    for i in demos:
        inte_census[i+'_ndp'] = inte_census[i+'_dmg']/inte_census['total_dmg'] / (inte_census[i+'_pop']/inte_census['total_pop'])
    inte_census = inte_census.fillna(0)
    #inte_census['nonwhite'+'_ndp'] = inte_census['nonwhite'+'_dmg']/inte_census['total_dmg'] / (inte_census['nonwhite'+'_pop']/inte_census['total_pop'])
    keep = ['tract_FIPS']#,'white_ndp','black_ndp','native_ndp','asian_ndp','pacific_ndp','other_ndp']
    for i in demos:
        keep.append(i+'_ndp')
        keep.append(i+'_dmg')
        keep.append(i+'_pop')
    #keep.append('nonwhite_ndp')
    inte_census_ = inte_census[keep]
#     root_out = r""
#     inte_census_.to_csv(root_out+"\\"+state+"_noise_dmg_equity_census.csv")
    
#     race_shp_ = race_shp[['TRACTCE','geometry']]  #merging for geometry
#     # ! #  need census tract geometry--NOT census block geometry
#     inte_census_ = pandas.merge(inte_census_,race_shp_,how='inner',on=['TRACTCE'])
    
    inte_census2 = pandas.melt(inte_census_,id_vars=['tract_FIPS'],var_name='metrics',value_name='values')
    inte_census2["Demographic"] = inte_census2.apply(lambda row: row['metrics'][:-4], axis=1)
    inte_census2['metrics'] = inte_census2.apply(lambda row: row['metrics'][-3:], axis=1)
    if c==0:
        census= pandas.DataFrame(columns=list(inte_census2.columns))
    census = pandas.concat([census,inte_census2])
    
    print("Total time:",round((time.time()-t0)/60,1),"mins")


#This block must be run before equity is called to set up the results storage
#National/State results will be stored in result
#County and Census Tract are stored separately
result = {'state':[], 'state_pop':[]}
for i in demos:
    result[i+'_pop'] = []
result['state_dmg'] = []
for i in demos:
    result[i+'_dmg'] = []
result['rural_dmg'] = []
result['urban_dmg'] = []
county= pandas.DataFrame()     
census = pandas.DataFrame()  


# Main running Section #
#calls the primary function 'equity' for each state
#alternatively it can be called for states individually
fails = []
for c in range(len(house_value['State'])):
    try:
        equity(c)
    except Exception as e:
        print("\n",house_value['State'][c], "failed")
        print(e)
        fails.append(house_value['State'][c])
    

results = pandas.DataFrame(result)

########***##################
#add US row
#first sum each column, store in list
#add list as row to results
###############***###########

us = [i for i in results.loc[:,results.columns!='state'].apply(lambda col: sum(col), axis=0)]
us.insert(0,'US')
results.loc[len(results.index)] = us
results['urban_pct'] = results['urban_dmg'] / (results['urban_dmg'] + results['rural_dmg'])
results['rural_pct'] = results['rural_dmg'] / (results['urban_dmg'] + results['rural_dmg'])

for i in demos:
    results[i+'_ndp'] = results[i+'_dmg']/results['state_dmg'] / (results[i+'_pop']/results['state_pop'])

keep = ['state','urban_pct','rural_pct']#,'white_ndp','black_ndp','native_ndp','asian_ndp','pacific_ndp','other_ndp']
for i in demos:
    keep.append(i+'_ndp')
results_ = results[keep]
results2 = pandas.melt(results_,id_vars=['state'],var_name='metrics',value_name='values')
results2["Demographic"] = results2.apply(lambda row: row['metrics'][:-4], axis=1)
results2['metrics'] = results2.apply(lambda row: row['metrics'][-3:], axis=1)

#may want to include date
results2.to_csv(root+"\\Equity_long.csv")
county.to_csv(root+"\\Equity_county.csv")
census.to_csv(root+"\\Equity_census.csv")