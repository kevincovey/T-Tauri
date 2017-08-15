def apVisit_Catalog_Output(filename,savefile):
    '''
    Notes:
    '''

    import functions
    import csv

    #Setting up Dictionary----------------------------------------------------------------------------|
    apVisit_Dict = {} #Not used right now. Will shift to it eventually.

    #Creating List------------------------------------------------------------------------------------|
    master = []
    with open(filename) as csvfile:

        reader = csv.reader(csvfile,delimiter=',')
        j = 0
        for row in reader:
            loc_id = int(row[1])
            twomass_id = '%s' % row[0]
            #print(loc_id,twomass_id,type(twomass_id))
            x = functions.apStar_to_apVisit(loc_id,twomass_id)
            for i in range(len(x)):
                plate = x[i][0]
                mjd = x[i][1]
                fiber = x[i][2]
                master.append((loc_id,twomass_id,plate,mjd,fiber))
                j += 1
                print('Appending %s to master' %j)
            
    with open(savefile,'w') as savefile:
        
        writer = csv.writer(savefile,delimiter = '\t')
        writer.writerow(('Location ID','2Mass ID','Plate','MJD','Fiber'))
        for i in range(len(master)):
            writer.writerow((master[i][0],master[i][1],
                            master[i][2],master[i][3],master[i][4]))
            print('Writing row %s' %i)

def find_nearest(array,value):
    from numpy import abs
    index = (abs(array-value)).argmin()
    return index

def apStar_to_apVisit(locid,twomassid):
    import apogee.tools.read as apread

    header = apread.apStar(locid,twomassid,ext=0,header=True)

    visits = header[1]['NVISITS']
    array = []
    for i in range(visits):
        x = i+1
        SFILE = 'SFILE%s' % x
        plate = header[1][SFILE][11:15]
        MJD = header[1][SFILE][16:21]
        fiber = header[1][SFILE][22:25]
        array.append((int(plate),int(MJD),fiber,int(visits)))
    return array,locid,twomassid

def Br_Equiv_Width_Apogee(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    from astropy.io import fits
    import functions

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    #Importing header via astropy--------------------------------------------------------------------------------------------------------------------|
    filename = functions.File_Path(plateid,MJD,fiber)
    main_header = fits.open(filename)

    #Barycentric Correction--------------------------------------------------------------------------------------------------------------------------|
    
    #vbcstring = 'BC' + str(1) Somehow need to figure out which visit this MJD applies to
    vbc = main_header[0].header['BC']
    vhelio = main_header[0].header['VHELIO']
    observed_wavelength,shift,rest_wavelength = functions.Barycentric_Correction(emission_line,vbc)
    
    #Equivalent Width Calculation--------------------------------------------------------------------------------------------------------------------|

    #Finding the centerline and checking that it matches the peak in the window
    centerline = find_nearest(wave,observed_wavelength) #Finds the closest element of wave for our observed peak
    centerline_check = Max_Flux_Check(wave,spec,centerline)
    
    '''if emission_line == 11:
        if spec[centerline] != spec[centerline_check[1]]:
            print('The centerline has changed from ' + str(centerline) + ' to ' + str(centerline_check[1]) + ' with a new flux of ' + str(centerline_check[0]) + 
            ' from ' + str(spec[centerline]) + '.')
            centerline = centerline_check[1]'''
    

    L1 = centerline - 240 # ~ 56 Angstroms
    L2 = centerline - 150 # ~ 35 Angstroms
    R1 = centerline + 150
    R2 = centerline + 240

    Fluxcontinuum = (np.sum(spec[L1:L2])+np.sum(spec[R1:R2])) / (len(spec[L1:L2])+len(spec[R1:R2]))
    EqW1 = 0

    if Fluxcontinuum == 0:

        EqW1 = 0
        EqW1_rounded = 0

    if Fluxcontinuum != 0:

        for i in range(L2,centerline):

            left_area = (wave[i+1]-wave[i])*(spec[i+1]-Fluxcontinuum)-(1./2.)*(wave[i+1]-wave[i])*(spec[i+1]-spec[i])
            EqW1 += left_area

        for i in range(centerline,R1):

            right_area = (wave[i+1]-wave[i])*(spec[i]-Fluxcontinuum)-(1./2.)*(wave[i+1]-wave[i])*(spec[i]-spec[i+1])
            EqW1 += right_area

        EqW_rounded = round(EqW1/Fluxcontinuum,5)
        EqW = EqW1/Fluxcontinuum
    
    return EqW,EqW_rounded,vbc,vhelio,Fluxcontinuum,centerline,shift
 
def Br_Equiv_Width_Plotter(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    import matplotlib.pyplot as plt
    import functions

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    #Values for plotter needed-----------------------------------------------------------------------------------------------------------------------|
    EqW,EqW_rounded,vbc,vhelio,Fluxcontinuum,centerline,shift = functions.Br_Equiv_Width(plateid,MJD,fiber,emission_line)

    #Plot averaged spectrum with EqW-----------------------------------------------------------------------------------------------------------------|
    title = str(plateid)+'-'+str(MJD)+'-'+str(fiber)+'-'+str(emission_line)
    fig,ax = plt.subplots(figsize=(16,8))
    plt.plot(wave+shift,spec,linewidth=2.5,label='Shifted')
    
    plt.axhline(y=Fluxcontinuum,ls='dashed',color='black')
    plt.axvline(x=wave[centerline]+shift,ls='dashed',color='r',label='Rest Emission')
    
    plt.legend(loc=1,prop={'size':18})
    plt.xlabel('Wavelength'+' '+'('+ r'$\AA$'+')', fontsize=24)
    plt.ylabel('Flux (erg s' + r'$^{-1}$'+' cm'+r'$^{-2}$' + r'$\AA^{-1}$'+')', fontsize=24)
    plt.suptitle(title,fontsize = 20)
    plt.xlim(wave[centerline]-40,wave[centerline]+40)
    plt.ylim(Fluxcontinuum-(1/2)*(spec[centerline]-Fluxcontinuum),Fluxcontinuum+2*(spec[centerline]-Fluxcontinuum))
    ax.tick_params(axis='both', labelsize=20)   
    #plt.show()
def Unshifted_Plotter(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    import matplotlib.pyplot as plt
    import functions

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    #Values for plotter needed-----------------------------------------------------------------------------------------------------------------------|
    EqW,EqW_rounded,vbc,vhelio,Fluxcontinuum,centerline,shift = functions.Br_Equiv_Width(plateid,MJD,fiber,emission_line)

    #Plot averaged spectrum with EqW-----------------------------------------------------------------------------------------------------------------|
    title = str(plateid)+'-'+str(MJD)+'-'+str(fiber)+'-'+str(emission_line)
    fig,ax = plt.subplots(figsize=(16,8))
    plt.plot(wave,spec,linewidth=2.5,label='Not Shifted')
    
    plt.axhline(y=Fluxcontinuum,ls='dashed',color='black')
    plt.axvline(x=wave[centerline]+shift,ls='dashed',color='r',label='Rest Emission')
    
    plt.legend(loc=1,prop={'size':18})
    plt.xlabel('Wavelength'+' '+'('+ r'$\AA$'+')', fontsize=24)
    plt.ylabel('Flux (erg s' + r'$^{-1}$'+' cm'+r'$^{-2}$' + r'$\AA^{-1}$'+')', fontsize=24)
    plt.suptitle(title,fontsize = 20)
    plt.xlim(wave[centerline]-40,wave[centerline]+40)
    plt.ylim(Fluxcontinuum-(1/2)*(spec[centerline]-Fluxcontinuum),Fluxcontinuum+2*(spec[centerline]-Fluxcontinuum))
    ax.tick_params(axis='both', labelsize=20)   
    #plt.show()
def Skyline_Plotter(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import functions

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|
    
    spec = apread.apVisit(plateid,MJD,fiber,ext=5,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    #Values for plotter needed-----------------------------------------------------------------------------------------------------------------------|
    EqW,EqW_rounded,vbc,vhelio,Fluxcontinuum,centerline,shift = functions.Br_Equiv_Width(plateid,MJD,fiber,emission_line)

    #Plot averaged spectrum with EqW-----------------------------------------------------------------------------------------------------------------|
    title = str(plateid)+'-'+str(MJD)+'-'+str(fiber)+'-'+str(emission_line)
    fig,ax = plt.subplots(figsize=(16,8))
    plt.plot(wave,spec,linewidth=2.5,label='Shifted')
    
    plt.axhline(y=Fluxcontinuum,ls='dashed',color='black')
    plt.axvline(x=wave[centerline]+shift,ls='dashed',color='r',label='Rest Emission')
    
    
    plt.xlabel('Wavelength'+' '+'('+ r'$\AA$'+')', fontsize=24)
    plt.ylabel('Flux (erg s' + r'$^{-1}$'+' cm'+r'$^{-2}$' + r'$\AA^{-1}$'+')', fontsize=24)
    plt.suptitle(title,fontsize = 20)
    plt.xlim(wave[centerline]-40,wave[centerline]+40)
    plt.ylim(0,3*Fluxcontinuum)
    ax.tick_params(axis='both', labelsize=20)   
    #plt.show()

def Barycentric_Correction(emission_line,vbc):

    #Constants---------------------------------------------------------------------------------------------------------------------------------------|
    
    n = (float(emission_line))**2 #Beginning electron level
    c = 299792.458 #Speed of light (km/s)
    rydberg_inf =  1.0973731568539*(10**7) #Rydberg constant (m^-1)
    electron = 9.10938356*(10**-31) #Mass of electron (kg)
    nucleus = 1.672621898*(10**-27) #Mass of hydrogen nucleus (kg)
    
    #Equations---------------------------------------------------------------------------------------------------------------------------------------|
    
    rydberg_red = rydberg_inf/(1+(electron/nucleus)) #Reduced mass Rydberg constant (m^-1)
    rest_wavelength1 = rydberg_red*((1./16.)-(1./n)) #(m^-1)
    rest_wavelength = 1/rest_wavelength1 #wavelength (m)
    observed_wavelength1 = rest_wavelength*(1-(vbc/c)) #Finding the location of the peak of the observed spectrum
    observed_wavelength = observed_wavelength1*(10**10)
    shift = (rest_wavelength-observed_wavelength1)*(10**10) #Finding the difference between rest and observed wavelengths

    #Returns-----------------------------------------------------------------------------------------------------------------------------------------|

    return observed_wavelength,shift,rest_wavelength

def Max_Flux_Check(x_axis,y_axis,centerline):
    import functions
    c = 299792.458
    v_window = 100
    right_shift = x_axis[centerline]*(1+(v_window/c))
    left_shift = x_axis[centerline]*(1-(v_window/c))
    leftwindow = functions.find_nearest(x_axis,left_shift)
    rightwindow = functions.find_nearest(x_axis,right_shift)
    y_max = max(y_axis[leftwindow:rightwindow])
    z = y_axis.tolist().index(y_max)
    return y_max,z


def Brackett_Ratios(plateid,mjd,fiber):
    
    import csv
    import os
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import numpy as np

    
    #Header------------------------------------------------------------------------------------------------------------------------------------------|

    filename = File_Path(plateid,mjd,fiber)
    main_header = fits.open(filename)
    
    loc_id = main_header[0].header['LOCID']
    twomass_id = main_header[0].header['OBJID']

    #Reading in the Visits file----------------------------------------------------------------------------------------------------------------------|
    

    fileDir = os.path.dirname(os.path.realpath('__file__'))
    filename = os.path.join(fileDir, '../Data/Average Visits.csv')
    visits = os.path.abspath(os.path.realpath(filename))
    
    br_num = np.asarray([11,12,13,14,15,16,17,18,19,20])
    br_value=np.zeros(10)

    with open(visits) as csvfile:

        reader = csv.DictReader(csvfile,delimiter='\t')

        for row in reader:
            if int(row['Location ID'])==loc_id and row['2Mass ID']==twomass_id:
                for i in range(10):
                    num = 11 + i
                    br = 'Br' + str(num) + ' Avg EqW'
                    br_value[i] = float(row[br])

    #Plotting---------------------------------------------------------------------------------------------------------------------------------------|    
    
    fig,ax=plt.subplots(figsize=(16,8))
    ax.tick_params(axis='both', labelsize=20)
    
    plt.plot(br_num,br_value/br_value[0])
    plt.ylabel('Br n>11 / Br 11',fontsize=24)
    plt.xlabel('n',fontsize=24)
    plt.show()
    

def File_Path(plateid,mjd,fiber):

    import os

    #Creating file path------------------------------------------------------------------------------------------------------------------------------|
    '''
    Notes:
        - Need to rework this to include Location ID from master list
    '''
    server = '/Volumes/CoveyData/APOGEE_Spectra/python_DR13/dr13/apogee/spectro/redux/r6/apo25m/'
    #server = '/Users/ballanr/Desktop/SummerResearch/dr13/apogee/spectro/redux/r6/apo25m/'
    plate = str(plateid)
    MJD = str(mjd)
    fiber_num = str(fiber)
    therest = 'apVisit-r6-'
    dashes = '-'
    f = '.fits'

    endname = therest + plate + dashes + MJD + dashes + fiber_num + f
    filename = os.path.join(server,plate,MJD,endname)

    return filename
    
def Balmer_Decrement_Plot():

    import matplotlib.pyplot as plt
    import numpy as np

    #Importing Decrement File------------------------------------------------------------------------------------------------------------------------|
    data = np.loadtxt('/Users/ballanr/Desktop/Fwd__Bracket_Decrement/x74.txt',delimiter = '\t',unpack=True)
    print(data[0])
    T = [3750,5000,7500,8750,10000,12500,15000]
    for i in range(len(T)):
        plt.plot(np.linspace(8,12.4,num=23),data[i],label = str(T[i]) + ' K')
        plt.scatter(np.linspace(8,12.4,num=23),data[i])
    plt.xticks(np.arange(8,12.6,0.2))
    for k in range(5):
        x = 8 + k
        plt.axvline(x,ls = 'dashed',linewidth = 0.5,color = 'black')
        g = 0.1 + (0.1*k)
        plt.axhline(g,ls = 'dashed', linewidth = 0.5, color = 'black')
    plt.xlabel('Log $n_e$')
    plt.ylabel('Transition Probabilities')
    plt.legend()
    plt.show()
    
def Probs_Calc(n,T):

    import numpy as np

    #Constants------------------------------------------------------------------------------------------------------------------------|
    k = 8.6173303*(10**(-5)) #eV/K
   
    #Calculations---------------------------------------------------------------------------------------------------------------------|
    n_i = 4**2
    n_k = n**2
    E_1 = -13.6/n_i
    E_2 = -13.6/n_k
    exponent = -(E_2 - E_1)/(k*T)
    boltz = 2*(n_k/n_i)*np.exp(exponent)

    return boltz

def Probs_Plots():
    import numpy as np
    import functions
    import matplotlib.pyplot as plt
    
    T = [3750,5000,7500,8750,10000,12500,15000]
    for j in range(len(T)):
        xx = []
        yy = []
        for i in range(20):
            x = 1 + i
            y = functions.Probs_Calc(x,T[j])
            y = np.log(y)
            #plt.scatter(x,y)
            xx.append(x)
            yy.append(y)
        plt.plot(xx,yy,label = str(T[j]) + ' K')
        plt.scatter(xx,yy)
        plt.xlabel('N')
        plt.ylabel('Probability')
    plt.legend(bbox_to_anchor=(1.25,1))
    #plt.xlim(2,20)
    plt.xticks(np.arange(2,21,1))
    #plt.ylim(-2.75,1)
    plt.show()

def Saha(n,T,density):

    '''
    This isn't really useful right now...this looks at the ionization states,
    not the transition probabilities...
    '''

    import numpy as np

    k_J = 1.380648*(10**(-23)) #Joules
    k_e = 8.617*(10**(-5)) # eVs
    
    n_i = n**2
    n_k = (n+1)**2
    n_e = density
    h = 6.626*(10**(-34))
    m_e = 9.109*(10**(-31))
    E_1 = -13.6/n_i
    E_2 = -13.6/n_k
    E = E_2 - E_1

    part1 = (1/n_e)
    #print(part1)
    part2 = (2*(n_k/n_i))
    part2 = 1
    #print(part2)
    part3 = (((2*np.pi*m_e*k_J*T)/(h**2))**(3/2))
    #print(part3)
    part4 = (np.exp(E_1/(k_e*T)))
    #print(part4)
    
    #A = np.log(part1*part2*part3*part4)
    A = part1*part2*part3*part4
    
    return A
    
    
def Boltzmann(n_upper,n_lower,T):
    
    import numpy as np

    k_e = 8.617*(10**(-5)) # eVs
    

    g_b = 2*(n_upper**2)
    g_a = 2*((n_lower)**2)
    del_E = -3.4+13.6
    E_b = -13.6/(n_upper**2)
    E_a = -13.6/(n_lower**2)

    B = (g_b / g_a)*np.exp(-(E_b - E_a)/(k_e * T))

    return B

def Saha_Boltzmann(n_upper,n_lower,ion,T,density):
    
    import functions
    import numpy as np

    x = functions.Saha(ion,T,density)
    y = functions.Boltzmann(n_upper,n_lower,T)
    
    SB = (y/(1+y))*(1/(1+x))

    return np.log(SB)

def SB_Plotter():

    import functions
    import numpy as np
    import matplotlib.pyplot as plt

    T = [3750,5000,7500,8750,10000,12500,15000]
    gg = [11,12,13,14,15,16,17,18,19,20]
    densities = []
    for i in range(23):
        y = 8 + (0.2*i)
        x = (100**3)*(np.exp(y))
        densities.append(x)
    #Plotter 1
    '''
    yy = []
    for i in range(len(T)):
        y = np.log(functions.Saha_Boltzmann(11,4,1,T[i]))
        yy.append(y)
    plt.plot(T,yy)
    plt.scatter(T,yy)
    '''
    #Plotter 2

    for i in range(len(T)):
        
        for j in range(len(densities)):
            yy = []
            for k in range(len(gg)):
                y = functions.Saha_Boltzmann(gg[k],4,1,T[i],densities[j])
                yy.append(y)
            plt.plot(gg,yy,label='Temp ' + str(T[i])+'density '+str(densities[j]))
            plt.scatter(gg,yy)
    plt.legend(bbox_to_anchor=(1,1))
    plt.show()

def SB_CSV(savefile):

    import csv
    import functions
    import numpy as np
    import matplotlib.pyplot as plt
    
    with open(savefile,'w') as savefile:

        densities = []
        for i in range(23):
            y = 8 + (0.2*i)
            x = (100**3)*(np.exp(y))
            densities.append(x)
        
        T = [3750,5000,7500,8750,10000,12500,15000]
        brackett = [11,12,13,14,15,16,17,18,19,20]

        names = ['Densities','3750','5000','7500','8750','10000','12500','15000']
        writer = csv.DictWriter(savefile,delimiter = '\t',fieldnames = names)
        writer.writeheader()
        b=[]
        for i in range(len(densities)):
            a = []
            for k in range(len(T)):
                y = functions.Saha_Boltzmann(11,4,1,T[k],densities[i])
                a.append(y)
            #writer.writerow({'Densities':densities[i],'3750':a[0],'5000':a[1],'7500':a[2],'8750':a[3],'10000':a[4],'12500':a[5],'15000':a[6]})
            b.append(a)
        for i in range(len(b)):
            for k in range(7):
                b[i][k] = b[i][k]/b[12][3]
        for i in range(len(densities)):
            writer.writerow({'Densities':densities[i],'3750':b[i][0],'5000':b[i][1],'7500':b[i][2],'8750':b[i][3],'10000':b[i][4],'12500':b[i][5],'15000':b[i][6]})

        for i in range(23):
            plt.plot(T,b[i])
            plt.scatter(T,b[i])
        plt.legend(bbox_to_anchor=(1,1))
        plt.show()

def hundred_plotter(filename):

    import csv
    import functions
    import matplotlib.pyplot as plt

    with open(filename) as csvfile:

        reader = csv.reader(csvfile,delimiter = '\t')
        i = 0
        x = []
        for row in reader:
            if i < 102:
                x.append(row)
                i += 1
            else:
                break
    #return x
    print(len(x))
    for k in range(len(x)):
        if k != 0:
            line = 11
            plate = int(x[k][2])
            MJD = int(x[k][3])
            if len(x[k][4]) == 3:
                fiber = x[k][4]
            if len(x[k][4]) == 2:
                fiber = '0'+x[k][4]
            if len(x[k][4]) == 1:
                fiber = '00'+x[k][4]
            y = functions.Br_Equiv_Width_Plotter(plate,MJD,fiber,line)
            print(fiber)
            title = x[k][2]+'-'+x[k][3]+'-'+fiber
            plt.savefig('/Users/ballanr/Desktop/Test/'+title+'-'+str(line)+'.png')
            plt.close()

def twomass_uniqueness(filename,savefile):
    import csv
    import functions
    import math

    uniques = [[str(4614),'2M05412754-0646074']] #Setting up array to check for unique 2M IDs
    nonuniques = [] #The array for storing non-unqiue 2M IDs
    
    g = 1
    gg = 579468

    with open(filename) as csvfile:

        reader = csv.reader(csvfile,delimiter='\t') 
        next(reader, None)  # skip the headers
        for row in reader:
            locid = row[0]
            twomassid = row[1]
            x = functions.Uniqueness_Filter(uniques, locid, twomassid) #Checks the location ID and 2M ID against the growing list of uniques
            
            if x == 1:
                uniques.append((row[0],row[1]))
                #print('Unique visit %s' %g)
                g += 1

            elif x == 0:
                uniques.append((row[0],row[1]))
                nonuniques.append((row[0],row[1]))
                #print('Non-unique visit %s' %g)
                g += 1
            
            percent = (g/gg)*100
            numbertimes = math.floor(percent)

            hashes = '#' * int(numbertimes)
            
            spaces = ' '* (100 - int(numbertimes))

            print('Sorting through files: %.5f done ['%percent + hashes + str(spaces) + ']',end='\r',flush = True)

    with open(savefile,'w') as savefile:
        
        writer = csv.writer(savefile,delimiter = '\t')
        
        writer.writerow(('Location ID','2Mass ID'))
        
        for i in range(len(nonuniques)):
            if i != 0:
                percent = (i/len(nonuniques))*100
                numbertimes = math.floor(percent)
                hashes = '#' * int(numbertimes)
                spaces = ' '* (100 - int(numbertimes))

                writer.writerow((nonuniques[i][1],nonuniques[i][0]))

                print('Writing rows: %.5f done ['%percent + hashes + str(spaces) + ']',end='\r',flush = True)
    
    return uniques,nonuniques

def Uniqueness_Filter(array, locid, twomassid):
    
    state = 0
    
    for a,b in array:

        if a != locid and b != twomassid:
            state = 1
        
        elif a != locid and b == twomassid:
            state = 0
    
    return state



def apVisit_Updated_Catalog(infile,rowstart):

    import functions
    import pandas as pd
    from astropy.io import fits as fits
    import numpy as np
    import itertools

    x = pd.read_csv(infile,delimiter = '\t')
    problems = []
    cols = ['Location ID','2Mass ID', 'Plate ID','MJD','Fiber','confidence1','confidence2','Mean','arearatios','SNR','Br 11 EqW']#,'Br 12 EqW','Br 13 EqW','Br 14 EqW','Br 15 EqW','Br 16 EqW','Br 17 EqW',
            #'Br 18 EqW','Br 19 EqW','Br 20 EqW']
    
    df = pd.DataFrame(columns = cols)
    g = 0
    #rowstart = 10000
    rowend =  rowstart + 79468
    for index,row in itertools.islice(x.iterrows(),rowstart,rowend):
        try:
            g+=1
            loc = row['Location ID']
            twomass = row['2Mass ID']
            print(str(g))
            plateid = row['Plate']
            mjd = row['MJD']
            if len(str(row['Fiber'])) == 3:
                fiber = str(row['Fiber'])
            elif len(str(row['Fiber'])) == 2:
                fiber = '0' + str(row['Fiber']) 
            else:
                fiber = '00' + str(row['Fiber'])
                    
            fitsfile = functions.File_Path(plateid,mjd,fiber)

            #this passes a filepath that astropy can open with fits, circumventing apogee entirely...
            
            openfile = fits.open(fitsfile)
            snr = openfile[0].header['SNR']
            vbc = openfile[0].header['BC']
            c = 299792.458
            lamshift = 1 + (vbc/c)
            
            fspec = openfile[1]
            fwave = openfile[4]
            wave = []
            flux = []
            
            for i in range(len(fwave.data[2])):
                wave.append(fwave.data[2][-i-1])
                flux.append(fspec.data[2][-i-1])
            for i in range(len(fwave.data[1])):
                wave.append(fwave.data[1][-i-1])
                flux.append(fspec.data[1][-i-1])
            for i in range(len(fwave.data[0])):
                wave.append(fwave.data[0][-i-1])
                flux.append(fspec.data[0][-i-1])
            
            openfile.close()
            
            newflux = functions.OH_Skylines_Remover(wave,flux)
            
            #now we run equiv width calc
            lines = [11,12,13,14,15,16,17,18,19,20]
            EqW = []
            equiv_check = 1000
            
            #just for line 11
            equiv_width,fcontinuum,shift,rest_wavelength = functions.Br_EqW(wave,newflux,11,vbc)
            equiv_check = equiv_width
            #for all the lines
            '''
            for i in range(10):
                
                equiv_width,fcontinuum,shift,rest_wavelength = functions.Br_EqW(wave,newflux,lines[i],vbc)
                EqW.append(equiv_width)
                
                if i == 0:
                    equiv_check = equiv_width
            '''
                
            wave = np.asarray(wave) * lamshift #check which direction shift is going
            
            if equiv_check > 0:
                
                confidence1,confidence2,Mean,arearatios = functions.Confidence_Level(wave,flux,snr)
                #confidence3 = confidence1/confidence2

                data = [int(loc),twomass,int(plateid),int(mjd),fiber,confidence1,confidence2,Mean,arearatios,snr,equiv_width]#EqW[0],EqW[1],EqW[2],EqW[3],EqW[4],EqW[5],EqW[6],EqW[7],EqW[8],EqW[9]]

                df.loc[len(df)+1] = data

                df1 = pd.DataFrame(wave,columns=['Wavelength'])
                df1['Flux'] = newflux
                filename = str(plateid) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
                df1.to_csv('/Users/ballanr/Desktop/File Outputs/Wave and Flux/'+filename,index=False)
                df1 = df1.iloc[0:0]

        except KeyError:
            print('Row '+str(g)+' has no BC value...')
            problems.append((loc,twomass,'KeyError'))
            openfile.close()
        
        except FileNotFoundError:
            print('Row '+str(g)+' doesn\'t exist...')
            problems.append((loc,twomass,'FileNotFound'))

    df.to_csv('/Users/ballanr/Desktop/File Outputs/'+str(rowstart)+'-'+str(rowend)+' End.csv',index=False)
    df = df.iloc[0:0]

    probs = pd.DataFrame(problems,columns = ['Location ID','2Mass ID','Problem Type'])
    probs.to_csv('/Users/ballanr/Desktop/File Outputs/'+str(rowstart)+'-'+str(rowend)+' Problems.csv',index=False)

def Br_EqW(wave,spec,line,vbc):
    
    import functions
    import numpy as np
    
    #wave = np.asarray(wave)
    observed_wavelength,shift,rest_wavelength = functions.Barycentric_Correction(line,vbc)
    rest = rest_wavelength*(10**10)
    centerline = functions.find_nearest(wave,rest)
    
    L1 = centerline - 301 # ~ 27.42 Angstroms
    L2 = centerline - 150 # ~ 17.21 Angstroms
    R1 = centerline + 150
    R2 = centerline + 301

    Fluxcontinuum = (np.sum(spec[L1:L2])+np.sum(spec[R1:R2])) / (len(spec[L1:L2])+len(spec[R1:R2]))
    EqW1 = 0

    if Fluxcontinuum == 0:

        EqW1 = 0
        EqW1_rounded = 0
        equivs = 0
        
    if Fluxcontinuum != 0:

        for i in range(L2,R1):

            trapezoid = (0.5)*(wave[i+1] - wave[i])*(spec[i+1] + spec[i] - (2*Fluxcontinuum))
            EqW1 += trapezoid
        
    #EqW_rounded = round(EqW1/Fluxcontinuum,5)
        equivs = EqW1/Fluxcontinuum
    
    return equivs,Fluxcontinuum,shift,rest_wavelength


def Confidence_Level(wave,flux,snr,restwave):

    import numpy as np
    import functions

    center = functions.find_nearest(wave,restwave)
    L1 = center - 301 # ~ 56 Angstroms
    L2 = center - 150 # ~ 35 Angstroms
    R1 = center + 150
    R2 = center + 301

    #method 1 
    leftwindow = np.mean(flux[L1:L2])
    rightwindow = np.mean(flux[R1:R2])
    cmean = (0.5)*(leftwindow + rightwindow)
    linemean = np.mean(flux[L2:R1])
    lefterr = np.std(flux[L1:L2])
    righterr = np.std(flux[R1:R2])
    cerr = (0.5)*(lefterr + righterr)
    confidence1 = (linemean - cmean)/cerr
    
    #method 2
    n_l = len(wave[L2:R1])
    n_c = len(wave[L1:L2])+len(wave[R1:R2])
    l = linemean
    c = cmean
    dellam = wave[L1+1]-wave[L1]
    r = l/c
    top = n_l*(dellam)**2
    bottom = snr**2
    sig = (top/bottom)*(r/n_c)*(r+n_c)
    ebbets = np.sqrt(sig)
    
    #method 4
    abovex = []
    abovey = []
    belowx = []
    belowy = []

    for k in range(R1-L2):
        x = wave[L2 + k]
        y = flux[L2 + k]
        gy = y - cmean
        
        if gy > 0:
            abovex.append(x)
            abovey.append(y)
        else:
            belowx.append(x)
            belowy.append(y)

    aboveArea = 0
    belowArea = 0

    for i in range(len(abovex)-1):

                trapezoid = (0.5)*(abovex[i+1] - abovex[i])*(abovey[i+1] + abovey[i] - (2*cmean))
                aboveArea += trapezoid

    for j in range(len(belowx)-1):

                trapezoid1 = (0.5)*(belowx[j+1] - belowx[j])*(-belowy[j+1] - belowy[j] + (2*cmean))
                belowArea += trapezoid1

    if belowArea != 0:
        funratio = (aboveArea / belowArea)
    else: 
        funratio = aboveArea
          
    return confidence1,ebbets,funratio

def skylines_cleaner(wave,flux):
    
    import functions
    import numpy as np

    windows = []
    wave = np.asarray(wave)

    removal_beginning = [15185.0,15240.0,15286.0,15330.0,15393.6,15429.9,15460.8,15473.3,15499.6,15508.0,15516.0,
                        15538.5,15545.0,15568.9,15595.9,15652.9,15701.3,15758.9,15861.5,15867.1,15874.3,15888.8,
                        15971.1,16028.8,16078.2,16126.4,16193.1,16233.4,16269.5,16278.8,16301.1,16314.8,16339.1,
                        16349.3,16356.3,16358.6,16387.1,16413.5,16475.1,16476.9,16478.3,16501.4,16527.4,16530.4,
                        16542.7,16552.9,16558.7,16585.4,16609.4,16611.6,16654.1,16658.2,16666.4,16688.2,16691.2,
                        16701.7,16707.7,16711.8,16717.7,16723.1,16725.1,16730.9,16753.4,16755.9,16761.5,16764.8,
                        16839.5,16852.8,16877.2,16883.0,16888.8,16890.3,16899.6,16902.4,16907.1,16908.6,16909.9,16913.7]

    removal_ending = [15189.0,15242.6,15289.0,15334.7,15396.6,15433.9,15463.8,15475.3,15502.3,15511.0,
                    15518.8,15543.4,15547.6,15571.4,15599.3,15658.9,15704.3,15761.8,15863.6,15871.0,
                    15877.0,15892.7,15974.4,16033.1,16081.6,16131.0,16196.5,16237.9,16271.4,16280.7,
                    16304.7,16318.7,16343.0,16353.7,16357.2,16361.7,16390.3,16416.2,16476.3,16477.9,
                    16479.7,16503.7,16528.9,16531.1,16543.5,16555.1,16560.2,16587.5,16610.7,16612.8,
                    16655.0,16661.0,16667.3,16690.4,16693.9,16703.8,16710.5,16712.5,16719.7,16724.5,
                    16725.9,16734.9,16755.2,16756.8,16762.5,16765.9,16841.6,16853.3,16877.8,16884.5,
                    16889.5,16891.2,16900.2,16905.5,16907.7,16909.3,16910.6,16914.3]

    for i in range(len(removal_beginning)):

        lwindow = functions.find_nearest(wave,removal_beginning[i])
        rwindow = functions.find_nearest(wave,removal_ending[i])

        windows.append((lwindow,rwindow))

    for i in range(len(removal_beginning)):

        lwindowelement = int(windows[i][0])
        rwindowelement = int(windows[i][1])
        leftwindow = wave[windows[i][0]]
        rightwindow = wave[windows[i][1]]
        leftflux = flux[windows[i][0]]
        rightflux = flux[windows[i][1]]

        slope = (rightflux - leftflux) / (rightwindow - leftwindow)

        for k in range(rwindowelement - lwindowelement):

            fluxvalue = slope*(wave[lwindowelement+k] - leftwindow) + leftflux
            flux[lwindowelement+k] = fluxvalue

    return flux


def DR15_Brackett_Catalog():

    import functions
    import pandas as pd
    from astropy.io import fits as fits
    import numpy as np
    import itertools
    from PyAstronomy.pyasl import helcorr as helcorr

    dr15table = pd.read_csv('/Users/ballanr/Desktop/File Outputs/forMdot_analysis.csv')
    problems = []
    cols = ['Location ID','2Mass ID', 'Plate ID','MJD','Fiber','S/R','Overall Confidence','Br11 Error',
            'Br12 Error','Br13 Error','Br14 Error','Br15 Error','Br16 Error','Br17 Error','Br18 Error',
            'Br19 Error','Br20 Error','Br11 EqW','Br12 EqW','Br13 EqW','Br14 EqW','Br15 EqW','Br16 EqW',
            'Br17 EqW','Br18 EqW','Br19 EqW','Br20 EqW']

    rowstart = 1000

    df = pd.DataFrame(columns = cols)
    g = 0

    for index,row in itertools.islice(dr15table.iterrows(),rowstart,None):
        try:
            g+=1
            loc = row['LOCID']
            twomass = row['2MASS_ID']
            print(str(g))
            plateid = row['PLATE']
            mjd = row['MJD']
            snr = row['S/N']
            if len(str(row['FIB'])) == 3:
                fiber = str(row['FIB'])
            elif len(str(row['FIB'])) == 2:
                fiber = '0' + str(row['FIB']) 
            else:
                fiber = '00' + str(row['FIB'])

            fitsfile = row['path-to-file-on-CoveyData']

            #this passes a filepath that astropy can open with fits, circumventing apogee entirely...
            
            openfile = fits.open(fitsfile)

            header = openfile[0].header

            ra = header['RA']
            dec = header['DEC'] 
            jd = header['JD-MID']
            telescope = header['TELESCOP']
            
            if telescope == 'apo25m':
                height = 2788
                longg = -105.4913
                lat = 36.4649
            else:
                height = 2380
                longg = -70.413336
                lat = -29.05256

            vbc,hjd = helcorr(longg,lat,height,ra,dec,jd)

            c = 299792.458
            lamshift = 1 + (vbc/c)
            
            fspec = openfile[1]
            ferr = openfile[2]
            fwave = openfile[4]
            wave = []
            flux = []
            error = []
            
            for i in range(len(fwave.data[2])):
                wave.append(fwave.data[2][-i-1])
                flux.append(fspec.data[2][-i-1])
                error.append(ferr.data[2][-i-1])
            for j in range(len(fwave.data[1])):
                wave.append(fwave.data[1][-j-1])
                flux.append(fspec.data[1][-j-1])
                error.append(ferr.data[2][-j-1])
            for k in range(len(fwave.data[0])):
                wave.append(fwave.data[0][-k-1])
                flux.append(fspec.data[0][-k-1])
                error.append(ferr.data[2][-k-1])
            
            openfile.close()
            
            newflux = functions.skylines_cleaner(wave,flux)
            
            #now we run equiv width calc
            lines = [11,12,13,14,15,16,17,18,19,20]
            EqW = []
            restwaves = []
            econfidences = []
            confidences = []
            equivs_error = []
            equiv_check = 1000
            
            wave = np.asarray(wave) * lamshift #check which direction shift is going

            for i in range(10):
                
                equiv_width,fcontinuum,shift,rest_wavelength = functions.Br_EqW(wave,newflux,lines[i],vbc)
                rest = rest_wavelength*(10**10)
                EqW.append(equiv_width)
                restwaves.append(rest)
                
                dEqW = functions.Br_Error(wave,newflux,error,lines[i])
                equivs_error.append(dEqW)
                
                if i == 0:
                    equiv_check = equiv_width
            
            if equiv_check > 0:
                
                for k in range(len(restwaves)):

                    Confidence,ebbet,arearatios = functions.Confidence_Level(wave,flux,snr,restwaves[k])
                    #econfidences.append(ebbet)
                    confidences.append(Confidence)

                data = [int(loc),twomass,int(plateid),int(mjd),fiber,snr,confidences[0],equivs_error[0],equivs_error[1],
                        equivs_error[2],equivs_error[3],equivs_error[4],equivs_error[5],equivs_error[6],equivs_error[7],
                        equivs_error[8],equivs_error[9],EqW[0],EqW[1],EqW[2],EqW[3],EqW[4],EqW[5],EqW[6],EqW[7],EqW[8],EqW[9]]

                df.loc[len(df)+1] = data

                df1 = pd.DataFrame(wave,columns=['Wavelength'])
                df1['Flux'] = newflux
                df1['Error'] = error
                filename = str(plateid) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
                df1.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Wave and Flux/'+filename,index=False)
                df1 = df1.iloc[0:0]                        

        except KeyError:
            print('Row '+str(g)+' has no BC value...')
            problems.append((loc,twomass,plateid,mjd,fiber,'KeyError'))
            openfile.close()
        
        except FileNotFoundError:
            print('Row '+str(g)+' doesn\'t exist...')
            problems.append((loc,twomass,plateid,mjd,fiber,'FileNotFound'))


    df.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Catalog-test.csv',index=False)
    df = df.iloc[0:0]

    probs = pd.DataFrame(problems,columns = ['Location ID','2Mass ID','Plate ID','MJD','Fiber','Problem Type'])
    probs.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Catalog Problems-test.csv',index=False)

def Br_Error(wave,flux,err,line):
    
    import functions
    import numpy as np

    emission,nah,nahh = functions.Barycentric_Correction(line,0)
    center = functions.find_nearest(wave,emission)
    
    wave1 = np.asarray(wave[center-301:center+301])
    flux1 = np.asarray(flux[center-301:center+301])
    err1 = np.asarray(err[center-301:center+301])
    
    L1 = 0
    L2 = 151
    R1 = 451
    R2 = 602
    
    lcheck = np.median(err1[L1:L2])
    for k in range(L1,L2):
        if err1[k] > 10*lcheck:
            err1[k] = 10*lcheck


    rcheck = np.median(err1[R1:R2])
    for k in range(R1,R2):
        if err1[k] > 10*rcheck:
            err1[k] = 10*rcheck


    
    Fc = (np.sum(flux1[L1:L2]) + np.sum(flux1[R1:R2])) / (len(flux1[L1:L2]) + len(flux1[R1:R2]))
    dFc = (np.sum(err1[L1:L2]) + np.sum(err1[R1:R2])) / (len(err1[L1:L2]) + len(err1[R1:R2]))
    

    fluxavg = 0 
    
    for i in range(L2,R1):

        fluxavg += flux1[i] - Fc
        
        #A += (0.5)*(wave1[i+1] - wave1[i])*(flux1[i+1] + flux1[i] - (2*Fc))
    
    N = len(wave1[L2:R1])
    fluxavg = fluxavg / N
    del_lam = wave1[2] - wave1[1]    
    A = N*del_lam*fluxavg

    dA = 0
    for k in range(L2,R1):

        #dA += (0.5)*(wave1[k+1] - wave1[k])*(err1[k+1] + err1[k] + (2*dFc))

        dA += err[k]**2

    dA = del_lam * np.sqrt(dA)

    dEqW = (1/(Fc**2)) * np.sqrt( (Fc*dA)**2 + (A*dFc)**2 )
    #EqW = A / Fc

    return dEqW

def Brackett_Ratios_Updated(plate,mjd,fiber):

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import itertools
    
    #finding the equivalent widths and errors for a given visit
    filefile = pd.read_csv('/Users/ballanr/Desktop/File Outputs/DR15/Catalog-test.csv')

    if len(str(fiber)) == 3:
        Fiber = str(fiber)
    elif len(str(fiber)) == 2:
        Fiber = '0' + str(fiber) 
    else:
        Fiber = '00' + str(fiber)
    equivs = []
    errors = []
    for index,row in itertools.islice(filefile.iterrows(),0,None):
        plateid = int(row['Plate ID'])
        MJD = int(row['MJD'])
        ffiber = str(row['Fiber'])
        
        if len(str(ffiber)) == 3:
            ffiber = str(ffiber)
        elif len(str(ffiber)) == 2:
            ffiber = '0' + str(ffiber) 
        else:
            ffiber = '00' + str(ffiber)

        if plateid == plate and MJD == mjd and ffiber == Fiber:

            for j in range(10):
                number = 11+j
                errorstring = 'Br' + str(number) + ' Error'
                equivstring = 'Br' + str(number) + ' EqW'
                hh = row[errorstring]
                gg = row[equivstring]
                equivs.append(gg)
                errors.append(hh)


    equivs = np.asarray(equivs)
    errors = np.asarray(errors)
    equivs = equivs/equivs[0]

    temps = ['3750 K','5000 K','7500 K','8750 K','10000 K','12500 K','15000 K']

    #loop through densities
    densityrange = np.arange(8,12.6,0.2)
    for i in range(len(densityrange)):
        
        density = densityrange[i]
        
        filename = '/Users/ballanr/Desktop/File Outputs/Brackett Decrements/Density Files/Density ' + str(8+(i*0.2)) + ' Ratios.csv'
        openfile = pd.read_csv(filename)

        y = []

        for j in range(len(temps)):

            y = openfile[temps[j]]
                    
            y = np.asarray(y)
            y = y/ y[0]
            plt.figure(figsize=(20,10))
            plt.plot(np.arange(11,21,1),y,color='purple',label = temps[j])
            plt.scatter(np.arange(11,21,1),equivs,color='purple',label='_nolabel_')
            plt.errorbar(np.arange(11,21,1),equivs,errors,ecolor='red',fmt='--o',elinewidth = 1,capsize = 5,label='Br Emission Lines')
            plt.xticks(np.arange(11,21,1))
            plt.legend(bbox_to_anchor=(1,1),fontsize=18)
            plt.grid(True,color='grey',ls='dashed',alpha=0.7)
            plt.title('Density ' + str(density))
            plt.xlabel('Brackett Lines',fontsize = 18)
            plt.ylabel('Ratios',fontsize = 18)
            #plt.show()
            savestring = '/Users/ballanr/Desktop/File Outputs/Brackett Decrements/test/D' + str(density) + 'T' + str(temps[j]) + '.pdf'
            plt.savefig(savestring,bbox_inches='tight',dpi=300)
            plt.close()