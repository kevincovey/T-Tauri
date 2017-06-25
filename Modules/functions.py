def apVisit_Catalog_Output(filename,savefile):
    '''
    Notes:
    '''

    import functions
    import csv

    #Setting up Dictionary----------------------------------------------------------------------------|
    apVisit_Dict = {}

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


def Catalog_Update():
    #add stuff here to update the catalog, reading and writing to a new file instead of running it again
    return blah

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
    return array

def Br_Equiv_Width(plateid,MJD,fiber,emission_line):
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
    observed_wavelength,shift,rest_wavelength = Barycentric_Correction(emission_line,vbc)
    
    #Equivalent Width Calculation--------------------------------------------------------------------------------------------------------------------|

    #Finding the centerline and checking that it matches the peak in the window
    centerline = find_nearest(wave,observed_wavelength) #Finds the closest element of wave for our observed peak
    centerline_check = Max_Flux_Check(wave,spec,centerline)

    if spec[centerline] != spec[centerline_check[1]]:
        print('The centerline has changed from ' + str(centerline) + ' to ' + str(centerline_check[1]) + ' with a new flux of ' + str(centerline_check[0]) + 
        ' from ' + str(spec[centerline]) + '.')
        centerline = centerline_check[1]
    

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
    
    return EqW,EqW_rounded,vbc,Fluxcontinuum,centerline,shift
 
def Br_Equiv_Width_Plotter(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import functions

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    #Values for plotter needed-----------------------------------------------------------------------------------------------------------------------|
    EqW,EqW_rounded,vbc,Fluxcontinuum,centerline,shift = functions.Br_Equiv_Width(plateid,MJD,fiber,emission_line)

    #Plot averaged spectrum with EqW-----------------------------------------------------------------------------------------------------------------|
    
    fig,ax = plt.subplots(figsize=(16,8))
    plt.plot(wave+shift,spec,linewidth=2.5,label='Shifted')
    #plt.plot(Lambda,spec1,linewidth=2.5,label='Unshifted')
    plt.axhline(y=Fluxcontinuum,ls='dashed',color='black')
    plt.axvline(x=wave[centerline]+shift,ls='dashed',color='r',label='Rest Emission')
    #plt.axvline(calculated_point2,ls=':',color='r',label='Star Emission')
    plt.legend(loc=1,prop={'size':18})
    plt.xlabel('Wavelength'+' '+'('+ r'$\AA$'+')', fontsize=24)
    plt.ylabel('Flux (erg s' + r'$^{-1}$'+' cm'+r'$^{-2}$' + r'$\AA^{-1}$'+')', fontsize=24)
    plt.xlim(wave[centerline]-40,wave[centerline]+40)
    plt.ylim(Fluxcontinuum-(1/2)*(spec[centerline]-Fluxcontinuum),Fluxcontinuum+2*(spec[centerline]-Fluxcontinuum))
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
    v_window = 500
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
    server = '/Volumes/CoveyData-1/APOGEE_Spectra/python_DR13/dr13/apogee/spectro/redux/r6/apo25m/'
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