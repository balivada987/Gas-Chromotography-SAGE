


import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys
import FastHadamardTransform as fht
import DataProcessing as dpr
import Instrumentation as ins
import CircRepTesting
from time import clock
import os

def calc_max_bolt_pdf(x, start, mode, height):
    """
    This method calculated the PDF of a Maxwell-Boltzmann (max-bolt)
    distribution for a given vector of x values. The equation for the
    Maxwell-Boltzmann distribution is defined by:

    f(x) = h * (sqrt(2/pi) * (x^2*exp(-x^2/(2*a^2)) / a^3))
    where:
    a = mode/sqrt(2)
    h = (height*exp(1)*a*sqrt(pi)/(2*sqrt(2))

    note: exp(1) = e

    :param x: The values at which the max-bolt pdf is being evaluated
    :param start: The x position that the max-bolt pdf peak first appears
    :param mode: The x position that has the heightest value of the max-bolt
        pdf
    :param height: The value of the max-bolt pdf at the mode
    :return: A vector of the max-bolt pdf values evaluated at the associated
        x-values
    """
    x = np.array(x)
    start = np.array(start)
    mode = np.array(mode)-start
    height = np.array(height)
    xcalc = x-start
    a = mode/np.sqrt(2.)
    h = (height*np.exp(1)*a*np.sqrt(np.pi))/(2*np.sqrt(2))
    f = h*(np.sqrt(2/np.pi)*(xcalc**2*np.exp(-xcalc**2/(2*a**2))/a**3))
    f[xcalc <= 0] = 0

    return f
#------------------------------------------------------Generates Barcode---------------------------------------
def generate_sample_barcode(ppoly, initialvals, n=None):
    """
    This method creates a barcode similar to the ones used in

    O. Trapp, Boosting the Throughput of Separation Techniques by
    "Multiplexing", Angewandte Chemie Int. Ed., 46, 5609-5613 (2007)

    These barcodes can be used to determine the sampling order for HTGC and
    are generated through the use of a shift register.

    :param ppoly: The primitive polynomial that is the characteristics
        circuit polynomial for the shift register. The primitive poly should
        be given as a vector of the powers whose coeffients are 1.
        i.e. P(x) = 1 + x + x^4 ==> ppoly = [0, 1, 4]
    :param initialvals: The initial values populating the shift register.
        All values should always be 0 or 1 and will populate the nodes
        from the closest to the return of the shift register to the furthest
        from the return in the shift register.
    :param n: The number of values the shift register should generate, if n is
        None then only the values from the first cycle of the shift register
        will be calculated.
    :return: The values of the codelet that the shift register generated.
    """
    # This will ensure that the polynomial used does not have repeats and is
    # sorted from the smallest order to largest order

    poly = sorted(set(ppoly))

    # 0 will always be an element in the polynomial, so we will check to
    # ensure that it is included, if not we will add it to poly.
    if 0 not in poly:
        poly.insert(0, 0)

    if (
            len(initialvals) != poly[-1] or
            set(initialvals) not in ({0, 1}, {0}, {1})
    ):
        # This means that the initial values inputted were not valid
        return None

    incompcycle = True
    vals = initialvals
    brcd = list()

    while incompcycle:
        brcd += [vals[0]]
        endval = sum([vals[j] for j in poly[:-1]]) % 2
        vals = vals[1:]
        vals += [endval]
        if vals == initialvals:
            incompcycle = False
    if n is not None:
        mult = int(n/len(brcd))+1
        slcbrcd = brcd*mult
        brcd = slcbrcd[:n+1]
    return brcd

#-----------------------------------------determines minimum and maximum values on important sampling parameters
def determine_sampling_parameters(
        bcodelen, timeint, maxret, sampfreq, numsample
):
    """
    This method will determine some of the minimum and maximum values on
    important sampling parameters to use high-throughput multiplexed gas
    chromatography (htMPGC). These parameters are total time duration that
    will be needed to perform htMPGC (tmax), the maximum number of
    analytes, numanalyte, that can be unambiguously determined in numsample
    samples (also referred to here as channels), and the minimum number of
    repetitive injections that each sample needs to be injected, rmin. These
    three parameters are defined by:

    tmax = bcodelen*timeint+maxret
    numanalyte = ((bcodelen+(maxret/timeint))/numsample)*sampfreq
    rmin = bcodelen/(2*numsample)

    These equations can be found in the supplemental material of:
    O. Trapp, Boosting the Throughput of Separation Techniques by
    "Multiplexing", Angewandte Chemie Int. Ed., 46, 5609-5613 (2007)

    :param bcodelen: The length of the barcode used for sampling, the barcode
        can be created using the generate_sample_barcode method. Note that
        if left to the default n for generate_sample_barcode this will be the
        number of intervals needed to complete 1 cycle of the barcode, and not
        the overall number of intervals that exists.
    :param timeint: The time bin interval between injections.
    :param maxret: The maximum retention time of any of the analytes (species)
        that are being analyzed.
    :param sampfreq: The frequency of data acquisition
    :param numsample: The number of samples (channels) that are being analyzed.
    :return: A vector containing the structure [tmax, numanalyte, rmin]
    """
    tmax = bcodelen*timeint+maxret
    numanalyte = ((bcodelen+(maxret/timeint))/numsample)*sampfreq
    rmin = bcodelen/(2*numsample)
    return [tmax, numanalyte, rmin]
#-------------------------------------Returns a modified version of barcode------------------------



#-------------------Generates convolute chromatrogram using bar codelet----------------------------------------------------------
def generate_convoluted_chromatogram(
        barcodelet, sampleparas, samplingdelay, samplingtimes
):
    """
    This method simulates the convoluted chromatogram from a high-throughput
    multiplexed GC experiment. This simulation assumes boltzmann distribution
    shaped elution peaks for each species.

    :param barcodelet: A modified version of the barcode with the sample
        number that should be sampled. The barcodelet can be generated with
        the define_codelets method.
    :param sampleparas: The sample parameters that should be used to model
        each sample's elution profiles. This should be a N by M by 3 list,
        where N is the number of unique samples that are being injected, M is
        the number of species in each sample (needs to be the same for all N
        sample), and 3 parameters that are in each individual M sub-list are
        (in order) the start, mode, and height parameters of the
        calc_max_bolt_pdf method.
    :param samplingdelay: The total delay (injection time plus delay between
        injections) between injecting two consecutive samples. ***bin int.***
    :param samplingtimes: The times at which each sample is taken by the
        detector. Aka the time values that each data point was taken during
        a GC experiment.
    :return: A vector of the intensity values for the convoluted
        chromatogram.
    """
    totdelaytime = 0
    conchrom = np.array([0]*len(samplingtimes))
    for i in barcodelet:
        if i == 0:
            totdelaytime += samplingdelay
        else:
            smplchroms = [
                calc_max_bolt_pdf(
                    samplingtimes, smpl[0]+totdelaytime,
                    smpl[1]+totdelaytime, smpl[2]
                ) for smpl in sampleparas[i-1]
                ]
            conchrom = conchrom + np.sum(smplchroms, axis=0)
            totdelaytime += samplingdelay
    return conchrom

#-------------------------------plot convoluted chromatogram--------------------------------------------
def plot_convoluted_chromatogram(samplingtimes, convolutedchromatogram):

    plt.figure(1)
    plt.plot(samplingtimes, convolutedchromatogram, "-")
    plt.title("Simulated Data: Single Injection")
    plt.xlabel("Time (min)")
    plt.ylabel("Relative Intensity (%)")

#------------------------splits data produced by convoluted chromatogram
def split_data_set(data, bins, mode):
    """
    This function splits the data produced by generate_convoluted_chromatogram
    into a specified number of bins. This allows the system of equations to be
    solved using a Hadamard matrix.

    :param data: The data set for the convoluted chromatogram.
    :param bins: The number of bins to divide the data into.
    """
    data_vector = []
    split_data = np.array_split(data,bins)
    for n in range(bins):
        if mode == "avg":
            data_avg = sum(split_data[n])/len(split_data[n])
            data_vector.append(data_avg)
        elif mode == "sum":
            data_vector.append(sum(split_data[n]))
        else:
            print("Please enter avg or sum for mode")
            return None
    return data_vector
#-------------- -------------transforms raw chromatogram into circular representation
def circ_rep(barcode, bin_int, sam_freq, cchrom):
    """
    This function transforms the raw chromatogram into its circular
    representation. This means that data points which exceed the sequence
    length divided by the sampling rate are added to the beginning of the
    chromatogram. A Hadamard transformation is then carried out on the
    resulting chromatogram.

    :param barcode: The overall injection sequence, regardless of the
    barcodelets or how many different samples were injected. This is used to
    determine the last data point, after which values must be added to the
    beginning of the chromatogram.
    :param bin_int: The time bin interval, or the combined injection time and
    wait time between injections (i.e. 6 second injection + 54 second wait
    is 60 second bin time).
    :param sam_freq: The sampling frequency, AKA oversampling factor. How often
    the data points are collected. (Technically a time span, not a frequency)
    :param cchrom: The raw data, or convoluted chromatogram.
    :return: A 1D array containing the circular representation of the data.
    """
    last_pt = int(len(barcode) * bin_int / sam_freq)
    data = list(cchrom[:last_pt])
    leftovers = list(cchrom[last_pt:])
    for i in range(len(leftovers)):
        j = i % len(data)
        data[j] = data[j]+leftovers[i]
    return data #np.array(data)


if __name__ == "__main__":
    """
    N2 Func: y = 7.04861*10^(-6)*x
    CH4 Func: y = 2.36348*10^(-6)*x
    y is the volume percent
    x is the peak area
    """
    start_time = clock()

    barcode = generate_sample_barcode([0, 1, 6], [1]*6, n=None)
    #print(len(barcode)), "<-- len(barcode)"
    #print(barcode, "<-- barcode")
    sam_freq = 0.00067   #actual (experimental) value is 0.00067 sec (min?)
    tesnumsample = 1
    tesbarcodelet = define_codelets(barcode, tesnumsample)
    smplparas = [
        #[[3.416, 3.615, 9664], [10.78, 11.265, 2509]],
        [[3.400, 3.620, 10640], [10.766, 11.273, 2460]]
    ]
    binint = 1
    simtimes = np.arange(0, len(barcode)*binint+11.273+5, sam_freq)
    #print len(simtimes)

    cchrom = generate_convoluted_chromatogram(
        tesbarcodelet, smplparas, binint, simtimes)

    # normalize data to 1 (before binning)
    data_max = max(cchrom)
    for i in range(len(cchrom)):
        cchrom[i] = cchrom[i]/data_max

    # write barcode to file for easier input to AHK
    #ins.write_for_ahk(barcode, "barcode")

    plot_convoluted_chromatogram(simtimes, cchrom)

    #data = circ_rep(barcode, binint, sam_freq, cchrom)
    #print len(data)

    """
    bs_times = range(len(data))
    plt.figure(9)
    plt.plot(bs_times, data, "-")
    plt.title("Data after circ_rep and before binning")
    """
#-----------------------------------------------------------------------------------------------------------------------
    #data_vector = split_data_set(data, len(barcode), "avg")

    """
    # normalize data to 1 (after binning)
    data_max = max(data_vector)
    for i in range(len(data_vector)):
        data_vector[i] = data_vector[i]/data_max
    """
    """
    # carry out Fast Hadamard Transform on simulated data
    simplex = fht.create_simplex_matrix(barcode)
    pi1 = fht.pi_one(simplex)
    pi2 = fht.pi_two(simplex, pi1)
    theta = fht.theta(pi1, pi2, data_vector)
    overview = fht.FHT(theta)

    # get array of times corresponding to bins instead of all data points
    times_array = np.array_split(simtimes[:int((len(barcode)*binint)/sam_freq)],
                  len(barcode))
    times = np.array([np.mean(i) for i in times_array])
    """
    """
    plt.figure(2)
    plt.plot(times, overview[::-1], "-") #times+binint
    plt.title("Simulated Overview Chromatogram")
    plt.xlabel("Time (min.)")
    plt.ylabel("Relative Intensity (%)")
    """
    """
    data_array = np.zeros((len(times),2))
    data_array[:,0] = times
    data_array[:,1] = overview[::-1].flatten()
    """
    # below this line is experimental


    #print "Working directory: ", os.getcwd()
    # GC data file must be located in the working directory
    raw_data = np.loadtxt("20160610_20%CH4_80%N2_32inj_1.txt", skiprows=147)
    exptimes = raw_data[:,0]
    expintensities =  raw_data[:,1]

    exp_max = max(expintensities)
    for i in range(len(expintensities)):
        expintensities[i] = expintensities[i]/exp_max

    plt.figure(3)
    plt.plot(exptimes, expintensities, "-")
    plt.title("Experimental Raw Data")

    # convert experimental data to circ. representation, then split into bins
    intensities = circ_rep(barcode, binint, sam_freq, expintensities)
    intensities = split_data_set(intensities, len(barcode), "avg")
    '''

    """
    # carry out Fast Hadamard Transform on experimental data
    theta2 = fht.theta(pi1, pi2, intensities)
    overview2 = fht.FHT(theta2)

    exp_data = np.zeros((len(times),2))
    exp_data[:,0] = times
    exp_data[:,1] = overview2[::-1].flatten()

    plt.figure(4)
    plt.plot(times+binint, overview2[::-1], "-")
    plt.title("Experimental Overview Chromatogram")
    """
    # overlapped plots to compare simulated vs. experimental results
    plt.figure(5)
    plt.plot(simtimes, cchrom, label="Simulated")
    plt.plot(exptimes, expintensities, label="Experimental")
    plt.title("Simulated vs. Experimental Raw Data")
    leg1 = plt.legend()
    for item in leg1.legendHandles:
        item.set_linewidth(3.0)
    plt.xlabel("Time (min.)")
    plt.ylabel("Relative Intensity (%)")
'''
    """
    plt.figure(6)
    plt.plot(times, overview[::-1], "-", label="Simulated")
    plt.plot(times, overview2[::-1], "-", label="Experimental")
    plt.title("Simulated vs. Experimental Overview Chromatograms")
    leg2 = plt.legend()
    for item in leg2.legendHandles:
        item.set_linewidth(3.0)
    plt.xlabel("Time (min.)")
    plt.ylabel("Relative Intensity (%)")


    # peak_params is an array containing:
    # 1) peak heights
    # 2) retention times
    # 3) peak indexes
    # 4) start times

    peak_params = dpr.peak_params(exp_data, 2, 1000)
    print peak_params, "<-- peak_params"

    # spec_params contains (in order) the start time, retention time and
    # height of each peak. This information is taken from peak_params

    spec_params = []
    for i in range(len(peak_params[0,:])):
        spec_params.append([peak_params[3,i],peak_params[1,i],peak_params[0,i]])
    print spec_params, "<-- spec_params"

    # fill_chroms is a list containing convoluted chromatograms obtained by
    # using generate_convoluted_chromatogram with the parameters in spec_params

    fill_chroms = []
    for i in spec_params:
        chrom = generate_convoluted_chromatogram(
        tesbarcodelet, [[i]], binint, times)
        fill_chroms.append(chrom)

    fill_matrix = np.zeros((len(times),len(spec_params)))
    for i in range(len(spec_params)):
        fill_matrix[:,i] = fill_chroms[i]


    #normalize to 1 per Trapp paper:
    fill_matrix_max = np.amax(fill_matrix)

    for i in range(len(fill_chroms)): #***this can be improved***
        for j in range(len(times)):
            fill_matrix[j,i] = fill_matrix[j,i]/fill_matrix_max

    print fill_matrix.shape
    #print intensities.shape

    soln = np.linalg.lstsq(fill_matrix, intensities)
    print soln
    """
    # above this line is experimental
    # below this line is simulated
    """
    # sim_peak_params contains (in order) from the simulated overview chrom.:
    # 1) peak heights
    # 2) retention times
    # 3) peak indexes
    # 4) start times
    sim_peak_params = dpr.peak_params(data_array, 2, 0.1, binint)
    #print(sim_peak_params, "<-- sim_peak_params")

    # sim_spec_params contains (in order) the start time, retention time and
    # height of each peak. This information is taken from sim_peak_params
    sim_spec_params = []
    for i in range(len(sim_peak_params[0,:])):
        sim_spec_params.append(
            [sim_peak_params[3,i],sim_peak_params[1,i],sim_peak_params[0,i]])
    print sim_spec_params, "<-- sim_spec_params"

    # sim_fill_chroms is a list containing convoluted chromatograms obtained by
    # using generate_convoluted_chromatogram with the parameters in
    # sim_spec_params
    sim_fill_chroms = []
    for i in sim_spec_params:
        chrom = generate_convoluted_chromatogram(
            tesbarcodelet, [[i]], binint, simtimes)
        sim_fill_chroms.append(chrom)

    # normalize recalculated raw data to 1
    recalc_data_max = max((max(i) for i in sim_fill_chroms))
    sim_fill_chroms = sim_fill_chroms / recalc_data_max
    """
    """
    plt.figure(10)
    plt.plot(simtimes, sim_fill_chroms[0], "-")
    plt.title("1st col of filling matrix")

    plt.figure(11)
    plt.plot(simtimes, sim_fill_chroms[1], "-")
    plt.title("2nd col of filling matrix")
    """
    """
    plt.figure(12)
    plt.plot(simtimes, cchrom, "-", label="Original")
    plt.plot(simtimes, sim_fill_chroms[0], "-", label="Recalc. N2")
    plt.plot(simtimes, sim_fill_chroms[1], "-", label="Recalc. CH4")
    plt.title("Original Raw Data and Split (not binned) Recalculated Raw Data")
    leg = plt.legend()
    for item in leg.legendHandles:
        item.set_linewidth(3.0)
    """
    """
    sim_fill_matrix = np.zeros((len(times),len(sim_spec_params)))
    for i in range(len(sim_spec_params)):
        sim_fill_matrix[:,i] = sim_fill_chroms[i]
    #print sim_fill_matrix, "<-- sim_fill_matrix"

    # use the built-in scipy least-square solver, instead of a Gauss-Jordan
    # algorithm, to solve for the "concentration" matrix
    sim_soln = np.linalg.lstsq(sim_fill_matrix, data_vector)
    #print(sim_soln[0:2], "<-- sim_soln")
    conc_vector = np.array(sim_soln[0])[np.newaxis].T
    #print(conc_vector, "<-- conc_vector")


    # multiply the "conc." vector by the filling matrix to obtain what should
    # be the raw data: check the error (difference) between them
    hopefully_raw_data = np.dot(sim_fill_matrix, conc_vector)
    binned_error = []
    for i in range(len(data_vector)):
        binned_error.append(data_vector[i] - hopefully_raw_data[i])
    """
    """
    plt.figure(7)
    plt.plot(times, data_vector, "-", label="Original")
    plt.plot(times, hopefully_raw_data, "-", label="Recalculated")
    plt.plot(times, binned_error, "-", label="Error")
    plt.title("Binned Raw Data and Recalculated (binned) Raw Data")
    leg = plt.legend()
    for item in leg.legendHandles:
        item.set_linewidth(3.0)
    """
    """
    # generate new "raw data" using parameters from overview chromatogram
    # scale peak heights by "concentration" vector
    for i in range(len(sim_spec_params)):
        sim_spec_params[i][2] = sim_spec_params[i][2]*sim_soln[0][i]
    #print(sim_spec_params, "<-- sim_spec_params")
    """
    """
    # generate new convoluted chromatogram from overview, normalize to 1
    tot_raw_data = generate_convoluted_chromatogram(
        tesbarcodelet, [sim_spec_params], binint, simtimes)
    tot_raw_data = tot_raw_data/np.max(tot_raw_data)
    """
    """
    # find error between simulated and recalculated raw data
    error = []
    for i in range(len(cchrom)):
        error.append(cchrom[i] - tot_raw_data[i])
    """
    """
    # plot original convoluted chromatogram (before anything) and new
    # convoluted chromatogram recalculated from parameters from overview
    plt.figure(13)
    plt.plot(simtimes, cchrom, "-", label="Original")
    plt.plot(simtimes, tot_raw_data, "-", label="Recalculated")
    plt.title("Original Raw Data and Recalculated Raw Data")
    leg = plt.legend()
    for item in leg.legendHandles:
        item.set_linewidth(3.0)
    plt.xlabel("Time (min.)")
    plt.ylabel("Relative Intensity (%)")
    """
    """
    # plot data after circ rep and binning; this is the data input to the FHT
    plt.figure(8)
    plt.plot(times, data_vector, "-")
    plt.title("Data after circ_rep and binning")
    """
    print("---- %s seconds ----" % (clock() - start_time))
