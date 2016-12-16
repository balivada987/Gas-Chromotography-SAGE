# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 09:19:05 2016


"""
import numpy as np
import matplotlib.pyplot as plt
import time


def local_maxes(data, eps):
    """
    This function finds the local maxima (peak heights) of a set of data.
    The value "eps" can be adjusted to include or exclude
    peaks of a given height.

    This function is not used in the main code (MultiPlex GC Simulation.py)
    because the function peak_params (below) is easier to use. peak_params
    can be used because GC has a constant zero baseline, so it will not work
    for other spectra (XRD, etc.)

    :param data: The data set being considered. This should be an n x 2
    array, where the 1st column contains the independent variable (time, etc.)
    and the 2nd column contains the dependent variable (intensity, height,
    etc.)
    :param eps: The value below which a peak will not be counted as a
    maximum.
    :return: A  vector containing the times at which the maxima occur, as well
    as the magnitudes of the maxima themselves.
    """
    times = data[:, 0]
    data = data[:, 1]
    # eps = 0.1*(max(data) - min(data))
    print(eps, "<-- eps")
    derivs = []
    step = 2
    for i in range(len(data) - step):
        deriv = (data[i + step] - data[i]) / (times[i + step] - times[i])
        if deriv < (eps * -1):
            derivs.append(-1)
        elif abs(deriv) <= eps:
            derivs.append(0)
        else:
            derivs.append(1)
    print(derivs, "<-- derivs")
    peak_indexes = []
    for i in range(len(derivs)):
        if derivs[i] == 1:
            if derivs[i + 1] == -1:
                peak_indexes.append(i + 1)
            elif derivs[i + 1] == 0:
                n = 0
                while derivs[i + n + 1] == 0 and (i + n + 1) != (len(derivs) - 1):
                    n += 1
                if derivs[i + n + 1] == -1:
                    peak_indexes.append(i + n + 1)
    print(peak_indexes, "<-- peak_indexes")
    peaks = []
    peak_times = []
    data = data.tolist()
    for i in peak_indexes:
        loc_max = max(data[(i - step):(i + step)])
        peaks.append(loc_max)
        max_index = data.index(loc_max)
        max_time = times[max_index]
        peak_times.append(max_time)
    print(peaks, "<-- peaks")
    print(peak_times, "<-- peak_times")
    return [peak_times, peaks]


def peak_params(data, step, cutoff, binint):
    """
    This function finds the local maxima (peak heights) of a set of data.
    All peaks whose heights are above the cutoff value are counted.

    :param data: The data set being considered. This should be an n x 2
    array, where the 1st column contains the independent variable (time, etc.)
    and the 2nd column contains the dependent variable (intensity, height,
    etc.)
    :param step: The interval over which the integral is calculated. This
    should be fairly large to minimize the contribution of noise, but cannot
    exceed the first peak index or the program will attempt to integrate
    over (nonexistent) negative time values.
    :param cutoff: The cutoff value for peak detection. If the peak height is
    above this value, it will be counted. Otherwise, it will be ignored.
    :return: An array containing (in order) the peak heights, the
    corresponding times, the indexes at which the peaks occur, and the start
    times of the peaks.
    """
    times = data[:, 0]
    data = data[:, 1]
    peak_heights = []
    peak_times = []
    peak_indexes = []
    for i in range(len(data)):
        if data[i] > cutoff and data[i] >= data[i + 1] and data[i] >= data[i - 1]:
            peak_heights.append(data[i])
            peak_times.append(times[i] + binint)  # add binint because it is also added
            # in the main function to obtain correct retention times
            peak_indexes.append(i)
    """
    print peak_heights, "<-- peak_heights"
    print peak_times, "<-- peak_times"
    print peak_indexes, "<-- peak_indexes"
    """
    integral = []
    for i in peak_indexes:
        if step <= i:
            integral.append(np.trapz(data[(i - step):(i + step)]))
    # print integral, "<-- integral"
    start_times = []
    for i in range(len(peak_indexes)):
        start_times.append(peak_times[i] - (4 * integral[i]) / (peak_heights[i] * np.exp(1) * np.sqrt(np.pi)))
    ret_array = np.zeros((4, len(peak_heights)))
    ret_array[0] = peak_heights
    ret_array[1] = peak_times
    ret_array[2] = peak_indexes
    ret_array[3] = start_times
    return ret_array


if __name__ == "__main__":
    start_time = time.clock()
    """
    times = range(13)
    data = [0,1,2,1,0,1,2,3,3,3,2,1,0]
    data_array = np.zeros((13,2))
    data_array[:,0] = times
    data_array[:,1] = data
    print data_array

    loc_maxes = local_maxes(data_array)
    print loc_maxes
    """
    times = range(6)
    data = [0, 0, 1, 2, 1, 0]
    integral = np.trapz(data[3:6])
    print(integral)

    print("---- %s seconds ----" % (time.clock() - start_time))