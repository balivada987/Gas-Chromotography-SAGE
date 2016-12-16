# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 13:33:12 2016


"""
import numpy as np
import matplotlib.pyplot as plt


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
        j = i % len(barcode)
        data[j] = data[j] + leftovers[i]
    return np.array(data)


if __name__ == "__main__":
    data = [1, 2, 1, 2, 2, 1, 2, 5]
    times = range(8)
    plt.figure(1)
    plt.plot(times, data, "-")

    brcd = [1, 1, 1, 1]
    bin_ = 1
    freq = 1
    times2 = range(len(brcd))

    circ = circ_rep(brcd, bin_, freq, data)
    plt.figure(2)
    plt.plot(times2, circ, "-")