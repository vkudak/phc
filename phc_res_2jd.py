import sys
import os
import numpy as np
import math
from datetime import datetime, timedelta


def get_julian_datetime(date):
    """
    Convert a datetime object into julian float.
    Args:
        date: datetime-object of date in question

    Returns: float - Julian calculated datetime.
    Raises:
        TypeError : Incorrect parameter type
        ValueError: Date out of range of equation
    """

    # Ensure correct format
    if not isinstance(date, datetime):
        raise TypeError('Invalid type for parameter "date" - expecting datetime')
    elif date.year < 1801 or date.year > 2099:
        raise ValueError('Datetime must be between year 1801 and 2099')

    # Perform the calculation
    julian_datetime = 367 * date.year - int((7 * (date.year + int((date.month + 9) / 12.0))) / 4.0) + int(
        (275 * date.month) / 9.0) + date.day + 1721013.5 + (
                          date.hour + date.minute / 60.0 + date.second / math.pow(60,
                                                                                  2)) / 24.0 - 0.5 * math.copysign(
        1, 100 * date.year + date.month - 190002.5) + 0.5

    return julian_datetime


def read_header(fname):
    header_list = []
    with open(fname) as f:
        for line in f:
            if len(line) < 100 or line[:4] == "  UT":
                header_list.append(line)
    # print(header_list)
    return header_list


def read_data(fname, header_len):
    date_time = []
    with open(fname) as f:
        line1 = f.readline()
        date, t = line1.split()
        # date = datetime.strptime(date, "%Y-%m-%d")
    date = date.strip()

    time = np.genfromtxt(fname, unpack=True, skip_header=header_len, usecols=(0,), dtype=None, encoding="utf-8")
    ImpB_fon, ImpV_fon, FonB, FonV, mB, mV, Az, El, Rg = \
        np.genfromtxt(filename, skip_header=header_len,
                      usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9,),
                      unpack=True)

    for i in range(0, len(time)):
        date_time0 = datetime.strptime(date + ' ' + time[i], "%Y-%m-%d %H:%M:%S.%f")
        if i > 0:
            # check if we have next day
            if date_time[i-1] > date_time0:
                date_time0 = date_time0 + timedelta(days=1)
        date_time.append(date_time0)

    return date_time, ImpB_fon, ImpV_fon, FonB, FonV, mB, mV, Az, El, Rg


def save_jd_data(filename_jd, header,
                 date_time, ImpB_fon, ImpV_fon, FonB, FonV, mB, mV, Az, El, Rg):

    jd = [get_julian_datetime(x) for x in date_time]
    with open(filename_jd, "w") as f_jd:
        for line in header[:-1]:
            f_jd.write(line)

        f_jd.write(
            "     JD             ImpB-fon     ImpV-fon         FonB       FonV         mB      mV     Az(deg)  El(deg)  Rg(Mm)\n")

        for i in range(0, len(date_time)):
            f_jd.write(
                "{:14.8f}  ".format(jd[i]) + " {:10.3f}".format(ImpB_fon[i]) + "  {:10.3f}".format(ImpV_fon[i]) +
                "      {:8.3f}".format(FonB[i]) + "    {:8.3f}".format(FonV[i]) +
                "    {:6.3f}".format(mB[i]) + "  {:6.3f}".format(mV[i]) +
                "   {:8.3f}".format(Az[i]) + "{:8.3f}".format(El[i]) + "  {:8.3f}".format(Rg[i]) + "\n"
            )

        # {'{:13.4f}'.format(flux[i])}      {'{:8.4f}'.format(flux_err[i])}


if __name__ == "__main__":
    filename = sys.argv[1]
    wd = os.path.dirname(filename)
    base = os.path.basename(filename)
    fname, ext = os.path.splitext(base)

    header = read_header(filename)
    date_time, ImpB_fon, ImpV_fon, FonB, FonV, mB, mV, Az, El, Rg = read_data(filename, len(header))

    f_jd_name = os.path.join(wd, fname + "_jd" + ext)
    # print(f_jd_name)

    save_jd_data(f_jd_name, header, date_time, ImpB_fon, ImpV_fon, FonB, FonV, mB, mV, Az, El, Rg)
