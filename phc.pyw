#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is a main unit of Photometry curves standardization program
Dependences: wx, matplotlib, ephem, numpy
"""
import wx
import wx.xrc as xrc

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import os
import locale
import ephem
import rw  # additional units
import pu  # additional units

import numpy as np
import math as m
from datetime import datetime, timedelta
from copy import deepcopy as copy
wx = wx  # just the trick :)

glist = []
mV = []
mB = []
RA, DEC = [], []
Avs = 0
Abs = 0
# SAT=rw.group()
El, Rg = [], []
tle = False
TLE_list = []
self_path = ''  # global path to pch.pyw


def mjd(utc_time):
    dt = utc_time - datetime(2000, 1, 1, 12, 0)
    jd = (dt.days + (dt.seconds + dt.microseconds / 1000000.0) / (24 * 3600.0) + 51544.5)
    return jd  # - 2400000.5     + 2451545


def readNPS():
    """ read NPS Catalog """
    global self_path
    if os.path.exists(self_path + '//NPS.cat'):
        f = open(self_path + '//NPS.cat', 'r')
        f.readline()
        global mV, mB
        global RA, DEC
        mB = []
        mV = []
        RA, DEC = [], []
        for line in f:
            s = line.split()
            mV.append(float(s[1]))
            mB.append(float(s[2]))
            RA.append(m.degrees(float(s[4]) / 15))
            DEC.append(m.degrees(float(s[5])))
        for i in range(0, 11):
            mB[i] = mB[i] + mV[i]
        f.close()
        return 1
    else:
        print 'Can''t find NPS Catalog NPS.cat'
        return 0


def Warn(parent, message, caption='Warning!'):
    dlg = wx.MessageDialog(parent, message, caption, wx.OK | wx.ICON_WARNING)
    dlg.ShowModal()
    dlg.Destroy()


class MyApp(wx.App):
    def OnInit(self):
        if os.path.exists("phc.xrc"):
            self.res = xrc.XmlResource("phc.xrc")

            self.frame = self.res.LoadFrame(None, 'MyFrame')
            self.list_box = xrc.XRCCTRL(self.frame, "list_box_1")
            self.notebook = xrc.XRCCTRL(self.frame, "Notebook")
            self.StatusBar = xrc.XRCCTRL(self.frame, "MFrame_statusbar")

            # NPS controls
            self.lb_nps = xrc.XRCCTRL(self.frame, "lb_nps")

            self.lb_nps.SetString(0, "NPS-0=3-4")
            self.lb_nps.SetString(1, "NPS-1=5-6")
            self.lb_nps.SetString(2, "NPS-2=7-8")
            self.lb_nps.SetString(3, "NPS-3=9-10")
            self.lb_nps.SetString(4, "NPS-4=11-12")
            self.lb_nps.SetString(5, "NPS-5=13-14")
            self.lb_nps.SetString(6, "NPS-6=15-16")
            self.lb_nps.SetString(7, "NPS-7=17-18")
            # self.lb_nps.SetString(8, "NPS-8=19-20")
            # self.lb_nps.SetString(9, "NPS-9=21-22")
            # self.lb_nps.SetString(10, "NPS-10=23-24")

            self.lb_nps_res = xrc.XRCCTRL(self.frame, "lb_nps_res")
            self.tc_cb = xrc.XRCCTRL(self.frame, "tc_cb")
            self.tc_cv = xrc.XRCCTRL(self.frame, "tc_cv")
            self.rb = xrc.XRCCTRL(self.frame, "radio_box_nps")

            # Standardization controls
            self.sc_sat_gr = xrc.XRCCTRL(self.frame, "spin_ctrl_1")
            self.sc_sfon_gr = xrc.XRCCTRL(self.frame, "spin_ctrl_2")
            self.chb_mo = xrc.XRCCTRL(self.frame, "checkbox_1")
            self.chb_inst = xrc.XRCCTRL(self.frame, "checkbox_2")
            self.tc_kb = xrc.XRCCTRL(self.frame, "text_ctrl_4")
            self.tc_kv = xrc.XRCCTRL(self.frame, "text_ctrl_5")
            self.tc_cospar = xrc.XRCCTRL(self.frame, "text_ctrl_6")
            self.tc_norad = xrc.XRCCTRL(self.frame, "text_ctrl_7")
            self.tc_name = xrc.XRCCTRL(self.frame, "text_ctrl_8")
            self.cb_intr = xrc.XRCCTRL(self.frame, "combo_box_1")
            self.btn_eph = xrc.XRCCTRL(self.frame, "button_2")
            # self.btn_phc=xrc.XRCCTRL(self.frame,"btn_phot_calc")

            # Binds
            self.frame.Bind(wx.EVT_MENU, self.OnOpen, id=xrc.XRCID('Open'))
            self.frame.Bind(wx.EVT_MENU, self.OnAddNPS, id=xrc.XRCID('Add_nps'))
            self.frame.Bind(wx.EVT_MENU, self.OnQuit, id=xrc.XRCID('Exit'))
            self.frame.Bind(wx.EVT_LISTBOX_DCLICK, self.OnGraph, id=xrc.XRCID('list_box_1'))
            self.frame.Bind(wx.EVT_BUTTON, self.OnNPS_calc, id=xrc.XRCID('btn_nps_calc'))
            self.frame.Bind(wx.EVT_LISTBOX_DCLICK, self.OnListDk, id=xrc.XRCID('lb_nps'))
            self.frame.Bind(wx.EVT_CHECKBOX, self.On_m0_ch, id=xrc.XRCID('checkbox_1'))
            self.frame.Bind(wx.EVT_CHECKBOX, self.On_Ninst_ch, id=xrc.XRCID('checkbox_2'))
            self.frame.Bind(wx.EVT_CHECKBOX, self.On_NSfon, id=xrc.XRCID('checkbox_3'))
            self.frame.Bind(wx.EVT_BUTTON, self.OnNPS_graph, id=xrc.XRCID('btn_nps_graph'))

            self.Bind(wx.EVT_BUTTON, self.OnStand, id=xrc.XRCID("btn_phot_calc"))
            self.Bind(wx.EVT_BUTTON, self.OnEph, id=xrc.XRCID("button_2"))
            self.Bind(wx.EVT_KEY_DOWN, self.OnKeyLb, id=xrc.XRCID('lb_nps'))
            # self.Bind(wx.EVT_CHAR, self.OnKeyFrame)
            ######

            self.frame.Size = (500, 500)
            self.frame.Show()
            global self_path
            self_path = os.path.dirname(os.path.abspath(__file__))
            wx.CallAfter(self.list_box.SetFocus)
        else:
            print "File phc.xrc don't find"
        return True

    def OnKeyLb(self, evt):
        x = evt.GetKeyCode()
        if x == wx.WXK_RETURN:
            self.OnListDk(self)
        elif x == wx.WXK_DOWN:
            item = self.lb_nps.GetSelection()
            self.lb_nps.SetSelection(item + 1)
        elif x == wx.WXK_UP:
            item = self.lb_nps.GetSelection()
            if item > 0:
                self.lb_nps.SetSelection(item - 1)
        elif x == 6:
            self.OnOpen(self)

    def On_m0_ch(self, evt):
        print self.chb_mo.GetValue()
        # if self.chb_mo.GetValue() == True:
        #     self.btn_eph.Disable()
        # else:
        #     self.btn_eph.Enable()

    def On_Ninst_ch(self, evt):
        print self.chb_inst.GetValue()
        # if self.chb_inst.GetValue() == True:
        #     self.btn_eph.Disable()
        # else:
        #     self.btn_eph.Enable()

    def On_NSfon(self, evt):
        self.lb_nps.SetString(0, "NPS-0=2-3")
        self.lb_nps.SetString(1, "NPS-1=4-5")
        self.lb_nps.SetString(2, "NPS-2=6-7")
        self.lb_nps.SetString(3, "NPS-3=8-9")
        self.lb_nps.SetString(4, "NPS-4=10-11")
        self.lb_nps.SetString(5, "NPS-5=12-13")
        self.lb_nps.SetString(6, "NPS-6=14-15")
        self.lb_nps.SetString(7, "NPS-7=16-17")
        self.lb_nps.SetString(8, "NPS-8=18-19")

    def OnOpen(self, evt):
        self.list_box.Clear()
        self.tc_cospar.SetValue('')
        self.tc_norad.SetValue('')
        self.tc_name.SetValue('')
        wildcard = "TXT(*.txt)|*.txt;*.TXT"
        dlg = wx.FileDialog(self.frame, message="Choose File", defaultDir=os.getcwd(),
                            defaultFile='', wildcard=wildcard, style=wx.OPEN | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            for path in paths:
                global glist
                try:
                    glist = rw.read(path)
                    for N in range(1, len(glist)):
                        self.list_box.Append('Gr#' + str(N) + '-' + str(glist[N].c))
                except:
                    Warn(self.frame, "Wrong wile format probably")
            self.StatusBar.SetStatusText('FileName=' + os.path.basename(path), 0)
        dlg.Destroy()

    def OnAddNPS(self, evt):
        wildcard = "TXT(*.txt)|*.txt;*.TXT"
        dlg = wx.FileDialog(self.frame, message="Choose File", defaultDir=os.getcwd(),
                            defaultFile='', wildcard=wildcard, style=wx.OPEN | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            for path in paths:
                global glist
                NPS_list = rw.read(path)
                cc = len(glist) - 1
                for N in range(1, len(NPS_list)):
                    glist.append(NPS_list[N])
                    self.list_box.Append('NPS_Gr#' + str(N + cc) + '-' + str(NPS_list[N].c))
            self.StatusBar.SetStatusText(os.path.basename(path), 1)
        dlg.Destroy()

    def OnGraph(self, evt):
        plt.rcParams['figure.figsize'] = 12, 6
        if glist != []:
            item = self.list_box.GetSelection() + 1
            # plt.figure('Group #' + str(item))
            dt = glist[item].dt
            UT = glist[item].Time
            Data = glist[item].Date.strip()
            print UT
            j = 0
            Tt = []
            UTd = datetime.strptime(Data.replace(" ", "")+' '+UT, "%d.%m.%y %H:%M:%S")
            while j < glist[item].c:
                Tt.append(UTd+timedelta(seconds=j*dt))
                j = j + 1
            # print(Tt)
            timeFmt = DateFormatter("%H:%M:%S")
            locale.setlocale(locale.LC_NUMERIC, 'C')  # for graph
            if np.mean(glist[item].B) != 65535:
                # plt.plot(range(glist[item].c), glist[item].B, 'b-')
                plt.plot(Tt, glist[item].B, 'b-')
            if np.mean(glist[item].V) != 65535:
                # plt.plot(range(glist[item].c), glist[item].V, 'g-')
                plt.plot(Tt, glist[item].V, 'g-')
            ax = plt.gca()
            ax.xaxis.set_major_formatter(timeFmt)
            Mean_B = np.mean(glist[item].B)
            Mean_V = np.mean(glist[item].V)
            STD_B = np.std(glist[item].B, axis=0)
            STD_B = 100.0 * STD_B / Mean_B
            STD_V = np.std(glist[item].V, axis=0)
            STD_V = 100.0 * STD_V / Mean_V
            plt.title("Date=%s  UT=%s   dt=%s\n Mean_B=%5.3f  STD_B=%5.3f%%  Mean_V=%5.3f  STD_V=%5.3f%%"
                      % (Data, UT, str(dt), Mean_B, STD_B, Mean_V, STD_V))
            plt.show()

    def OnQuit(self, evt):
        self.Exit()

    def OnListDk(self, evt):
        selected = self.lb_nps.GetSelection()
        text = self.lb_nps.GetString(selected)
        renamed = wx.GetTextFromUser('Add groups', 'Add groups dialog', text)
        if renamed != '':
            self.lb_nps.SetString(selected, renamed)
        wx.CallAfter(self.lb_nps.SetFocus)

    def OnNPS_graph(self, evt):
        plt.show()

    def OnNPS_calc(self, evt):
        # print self.rb.GetStringSelection()
        if readNPS() == 0:
            Warn(self.frame, 'Can''t find NPS Catalog NPS.cat')
            return
            # exit()
        if len(glist) == 0:
            Warn(self.frame, "Open some file first")
            return
        self.lb_nps_res.Clear()
        Cb = float(self.tc_cb.GetValue())
        Cv = float(self.tc_cv.GetValue())
        nps_list = self.lb_nps.GetItems()
        NPS = range(11)
        for nps in nps_list:
            q = nps.split('=')
            if len(q) > 1:
                n = q[0].split('-')[1]
                NPS[int(n)] = q[1].split('-')
        # print NPS
        star = []
        # Ab = []
        # Av = []
        Zb, Zv = [], []
        Mza = []
        for i in range(0, 11):
            if len(NPS[i]) > 1:
                st = rw.star()
                ind1 = int(NPS[i][0])  # signal
                ind2 = int(NPS[i][1])  # fon
                st.mB = mB[i]
                st.mV = mV[i]
                st.bmv = st.mB - st.mV
                if self.rb.GetSelection() == 0:  # Graphical way
                    st.RA = RA[i]
                    st.DEC = DEC[i]
                    st.UT = glist[ind1].Time
                    st_date = glist[ind1].Date
                    # print st.RA/15, st.DEC, st.UT
                    station = ephem.Observer()
                    station.lat = '48.5635505'
                    station.long = '22.453751'
                    station.elevation = 231.1325
                    starr = ephem.FixedBody()
                    starr._ra = ephem.hours('%s' % st.RA)
                    starr._dec = ephem.degrees('%s' % st.DEC)
                    day, month, year = st_date.split('.')
                    y = int(year)
                    if y > 50:
                        y = y + 1900
                    else:
                        y = y + 2000
                    station.date = '%i/%s/%s %s' % (y, month.strip(), day.strip(), st.UT)
                    starr.compute(station)
                    # print m.degrees(starr.alt)
                    st.Mz = 1 / m.cos(m.radians(90 - m.degrees(starr.alt)))

                st.B = np.array(glist[ind1].B)
                st.V = np.array(glist[ind1].V)

                fonB = np.array(glist[ind2].B)
                st.ImpB = st.B.mean(axis=0) - fonB.mean(axis=0)
                st.RMSb = st.B.std(axis=0)

                fonV = np.array(glist[ind2].V)
                st.ImpV = st.V.mean(axis=0) - fonV.mean(axis=0)
                st.RMSv = st.V.std(axis=0)
                # print i, st.ImpV, st.RMSv

                if st.ImpB > 0:
                    Zb.append(st.mB - Cb * (st.bmv) + 2.5 * m.log10(st.ImpB))
                if st.ImpV > 0:
                    Zv.append(st.mV - Cv * (st.bmv) + 2.5 * m.log10(st.ImpV))
                Mza.append(st.Mz)
                star.append(st)
        global Abs
        global Avs
        Mzab = copy(Mza)

        if self.rb.GetSelection() == 0:  # Graphical way
            if len(Zv) > 0:
                print 'Zv=', Zv
                print 'Mza', Mza
                kv, Avs, res_max, ind = pu.lsqFit(Zv, Mza)
                while res_max > 0.1:
                    print "Residual", res_max
                    print 'delete V', ind
                    Zv = np.delete(Zv, ind)
                    Mza = np.delete(Mza, ind)
                    kv, Avs, res_max, ind = pu.lsqFit(Zv, Mza)
                self.lb_nps_res.Append('Av=%2.4f' % Avs)
                self.lb_nps_res.Append('Av_Residual=%2.4f' % res_max)
                self.lb_nps_res.Append('Kv=%2.4f' % abs(kv))
                # y=kx+A
                plt.plot(Mza, Zv, 'go')
                m_min = np.min(Mza)
                m_max = np.max(Mza)
                plt.plot([m_min, m_max], [kv * m_min + Avs, kv * m_max + Avs], 'g')

            if len(Zb) > 0:
                kb, Abs, res = pu.lsqFit(Zb, Mzab)
                while res > 0.1:
                    X = 0
                    for j in range(0, len(Zb)):
                        dz = abs(Zb[j] - (kb * Mzab[j] + Abs))
                        if dz > 0.1 and dz > X:
                            i = j
                            X = dz
                    Zb = np.delete(Zb, i)
                    print 'delete B', i
                    Mzab = np.delete(Mzab, i)
                    kb, Abs, res = pu.lsqFit(Zb, Mzab)
                self.lb_nps_res.Append('Ab=%2.4f' % Abs)
                self.lb_nps_res.Append('Ab_Residual=%2.4f' % res)
                self.lb_nps_res.Append('Kb=%2.4f' % abs(kb))
                plt.plot(Mzab, Zb, 'bo')
                m_min = np.min(Mzab)
                m_max = np.max(Mzab)
                plt.plot([m_min, m_max], [kb * m_min + Abs, kb * m_max + Abs], 'b')

            if (len(Zv) < 1 and Avs != 0) or (len(Zb) < 1 and Abs != 0):
                wx.MessageBox('You are doing something wrong - No Ab (or Av) coefficient calculated', 'Warning!!!', 1)

        if self.rb.GetSelection() == 1:  # Standard way
            if len(Zv) > 0:
                Zv = pu.RMS_del(Zv, 0.1)
                Avs = np.mean(Zv)
                res = np.std(Zv)
                self.lb_nps_res.Append('Av=%2.4f' % Avs)
                self.lb_nps_res.Append('Av_StdDev=%2.4f' % res)
            if len(Zb) > 0:
                Zb = pu.RMS_del(Zb, 0.1)
                Avs = np.mean(Zb)
                res = np.std(Zb)
                self.lb_nps_res.Append('Ab=%2.4f' % Abs)
                self.lb_nps_res.Append('Ab_StdDev=%2.4f' % res)

    def OnEph(self, evt):
        global tle
        tle = False
        El = np.array([])
        Rg = np.array([])
        wildcard = "TLE(*.tle,*.txt)|*.tle;*.TLE;*.txt;*.TXT| S3(*.s3)|*.s3"
        dlg = wx.FileDialog(
            self.frame, message="Choose ephemeris File or TLE File",
            defaultDir=os.getcwd(),
            defaultFile='',
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
        )
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            for path in paths:
                print path
                ephf = open(path, 'r')
                ext = os.path.splitext(path)
                ext = ext[1]
                if ext in [".s3", ".S3"]:
                    ephf.readline()
                    for line in ephf:
                        if line != '':
                            r = line.split()
                            El = np.append(El, float(r[2]))
                            Rg = np.append(Rg, float(r[3]))
                elif ext in [".tle", ".TLE", ".txt", ".TXT"]:
                    i = 1
                    tle1 = []
                    for l in ephf:
                        if l != '':
                            tle1.append(l.strip())
                            i = i + 1
                            if i > 3:
                                TLE_list.append(tle1)
                                tle1 = []
                                i = 1
                    tle = True
                ephf.close()
        dlg.Destroy()

    def calc_from_tle(self, count, date, time, dt, COSPAR, NORAD, NAME):
        tle = []
        for i in range(0, len(TLE_list)):
            l2 = TLE_list[i][1].split()
            cosp = l2[2]
            nor = l2[1]
            name = TLE_list[i][0]
            if cosp == COSPAR:
                c = l2[2]
                no = l2[1]
                n = TLE_list[i][0]
                tle = TLE_list[i]
        if tle == []:
            for i in range(0, len(TLE_list)):
                l2 = TLE_list[i][1].split()
                nor = l2[1]
                if nor == NORAD:
                    c = l2[2]
                    no = l2[1]
                    n = TLE_list[i][0]
                    tle = TLE_list[i]
        if tle == []:
            for i in range(0, len(TLE_list)):
                l2 = TLE_list[i][1].split()
                name = TLE_list[i][0]
                if name == NAME:
                    c = l2[2]
                    no = l2[1]
                    n = TLE_list[i][0]
                    tle = TLE_list[i]

        # Calculating  El, Rg
        station = ephem.Observer()
        station.lat = '48.5635505'
        station.long = '22.453751'
        station.elevation = 231.1325

        sat = ephem.readtle(tle[0], tle[1], tle[2])
        try:
            station.date = datetime.strptime(date.strip().replace(" ", "") + ' ' + time, "%d.%m.%y %H:%M:%S.%f")
        except Exception:
            print "Error. Trying Options #2 - time without seconds fractions"
            station.date = datetime.strptime(date.strip().replace(" ", "") + ' ' + time, "%d.%m.%y %H:%M:%S")
        el = []
        rg = []
        print station.date
        j = 1
        while j < count:
            sat.compute(station)
            el.append(m.degrees(sat.alt))
            rg.append(sat.range / 1000.0)  # km
            station.date = ephem.Date(station.date + dt / 3600 / 24)
            # print station.date, sat.alt, sat.range
            j = j + 1
        no = no[:-1]  # no U
        return el, rg, n, no, c, tle

    def OnStand(self, evt):
        global glist
        MAX_M = 15
        MAX_N = 65535
        Kb = float(self.tc_kb.GetValue())
        Kv = float(self.tc_kv.GetValue())
        COSPAR = self.tc_cospar.GetValue()
        NORAD = self.tc_norad.GetValue()
        NAME = self.tc_name.GetValue()
        if NAME == '' and NORAD == '' and COSPAR == '':
            Warn(self.frame, 'No satellite NAME or NUMBER entered!')
            return  # exit()
        if int(self.sc_sfon_gr.GetValue()) > 0:
            sn = self.sc_sat_gr.GetValue()
            fn = self.sc_sfon_gr.GetValue()
            SAT = copy(glist[sn])
            if glist[sn].c == glist[fn].c:
                SAT.B = SAT.B - glist[fn].B
                SAT.V = SAT.V - glist[fn].V
            else:  # if SAT group > fon group
                FON_B = copy(glist[fn].B)
                FON_V = copy(glist[fn].V)
                count = copy(glist[fn].c)
                ###B
                FON_B2 = pu.interp(FON_B, count, SAT.c)  # - 1
                FON_B2 = np.array(FON_B2)
                ###V
                FON_V2 = pu.interp(FON_V, count, SAT.c)  # - 1
                FON_V2 = np.array(FON_V2)
                print "Interpolation DONE!!!"
                print 'Count=', len(FON_V2)
                while len(FON_V2) > SAT.c:
                    FON_V2 = np.delete(FON_V2, -1)
                    FON_B2 = np.delete(FON_B2, -1)
                    print 'Count=', len(FON_V2)
                '''
                locale.setlocale(locale.LC_NUMERIC, 'C')#for graph
                plt.figure()
                plt.plot(range(count),FON_V,'x',range(SAT.c),FON_V2,'b')
                plt.title('Cubic-spline interpolation')
                plt.show()'''
                #
                FON_B2 = pu.RMS_mean(FON_B2, MAX_M)
                FON_V2 = pu.RMS_mean(FON_V2, MAX_M)
                SAT.B = SAT.B - FON_B2
                SAT.V = SAT.V - FON_V2
        else:  # No FON
            sn = self.sc_sat_gr.GetValue()
            SAT = copy(glist[sn])
        # Standardization
        global Abs
        global Avs
        if self.chb_inst.GetValue() is False:
            # m0
            for i in range(SAT.c):
                if SAT.B[i] == 0 or SAT.B[i] < 0.1 or SAT.B[i] == 65535:
                    SAT.B[i] = MAX_M  # **************
                else:
                    SAT.B[i] = -2.5 * m.log10(SAT.B[i]) + Abs

                if SAT.V[i] == 0 or SAT.V[i] < 0.1 or SAT.V[i] == 65535:
                    SAT.V[i] = MAX_M  # **************
                else:
                    SAT.V[i] = -2.5 * m.log10(SAT.V[i]) + Avs
        # m_z and m_ro#
        if (self.chb_mo.GetValue() == False) and (self.chb_inst.GetValue() == False):  # if ephemeris is present
            global tle
            print "tle file=", tle
            if tle:
                El, Rg, name, nor, cosp, tle_lines = self.calc_from_tle(SAT.c + 1,
                                                                        SAT.Date.strip(),
                                                                        SAT.Time,
                                                                        SAT.dt,
                                                                        COSPAR,
                                                                        NORAD,
                                                                        NAME)
                self.tc_cospar.SetValue(cosp)
                self.tc_norad.SetValue(nor)
                self.tc_name.SetValue(name)
                COSPAR = cosp
                NORAD = nor
                NAME = name
            for i in range(SAT.c):
                # print i ##m.radians
                if self.rb.GetSelection == 1:  # Mean A coefficient
                    mz = 1 / m.cos(m.radians(42)) - 1 / m.cos(m.radians(90 - El[i]))
                else:  # Graphical way
                    mz = 1 / m.cos(m.radians(90 - El[i]))
                mzb = Kb * mz
                mzv = Kv * mz
                mr = -5 * m.log10(Rg[i] / 1000.0)
                if SAT.B[i] != MAX_M:
                    SAT.B[i] = SAT.B[i] - mzb + mr  # + or-  Mz ???
                if SAT.V[i] != MAX_M:
                    SAT.V[i] = SAT.V[i] - mzv + mr
        # Graphic!!!
        if len(SAT.B) > 0:
            SAT.B = np.array(SAT.B)
            minB = SAT.B.min(axis=0)
            maxB = SAT.B.max(axis=0)
            if (maxB == 65535) and (SAT.B.mean(axis=0) == 65535):
                maxB = 0
        if len(SAT.V) > 0:
            SAT.V = np.array(SAT.V)
            minV = SAT.V.min(axis=0)
            maxV = SAT.V.max(axis=0)
            if (maxV == 65535) and (SAT.V.mean(axis=0) == 65535):
                maxV = 0
        # print 'MIN=',minB, minV
        # print 'MAX=',maxB, maxV
        # NENT is added to get time on xaxis and graph title
        jj = 0
        Tt = []
        timeFmt = DateFormatter("%H:%M:%S")
        UTd = datetime.strptime(SAT.Date.strip().replace(" ", "")+' '+SAT.Time, "%d.%m.%y %H:%M:%S")
        print UTd
        while jj < SAT.c:
            Tt.append(UTd+timedelta(seconds=jj*SAT.dt))
            jj = jj + 1
        # ##########################################################################
        plt.rcParams['figure.figsize'] = 12, 6
        Tmin = min(Tt)
        Tmax = max(Tt)
        if np.mean(SAT.B) == MAX_M:
            minB = minV
            maxB = maxV
        if np.mean(SAT.V) == MAX_M:
            minV = minB
            maxV = maxB
        if self.chb_inst.GetValue() is False:
            plt.axis([Tmin, Tmax, max(maxB, maxV), min(minB, minV)])
        else:
            plt.axis([Tmin, Tmax, min(minB, minV), max(maxB, maxV)])
        locale.setlocale(locale.LC_NUMERIC, 'C')  # for graph

        if np.mean(SAT.B) not in [MAX_M, MAX_N]:  # 65535:
            plt.plot(Tt, SAT.B, 'b-', label="B")
        if np.mean(SAT.V) not in [MAX_M, MAX_N]:  # 65535:
            plt.plot(Tt, SAT.V, 'g-', label="V")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(timeFmt)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title("Satellite=%s \n Date=%s  UT=%s   dt=%s" % (NAME, SAT.Date.strip(), SAT.Time, str(SAT.dt)))
        plt.show()

        if NORAD == '':
            # TryFindNames(NAME) # write!!!!!!!!!!!!!!!!!!
            tle = []
            if tle == []:
                for i in range(0, len(TLE_list)):
                    l2 = TLE_list[i][1].split()
                    name = TLE_list[i][0]
                    if name == NAME:
                        COSPAR = l2[2]
                        NORAD = l2[1][:-1]
                        NAME = TLE_list[i][0]
                        tle = TLE_list[i]

        dlg = wx.MessageDialog(self.frame, 'Save file?', 'Save dialog', wx.YES | wx.NO)
        if dlg.ShowModal() == wx.ID_YES:
            ''' SAVE RESULT'''
            # print "PATH=", os.getcwd() + '/' + NAME
            wildcard = "PHC(*.phc)|*.phc| UMOSS(*.umoss)|*.umoss| PH1(*.ph1)|*.ph1"
            if not os.path.exists(os.getcwd() + '/' + NAME):
                os.makedirs(os.getcwd() + '/' + NAME)
            dlgS = wx.FileDialog(
                self.frame, message="Save File in different formats",
                defaultDir=os.getcwd() + '/' + NAME,
                defaultFile=NAME + '.phc',
                wildcard=wildcard,
                style=wx.SAVE | wx.CHANGE_DIR
            )
            if dlgS.ShowModal() == wx.ID_OK:
                ffpath = dlgS.GetPaths()[0]
                if dlgS.GetFilterIndex() == 0 and not ffpath.endswith(".phc"):
                    ffpath += ".phc"
                if dlgS.GetFilterIndex() == 1 and not ffpath.endswith(".umoss"):
                    ffpath += ".umoss"
                if dlgS.GetFilterIndex() == 2 and not ffpath.endswith(".ph1"):
                    ffpath += ".ph1"
                ext = os.path.splitext(ffpath)
                ext = ext[1]
                # if NORAD == '':
                #     #TryFindNames(NAME) # write!!!!!!!!!!!!!!!!!!
                #     tle = []
                #     if tle == []:
                #         for i in range(0, len(TLE_list)):
                #             l2 = TLE_list[i][1].split()
                #             name = TLE_list[i][0]
                #             if name == NAME:
                #                 COSPAR = l2[2]
                #                 NORAD = l2[1]
                #                 NAME = TLE_list[i][0]
                #                 tle = TLE_list[i]
                if ext == '.umoss':
                    p = os.path.dirname(ffpath)
                    d0 = datetime.strptime(SAT.Date.strip().replace(" ", "") +
                                           ' ' + SAT.Time, "%d.%m.%y %H:%M:%S")
                    d1 = d0.strftime("%Y.%m.%d - %H:%M:%S")
                    d2 = d0 + timedelta(seconds=SAT.c * SAT.dt)
                    d2 = d2.strftime("%Y.%m.%d - %H:%M:%S")
                    # B chanel
                    ff1 = p + '//' + NORAD + '_' + d0.strftime("%y%m%d_%H%M") + 'B.U11'
                    ff2 = p + '//' + NORAD + '_' + d0.strftime("%y%m%d_%H%M") + 'V.U11'
                    f1 = open(ff1, 'w')
                    f2 = open(ff2, 'w')
                    f1.write('PH3 U11 B 80 0.005\n')
                    f2.write('PH3 U11 V 80 0.005\n')
                    cy = COSPAR[0:2]
                    print 'cy=', cy
                    if int(cy) > 50:
                        cy = '19' + cy
                    else:
                        cy = '20' + cy
                    cn = COSPAR[2:5]
                    ct = COSPAR[5:]
                    f1.write("%s-%s%-3s  %s  %s\n" % (cy, cn, ct, NORAD, NAME))
                    f2.write("%s-%s%-3s  %s  %s\n" % (cy, cn, ct, NORAD, NAME))
                    f1.write(d1 + '\n')
                    f1.write(d2 + '\n')
                    f2.write(d1 + '\n')
                    f2.write(d2 + '\n')
                    if tle:
                        f1.write('%\n%' + tle_lines[0] + '\n')
                        f1.write('%' + tle_lines[1] + '\n')
                        f1.write('%' + tle_lines[2] + '\n%\n')
                        f2.write('%\n%' + tle_lines[0] + '\n')
                        f2.write('%' + tle_lines[1] + '\n')
                        f2.write('%' + tle_lines[2] + '\n%\n')
                    else:
                        f1.write('% \n')
                        f2.write('% \n')
                    f1.write('Time (MJD)        Mst\n')
                    f2.write('Time (MJD)        Mst\n')
                    for i in range(SAT.c):
                        tt = i * SAT.dt
                        T = d0 + timedelta(seconds=tt)
                        f1.write('%5.8f  %2.3f\n' % (mjd(T), SAT.B[i]))
                        f2.write('%5.8f  %2.3f\n' % (mjd(T), SAT.V[i]))
                    f1.write('END')
                    f2.write('END')
                    f1.close()
                    f2.close()
                if ext == '.ph1':
                    p = os.path.dirname(ffpath)
                    d0 = datetime.strptime(SAT.Date.strip().replace(" ", "") +
                                           ' ' + SAT.Time, "%d.%m.%y %H:%M:%S")
                    d1 = d0.strftime("%Y.%m.%d - %H:%M:%S")
                    d2 = d0 + timedelta(seconds=SAT.c * SAT.dt)
                    d2 = d2.strftime("%Y.%m.%d - %H:%M:%S")
                    # B chanel
                    ff1 = p + '//' + NORAD + '_' + d0.strftime("%y%m%d_%H%M") + 'B.U11'
                    ff2 = p + '//' + NORAD + '_' + d0.strftime("%y%m%d_%H%M") + 'V.U11'
                    f1 = open(ff1, 'w')
                    f2 = open(ff2, 'w')
                    f1.write('PH1 U11 B 80 0.005\n')
                    f2.write('PH1 U11 V 80 0.005\n')
                    if tle:
                        cy = COSPAR[0:2]
                        if int(cy) > 50:
                            cy = '19' + cy
                        else:
                            cy = '20' + cy
                        cn = COSPAR[2:5]
                        ct = COSPAR[5:]
                        f1.write("%s-%s%-3s  %s  %s\n" % (cy, cn, ct, NORAD, NAME))
                        f2.write("%s-%s%-3s  %s  %s\n" % (cy, cn, ct, NORAD, NAME))
                    else:
                        f1.write("%s  %s\n" % (NORAD, NAME))
                        f2.write("%s  %s\n" % (NORAD, NAME))
                    f1.write(d1 + '\n')
                    f1.write(d2 + '\n')
                    f2.write(d1 + '\n')
                    f2.write(d2 + '\n')
                    f1.write('% \n')
                    f2.write('% \n')
                    f1.write('Time (MJD)        N_inst\n')
                    f2.write('Time (MJD)        N_inst\n')
                    for i in range(SAT.c):
                        tt = i * SAT.dt
                        T = d0 + timedelta(seconds=tt)
                        f1.write('%5.8f  %2.3f\n' % (mjd(T), SAT.B[i]))
                        f2.write('%5.8f  %2.3f\n' % (mjd(T), SAT.V[i]))
                    f1.write('END')
                    f2.write('END')
                    f1.close()
                    f2.close()
                else:  # UzhNU format
                    print 'Saving in LKD UzhNU format'
                    f = open(ffpath, 'w')
                    T0 = datetime.strptime(SAT.Date.strip().replace(" ", "")
                                           + ' ' + SAT.Time, '%d.%m.%y %H:%M:%S')
                    f.write(str(T0) + '\n')
                    f.write(str(T0 + timedelta(seconds=SAT.c * SAT.dt)) + '\n')
                    f.write('dt=' + str(SAT.dt) + '\n')
                    f.write('COSPAR ID=' + COSPAR + '\n')
                    f.write('NORAD ID=' + NORAD + '\n')
                    f.write('NAME=' + NAME + '\n')
                    if self.chb_mo.GetValue():
                        f.write("NO standardization for mZ an Range... only m0 magnitude is given!!! \n")
                    if self.chb_inst.GetValue():
                        f.write("NO magnitudes... only instrument count is given!!! \n")

                    if not self.chb_inst.GetValue():
                        f.write('  UT TIME            mB      mV      El(deg) Rg(Mm) \n')
                    else:
                        f.write('  UT TIME           N_B     N_V      El(deg) Rg(Mm) \n')
                    for i in range(SAT.c):
                        tt = i * SAT.dt
                        T = T0 + timedelta(seconds=tt)
                        f.write(T.strftime('%H:%M:%S.%f'))
                        f.write(str("  %8.3f" % SAT.B[i]))
                        f.write(str("%8.3f" % SAT.V[i]))
                        if (not self.chb_mo.GetValue()) and (not self.chb_inst.GetValue()):
                            f.write("    %5.3f  %5.3f \n" % (El[i], Rg[i]))
                        else:
                            f.write('\n')
                    f.close()
            dlgS.Destroy()
        dlg.Destroy()


if __name__ == "__main__":
    app = MyApp(False)
    app.MainLoop()
