#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import numpy as np
import re, sqlite3, os

from gen_nubase import gen_nubase

class IID():
    '''
    A script for auto calibrating the ion identification result based on the input Schottky spectrum
    '''
    # CODATA 2018
    c = 299792458 # speed of light in m/s
    e = 1.602176634e-19 # elementary charge in C
    me = 5.48579909065e-4 # electron mass in u
    u2kg = 1.66053906660e-27 # amount of kg per 1 u
    MeV2u = 1.07354410233e-3 # amount of u per 1 MeV

    # time unit
    time_unit = {'Yy': 31536000*1e24, 'Zy': 31536000*1e21, 'Ey': 31536000*1e18, 'Py': 31536000*1e15, 'Ty': 31536000*1e12, 'Gy': 31536000*1e9, 'My': 31536000*1e6, 'ky': 31536000*1e3, 'y': 31536000, 'd': 86400, 'h': 3600, 'm': 60, 's': 1, 'ms': 1e-3, 'us': 1e-6, 'ns': 1e-9, 'ps': 1e-12, 'fs': 1e-15, 'as': 1e-18, 'zs': 1e-21, 'ys': 1e-24}

    def __init__(self, lppion, cen_freq, span, gamma_t, delta_Brho_over_Brho, gamma_setting, delta_v_over_v=1e-6, L_CSRe=128.8, nubase_update=False, verbose=False):
        '''
        extract all the secondary fragments and their respective yields calculated by LISE++
        (including Mass, Half-life, Yield of all the fragments)
        lppion:     LISE++ output file to be loaded
        cen_freq:   center frequency of the spectrum in MHz
        span:       span of the spectrum in kHz
        L_CSRe:     circumference of CSRe in m, default value 128.8
        gamma_t:                the gamma_t value for isochronous mode
        delta_Brho_over_Brho:   the ΔBrho/Brho for isochronous mode, [%]
        gamma_setting:          the velocity (gamma) for the e-cooler setting
        delta_v_over_v:         the Δv/v for the e-cooler setting, represent the e-cooler capacity of cooling
        '''
        self.cen_freq = cen_freq # MHz
        self.span = span # kHz
        self.gamma_t = gamma_t
        self.delta_Brho_over_Brho = delta_Brho_over_Brho # %
        self.gamma_setting = gamma_setting
        self.delta_v_over_v = delta_v_over_v # %
        self.L_CSRe = L_CSRe # m
        self.verbose = verbose
        # check database ionic_data.db exists, if not create one
        if nubase_update or ((not nubase_update) and (not os.path.exists("./ionic_data.db"))):
            gen_nubase()
        self.conn = sqlite3.connect("./ionic_data.db")
        self.cur = self.conn.cursor()
        # check table ioncidata exists
        if self.cur.execute("SELECT count(name) FROM sqlite_master WHERE type='table' and name='IONICDATA'").fetchone()[0] == 0:
            gen_nubase()

        self.cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
        table_name = [(part[0],) for part in self.cur.fetchall() if (part[0] != 'IONICDATA')]
        if len(table_name) > 0:
            self.cur.executescript(';'.join(["DROP TABLE IF EXISTS %s" %i for i in table_name]))
            self.conn.commit()
        self.cur.execute('''CREATE TABLE IF NOT EXISTS LPPDATA
                (A          INT         NOT NULL,
                ElEMENT     CHAR(2)     NOT NULL,
                Q           INT         NOT NULL,
                ION         TEXT        NOT NULL,
                YIELD       REAL);''')
        self.cur.execute("DELETE FROM LPPDATA")
        self.conn.commit()
        with open(lppion, encoding='latin-1') as lpp:
            while True:
                line = lpp.readline().strip()
                if line == "[D4_DipoleSettings]":
                    self.Brho = float(lpp.readline().strip().split()[2]) # Tm
                elif line == "[Calculations]":
                    break
            for line in lpp:
                segment = line.strip().split(',')[0].split()
                A, element, Q = re.split("([A-Z][a-z]?)", segment[0]+segment[1][:-1])
                self.cur.execute("INSERT INTO LPPDATA (A,ELEMENT,Q,ION,YIELD) VALUES (?,?,?,?,?)", (A, element, Q, ''.join([A,element,Q]), segment[-1][1:]))
            self.conn.commit()
        # reset yield for both fission (including pps for AFhihg, AFmid and AFlow) and PF processing
        result = self.cur.execute("SELECT sum(YIELD), ION FROM LPPDATA GROUP BY ION").fetchall()
        self.cur.executemany("UPDATE LPPDATA SET YIELD=? WHERE ION=?", result)
        self.conn.commit()
        self.cur.execute("CREATE TABLE TEMPTABLE as SELECT DISTINCT * FROM LPPDATA")
        self.cur.execute("DROP TABLE LPPDATA")
        self.cur.execute("ALTER TABLE TEMPTABLE RENAME TO LPPDATA")
        self.conn.commit()
        self.calc_isochronous_peak()

    def prepare_result(self):
        self.cur.execute("SELECT count(name) FROM sqlite_master WHERE type='table' AND name='OBSERVEDION'")
        if self.cur.fetchone()[0] == 1:
            self.cur.execute("DROP TABLE OBSERVEDION")
            self.conn.commit()
        self.cur.execute('''CREATE TABLE IF NOT EXISTS OBSERVEDION
                (ID         INT,
                A          INT          NOT NULL,
                ElEMENT     CHAR(2)     NOT NULL,
                Q           INT         NOT NULL,
                Z           INT         NOT NULL,
                N           INT         NOT NULL,
                MASS        DOUBLE,
                SOURCE      TEXT,                
                ION         TEXT        NOT NULL,
                YIELD       REAL,     
                TYPE        TEXT,
                ISOMERIC    CHAR(1),
                HALFLIFE    TEXT,
                PRIMARY KEY (ID));''')
        self.cur.execute("DELETE FROM OBSERVEDION")
        self.cur.execute("INSERT INTO OBSERVEDION(A,ELEMENT,Q,Z,N,MASS,SOURCE,ION,YIELD,TYPE,ISOMERIC,HALFLIFE) \
                SELECT LPPDATA.A, LPPDATA.ELEMENT, LPPDATA.Q, IONICDATA.Z, IONICDATA.N, IONICDATA.MASS, IONICDATA.SOURCE, LPPDATA.ION, LPPDATA.YIELD, IONICDATA.TYPE, IONICDATA.ISOMERIC, IONICDATA.HALFLIFE \
                FROM IONICDATA \
                INNER JOIN LPPDATA ON IONICDATA.Q=LPPDATA.Q AND IONICDATA.ELEMENT=LPPDATA.ELEMENT AND IONICDATA.A=LPPDATA.A")
        # reset the yield of the isometric_state
        result = self.cur.execute("SELECT YIELD, ION, ISOMERIC FROM OBSERVEDION WHERE ISOMERIC!=0").fetchall()
        ## yield(isometric_state) = yield(bare) * 10**(-isometric_state)
        #re_set = [(item[0]*10**(-int(item[2])), item[1], item[2]) for item in result]
        # yield(isometric_state) = yield(bare) * 1/2
        re_set = [(item[0]*0.5, item[1], item[2]) for item in result]
        self.cur.executemany("UPDATE OBSERVEDION SET YIELD=? WHERE ION=? AND ISOMERIC=?", re_set)
        self.conn.commit()

    def calc_isochronous_peak(self):
        '''
        calculate peak locations of the Schottky signals from secondary fragments visible in the pre-defined frequency range (isochronous mode)
        gamma_t:                the gamma_t value for isochronous mode
        delta_Brho_over_Brho:   the ΔBrho/Brho for isochronous mode, [%]
        '''
        self.prepare_result()
        # initial table isochronousion
        self.cur.execute("SELECT count(name) FROM sqlite_master WHERE type='table' AND name='ISOCHRONOUSION'")
        if self.cur.fetchone()[0] == 1:
            self.cur.execute("DROP TABLE ISOCHRONOUSION")
            self.conn.commit()
        self.cur.execute('''CREATE TABLE IF NOT EXISTS ISOCHRONOUSION
                (ION         TEXT        NOT NULL,
                ISOMERIC    CHAR(1),
                MASS        DOUBLE,
                GAMMA       REAL,
                SOURCE      TEXT,                
                YIELD       REAL,
                PEAKLOC     REAL,
                PEAKWIDTH   REAL,
                PEAKHEIGHT  REAL,
                HARMONIC    INT,
                REVFREQ     REAL,     
                TYPE        TEXT,
                HALFLIFE    TEXT,
                WEIGHT      DOUBEL);''')
        self.conn.commit()
        self.cur.execute("SELECT MASS, Q, YIELD, ION, ISOMERIC FROM OBSERVEDION")
        mass, Q, ion_yield, ion, isometric_state = np.array(self.cur.fetchall()).T
        mass, Q, ion_yield = mass.astype(np.float64), Q.astype(np.float64), ion_yield.astype(np.float64)
        gamma_beta = self.Brho / mass * Q / self.c / self.u2kg * self.e
        beta = gamma_beta / np.sqrt(1 + gamma_beta**2)
        gamma = 1 / np.sqrt(1 - beta**2)
        #energy = (gamma - 1) / self.MeV2u # MeV/u
        rev_freq = beta * self.c / self.L_CSRe / 1e6 # MHz
        weight = ion_yield * Q**2 * rev_freq**2
        lower_freq, upper_freq = self.cen_freq - self.span/2e3, self.cen_freq + self.span/2e3 # MHz, MHz
        i = 0
        while (i < len(ion)):
            temp = self.cur.execute("SELECT ION,ISOMERIC,MASS,SOURCE,YIELD,TYPE,HALFLIFE FROM OBSERVEDION WHERE ION=? AND ISOMERIC=?", (ion[i],isometric_state[i])).fetchone()
            # drop the ion of life time < 10 ms
            if self.life_transform(temp[6]) <= 10 * self.time_unit['ms']:
                i += 1
                continue
            harmonics = np.arange(np.ceil(lower_freq/rev_freq[i]), np.floor(upper_freq/rev_freq[i])+1).astype(int)
            peak_width = np.abs(1 / gamma[i]**2 - 1 / self.gamma_t**2) * self.delta_Brho_over_Brho *1e-2 * rev_freq[i] * harmonics if gamma[i] != self.gamma_t else  np.nonzero(np.abs(1 / gamma**2 - 1 / self.gamma_t**2)).min()*0.5 * self.delta_Brho_over_Brho *1e-2 * rev_freq[i] * harmonics
            peak_height = weight[i] / peak_width
            # filter harmonics
            if len(harmonics) == 1:
                self.cur.execute("INSERT INTO ISOCHRONOUSION(ION,ISOMERIC,MASS,SOURCE,YIELD,TYPE,HALFLIFE,GAMMA,WEIGHT,PEAKLOC,PEAKWIDTH,PEAKHEIGHT,REVFREQ,HARMONIC) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (*temp, gamma[i], weight[i], (harmonics[-1]*rev_freq[i]-self.cen_freq)*1e3, peak_width[-1], peak_height[-1], rev_freq[i], int(harmonics[-1])))
            elif len(harmonics) > 1:
                re_set = [(*temp, gamma[i], weight[i], (harmonics[j]*rev_freq[i]-self.cen_freq)*1e3, peak_width[j], peak_height[j], rev_freq[i], int(harmonics[j])) for j in range(len(harmonics))]
                self.cur.executemany("INSERT INTO ISOCHRONOUSION(ION,ISOMERIC,MASS,SOURCE,YIELD,TYPE,HALFLIFE,GAMMA,WEIGHT,PEAKLOC,PEAKWIDTH,PEAKHEIGHT,REVFREQ,HARMONIC) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", re_set)
            else:
                pass
            i += 1
        self.conn.commit()
                
    def calc_ecooler_peak(self):
        '''
        calculate peak locations of the Schottky signals from secondary fragments visible in the pre-defined frequency range (ecooler on)
        gamma_t:                the gamma_t value for the ring mode
        delta_Brho_over_Brho:   the δBrho/Brho for the ring mode
        gamma_setting:          the velocity (gamma) for the e-cooler setting
        delta_v_over_v:         the Δv/v for the e-cooler setting, represent the e-cooler capacity of cooling
        '''
        self.cur.execute("SELECT count(name) FROM sqlite_master WHERE type='table' AND name='ECOOLERION'")
        if self.cur.fetchone()[0] == 1:
            self.cur.execute("DROP TABLE ECOOLERION")
            self.conn.commit()
        self.cur.execute('''CREATE TABLE IF NOT EXISTS ECOOLERION
                (ION         TEXT        NOT NULL,
                ISOMERIC    CHAR(1),
                MASS        DOUBLE,
                PSEUDOGAMMA REAL,      
                SOURCE      TEXT,
                YIELD       REAL,
                PEAKLOC     REAL,
                PEAKWIDTH   REAL,
                PEAKHEIGHT  REAL,
                HARMONIC    INT,
                REVFREQ     REAL,     
                TYPE        TEXT,
                HALFLIFE    TEXT,
                WEIGHT      DOUBEL);''')
        self.conn.commit()
        self.cur.execute("SELECT ION, ISOMERIC FROM OBSERVEDION")
        ion, isometric_state = np.array(self.cur.fetchall()).T
        beta = np.sqrt(1 - 1/self.gamma_setting**2)
        upper_mass_over_Q = self.Brho * (1 + self.delta_Brho_over_Brho * 0.5e-2) / self.gamma_setting / beta / self.c / self.u2kg * self.e
        lower_mass_over_Q = self.Brho * (1 - self.delta_Brho_over_Brho * 0.5e-2) / self.gamma_setting / beta / self.c / self.u2kg * self.e
        lower_freq, upper_freq = self.cen_freq - self.span/2e3, self.cen_freq + self.span/2e3 # MHz, MHz
        # load the result into tale ecoolerion
        re_set = []
        i = 0
        while (i < len(ion)):
            temp = self.cur.execute("SELECT ION,ISOMERIC,MASS,SOURCE,YIELD,TYPE,HALFLIFE,Q FROM OBSERVEDION WHERE ION=? AND ISOMERIC=?", (ion[i],isometric_state[i])).fetchone()
            # drop the ion of life time < 1 s
            if self.life_transform(temp[6]) <= 1 * self.time_unit['s']:
                i += 1
                continue
            ion_yield, mass, Q = float(temp[4]), float(temp[2]), float(temp[-1])
            # filter mass / Q
            if mass / Q <= upper_mass_over_Q and mass / Q >= lower_mass_over_Q:
                ion_Brho = mass / Q * self.gamma_setting * beta * self.c * self.u2kg / self.e
                ion_L = self.L_CSRe * (1 + (ion_Brho-self.Brho)/self.Brho/self.gamma_t**2)
                ion_rev_freq = beta * self.c / ion_L *1e-6 # MHz
                ion_weight = ion_yield * Q**2 * ion_rev_freq
                harmonics = np.arange(np.ceil(lower_freq/ion_rev_freq), np.floor(upper_freq/ion_rev_freq)+1).astype(int)
                peak_width = np.abs(1 - self.gamma_setting**2 / self.gamma_t**2) * self.delta_v_over_v * ion_rev_freq * harmonics if self.gamma_setting != self.gamma_t else 1e-5 * self.delta_v_over_v * ion_rev_freq * harmonics
                peak_height = ion_weight / peak_width
                ion_pseudo_gamma_beta = self.Brho / mass * Q /self.c / self.u2kg * self.e
                ion_pseudo_beta = ion_pseudo_gamma_beta / np.sqrt(1 + ion_pseudo_gamma_beta**2)
                ion_pseudo_gamma = 1 / np.sqrt(1 - ion_pseudo_beta**2)
                # filter harmonics
                if len(harmonics) == 1:
                    self.cur.execute("INSERT INTO ECOOLERION(ION,ISOMERIC,MASS,SOURCE,YIELD,TYPE,HALFLIFE,WEIGHT,PEAKLOC,PEAKWIDTH,PEAKHEIGHT,REVFREQ,HARMONIC,PSEUDOGAMMA) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (*temp[:-1], ion_weight, (harmonics[-1]*ion_rev_freq-self.cen_freq)*1e3, peak_width[-1], peak_height[-1], ion_rev_freq, int(harmonics[-1]), ion_pseudo_gamma))
                elif len(harmonics) > 1:
                    re_set = [(*temp[:-1], ion_weight, (harmonics[j]*ion_rev_freq-self.cen_freq)*1e3, peak_width[j], peak_height[j], ion_rev_freq, int(harmonics[j]), ion_pseudo_gamma) for j in range(len(harmonics))]
                    self.cur.executemany("INSERT INTO ECOOLERION(ION,ISOMERIC,MASS,SOURCE,YIELD,TYPE,HALFLIFE,WEIGHT,PEAKLOC,PEAKWIDTH,PEAKHEIGHT,REVFREQ,HARMONIC,PSEUDOGAMMA) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", re_set)
                else:
                    pass
            i += 1
        self.conn.commit()

    def update_gamma_t(self, gamma_t, ec_on=False):
        '''
        using the specific gamma_t to update whole data
        '''
        self.gamma_t = gamma_t
        self.calc_isochronous_peak()
        if ec_on:
            self.calc_ecooler_peak()

    def update_delta_Brho_over_Brho(self, delta_Brho_over_Brho, ec_on=False):
        '''
        using the specific ΔΒρ/Βρ to update whole data
        '''
        self.delta_Brho_over_Brho = delta_Brho_over_Brho
        self.calc_isochronous_peak()
        if ec_on:
            self.calc_ecooler_peak()

    def update_cen_freq(self, cen_freq, ec_on=False):
        '''
        using the specific center frequency [MHz] to update whole data
        '''
        self.cen_freq = cen_freq # [MHz]
        self.calc_isochronous_peak()
        if ec_on:
            self.calc_ecooler_peak()

    def update_span(self, span, ec_on=False):
        '''
        using the specific span [kHz] to update whole data
        '''
        self.span = span # [kHz]
        self.calc_isochronous_peak()
        if ec_on:
            self.calc_ecooler_peak()

    def calibrate_Brho(self, Brho, ec_on=False):
        '''
        using the measured Brho with the identified ion to calibrate
        Brho:       the magnetic rigidity of the target ion in Tm
        '''
        self.Brho = Brho
        self.calc_isochronous_peak()
        if ec_on:
            self.calc_ecooler_peak()

    def calibrate_peak_loc(self, ion, isometric_state, peak_loc, harmonic):
        '''
        using the measured peak location with the identified ion to calibrate the magnetic rigidity of CSRe
        ion:        a string in the format of AElementQ, e.g., 3He2 
        isometric_state: a integer of isometric state, e.g., 0
        peak_loc:   peak location in kHz after deduction of the center frequency
        harmonic:   harmonic number
        '''
        mass, Q = self.cur.execute("SELECT MASS, Q FROM OBSERVEDION WHERE ION=? AND ISOMERIC=?", (ion, isometric_state)).fetchone()
        rev_freq = (self.cen_freq + peak_loc/1e3) / harmonic # MHz
        beta = rev_freq * self.L_CSRe / self.c * 1e6
        gamma = 1 / np.sqrt(1 - beta**2)
        self.Brho = gamma * beta * mass / Q * self.c * self.u2kg / self.e # Tm
        return self.Brho

    def calibrate_ecooler(self, gamma_setting, delta_v_over_v):
        '''
        using the specific gamma_setting, Δv/v with the identified ion for calibrate
        gamma_setting:          the velocity (gamma) for the e-cooler setting
        delta_v_over_v:         the Δv/v for the e-cooler setting, represent the e-cooler capacity of cooling, [%]
        '''
        self.gamma_setting = gamma_setting
        self.delta_v_over_v = delta_v_over_v
        self.calc_ecooler_peak()
    
    def life_transform(self, half_life):
        '''
        transform the halflife in SQLite Table to lifetime [s]
        '''
        temp = half_life.split() # input half life in str
        if len(temp) == 1:
            if temp[0] == 'stbl':
                return 1e64
            else:
                return -1
        else:
            if temp[-1].isalpha():
                try:
                    temp_0 = float(temp[0]) / np.log(2) * self.time_unit[temp[-1]]
                except:
                    temp_0 = float(''.join(re.split("(\d+)", temp[0])[1:])) / np.log(2) * self.time_unit[temp[-1]]
                return temp_0
            else:
                return -1


if __name__ == "__main__":
    iid = IID("./238U92.lpp", 243.5, 3000)
