#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import pandas as pd
import numpy as np
import os, re, sqlite3

def gen_nubase():
    # CODATA 2018
    me = 5.48579909065e-4 # electron mass in u
    keV2u = 1.07354410233e-6 # amount of u per 1 keV
    orb_e = {0: "bare", 1: "H-like", 2: "He-like", 3: "Li-like", 4: "Be-like", 5: "B-like", 6: "C-like", 7: "N-like", 8: "O-like", 9: "F-like", 10: "Ne-like", 11: 'Na-like', 12: 'Mg-like', 13: 'Al-like', 14: 'Si-like', 15: 'P-like', 16: 'S-like', 17: 'Cl-like', 18: 'Ar-like', 19: 'K-like', 20: 'Ca-like', 21: 'Sc-like', 22: 'Ti-like', 23: 'V-like', 24: 'Cr-like', 25: 'Mn-like', 26: 'Fe-like', 27: 'Co-like', 28: 'Ni-like', 29: 'Cu-like', 30: 'Zn-like', 31: 'Ga-like', 32: 'Ge-like', 33: 'As-like', 34: 'Se-like', 35: 'Br-like', 36: 'Kr-like', 37: 'Rb-like', 38: 'Sr-like', 39: 'Y-like', 40: 'Zr-like', 41: 'Nb-like', 42: 'Mo-like', 43: 'Tc-like', 44: 'Ru-like', 45: 'Rh-like', 46: 'Pd-like', 47: 'Ag-like', 48: 'Cd-like', 49: 'In-like', 50: 'Sn-like', 51: 'Sb-like', 52: 'Te-like', 53: 'I-like', 54: 'Xe-like', 55: 'Cs-like', 56: 'Ba-like', 57: 'La-like', 58: 'Ce-like', 59: 'Pr-like', 60: 'Nd-like', 61: 'Pm-like', 62: 'Sm-like', 63: 'Eu-like', 64: 'Gd-like', 65: 'Tb-like', 66: 'Dy-like', 67: 'Ho-like', 68: 'Er-like', 69: 'Tm-like', 70: 'Yb-like', 71: 'Lu-like', 72: 'Hf-like', 73: 'Ta-like', 74: 'W-like', 75: 'Re-like', 76: 'Os-like', 77: 'Ir-like', 78: 'Pt-like', 79: 'Au-like', 80: 'Hg-like', 81: 'Tl-like', 82: 'Pb-like', 83: 'Bi-like', 84: 'Po-like', 85: 'At-like', 86: 'Rn-like', 87: 'Fr-like', 88: 'Ra-like', 89: 'Ac-like', 90: 'Th-like', 91: 'Pa-like', 92: 'U-like', 93: 'Np-like', 94: 'Pu-like', 95: 'Am-like', 96: 'Cm-like', 97: 'Bk-like', 98: 'Cf-like', 99: 'Es-like', 100: 'Fm-like', 101: 'Md-like', 102: 'No-like', 103: 'Lr-like', 104: 'Rf-like', 105: 'Db-like', 106: 'Sg-like', 107: 'Bh-like', 108: 'Hs-like', 109: 'Mt-like', 110: 'Ds-like'} # ionic charge state related to its orbital electron count
    try:
        if os.path.exists("ionic_data.db"):
            os.remove("ionic_data.db")
        # create the database
        conn = sqlite3.connect('ionic_data.db')
        c = conn.cursor()
        # create the table for ionic_data 
        c.execute('''CREATE TABLE IONICDATA
                (A          INT         NOT NULL,
                ElEMENT     CHAR(2)     NOT NULL,
                Q           INT         NOT NULL,
                Z           INT         NOT NULL,
                N           INT         NOT NULL,
                ISOMERIC    CHAR(1)     NOT NULL,     
                TYPE        TEXT        NOT NULL,
                MASS        DOUBLE,
                MASSACC     REAL,
                SOURCE      TEXT,
                JPI         TEXT,
                HALFLIFE    TEXT,
                DECAYMODE   TEXT);''')
        bind_eng = pd.read_csv("./binding_energy.csv", comment='#')
        with open("./nubase2020.txt") as nubase:
            for _ in '_'*25:
                nubase.readline()
            for l in nubase:
                A, Z, Q, isomer_state, element, mass, jpi, mass_accuracy, decay_mode = int(l[:3]), int(l[4:7]), int(l[4:7]), l[7], re.split('(\d+)', l[11:16])[-1][:2], l[18:31].split(), ','.join(l[88:102].split()), l[31:42].split(), l[119:209]
                element = element[0] if element[1]==' ' else element
                stubs = l[69:80].split()
                half_life = stubs[0].rstrip('#') if len(stubs) > 0 else "n/a"
                half_life += ' ' + stubs[1] if (half_life[-1].isdigit() and len(stubs)>1) else ""
                if len(mass) == 0:
                    mass = " "
                    mass_accuracy = " "
                elif mass[0][-1] == "#":
                    mass, source = A + float(mass[0][:-1]) * keV2u, "estimated"
                    mass_accuracy = float(mass_accuracy[0][:-1]) * keV2u 
                else:
                    mass, source = A + float(mass[0]) * keV2u, "measured"
                    mass_accuracy = float(mass_accuracy[0]) * keV2u
                while True:
                    if mass == " ":
                        break
                    if Z == 0 and Q == 0:
                        _mass =  mass
                        pass
                    elif Z==Q and (len(bind_eng[(bind_eng['Element']==element.split(" ")[0])&(bind_eng['Q']==0)]) > 0):
                        #print("{:}: {:}".format(element.split(" ")[0], Q))
                        _mass = mass - Q*me + bind_eng[(bind_eng['Element']==element.split(" ")[0])&(bind_eng['Q']==0)]["Binding Energy"].values[0]/1e3*keV2u
                    elif Q>0 and len(bind_eng[(bind_eng['Element']==element.split(" ")[0])&(bind_eng['Q']==Q)]) > 0:
                        #print("{:}: {:}".format(element.split(" ")[0], Q))
                        _mass = mass - Q*me + bind_eng[(bind_eng['Element']==element.split(" ")[0])&(bind_eng['Q']==Q)]["Binding Energy"].values[0]/1e3*keV2u
                        jpi = ''
                        decay_mode = ''
                    else:
                        break
                    c.execute("INSERT INTO IONICDATA (A,ELEMENT,Q,Z,N,ISOMERIC,TYPE,MASS,MASSACC,SOURCE,HALFLIFE, JPI, DECAYMODE)\
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (A, element, Q, Z, A-Z, isomer_state, orb_e[int(Z-Q)], _mass, mass_accuracy, source, half_life, jpi, decay_mode))
                    Q -= 1
        conn.commit()
        conn.close()
    except FileNotFoundError:
        print("Error: cannot find the files of nubase2020 and binding energy!")

if __name__ == '__main__':
    gen_nubase()
