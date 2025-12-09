from lsst.daf.butler import Butler
import numpy as np
from astropy.table import Table
import treegp
from tqdm import tqdm
import os
from smatch.matcher import Matcher
import astropy.units as units
from lsst.geom import SpherePoint, degrees

import pickle

#import lsst.afw.cameraGeom as cameraGeom
#from lsst.obs.lsst import LsstCam

__all__ = ["getDataBFVisit"]


detectorsId = [
2,
3,
4,
5,
6,
7,
8,
9,
10,
11,
12,
13,
14,
15,
16,
17,
18,
21,
22,
23,
24,
25,
26,
28,
29,
31,
32,
33,
34,
35,
36,
37,
38,
39,
40,
41,
42,
43,
44,
45,
46,
47,
48,
49,
50,
51,
52,
53,
54,
55,
56,
57,
58,
59,
60,
61,
62,
63,
64,
66,
67,
69,
70,
71,
72,
73,
74,
75,
76,
77,
78,
79,
80,
81,
82,
83,
84,
85,
86,
87,
88,
89,
90,
91,
92,
93,
94,
95,
96,
97,
98,
99,
100,
101,
102,
103,
104,
105,
106,
107,
108,
109,
110,
111,
112,
113,
114,
115,
116,
117,
118,
119,
121,
124,
125,
126,
127,
128,
129,
130,
131,
132,
133,
134,
135,
136,
137,
138,
139,
140,
141,
142,
143,
144,
145,
146,
147,
148,
149,
150,
151,
152,
153,
154,
155,
156,
157,
159,
160,
162,
163,
164,
165,
166,
167,
170,
171,
172,
173,
174,
175,
176,
177,
178,
179,
180,
181,
182,
183,
184,
185,
186,
]




columns_name = [
    'slot_Shape_xx', 'slot_Shape_yy', 'slot_Shape_xy',
    'slot_PsfShape_xx', 'slot_PsfShape_xy', 'slot_PsfShape_yy',
    'coord_ra', 'coord_dec', 'slot_Centroid_x', 'slot_Centroid_y',
    'psf_max_value', 'calib_psf_candidate', 'slot_PsfFlux_mag',
]


MODEL = "PSFEx"

def getTractsTable(butler, table):

    ra = np.array(table["coord_ra"].to(units.deg))
    dec =  np.array(table["coord_dec"].to(units.deg))
    skymap = butler.get('skyMap', skymap='lsst_cells_v1')
    tracts = []
    
    for r, d in zip(ra, dec):
        sp = SpherePoint(r, d, degrees)
        tract = skymap.findTract(sp)
        tracts.append(tract.getId())
    tracts = list(set(tracts))
    return tracts

def getDataBFVisit(collection, repOut, visit, bandId, matchFGCM=False):

    dic = {}
    butler = Butler("/repo/main", collections=collection)
    if matchFGCM:
        butlerFGCM = Butler("/repo/main", collections="LSSTCam/runs/DRP/20250604_20250921/w_2025_39/DM-52645")

    dic.update({visit:{}})

    for DETECTOR in tqdm(detectorsId):
        try:
            finalized_src_table = butler.get("single_visit_star_unstandardized", visit=visit, instrument="LSSTCam", detector=DETECTOR, parameters={"columns": columns_name})
            
            table = finalized_src_table[finalized_src_table['calib_psf_candidate']]

            table['ixx_src'] = table['slot_Shape_xx']
            table['ixy_src'] = table['slot_Shape_xy']
            table['iyy_src'] = table['slot_Shape_yy']
            
            table['ixx_psf'] = table['slot_PsfShape_xx']
            table['ixy_psf'] = table['slot_PsfShape_xy']
            table['iyy_psf'] = table['slot_PsfShape_yy']
            
            table['T_src'] = table['ixx_src'] + table['iyy_src']
            table['e1_src'] = (table['ixx_src'] - table['iyy_src']) / table['T_src']
            table['e2_src'] = 2*table['ixy_src'] / table['T_src']
            
            table['T_psf'] = table['ixx_psf'] + table['iyy_psf']
            table['e1_psf'] = (table['ixx_psf'] - table['iyy_psf']) / table['T_psf']
            table['e2_psf'] = 2*table['ixy_psf'] / table['T_psf']

            # match with FGCM:
            if matchFGCM:

                tracts =  getTractsTable(butlerFGCM, table)

                fgcm_standard_star_dict = {}
                for tract in tracts:
                    fgcmTable = butlerFGCM.get("fgcm_standard_star", instrument="LSSTCam", skymap="lsst_cells_v1", tract=tract)
                    fgcm_standard_star_dict[tract] = fgcmTable

                fgcm_standard_star_cat = []

                for tract in fgcm_standard_star_dict:
                    astropy_fgcm = fgcm_standard_star_dict[tract]
                    table_fgcm = np.asarray(astropy_fgcm)
                    fgcm_standard_star_cat.append(table_fgcm)

                fgcm_standard_star_cat = np.concatenate(fgcm_standard_star_cat)
                fgcmCat = fgcm_standard_star_cat

                raSrc = np.array(table['coord_ra'].to(units.degree))
                decSrc = np.array(table['coord_dec'].to(units.degree))

                with Matcher(raSrc, decSrc) as matcher:
                    idx, idxSrcCat, idxColorCat, d = matcher.query_radius(
                        fgcm_standard_star_cat["ra"],
                        fgcm_standard_star_cat["dec"],
                        1. / 3600.0,
                        return_indices=True,
                    )

                u_fgcm = np.zeros_like(raSrc) * np.nan
                g_fgcm = np.zeros_like(raSrc) * np.nan
                r_fgcm = np.zeros_like(raSrc) * np.nan
                i_fgcm = np.zeros_like(raSrc) * np.nan
                z_fgcm = np.zeros_like(raSrc) * np.nan
                y_fgcm = np.zeros_like(raSrc) * np.nan

                for idSrc, idColor in zip(idxSrcCat, idxColorCat):
                    u_fgcm[idSrc] = fgcmCat[f'mag_u'][idColor]
                    g_fgcm[idSrc] = fgcmCat[f'mag_g'][idColor]
                    r_fgcm[idSrc] = fgcmCat[f'mag_r'][idColor]
                    i_fgcm[idSrc] = fgcmCat[f'mag_i'][idColor]
                    z_fgcm[idSrc] = fgcmCat[f'mag_z'][idColor]
                    y_fgcm[idSrc] = fgcmCat[f'mag_y'][idColor]
            else:
                u_fgcm = None
                g_fgcm = None
                r_fgcm = None
                i_fgcm = None
                z_fgcm = None
                y_fgcm = None


            dic[visit].update({
                DETECTOR: {
                    'ixx_src': np.array(table['ixx_src']),
                    'iyy_src': np.array(table['iyy_src']),
                    'ixy_src': np.array(table['ixy_src']),
                    'ixx_psf': np.array(table['ixx_psf']),
                    'iyy_psf': np.array(table['iyy_psf']),
                    'ixy_psf': np.array(table['ixy_psf']),
                    'dT_T': np.array((table['T_src'] - table['T_psf']) / table['T_src']),
                    'de1': np.array(table['e1_src'] - table['e1_psf']),
                    'de2': np.array(table['e2_src'] - table['e2_psf']),
                    'ra': np.array(table['coord_ra'].to(units.deg)),
                    'dec': np.array(table['coord_dec'].to(units.deg)),
                    'x': np.array(table['slot_Centroid_x']),
                    'y': np.array(table['slot_Centroid_y']),
                    'detector': DETECTOR,
                    'psf_max_value': np.array(table['psf_max_value']),
                    'u_fgcm': u_fgcm,
                    'g_fgcm': g_fgcm,
                    'r_fgcm': r_fgcm,
                    'i_fgcm': i_fgcm,
                    'z_fgcm': z_fgcm,
                    'y_fgcm': y_fgcm,
                    'bandId': bandId,
                }
            })

        except:
            dic[visit].update({
                DETECTOR: None,
            })
            print(f"{visit} {DETECTOR} is not working")

    f = open(os.path.join(repOut, f'dic_{MODEL}_{visit}.pkl'), 'wb')
    pickle.dump(dic, f)
    f.close()