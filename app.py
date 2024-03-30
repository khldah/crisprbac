from flask import Flask,  request, url_for, redirect, render_template, json
import mysql.connector as mysql
from string import *
import regex
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
import re
import math
import pickle
import numpy as np
import os
import subprocess
from subprocess import Popen, PIPE, DEVNULL

app = Flask(__name__)

# remove cache
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
app.config["CACHE_TYPE"] = "null"

# DEFINE route to main.html
@app.route('/')
def index():
    return render_template('main.html')


# DEFINE route to help.html
@app.route('/help', methods=['GET', 'POST'])
def help():
    return render_template('help.html')

# DEFINE route to help.html
@app.route('/sample_result', methods=['GET', 'POST'])
def sample_result():
    return render_template('sample_result.html')


# DEFINE route to result.html and process form
@app.route('/result', methods=['GET', 'POST'])
def result():
    if request.method == 'POST':
        
        view_result = ""
        
        job_name = request.form.get("jobname")
        target_genome = request.form.get("targetgenome")
        target_sequence = request.form.get("sequence")
        target_seq = target_sequence.replace('\n', '').replace(' ', '').replace('\r', '')
        pamseq = request.form.get('pamsequence')
        
        # find pam sequence in target sequence
        pam = target_seq.find(pamseq)
        
        # retrieve genome sequence from database
        db = mysql.connect(host="localhost", user="root",
            passwd="", db="crisprbac")
        cursor = db.cursor()
        cursor.execute("SELECT sequence FROM genome WHERE refseq_id = %s", [target_genome])
        genome_seq = cursor.fetchone() 
    
        
        # ############################################################
        
        # REVERSE COMPLEMENTS A GIVEN STRING #########################################
        def revcom(s):
            basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
            letters = list(s[::-1])
            letters = [basecomp[base] for base in letters]
            return ''.join(letters)

        # UNPICKLE MISMATCH SCORES AND PAM SCORES #####################################
        def get_mm_pam_scores():
            try:
                mm_scores = pickle.load(open('mismatch_score.pkl','rb'))
                pam_scores = pickle.load(open('pam_scores.pkl','rb'))
                return (mm_scores,pam_scores)
            except:
                raise Exception("Could not find file with mismatch scores or PAM scores")

        # CALCULATE CFD SCORE ########################################################
        def calc_cfd(grna2,offt2,pam_offt):
            mm_scores,pam_scores = get_mm_pam_scores()
            score = 1
            offt2 = offt2.replace('T','U')
            grna2 = grna2.replace('T','U')
            offt_list = list(offt2)
            grna_list = list(grna2)
            for i,ol in enumerate(offt_list):
                if grna_list[i] == ol:
                    score*=1
                else:
                    key = 'r'+grna_list[i]+':d'+revcom(ol)+','+str(i+1)
                    score*= mm_scores[key]
            score*=pam_scores[pam_offt]
            return (score)

        #  VIEW RESULT ################################################################
        
        view_result += "<h3>Job name: %s</h3>" % (job_name)

        view_result += '''

                <table id="grnadesign">
                    <thead>
                    <tr>
                        <th>gRNA sequence</th>
                        <th>Position </th>
                        <th>GC_20nt (%) </th>
                        <th>Tm_20nt (&deg;C) </th>
                        <th>Model 3 Score </th>
                        <th>Model 4 score</th>
                        <th>Off-target </th>
                        </tr>
                    </thead>
        '''

        while pam != -1:
            seq_npam = target_seq[:pam]
            ngg = seq_npam[-1]
            grna = seq_npam[-21:-1]
            grna2 = grna + ngg + pamseq
            grna3 = grna + ngg

            length = len(grna)
            if length == 20:
                view_result += '''<tr class="expandable">'''

                # Column candidate gRNA
                view_result += "<td>"
                view_result += "%.20s <mark>%s%s</mark>" % (grna, ngg, pamseq)
                view_result += "</td>"

                # Column position
                view_result += "<td> "
                for y in genome_seq:
                    grna_loc = re.search(grna, y)
                    
                    if grna_loc is None:
                        view_result += "NA"
                    else:
                        # view_result += '%d - %d' % (grna_loc.start(), grna_loc.end())
                        view_result += '%s - %s' % (str(grna_loc.start()), str(grna_loc.end()))
                        
                view_result += "</td>"

                # evaluate top 150 features #######################################
                nonseed = grna[0:9]  # from P1-P9
                seed = grna[9:]  # from P10-P20

                # gc content
                gc20nt = GC(grna)
                gcseed = GC(seed)
                gcnon = GC(nonseed)

                if 40 <= gc20nt <= 60:
                    gc20nt_bin = 1
                else:
                    gc20nt_bin = 0

                if 40 <= gcseed <= 60:
                    gcseed_bin = 1
                else:
                    gcseed_bin = 0

                if 45 <= gcnon <= 66.7:
                    gcnon_bin = 1
                else:
                    gcnon_bin = 0

                #Tm
                tm20nt = mt.Tm_NN(grna)
                tmseed = mt.Tm_NN(seed)
                tmnon = mt.Tm_NN(nonseed)

                if 50 <= tm20nt <= 60:
                    tm20nt_bin = 1
                else:
                    tm20nt_bin = 0

                if 25 <= tmseed <= 35:
                    tmseed_bin = 1
                else:
                    tmseed_bin = 0

                if 15 <= tmnon <= 25:
                    tmnon_bin = 1
                else:
                    tmnon_bin = 0

                # MFE
                p = Popen('ViennaRNA_Package/RNAfold.exe', shell=False, stdin=PIPE, stdout=PIPE, stderr=DEVNULL)
                ans = p.communicate(grna.encode())
                res = [i for i in ans if i] # remove 'NONE' in list

                out_mfe = ''.join(map(bytes.decode, res)) #convert bytes to string.
                r_out_mfe = out_mfe.replace("\r\n", "")
                mfe = float(r_out_mfe[-6:-1])

                if -1 <= mfe <= 1:
                    mfe_bin = 1
                else:
                    mfe_bin = 0

                # position specific mono-nt
                p1 = grna[0]
                p3 = grna[2]
                p4 = grna[3]
                p5 = grna[4]
                p7 = grna[6]
                p8 = grna[7]
                p9 = grna[8]
                p10 = grna[9]
                p11 = grna[10]
                p12 = grna[11]
                p13 = grna[12]
                p14 = grna[13]
                p15 = grna[14]
                p16 = grna[15]
                p17 = grna[16]
                p18 = grna[17]
                p19 = grna[18]
                p20 = grna[19]

                p1_a = p1_t = p1_g = p1_c = 0
                p3_a = p3_t = p3_g = p3_c = 0
                p4_a = p4_t = p4_g = p4_c = 0
                p5_a = p5_t = p5_g = p5_c = 0
                p7_a = p7_t = p7_g = p7_c = 0
                p8_a = p8_t = p8_g = p8_c = 0
                p9_a = p9_t = p9_g = p9_c = 0
                p10_a = p10_t = p10_g = p10_c = 0
                p11_a = p11_t = p11_g = p11_c = 0
                p12_a = p12_t = p12_g = p12_c = 0
                p13_a = p13_t = p13_g = p13_c = 0
                p14_a = p14_t = p14_g = p14_c = 0
                p16_a = p16_t = p16_g = p16_c = 0
                p17_a = p17_t = p17_g = p17_c = 0
                p18_a = p18_t = p18_g = p18_c = 0
                p19_a = p19_t = p19_g = p19_c = 0
                p20_a = p20_t = p20_g = p20_c = 0

                for z in p1:
                    if z == 'A':
                        p1_a += 1
                    elif z == 'T':
                        p1_t += 1
                    elif z == 'G':
                        p1_g += 1
                    elif z == 'C':
                        p1_c += 1

                for z in p7:
                    if z == 'A':
                        p7_a += 1
                    elif z == 'T':
                        p7_t += 1
                    elif z == 'G':
                        p7_g += 1
                    elif z == 'C':
                        p7_c += 1

                for z in p8:
                    if z == 'A':
                        p8_a += 1
                    elif z == 'T':
                        p8_t += 1
                    elif z == 'G':
                        p8_g += 1
                    elif z == 'C':
                        p8_c += 1

                for z in p9:
                    if z == 'A':
                        p9_a += 1
                    elif z == 'T':
                        p9_t += 1
                    elif z == 'G':
                        p9_g += 1
                    elif z == 'C':
                        p9_c += 1

                for z in p11:
                    if z == 'A':
                        p11_a += 1
                    elif z == 'T':
                        p11_t += 1
                    elif z == 'G':
                        p11_g += 1
                    elif z == 'C':
                        p11_c += 1

                for z in p12:
                    if z == 'A':
                        p12_a += 1
                    elif z == 'T':
                        p12_t += 1
                    elif z == 'G':
                        p12_g += 1
                    elif z == 'C':
                        p12_c += 1

                for z in p13:
                    if z == 'A':
                        p13_a += 1
                    elif z == 'T':
                        p13_t += 1
                    elif z == 'G':
                        p13_g += 1
                    elif z == 'C':
                        p13_c += 1

                for z in p14:
                    if z == 'A':
                        p14_a += 1
                    elif z == 'T':
                        p14_t += 1
                    elif z == 'G':
                        p14_g += 1
                    elif z == 'C':
                        p14_c += 1

                for z in p16:
                    if z == 'A':
                        p16_a += 1
                    elif z == 'T':
                        p16_t += 1
                    elif z == 'G':
                        p16_g += 1
                    elif z == 'C':
                        p16_c += 1

                for z in p17:
                    if z == 'A':
                        p17_a += 1
                    elif z == 'T':
                        p17_t += 1
                    elif z == 'G':
                        p17_g += 1
                    elif z == 'C':
                        p17_c += 1

                for z in p18:
                    if z == 'A':
                        p18_a += 1
                    elif z == 'T':
                        p18_t += 1
                    elif z == 'G':
                        p18_g += 1
                    elif z == 'C':
                        p18_c += 1

                for z in p19:
                    if z == 'A':
                        p19_a += 1
                    elif z == 'T':
                        p19_t += 1
                    elif z == 'G':
                        p19_g += 1
                    elif z == 'C':
                        p19_c += 1

                for z in p20:
                    if z == 'A':
                        p20_a += 1
                    elif z == 'T':
                        p20_t += 1
                    elif z == 'G':
                        p20_g += 1
                    elif z == 'C':
                        p20_c += 1

                # position specific dint
                p7_di = grna[6] + grna[7]
                p8_di = grna[7] + grna[8]
                p11_di = grna[10] + grna[11]
                p12_di = grna[11] + grna[12]
                p13_di = grna[12] + grna[13]
                p19_di = grna[18] + grna[19]

                p7_gt = p7_ga = p8_tt = p8_ga = 0
                p11_gc = p11_ga = 0
                p12_ca = p12_ct = 0
                p13_ac = p13_at = 0
                p19_ac = p19_ga = p19_cc = 0

                if p7_di == 'GT':
                    p7_gt += 1
                elif p7_di == 'GA':
                    p7_ga += 1 #unused
                    
                if p8_di == 'TT':
                    p8_tt += 1
                elif p8_di == 'GA':
                    p8_ga += 1 #unused

                if p11_di == 'GC':
                    p11_gc += 1
                elif p11_di == 'GA':
                    p11_ga += 1 #unused

                if p12_di == 'CA':
                    p12_ca += 1
                elif p12_di == 'CT':
                    p12_ct += 1 #unused

                if p13_di == 'AC':
                    p13_ac += 1
                elif p13_di =='AT':
                    p13_at += 1 #unused

                if p19_di == 'AC':
                    p19_ac += 1
                elif p19_di == 'GA':
                    p19_ga += 1
                elif p19_di == 'CC':
                    p19_cc += 1

                # non_specific mono-nt
                count_a = grna.count('A')
                count_t = grna.count('T')
                count_g = grna.count('G')
                count_c = grna.count('C')

                # non_specific di-nt
                count_aa = grna.count('AA')
                count_at = grna.count('AT')
                count_ag = grna.count('AG')
                count_ac = grna.count('AC')
                count_ta = grna.count('TA')
                count_tt = grna.count('TT')
                count_tg = grna.count('TG')
                count_tc = grna.count('TC')
                count_ga = grna.count('GA')
                count_gt = grna.count('GT')
                count_gg = grna.count('GG')
                count_gc = grna.count('GC')
                count_ca = grna.count('CA')
                count_ct = grna.count('CT')
                count_cg = grna.count('CG')
                count_cc = grna.count('CC')

                # non_specific tri-nt
                count_aaa = grna.count('AAA')
                count_aat = grna.count('AAT')
                count_aac = grna.count('AAC')
                count_ata = grna.count('ATA')
                count_att = grna.count('ATT')
                count_atg = grna.count('ATG')
                count_atc = grna.count('ATC')
                count_aga = grna.count('AGA')
                count_agt = grna.count('AGT')
                count_agc = grna.count('AGC')
                count_aca = grna.count('ACA')
                count_act = grna.count('ACT')
                count_acg = grna.count('ACG')
                count_acc = grna.count('ACC')
                count_taa = grna.count('TAA')
                count_tac = grna.count('TAC')
                count_tta = grna.count('TTA')
                count_ttt = grna.count('TTT')
                count_ttg = grna.count('TTG')
                count_ttc = grna.count('TTC')
                count_tga = grna.count('TGA')
                count_tgt = grna.count('TGT')
                count_tgg = grna.count('TGG')
                count_tgc = grna.count('TGC')
                count_tca = grna.count('TCA')
                count_tcg = grna.count('TCG')
                count_tcc = grna.count('TCC')
                count_gaa = grna.count('GAA')
                count_gat = grna.count('GAT')
                count_gac = grna.count('GAC')
                count_gta = grna.count('GTA')
                count_gtt = grna.count('GTT')
                count_gtg = grna.count('GTG')
                count_gtc = grna.count('GTC')
                count_ggt = grna.count('GGT')
                count_ggc = grna.count('GGC')
                count_gca = grna.count('GCA')
                count_gct = grna.count('GCT')
                count_gcg = grna.count('GCG')
                count_gcc = grna.count('GCC')
                count_caa = grna.count('CAA')
                count_cat = grna.count('CAT')
                count_cag = grna.count('CAG')
                count_cac = grna.count('CAC')
                count_ctt = grna.count('CTT')
                count_ctg = grna.count('CTG')
                count_ctc = grna.count('CTC')
                count_cga = grna.count('CGA')
                count_cgt = grna.count('CGT')
                count_cgg = grna.count('CGG')
                count_cgc = grna.count('CGC')
                count_cca = grna.count('CCA')
                count_cct = grna.count('CCT')
                count_ccg = grna.count('CCG')

                # print to column result ##############################################
                # Column gc content (entire 20nt)
                view_result += "<td>"
                view_result += str(gc20nt)
                view_result += "</td>"

                # Column gc content (seed region)
                # print ("<td>")
                # print ('%.2f' % gcseed)
                # print ('</td>')

                # Column Tm (entire 20nt)
                view_result += "<td>"
                tm20nt = '%0.2f' % (tm20nt)
                view_result += str(tm20nt)

                # Column Tm (seed)
                # print('</td><td>')
                # tmseed = '%0.2f' % (tmseed)
                # print (tmseed)



                # Column on-target score Model 3 (RF + SVM)########################
                view_result += '</td><td>'

                sum_weight_model_3 = (
                    (gc20nt_bin * 0.4087) + (gcseed_bin * 1.9847) + (gcnon_bin * (-1.8354)) +
                    (tm20nt_bin * (-4.015)) + (tmseed_bin * 1.1525) + (tmnon_bin * 4.0558) +
                    (mfe * (-1.4621)) +
                    
                    (p1_c * 0.0963) + (p1_t * 0.1211) +
                    (p7_a * (-0.0188)) + (p7_t * 0.03) + (p7_c * (-0.0551)) + (p7_g  * 0.0439) + 
                    (p8_a * (-0.1861)) + (p8_t * (-0.3744)) + (p8_c * (-0.1002)) +
                    (p9_a * (-0.216)) + (p9_g * 0.2833) + (p9_t * (-0.088)) + (p9_c * 0.0207) +
                    
                    (p11_g * 0.0315) + (p12_c * (-0.0979)) +
                    (p13_a * (-0.0775)) + 
                    (p14_c * 0.0381) + (p16_g * 0.3132) +
                    (p17_c * (-0.5202)) + 
                    (p18_g * 0.4536) + (p18_t * (-0.2247)) + 
                    (p19_a * 0.2323) + (p19_g * (0.0807)) + (p19_t * (-0.263)) + (p19_c * (-0.0499)) +
                    (p20_a * (-0.2915)) + (p20_c * 0.0501) + (p20_g * 0.6645) + (p20_t * (-0.4232)) +
                    
                    (p7_gt * 1.638) + (p8_tt * 0.8972) + 
                    (p11_gc * 0.3489) +
                    (p12_ca * 0.3614) + 
                    (p13_ac * 0.9264) + 
                    (p19_ac * 0.8322) + (p19_ga * (-0.2288)) +(p19_cc * 1.1217) +
                    
                    
                    (count_a * (-0.2809)) + (count_g * (0.9933)) + (count_t * (-0.0921)) + (count_c * (-0.5505)) +
                    
                    (count_aa * (-0.8667)) + (count_ac * (-1.3666)) + (count_ag * (-1.9701)) + (count_at * (-2.1051)) + 
                    (count_ca * 0.7552) + (count_cc * 1.128) + (count_cg * (-0.4438)) + (count_ct * (-1.4744)) + 
                    (count_ga * (-1.3539)) + (count_gc * (-0.0575)) + (count_gg * 0.0525) + (count_gt * (-0.2256)) + 
                    (count_ta * (-0.7657)) + (count_tc * 0.0308) + (count_tg * (-3.6773)) + (count_tt * 0.5699) + 
                    
                    (count_aaa * (-0.0128)) + (count_aac * 0.9087) + (count_aat * 0.2856) + 
                    (count_aca * (-0.4023)) + (count_acc * (-0.2993)) + (count_acg * (-0.8448)) + (count_act * 0.5283) + 
                    (count_aga * 0.6506) + (count_agc * 0.4801) + (count_agt * 0.3879) + 
                    (count_ata * 0.7234) + (count_atc * 1.1278) + (count_atg * 0.931) + (count_att * 0.5874) + 
                    
                    (count_caa * (-0.3416)) + (count_cac * 0.2187) + (count_cag * (-0.3321)) + (count_cat * (-0.1825)) + 
                    (count_cca * 0.7264) + (count_ccg * -0.554) + (count_cct * 1.3832) + 
                    (count_cga * 0.6721) + (count_cgc * 0.5575) + (count_cgg * 1.8554) + (count_cgt * 0.2597) + 
                    (count_ctc * 0.1707) + (count_ctg * 0.1932) + (count_ctt * 0.492) +
                    
                    (count_gaa * 0.586) + (count_gac * 1.2603) + (count_gat * 0.5458) +
                    (count_gca * 0.0542) + (count_gcc * (-1.1678)) + (count_gcg * (-1.1669)) + (count_gct * 0.5598) +
                    (count_ggc * (-0.3735)) + (count_ggt * (-0.2959)) + 
                    (count_gta * (-0.4283)) + (count_gtc * 0.3377) + (count_gtg * 0.316) + (count_gtt * 0.402) + 
                    
                    (count_taa * 0.0932) + (count_tac * 0.0187) + 
                    (count_tca * (-1.2672)) + (count_tcc * (-1.8338)) + (count_tcg * (-1.7346)) + 
                    (count_tga * 1.1554) + (count_tgc * 1.0066) + (count_tgg * 0.9546) + (count_tgt * 1.2865) + 
                    (count_tta * (-0.3283)) + (count_ttc * (-0.3375)) + (count_ttg * (-0.9357)) + (count_ttt * (-0.6319))
                )

                # gsj = intercept + sum_weight
                gsj_model_3 = 1.6938508 + sum_weight_model_3

                on_score_model_3 = 1 / (1 + math.exp(gsj_model_3))

                view_result += str("%.5f" % (on_score_model_3))
                
                
                
                # Column on-target score Model 4 (RF + LR) ########################
                view_result += '</td><td>'

                sum_weight_model_4 = (
                    (gc20nt_bin * 0.3472) + (gcseed_bin * 1.6664) + (gcnon_bin * (-1.502)) +
                    (tm20nt_bin * (-1.3501)) + (tmseed_bin * 0.5889) + (tmnon_bin * 3.1337) +
                    (mfe * (-1.4908)) +
                    
                    (p1_c * 0.1169) + (p1_t * 0.0616) +
                    (p7_a * 0.0169) + (p7_t * 0.0026) + (p7_c * 0.005) + (p7_g  * 0.0499) + 
                    (p8_a * (-0.1351)) + (p8_t * (-0.4211)) + (p8_c * (-0.0852)) +
                    (p9_a * -0.2177) + (p9_g * 0.3036) + (p9_t * (-0.0307)) + (p9_c * 0.0192) +
                    
                    (p11_g * 0.0335) + (p12_c * (-0.1369)) +
                    (p13_a * (-0.0825)) + 
                    (p14_c * 0.1137) + (p16_g * 0.3525) +
                    (p17_c * (-0.5522)) + 
                    (p18_g * 0.5108) + (p18_t * (-0.2983)) + 
                    (p19_a * 0.3024) + (p19_g * 0.0877) + (p19_t * (-0.2558)) + (p19_c * (-0.0598)) +
                    (p20_a * (-0.2515)) + (p20_c * 0.2015) + (p20_g * 0.651) + (p20_t * (-0.5266)) +
                    
                    (p7_gt * 1.9919) + (p8_tt * 1.2282) + 
                    (p11_gc * 0.4386) +
                    (p12_ca * 0.4939) + 
                    (p13_ac * 0.8174) + 
                    (p19_ac * 0.4998) + (p19_ga * (-0.3521)) +(p19_cc * 0.9542) +
                    
                    
                    (count_a * (-0.4992)) + (count_g * 0.7499) + (count_t * 0.2444) + (count_c * (-0.3489)) +
                    
                    (count_aa * (-0.3071)) + (count_ac * (-0.9096)) + (count_ag * (-1.4796)) + (count_at * (-1.0431)) + 
                    (count_ca * 0.7176) + (count_cc * 1.0943) + (count_cg * (-0.5988)) + (count_ct * (-1.0479)) + 
                    (count_ga * (-0.5958)) + (count_gc * (-0.2657)) + (count_gg * 0.0761) + (count_gt * 0.0455) + 
                    (count_ta * (-0.338)) + (count_tc * 0.1305) + (count_tg * (-2.1386)) + (count_tt * 1.0485) + 
                    
                    (count_aaa * 0.3709) + (count_aac * 0.9393) + (count_aat * 0.1566) + 
                    (count_aca * (-0.5274)) + (count_acc * (-0.291)) + (count_acg * (-0.849)) + (count_act * 0.4239) + 
                    (count_aga * 0.3422) + (count_agc * 0.4058) + (count_agt * 0.3406) + 
                    (count_ata * 0.605) + (count_atc * 0.9266) + (count_atg * 0.4657) + (count_att * 0.3162) + 
                    
                    (count_caa * (-0.3753)) + (count_cac * (0.2859)) + (count_cag * (-0.2684)) + (count_cat * (-0.2828)) + 
                    (count_cca * 0.854) + (count_ccg * (-0.3878)) + (count_cct * 1.4651) + 
                    (count_cga * 0.3747) + (count_cgc * 0.2135) + (count_cgg * 1.8149) + (count_cgt * 0.2219) + 
                    (count_ctc * (-0.0482)) + (count_ctg * (-0.083)) + (count_ctt * 0.3343) +
                    
                    (count_gaa * 0.3526) + (count_gac * 1.0839) + (count_gat * 0.4857) +
                    (count_gca * 0.4287) + (count_gcc * (-0.9941)) + (count_gcg * (-0.9354)) + (count_gct * 0.7785) +
                    (count_ggc * (-0.4888)) + (count_ggt * 0.0396) + 
                    (count_gta * (-0.5258)) + (count_gtc * 0.1349) + (count_gtg * (-0.0969)) + (count_gtt * 0.3619) + 
                    
                    (count_taa * (-0.1051)) + (count_tac * 0.125) + 
                    (count_tca * (-1.1465)) + (count_tcc * (-1.7627)) + (count_tcg * (-1.4781)) + 
                    (count_tga * 0.3007) + (count_tgc * 0.3551) + (count_tgg * 0.4658) + (count_tgt * 0.4415) + 
                    (count_tta * (-0.0813)) + (count_ttc * (-0.2823)) + (count_ttg * (-1.2895)) + (count_ttt * (-0.3301))
                )

                # gsj = intercept + sum_weight
                gsj_model_4 = 0.356733156 + sum_weight_model_4

                on_score_model_4 = 1 / (1 + math.exp(gsj_model_4))

                view_result += str("%.5f" % (on_score_model_4))
                
                

                # Column mismatch_no #############################
                view_result += "</td><td>"

                for y in genome_seq:
                    # offtarget = regex.findall(r"(" + grna+"){1<=s<=5}([atgcATGC][AG][G])", y)   # y = genome_seq

                    offtarget = regex.findall(r"(" + grna3+"){1<=s<=5}([G][G])", y)
                    
                    no_offtarget = len(offtarget)
                    view_result += str(no_offtarget)

                    view_result += '''<input type="button" value="+" class="more">
                        </td></tr>


                        <tr>
                            <td id="sub-td"></td>
                            <td colspan="6" style="text-align:left; font-weight:bold;" id="sub-td"><br>
                            <span style="color:#3A506B;">OFF-TARGET PREDICTION:</span>
                            </td>
                        </tr>

                        </tr>
                        <tr>
                            <td id="sub-td"></td>
                            <td colspan="2" id="sub-td" style="font-weight:bold; color:#303030;">Off-target</td>
                            <td id="sub-td" style="font-weight:bold; color:#303030;">Mismatches</td>
                            <td id="sub-td" style="font-weight:bold; color:#303030;">Start</td>
                            <td id="sub-td" style="font-weight:bold; color:#303030;">End</td>
                            <td id="sub-td" style="font-weight:bold; color:#303030;">CFD Score</td>
                        </tr>

                    '''

                    for x in offtarget:
                        offt = ''.join(x)
                        offt_loc = re.search(offt, y)

                        ngg_offt = offt[-3:]
                        offt2 = offt[:20]

                        match = []
                        m = []
                        mismatch = 0


                        # grna2 = grna + ngg + pamseq
                        ngg_grna = ngg + pamseq


                        for a, b in zip(grna, offt2):
                            if a == b:
                                match.append(b)
                            else:
                                match.append('<span style="color:red;">')
                                match.append(b)
                                match.append('</span>')
                                mismatch += 1

                        for a, b in zip(ngg_grna, ngg_offt):
                            if a == b:
                                m.append(b)
                            else:
                                m.append('<span style="color:red;">')
                                m.append(b)
                                m.append('</span>')
                                mismatch += 1

                        matched = "".join(match)
                        matched_ngg = "".join(m)

                        view_result += '''
                        <tr>
                            <td id="sub-td"></td>
                            <td colspan="2" id="sub-td" class="seq">
                            '''

                        view_result += str("%s <mark>%s</mark>" % (matched, matched_ngg))

                        view_result += '''</td>
                            <td id="sub-td">'''
                        view_result += str("%g" % (mismatch))

                        view_result += '''
                            </td>
                            <td id="sub-td">'''

                        view_result += str(offt_loc.start())

                        view_result += '''</td>
                            <td id="sub-td">'''

                        view_result += str(offt_loc.end())

                        view_result += '''</td>
                        <td id="sub-td">'''

                        mm_scores,pam_scores = get_mm_pam_scores()
                        m_grna = re.search('[^ATCG]',grna2)
                        m_off = re.search('[^ATCG]',offt)
                        if (m_grna is None) and (m_off is None):
                            pam_offt = offt[-2:]
                            cfd_score = calc_cfd(grna2,offt2,pam_offt)
                            view_result += str("%.6f" % (cfd_score))

                        view_result += '''
                        </td>
                        </tr>
                        '''
                    view_result +='''<tr><td colspan="7" id="sub-td"></td></tr>'''

                pam = target_seq.find("GG", pam+1)

            else:
                pam = target_seq.find("GG", pam+1)

        view_result += '''</table>'''
        
        
        
        
        
        # ###########################################################
        
        # return(view_result)

    return render_template('result.html', result=view_result)
    # return


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000,debug=True)
