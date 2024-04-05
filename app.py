from flask import Flask,  request, url_for, redirect, render_template, json
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
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from Bio.Blast import NCBIXML
from fuzzysearch import find_near_matches

app = Flask(__name__)

target_genomes = {
    'GCA_000006945.2': 'Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 ',
    'GCF_000018445.1': 'Acinobacter baumannii str. ACICU ',
    'GCF_000401555.1': 'Aeromonas hydrophila str. ML09-119 ',
    'GCF_000014805.1': 'Aeromonas hydrophila subsp. hydrophila str. ATCC 7966 ',
    'GCF_000196395.1': 'Aeromonas salmonicida subsp. salmonicida str. A449 ',
    'GCF_000204115.1': 'Aeromonas veronii str. B565 ',
    'GCF_000024505.1': 'Anaplasma centrale str. Israel ',
    'GCF_000020305.1': 'Anaplasma marginale str. Florida ',
    'GCF_000013125.1': 'Anaplasma phagocytophilum str. HZ ',
    'GCA_000007845.1': 'Bacillus anthracis str. Ames ',
    'GCA_000015445.1': 'Bacillus bacilliformis str. KC583 ',
    'GCA_000007825.1': 'Bacillus cereus str. ATCC 14579 ',
    'GCF_000009825.1': 'Bacillus clausii str. KSM-K16 ',
    'GCF_000011145.1': 'Bacillus halodurans str. C-125 ',
    'GCF_000046705.1': 'Bacillus henselae str. Houston-1 ',
    'GCF_000046685.1': 'Bacillus quintana str. Toulouse ',
    'GCF_000009045.1': 'Bacillus subtilis subsp. subtilis str. 168 ',
    'GCF_000008505.1': 'Bacillus thuringiensis serovar konkukian str. 97-27 ',
    'GCF_000196435.1': 'Bacillus tribocorum str. CIP 105476 ',
    'GCF_000018825.1': 'Bacillus weihenstephanensis str. KBAB4 ',
    'GCF_000070465.1': 'Bordetella avium str. 197N ',
    'GCF_000195675.1': 'Bordetella bronchiseptica str. RB50 ',
    'GCF_000195695.1': 'Bordetella parapertussis str. 12822 ',
    'GCF_000195715.1': 'Bordetella pertussis str. Tohama I ',
    'GCA_000008145.1': 'Brucella abortus bv. 1 str. 9-941 chromosome I ',
    'GCA_000018525.1': 'Brucella canis str. ATCC 23365 chromosome I ',
    'GCF_000054005.1': 'Brucella melitensis biovar Abortus str. 2308 chromosome I ',
    'GCF_000007125.1': 'Brucella melitensis bv. 1 str. 16M chromosome I ',
    'GCF_000022625.1': 'Brucella melitensis str. ATCC 23457 chromosome I ',
    'GCA_000016845.1': 'Brucella ovis str. ATCC 25840 chromosome I ',
    'GCA_000007505.1': 'Brucella suis str. 1330 chromosome I ',
    'GCF_000009485.1': 'Burkholderia cenocepacia str. J2315 chromosome I ',
    'GCA_000011705.1': 'Burkholderia mallei str. ATCC 23344 chromosome I ',
    'GCA_000011545.1': 'Burkholderia pseudomallei str. K96243 chromosome I ',
    'GCA_000012365.1': 'Burkholderia thailandensis str. E264 chromosome I ',
    'GCF_000015085.1': 'Campylobacter fetus subsp. fetus str. 82-40 ',
    'GCF_000011865.1': 'Campylobacter jejuni str. RM1221 ',
    'GCF_000017905.1': 'Campylobacter jejuni subsp. jejuni str. 81116 ',
    'GCF_000026025.1': 'Chlamydia abortus str. S26/3 ',
    'GCA_000007605.1': 'Chlamydia caviae str. GPIC ',
    'GCA_000009945.1': 'Chlamydia felis str. Fe/C-56 ',
    'GCA_000006685.1': 'Chlamydia muridarum str. Nigg ',
    'GCF_000011165.1': 'Chlamydia pneumoniae str. J138 ',
    'GCA_000068585.1': 'Chlamydia trachomatis str. 434/Bu ',
    'GCA_000012125.1': 'Chlamydia trachomatis str. A/HAR-13 ',
    'GCA_000008725.1': 'Chlamydia trachomatis str. D/UW-3/CX ',
    'GCF_000068525.2': 'Chlamydia trachomatis str. L2b/UCH-1/proctitis ',
    'GCF_000008765.1': 'Clostridium acetobutylicum str. ATCC 824 ',
    'GCF_000016965.1': 'Clostridium beijerinckii str. NCIMB 8052 ',
    'GCF_000063585.1': 'Clostridium botulinum A str. ATCC 3502 ',
    'GCF_000022765.1': 'Clostridium botulinum A2 str. Kyoto ',
    'GCA_000019545.1': 'Clostridium botulinum A3 str. Loch Maree ',
    'GCF_000019305.1': 'Clostridium botulinum B1 str. Okra ',
    'GCF_000020285.1': 'Clostridium botulinum E3 str. Alaska E43 ',
    'GCF_000017065.1': 'Clostridium botulinum F str. Langeland ',
    'GCF_000009205.2': 'Clostridium difficile str. 630 ',
    'GCF_000014125.1': 'Clostridium novyi str. NT ',
    'GCF_000013285.1': 'Clostridium perfringens str. ATCC 13124 ',
    'GCF_000007625.1': 'Clostridium tetani str. E88 ',
    'GCF_000015865.1': 'Clostridium thermocellum str. ATCC 27405 ',
    'GCF_000195815.1': 'Corynebacterium diphtheriae str. NCTC 13129 ',
    'GCF_000011305.1': 'Corynebacterium efficiens str. YS-314 ',
    'GCF_000011325.1': 'Corynebacterium glutamicum str. ATCC 13032 ',
    'GCF_000006605.1': 'Corynebacterium jeikeium str. K411 ',
    'GCF_000143705.2': 'Corynebacterium pseudotuberculosis str. FRC41 ',
    'GCF_000007765.2': 'Coxiella burnetii str. RSA 493 ',
    'GCF_000281195.1': 'Enterococcus faecalis str. D32 ',
    'GCF_000007785.1': 'Enterococcus faecalis str. V583 ',
    'GCF_000174395.2': 'Enterococcus faecium str. DO ',
    'GCF_000299255.1': 'Escherichia coli O104:H4 str. 2009EL-2050 ',
    'GCF_000008865.2': 'Escherichia coli O157:H7 str. Sakai ',
    'GCF_000285655.3': 'Escherichia coli O25b:H4-ST131 strain:EC958 ',
    'GCF_000091005.1': 'Escherichia coli O26:H11 str. 11368 ',
    'GCF_000027125.1': 'Escherichia coli O44:H18 str. 042 ',
    'GCF_000245515.1': 'Escherichia coli O55:H7 str. RM12579 ',
    'GCF_000227625.1': 'Escherichia coli O7:K1 str. CE10 ',
    'GCF_000183345.1': 'Escherichia coli O83:H1 str. NRG 857C ',
    'GCF_000013305.1': 'Escherichia coli str. 536 ',
    'GCF_000014845.1': 'Escherichia coli str. APEC O1 ',
    'GCF_000005845.2': 'Escherichia coli str. K-12 substr. MG1655 ',
    'GCF_000026305.1': 'Escherichia coli str. SMS-3-5 ',
    'GCF_000212715.2': 'Escherichia coli str. UMNK88 ',
    'GCF_000262205.1': 'Francisella noatunensis subsp. orientalis str. Toba 04 ',
    'GCF_000168775.2': 'Francisella tularensis subsp. holarctica str. FSC200 ',
    'GCF_000018925.1': 'Francisella tularensis subsp. mediasiatica str. FSC147 ',
    'GCF_000008985.1': 'Francisella tularensis subsp. tularensis str. SCHU S4 ',
    'GCF_000007945.1': 'Haemophilus ducreyi 35000HP ',
    'GCF_000027305.1': 'Haemophilus influenzae Rd KW20 ',
    'GCF_000009305.1': 'Helicobacter acinonychis str. Sheeba ',
    'GCF_000007905.1': 'Helicobacter hepaticus ATCC 51449 ',
    'GCF_000091985.1': 'Helicobacter mustelae 12198 ',
    'GCF_000008525.1': 'Helicobacter pylori 26695 ',
    'GCF_000009885.1': 'Klebsiella pneumoniae subsp. pneumoniae NTUH-K2044 ',
    'GCF_000019565.1': 'Klebsiella variicola strain 342 ',
    'GCF_000091785.1': 'Legionella longbeachae NSW150 ',
    'GCF_000008485.1': 'Legionella pneumophila subsp. pneumophila str. Philadelphia 1 ',
    'GCF_000195795.1': 'Listeria innocua Clip11262 ',
    'GCF_000252975.1': 'Listeria ivanovii subsp. ivanovii PAM 55 ',
    'GCF_000196035.1': 'Listeria monocytogenes EGD-e ',
    'GCF_000210795.2': 'Listeria monocytogenes serotype 7 str. SLCC2482 ',
    'GCF_000008285.1': 'Listeria monocytogenes str. 4b F2365 ',
    'GCF_000210815.2': 'Listeria monocytogenes strain SLCC2372, serotype 1/2c ',
    'GCF_000197755.2': 'Listeria monocytogenes strain SLCC2755, serotype 1/2b ',
    'GCF_000027145.1': 'Listeria seeligeri serovar 1/2b str. SLCC3954 ',
    'GCF_000060285.1': 'Listeria welshimeri serovar 6b str. SLCC5334 ',
    'GCF_000277775.2': 'Mycobacterium abscessus subsp. massiliense str. GO 06 ',
    'GCF_000253355.1': 'Mycobacterium africanum GM041182 ',
    'GCF_000007865.1': 'Mycobacterium avium subsp. paratuberculosis str. k10 ',
    'GCF_000234725.1': 'Mycobacterium bovis BCG str. Mexico ',
    'GCF_000328805.1': 'Mycobacterium canettii CIPT 140060008 ',
    'GCF_000184435.1': 'Mycobacterium gilvum Spyr1 ',
    'GCF_000298095.1': 'Mycobacterium indicus pranii MTCC 9506 ',
    'GCF_000195955.2': 'Mycobacterium tuberculosis H37Rv ',
    'GCF_000013925.1': 'Mycobacterium ulcerans Agy99 ',
    'GCF_000063605.1': 'Mycoplasma agalactiae PG2 ',
    'GCF_000020065.1': 'Mycoplasma arthritidis 158L3-1 ',
    'GCF_000012765.1': 'Mycoplasma capricolum subsp. capricolum ATCC 27343 ',
    'GCF_000026765.1': 'Mycoplasma conjunctivae HRC/581T ',
    'GCF_000148625.1': 'Mycoplasma fermentans JER ',
    'GCF_000092585.1': 'Mycoplasma gallisepticum str. R(low)',
    'GCF_000027325.1': 'Mycoplasma genitalium G37 ',
    'GCF_000085865.1': 'Mycoplasma hominis ATCC 23114 ',
    'GCF_000008405.1': 'Mycoplasma hyopneumoniae 232 ',
    'GCF_000008365.1': 'Mycoplasma mobile 163K ',
    'GCF_000011445.1': 'Mycoplasma mycoides subsp. mycoides SC str. PG1 ',
    'GCF_000011225.1': 'Mycoplasma penetrans HF-2 DNA ',
    'GCF_000027345.1': 'Mycoplasma pneumoniae M129 ',
    'GCF_000195875.1': 'Mycoplasma pulmonis UAB CTIP ',
    'GCF_000008245.1': 'Mycoplasma synoviae 53 ',
    'GCF_000020105.1': 'Neisseria gonorrhoeae NCCP11945 ',
    'GCF_000196295.1': 'Neisseria lactamica 020-06 ',
    'GCF_000008805.1': 'Neisseria meningitidis MC58 ',
    'GCF_000006765.1': 'Pseudomonas aeruginosa PAO1 ',
    'GCF_000026105.1': 'Pseudomonas entomophila str. L48 ',
    'GCF_000009225.2': 'Pseudomonas fluorescens SBW25 ',
    'GCF_000007565.2': 'Pseudomonas putida KT2440 ',
    'GCF_000012205.1': 'Pseudomonas savastanoi pv. phaseolicola 1448A ',
    'GCF_000013785.1': 'Pseudomonas stutzeri A1501 ',
    'GCF_000023005.1': 'Rickettsia africae ESF-5 ',
    'GCF_000018205.1': 'Rickettsia akari str. Hartford ',
    'GCF_000018245.1': 'Rickettsia bellii OSU 85-389 ',
    'GCF_000014345.1': 'Rickettsia canadensis str. McKiel ',
    'GCF_000007025.1': 'Rickettsia conorii str. Malish 7 ',
    'GCF_000016625.1': 'Rickettsia massiliae MTU5 ',
    'GCF_000021525.1': 'Rickettsia peacockii str. Rustic ',
    'GCF_000195735.1': 'Rickettsia prowazekii str. Madrid E ',
    'GCF_000017445.4': 'Rickettsia rickettsii str. Iowa ',
    'GCF_000008045.1': 'Rickettsia typhi str. Wilmington ',
    'GCF_000018625.1': 'Salmonella enterica subsp. arizonae serovar 62:z4,z23:-- ',
    'GCF_000020885.1': 'Salmonella enterica subsp. enterica serovar Agona str. SL483 ',
    'GCF_000008105.1': 'Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67 ',
    'GCF_000020925.1': 'Salmonella enterica subsp. enterica serovar Dublin str. CT_02021853 ',
    'GCF_000009505.1': 'Salmonella enterica subsp. enterica serovar Enteritidis str. P125109 ',
    'GCF_000009525.1': 'Salmonella enterica subsp. enterica serovar Gallinarum str. 287/91 ',
    'GCF_000020705.1': 'Salmonella enterica subsp. enterica serovar Heidelberg str. SL476 ',
    'GCF_000016045.1': 'Salmonella enterica subsp. enterica serovar Newport str. SL254 ',
    'GCF_000026565.1': 'Salmonella enterica subsp. enterica serovar Paratyphi A str. AKU_12601 ',
    'GCF_000011885.1': 'Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150 ',
    'GCF_000018705.1': 'Salmonella enterica subsp. enterica serovar Paratyphi B str. SPB7 ',
    'GCF_000018385.1': 'Salmonella enterica subsp. enterica serovar Paratyphi C strain RKS4594 ',
    'GCF_000020745.1': 'Salmonella enterica subsp. enterica serovar Schwarzengrund str. CVM19633 ',
    'GCF_000195995.1': 'Salmonella enterica subsp. enterica serovar Typhi str. CT18 ',
    'GCF_000007545.1': 'Salmonella enterica subsp. enterica serovar Typhi Ty2 ',
    'GCF_000006945.2': 'Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 ',
    'GCF_000012025.1': 'Shigella boydii Sb227 ',
    'GCF_000012005.1': 'Shigella dysenteriae Sd197 ',
    'GCF_000007405.1': 'Shigella flexneri 2a str. 2457T ',
    'GCF_000013585.1': 'Shigella flexneri 5 str. 8401 ',
    'GCF_000092525.1': 'Shigella sonnei Ss046 ',
    'GCF_000009005.1': 'Staphylococcus aureus RF122 ',
    'GCF_000010445.1': 'Staphylococcus aureus subsp. aureus Mu3 ',
    'GCF_000009665.1': 'Staphylococcus aureus subsp. aureus Mu50 ',
    'GCF_000011925.1': 'Staphylococcus epidermidis RP62A ',
    'GCF_000009865.1': 'Staphylococcus haemolyticus JCSC1435 ',
    'GCF_000010125.1': 'Staphylococcus saprophyticus subsp. saprophyticus ATCC 15305 ',
    'GCF_000012705.1': 'Streptococcus agalactiae A909 ',
    'GCF_000010705.1': 'Streptococcus dysgalactiae subsp. equisimilis GGS_124 ',
    'GCF_000026585.1': 'Streptococcus equi subsp. equi 4047 ',
    'GCF_000020765.1': 'Streptococcus equi subsp. zooepidemicus MGCS10565 ',
    'GCF_000017005.1': 'Streptococcus gordonii str. Challis substr. CH1 ',
    'GCF_000091645.1': 'Streptococcus mutans NN2025 ',
    'GCF_000007045.1': 'Streptococcus pneumoniae R6 ',
    'GCF_000011765.3': 'Streptococcus pyogenes MGAS5005 ',
    'GCF_000014205.1': 'Streptococcus sanguinis SK36 ',
    'GCF_000014305.1': 'Streptococcus suis 05ZYH33 ',
    'GCF_000014485.1': 'Streptococcus thermophilus LMD-9 ',
    'GCF_000006745.1': 'Vibrio cholerae O1 biovar eltor str. N16961 chromosome I ',
    'GCF_000011805.1': 'Vibrio fischeri ES114 chromosome I ',
    'GCF_000196095.1': 'Vibrio parahaemolyticus RIMD 2210633 ',
    'GCF_000039765.1': 'Vibrio vulnificus CMCP6 chromosome I ',
    'GCF_000009345.1': 'Yersinia enterocolitica subsp. enterocolitica 8081 ',
    'GCF_000192105.1': 'Yersinia enterocolitica subsp. palearctica 105.5R®',
    'GCF_000007885.1': 'Yersinia pestis biovar Microtus str. 91001 ',
    'GCF_000016945.1': 'Yersinia pseudotuberculosis IP 31758 ',
}

# remove cache
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
app.config["CACHE_TYPE"] = "null"

# DEFINE route to main.html
@app.route('/')
def index():
    return render_template('main.html', target_genomes=target_genomes)


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
        
        print("Receiving requests...")
        job_name = request.form.get("jobname")
        target_genome = request.form.get("targetgenome")
        target_sequence = request.form.get("sequence")
        # target_seq = target_sequence.replace('\n', '').replace(' ', '').replace('\r', '')
        pamseq = request.form.get('pamsequence')
        
        # # find pam sequence in target sequence
        # pam = target_seq.find(pamseq)
        
        print("Finding candidate gRNAs...")
        # find candidate gRNAs in target sequence (include PAM site)
        candidate_gRNAs = find_candidate_gRNAs(target_sequence, pamseq)         
        print(len(candidate_gRNAs))
        
        # retrieve genome sequence from database
        print(f"Retrieving target genome ({target_genome}) sequence...")
        target_genome_file = os.path.join('target_genome_file', target_genome + ".fna")
        genome_seq = load_target_genome(target_genome_file)
    
        
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

        # while pam != -1:
        for grna_pam in candidate_gRNAs:
            
            print(f"#########Running analysis for {grna_pam}...#########")
            
            # grna_pam 
            ngg = grna_pam[-3]
            grna = grna_pam[:20]
            grna2 = grna_pam
            grna3 = grna_pam[:21]

            length = len(grna)
            
            view_result += '''<tr class="expandable">'''

            # Column candidate gRNA
            view_result += "<td>"
            view_result += "%.20s <mark>%s%s</mark>" % (grna, ngg, pamseq)
            view_result += "</td>"

            # Column position
            view_result += "<td> "
            
            print(f"Finding start_pos and end_pos for grna: {grna}")
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

            print(f"Calculate GC content of grna: {grna}...")
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

            print(f"Calculate Tm of grna: {grna}...")
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

            print(f"Calculate MFE of grna: {grna}...")
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

            print(f"Analysing position specific mono-nt of grna: {grna}...")
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

            print(f"Analysing position specific di-nt of grna: {grna}...")
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

            print(f"Calculating on-target score of grna: {grna}...")
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

            print(f"Finding off-target of grna: {grna}...")
            # off_targets = find_off_targets(grna_pam, genome_seq, allowed_mismatches=3)
            
            find_offtargets = regex.findall(r"(" + grna3+"){1<=s<=5}([G][G])", genome_seq)
            offtargets = [sequence + pam for sequence, pam in find_offtargets]
            no_offtarget = len(offtargets)
            
            print(offtargets)
            print(no_offtarget)
            
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
            
            # for y in genome_seq:
            for offtarget in offtargets:
            
                # offtarget = regex.findall(r"(" + grna+"){1<=s<=5}([atgcATGC][AG][G])", y)   # y = genome_seq

                # offtarget = regex.findall(r"(" + grna3+"){1<=s<=5}([G][G])", y)
                # offtarget = regex.findall(r"(" + grna3+"){1<=s<=5}([G][G])", genome_seq)
                # no_offtarget = len(offtarget)
                
                # offtarget_sequence = y[1]
                # no_offtarget = y[2]
                # highlighted_offtarget_sequence = y[3]
                # offtarget_position = y[0]

                # for x in offtarget:
                     
                # offt = ''.join(x)
                # offt = x[0] + x[1]
                print(f"Running analysis for off-target: {offtarget}...")
                offt_loc = re.search(offtarget, genome_seq)

                ngg_offt = offtarget[-3:]
                offt2 = offtarget[:20]

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

                print(f"Calculate off-target score for: {offtarget}...")
                
                mm_scores,pam_scores = get_mm_pam_scores()
                m_grna = re.search('[^ATCG]',grna2)
                m_off = re.search('[^ATCG]',offtarget)
                if (m_grna is None) and (m_off is None):
                    pam_offt = offtarget[-2:]
                    cfd_score = calc_cfd(grna2,offt2,pam_offt)
                    view_result += str("%.6f" % (cfd_score))

                view_result += '''
                </td>
                </tr>
                '''
                view_result +='''<tr><td colspan="7" id="sub-td"></td></tr>'''


            
        view_result += '''</table>'''
        
        
        
        
        
        # ###########################################################
        
        # return(view_result)

    return render_template('result.html', result=view_result)
    # return


def load_target_genome(target_genome_file):
    for record in SeqIO.parse(target_genome_file, "fasta"):
        return str(record.seq)

def find_candidate_gRNAs(target_sequence, pam_sequence):
    candidate_grnas = []
    candidate_grnas_pam = []

    target_seq_file = StringIO(target_sequence)

    # Parse the target sequence from the provided string
    for record in SeqIO.parse(target_seq_file, "fasta"):
        seq = str(record.seq)

        for i in range(len(seq) - 23):  # Assuming gRNA length is 20 (not including PAM)
            grna_sequence = seq[i:i+20]  # Extract potential gRNA sequence
            pam_candidate = seq[i+21:i+23]  # Extract candidate PAM sequence
            
            grna_and_pam = f"{grna_sequence}{pam_candidate}"

            # Check if the candidate sequence is a valid gRNA (ending with NGG)
            if pam_candidate == pam_sequence:
                candidate_grnas.append(grna_sequence)
                candidate_grnas_pam.append(grna_and_pam)
            
    return candidate_grnas_pam

def find_off_targets(grna_pam, genome_seq, allowed_mismatches):
    
    off_targets=[]
    grna_length = len(grna_pam)
    
    for i in range(len(genome_seq) - grna_length):
        window_sequence = genome_seq[i:i+grna_length]
        near_matches = find_near_matches(grna_pam, window_sequence, max_l_dist=allowed_mismatches)
        
        for offtarget in near_matches:
            print(f"Analysing off-target: {offtarget}...")
            
            match_position = i + offtarget.start
            match_sequence = genome_seq[match_position:match_position+grna_length]
            num_mismatches = offtarget.dist
            
            if match_sequence[-2:] == 'GG':
                print(f"Highlighting mismatches: {match_sequence}...")
                
                highlighted_sequence = ''
                for j, char in enumerate(match_sequence):
                    if grna_pam[j] != char:
                        highlighted_sequence += f'<span style="color: red;">{char}</span>'
                    else:
                        highlighted_sequence += char
            
            
                off_targets.append((match_position, match_sequence, num_mismatches, highlighted_sequence))
    
    return off_targets


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000,debug=True)
