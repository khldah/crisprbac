#!C:\Python27\python.exe
#from Bio.Seq import Seq
from string import *
import cgi, cgitb

form = cgi.FieldStorage()

target_seq = form.getvalue('sequence')
pam=find(target_seq,"GG")

print "Content-type:text/html\r\n\r\n"

print '''
<!doctype html>
<html>
    <head>
        <title>CRISPR-Bac</title>
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
        <meta name="keywords" content="CRISPR, gRNA design, guide sequence, CRISPR Cas9, genome editing, genome engineering, off-target" />

        <link href="https://fonts.googleapis.com/css?family=Righteous" rel="stylesheet">
        <link href="https://fonts.googleapis.com/css?family=Roboto:300,400,500,700&display=swap" rel="stylesheet">
        <link href="https://fonts.googleapis.com/css?family=Playfair+Display:600&display=swap" rel="stylesheet">
        <link href="https://fonts.googleapis.com/css?family=Quicksand:700&display=swap" rel="stylesheet">
        <link href="./css/style.css" rel="stylesheet" type="text/css" media="all" />

        <div class="navbar">
          <ul>
            <li><a href="contact.html">Contact</a></li>
            <li><a href="manual.html">Help</a></li>
            <li><a href="about.html">About</a></li>
            <li><a href="index.html">Design</a></li>
          </ul>
        </div>
    </head>
    <body>
        <div class="main">
            <h1>CRISPR<span style="color:#E83F6F;">BAC</span></h1>
            <h2>CRISPR guide RNA design for bacteria</h2>
                <table>
                    <tr>
                        <th>gRNA sequence</th>
                        <th>On-target score</th>
                        <th>Number of mismatch(es)</th>
                    </tr>
'''

while pam != -1:
    seq_npam = target_seq[:pam]
    grna = seq_npam[-21:-1]


    length=len(grna)
    if length == 20:
        print "<tr>"
        print "<td class='seq'>"
        print "%.20s" % (grna)
        print "</td>"
        print "<td>khalidah</td>"
        print "<td>2 mismatches</td>"
        print "</tr>"

        pam=find(target_seq,"GG", pam+1)

    else:
        pam=find(target_seq,"GG", pam+1)

print '''
    </table>
    </div>

</body>
</html>
'''
