# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:06:00 2020

@author: Jixiao He
"""

import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline

uniprot = ['P24941', 'O14733', 'P54646', 'Q07912', 'P28482', 'P35557', 'P43403', 'O14965',
           'P49674', 'P06493', 'Q8IVH8', 'P14618', 'O43781', 'P31749', 'P31751', 'Q9NWZ3',
           'Q13882', 'P27707', 'P51617', 'Q13131', 'O96017', 'P22392', 'Q00535', 'P00519',
           'Q96KB5', 'Q14680', 'Q13153', 'P41743', 'P27361', 'O96013', 'Q9H8S9', 'O43318',
           'O60825', 'P15531', 'Q9UQB9', 'P78356', 'P51955', 'P60891', 'Q99558', 'Q13976',
           'P00558', 'Q16566', 'Q08881', 'P49841', 'P43405', 'P53350', 'Q9BUB5', 'Q04759',
           'Q9NQU5', 'Q06187', 'Q15118', 'Q05397', 'Q16539', 'Q16875', 'P07332', 'Q7KZI7',
           'P51812', 'Q9UK32', 'P33981', 'O43353', 'Q15759', 'O75460', 'P36888', 'Q16512',
           'Q92918', 'P07949', 'Q96GD4', 'Q9NYV4', 'Q13546', 'O95747', 'P41279', 'O43683',
           'O15530', 'P35790', 'P52333', 'O00141', 'Q15208', 'O75385', 'Q9NYL2', 'P50750',
           'O43293', 'O14920', 'P00533']          
for i in uniprot:
    #BLAST
    fasta_string = open("D:/subject/Python script/evol/fasta/" + i + ".fasta").read()
    result_handle_P = NCBIWWW.qblast("blastp", "swissprot", fasta_string, descriptions = 1000, alignments = 1000,
                                     hitlist_size = 1000)
    
    #output blast results
    save_file = open("D:/subject/Python script/evol/blast/swissprot/"+ i + "_blast.xml", "w")  #single BLAST result 
    save_file.write(result_handle_P.read())
    save_file.close()
    result_handle_P.close()
    
    result_handle = open("D:/subject/Python script/evol/blast/swissprot/"+ i + "_blast.xml")
    blast_record = NCBIXML.read(result_handle)
    
    #write aligment resutls
    for alignment in blast_record.alignments: 
        # Each alignment is a Blast.Record.Alignment
        for hsp in alignment.hsps:
            txt_file = open("D:/subject/Python script/evol/blast_aln/swissprot/" + i + "_blast.txt", "a", encoding="utf-8")  #open file
            txt_file.write('>' + alignment.title)
            txt_file.write("\n")
            txt_file.write(hsp.sbjct)
            txt_file.write("\n")
            txt_file.close()
            
    os.chdir("D:/subject/Python script/evol/clustal_omega/swissprot")
    
    #Omega
    in_file = i + "_blast.txt"
    out_file = i + "_omega_aln.fasta"
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
    txt_file = open("omega.txt", "a", encoding="utf-8")
    txt_file.write(str(clustalomega_cline))
    txt_file.write("\n")
    txt_file.close()

        
        
        
        
        
        
        
        