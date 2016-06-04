#!/usr/bin/env python3

################################################################################
#  Soubor:    get_tRNA.py                                                      #
#  Autor:     Radim KUBIŠ, xkubis03                                            #
#  Vytvořeno: 18. dubna 2014                                                   #
#                                                                              #
#  Projekt do předmětu Bioinformatika (BIF).                                   #
#                                                                              #
#  Skript vyhledává sekvence tRNA genů v genomu. Jediným POVINNÝM parametrem   #
#  skriptu je název FASTA souboru, který obsahuje genom.                       #
#                                                                              #
#  Výstupem je multi-FASTA formát obsahující seznam nalezených tRNA genů.      #
#                                                                              #
################################################################################

# Import systémové knihovny
import sys
# Import knihovny pro regulární výrazy
import re

# Třída pro jeden záznam tRNA genu
class tRNA():
    #------------------------------------------------------------#
    # Konstruktor záznamu tRNA genu v genomu                     #
    #                                                            #
    # match     - nález regulárního výrazu                       #
    # direction - směr vlákna nálezu                             #
    # position  - pozice začátku genu od začátku přímého vlákna  #
    # id        - pořadové číslo genu v genomu                   #
    #------------------------------------------------------------#
    def __init__(self, match, direction, position, id):
        # Uložení směru vlákna
        self.direction = direction
        # Uložení aminokyseliny podle antikodonu
        self.acid = getAcid(match.group(1))
        # Uložení antikodonu
        self.anticodon = match.group(1)
        # Uložení délky genu
        self.length = len(match.group(0))
        # Uložení pozice začátku genu
        self.position = position
        # Uložení pořadového čísla genu
        self.id = id
        # Uložení bází genu
        self.bases = match.group(0)

    #------------------------------------------------------------#
    # Funkce pro tisk FASTA záznamu tRNA genu                    #
    #------------------------------------------------------------#
    def printFASTA(self):
        # Tisk hlavičky FASTA záznamu
        print(">tRNA_%(id)s|%(aa)s%(ac)s|%(pos)s|%(len)s|%(dir)s" % {
            'id': self.id,
            'aa': self.acid,
            'ac': self.anticodon,
            'pos': self.position,
            'len': self.length,
            'dir': self.direction
        })
        # Tisk řádků s bázemi ve FASTA formátu
        print(self.bases[:60])
        print(self.bases[60:])

# TABULKA KODONŮ
# Seznam bází
bases = ['T', 'C', 'A', 'G']
# Vytvoření seznam kodonů podle bází
codons = [i+j+k for i in bases for j in bases for k in bases]
# Pořadí aminokyselin v seznamu kodonů
acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
# Vytvoření tabulky kodonů s aminokyselinami
table = dict(zip(codons, acids))
# Slovník zkratek aminokyselin
acidDictionary = {
    'F': 'Phe', 'L': 'Leu', 'I': 'Ile', 'M': 'Met', 'V': 'Val', 'S': 'Ser',
    'P': 'Pro', 'T': 'Thr', 'A': 'Ala', 'Y': 'Tyr', '*': 'Stp', 'H': 'His',
    'Q': 'Gln', 'N': 'Asn', 'K': 'Lys', 'D': 'Asp', 'E': 'Glu', 'C': 'Cys',
    'W': 'Trp', 'R': 'Arg', 'G': 'Gly'
}
# Skupina antikodonů, které kódují aminokyseliny (bez stop kodonů)
anticodonGroup = "(TTT|TTC|TTG|TCT|TCC|TCG|TAT|TAC|TAA|TAG|TGT|TGC|TGA|TGG|CTT|CTC|CTG|CCT|CCC|CCA|CCG|CAT|CAC|CAA|CAG|CGT|CGC|CGA|CGG|ATT|ATC|ATA|ATG|ACT|ACC|ACA|ACG|AAT|AAC|AAA|AAG|AGT|AGC|AGA|AGG|GTT|GTC|GTA|GTG|GCT|GCC|GCA|GCG|GAT|GAC|GAA|GAG|GGT|GGC|GGA|GGG)"

#------------------------------------------------------------------------------#
#  Funkce vracející komplementární řetězec k zadanému řetězci bází             #
#                                                                              #
#  string - řetězec, ke kterému funkce vrací řetězec komplementární            #
#------------------------------------------------------------------------------#
def getComplementaryString(string):
    # Reverzace řetězce
    reverse = string[::-1]
    # Komplementarita řetězce ATGC->1234
    numbers = reverse.replace('A', '1').replace('T', '2').replace('G', '3').replace('C', '4')
    # Komplementarita řetězce 1234->TACG
    complement = numbers.replace('1', 'T').replace('2', 'A').replace('3', 'C').replace('4', 'G')
    # Navrácení komplementu
    return complement

#------------------------------------------------------------------------------#
#  Funkce vracející zkratku aminokyseliny na základě antikodonu                #
#                                                                              #
#  anticodon - antikodon, pro který funkce hledá aminokyselinu                 #
#------------------------------------------------------------------------------#
def getAcid(anticodon):
    # Získání komplementárního řetězce k antikodonu
    index = getComplementaryString(anticodon)
    # Navrácení zkratky aminokyseliny z tabulky
    return acidDictionary[table[index]]

# KONTROLA VSTUPNÍHO PARAMETRU
# Pokud není zadán povinný parametr skriptu
if len(sys.argv) < 2:
    # Tisk chyby
    sys.stderr.write("ERROR: Required argument <input_file> is not set\n")
    # Vyprázdnění výstupního bufferu standardního chybového výstupu
    sys.stderr.flush()
    # Tisk nápovědy použití skriptu
    print("Usage: %s <input_file>\n" % sys.argv[0])
    print("    <input_file> - FASTA file with genome sequence (required)")
    # Ukončení skriptu s chybou
    sys.exit(1)

# Proměnná pro přímý řetězec
directString = ""
# Proměnná pro komplementární řetězec
complementaryString = ""

try:
    # Otevření vstupního souboru
    inputFile = open(sys.argv[1], "r")
    # Přeskočení hlavičky vstupního souboru
    inputFile.readline()
    # Převod z FASTA formátu na jednořádkový řetězec
    for line in inputFile.readlines():
        directString += line.strip()
    # Uzavření vstupního souboru
    inputFile.close()
except FileNotFoundError:
    # Pokud vstupní soubor neexistuje,
    # tisk chyby
    sys.stderr.write("ERROR: File '%s' does not exist\n" % sys.argv[1])
    # Ukončení skriptu s chybou
    sys.exit(1)

# Získání komplementárního řetězce
complementaryString = getComplementaryString(directString)

# Regulární výraz pro hledání tRNA genů
expression = re.compile("[ATGC]{12,13}A[AG][ATGC]{1,3}G[ATGC]{11,14}[ACT]T"+anticodonGroup+"[AG][ATGC]{11,31}GTTC[AG]A[ATGC][TC]C[ATGC]{12}CCA")

# Seznam nalezených tRNA genů
tRNAs = []
# Pořadí nalezeného tRNA genu, číslování od jedničky
id = 1
# Uložení délky celého genomu
genomeLength = len(directString)

# Procházení výskytů regulárního výrazu na přímém vlákně
for match in expression.finditer(directString):
    # Uložení nalezeného tRNA genu do seznamu
    tRNAs.append(tRNA(match, "+", (match.start() + 1), id))
    # Inkrementace pořadového čísla
    id += 1

# Procházení výskytů regulárního výrazu na komplementárním vlákně
for match in expression.finditer(complementaryString):
    # Uložení nalezeného tRNA genu do seznamu
    tRNAs.append(tRNA(match, "-", (genomeLength - match.start()), id))
    # Inkrementace pořadového čísla
    id += 1

# Procházení uložených tRNA genů
for item in tRNAs:
    # Tisk FASTA formátu tRNA genu na výstup
    item.printFASTA()

# Ukončení skriptu bez chyby
sys.exit(0)
