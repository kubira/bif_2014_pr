#!/usr/bin/env python3

################################################################################
#  Soubor:    compare_tRNA.py                                                  #
#  Autor:     Radim KUBIŠ, xkubis03                                            #
#  Vytvořeno: 27. dubna 2014                                                   #
#                                                                              #
#  Projekt do předmětu Bioinformatika (BIF).                                   #
#                                                                              #
#  Skript  porovnává   multi-FASTA  soubor  tRNA   genů  nalezených  skriptem  #
#  get_tRNA.py  s  multi-FASTA souborem známých tRNA  genů. Při spuštění jsou  #
#  očekávány dva POVINNÉ parametry - názvy obou multi-FASTA souborů.           #
#                                                                              #
#  Výstupem jsou dva seznamy hlaviček FASTA záznamů. První obsahuje  hlavičky  #
#  nalezených genů, které nejsou součástí známých tRNA struktur. Druhý seznam  #
#  obsahuje hlavičky  známých tRNA genů,  které  skript get_tRNA.py nenalezl.  #
#  Seznamy jsou odděleny prázdným řádkem.                                      #
#                                                                              #
################################################################################

# Import systémového modulu
import sys
# Import modulu regulárních výrazů
import re

# Třída pro FASTA záznam
class FASTAItem():
    #-------------------------------------#
    # Inicializační funkce FASTA záznamu  #
    #                                     #
    # head     - FASTA hlavička           #
    # sequence - FASTA sekvence           #
    #-------------------------------------#
    def __init__(self, head, sequence):
        # Uložení hlavičky
        self.head = head
        # Uložení sekvence
        self.sequence = sequence

# KONTROLA VSTUPNÍCH PARAMETRŮ
# Pokud nejsou zadány povinné parametry skriptu
if len(sys.argv) < 3:
    # Pokud nebyl zadán žádný soubor
    if len(sys.argv) < 2:
        # Tisk chyby uživatelského souboru
        sys.stderr.write("ERROR: Required argument <your_file> is not set\n")
    # Tisk chyby souboru z databáze
    sys.stderr.write("ERROR: Required argument <database_file> is not set\n")
    # Vyprázdnění výstupního bufferu standardního chybového výstupu
    sys.stderr.flush()
    # Tisk nápovědy použití skriptu
    print("Usage: %s <your_file> <database_file>\n" % sys.argv[0])
    print("    <your_file>     - multi-FASTA file from get_tRNA.py (required)")
    print("    <database_file> - multi-FASTA file from tRNA database (required)")
    # Ukončení skriptu s chybou
    sys.exit(1)

# Seznam pro nalezená tRNA
tRNAUser = []
# Seznam pro databázová tRNA
tRNADatabase = []

try:
    # Otevření vstupního souboru s nalezenými tRNA
    userFile = open(sys.argv[1], "r")
    # Načtení celého souboru a rozdělení podle konce řádku do seznamu
    lines = userFile.read().split('\n')
    
    # Nastavení indexu do seznamu řádků na první pozici
    index = 0
    # Dokud nejsem na konci seznamu
    while index < len(lines):
        # Přeskočení prázdných řádků mezi jednotlivými záznamy
        if lines[index] == '':
            # Posun na další řádek
            index += 1
            # Další krok cyklu
            continue
        # Vytvoření nového objektu FASTA záznamu
        tRNAUser.append(FASTAItem(lines[index], lines[index+1]+lines[index+2]))
        # Posunutí pozice v seznamu řádků o 3 řádky (jeden záznam)
        index += 3
    
    # Uzavření vstupního souboru s nalezenými tRNA
    userFile.close()
except FileNotFoundError:
    # Pokud vstupní soubor s nalezenými tRNA neexistuje,
    # tisk chyby
    sys.stderr.write("ERROR: Your file '%s' does not exist\n" % sys.argv[1])
    # Ukončení skriptu s chybou
    sys.exit(1)

try:
    # Otevření vstupního souboru se známými tRNA
    databaseFile = open(sys.argv[2], "r")
    # Načtení celého souboru a rozdělení podle konce řádku do seznamu
    lines = databaseFile.read().split('\n')
    
    # Nastavení indexu do seznamu řádků na první pozici
    index = 0
    # Dokud nejsem na konci seznamu
    while index < len(lines):
        # Přeskočení prázdných řádků mezi jednotlivými záznamy
        if lines[index] == '':
            # Posun na další řádek
            index += 1
            # Další krok cyklu
            continue
        # Vytvoření nového objektu FASTA záznamu
        tRNADatabase.append(FASTAItem(lines[index], lines[index+1]+lines[index+2]))
        # Nastavení databázového tRNA jako nenalezeného
        tRNADatabase[-1].found = 0
        # Posunutí pozice v seznamu řádků o 3 řádky (jeden záznam)
        index += 3
    
    # Uzavření vstupního souboru se známými tRNA
    databaseFile.close()
except FileNotFoundError:
    # Pokud vstupní soubor se známými tRNA neexistuje,
    # tisk chyby
    sys.stderr.write("ERROR: Database file '%s' does not exist\n" % sys.argv[2])
    # Ukončení skriptu s chybou
    sys.exit(1)

# Procházení databázových tRNA
for dbItem in tRNADatabase:
    # Procházení uživatelských tRNA, která zatím nebyla nalezena v databázových
    # a jejich porovnání s aktuálním databázovým tRNA
    for usrItem in tRNAUser:
        # Porovnání aktuálního databázového tRNA s aktuálním uživatelským tRNA
        if re.match(".*"+dbItem.sequence+"\Z", usrItem.sequence):
            # Pokud aktuální databázový tRNA souhlasí s aktuálním uživatelským
            
            # Výpočet procentuálního překryvu obou tRNA
            podil = len(dbItem.sequence)/len(usrItem.sequence)
            # Pokud je překryv alespoň 80%
            if podil >= 0.80:
                # Vymazání uživatelského tRNA ze seznamu
                tRNAUser.remove(usrItem)
                # Označení databázového tRNA jako nalezeného
                dbItem.found = 1
                # Konec procházení uživatelských tRNA
                break

# Výpis zbylých/nekorespondujících uživatelských tRNA
for item in tRNAUser:
    # Výpis jednoho zbylého uživatelského tRNA
    print(item.head[1:])

# Tisk prázdného řádku
print()

# Výpis nenalezených databázových tRNA
for item in tRNADatabase:
    # Pokud nebylo tRNA nalezeno
    if item.found == 0:
        # Výpis hlavičky nenalezeného databázového tRNA
        print(item.head)

# Ukončení skriptu bez chyby
sys.exit(0)
