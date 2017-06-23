import os

PATH_NAME = "./ligands_data_set"
DB_NAME = "db.mol2"
STATISTICS_NAME = "statistics.txt"

MOLECULE_BEGIN = "@<TRIPOS>MOLECULE\n"
MOLECULE_ATOMS = "@<TRIPOS>ATOM\n"
MOLECULE_BONDS = "@<TRIPOS>BOND\n"

def existence(a,b):

    s = set(b)

    for i,x in enumerate(a):
        if x in s:
            return True

    return False


def createDB():

    numbers =[]
    characters = []

    for n in range(0,9):
        numbers += str(n)

    for c in range(65,91):
        characters += chr(c)

    for c in range(97, 123):
        characters += chr(c)

    alfanumeric = numbers + characters

    n_file = 0

    max_molecules = 0
    max_atoms = 0
    max_bonds = 0

    #total (for all the database)
    t_mols = 0
    t_atoms = 0
    t_bonds = 0


    DB_file = open(DB_NAME, 'w')
    statistics_file = open(STATISTICS_NAME, 'w')


    directory = os.fsencode(PATH_NAME)

    for file in os.listdir(PATH_NAME):

        filename = PATH_NAME + '/' + os.fsdecode(file)

        if filename.endswith(".mol2"):
            
            print ("Opening the: <" + filename + ">") 
            file = open(filename, 'r')
            
            n_file += 1

            #total_file
            t_f_mols = 0
            t_f_atoms = 0
            t_f_bonds = 0

            #number (for each molecule)
            n_atoms = 0
            n_bonds = 0

            count_atoms = False
            count_bonds = False

            #sto leggendo un solo file alla volta
            for line in file:

                if not existence(line, alfanumeric):
                    continue
                
                if count_atoms and existence(line, numbers):
                    t_atoms += 1
                    t_f_atoms += 1
                    n_atoms += 1

                elif count_bonds and existence(line, numbers):
                    t_bonds += 1
                    t_f_bonds += 1
                    n_bonds += 1

                if MOLECULE_BEGIN in line:
                    t_mols += 1
                    t_f_mols += 1
                    count_bonds = False

                    if n_atoms > max_atoms:
                        max_atoms = n_atoms
                    if n_bonds > max_bonds:
                        max_bonds = n_bonds

                    n_atoms = 0
                    n_bonds = 0

                elif MOLECULE_ATOMS in line:
                    count_atoms = True

                elif MOLECULE_BONDS in line:
                    count_atoms = False
                    count_bonds = True

                DB_file.write(line)

            str_res_file = filename + " has: " + str(t_f_mols) + " molecules, " + str(t_f_atoms) + " atoms, " + str(t_f_bonds) + " bonds.\n"

            print (str_res_file)
            statistics_file.write(str_res_file)

            file.close()

            if n_atoms > max_atoms:
                max_atoms = n_atoms
            if n_bonds > max_bonds:
                max_bonds = n_bonds
            if t_f_mols > max_molecules:
                max_molecules = t_f_mols
            
            t_f_mols = 0

        else:
            print ("Not a .mol2 file\n")

    DB_file.close()

    print ("I read " + str(n_file) + " files.")

    str_DB = "\nThe database contains: " + str(n_file) + " files, " + str(t_mols) + " molecules, " + str(t_atoms) + " atoms, " + str(t_bonds) + " bonds."
    print (str_DB)
    statistics_file.write(str_DB)
    
    string_max_mols = "The max number of molecules in a file is " + str(max_molecules)
    string_max_atoms = "The max number of atoms in a file is " + str(max_atoms)
    string_max_bonds = "The max number of bonds in a file is " + str(max_bonds)
    print (string_max_mols)
    print (string_max_atoms)
    print (string_max_bonds)

    statistics_file.write(string_max_mols)
    statistics_file.write('\n')
    statistics_file.write(string_max_atoms)
    statistics_file.write('\n')
    statistics_file.write(string_max_bonds)
    statistics_file.write('\n')

    statistics_file.close()


if __name__ == "__main__":
    createDB()
    print ("\n\nDB created.")
    print ("Adios!")
