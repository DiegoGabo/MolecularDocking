import subprocess
from tempfile import TemporaryFile
from xlwt import Workbook

DB_DIMENSION = 3961 #settare a 3961
ripetition = 5

delimitation = '-'*20
exetime_str = "Execution time: "
thr_str = "Throughput: "

M = "C++"
OC = "OPENCL-CPU"
OG = "OPENCL-GPU"


def run_program(name, n, l, ti_sh, th_sh):

    program = name
    l.write("I'm running {program} with {n_mols} molecules.\n".format(program=program, n_mols=str(n)))

    exetime_l = []
    thr_l = []

    for i in range(ripetition):
        with open("output.txt", 'w') as output:
            #disattivare per l'esecuzione
            # subprocess.call(["./test"], stdout=output, universal_newlines=True)
            
            #attivare per l'esecuzione
            r = ""
            n_opt = "--n"
            n_opt_value = str(n)
            dev_opt = "--d"
            dev_opt_value = ""
            
            if program == M:
                r = "./main"
                subprocess.call([r, n_opt, n_opt_value], stdout=output, stderr=error, universal_newlines=True)
            elif program == OC:
                r = "./opencl"
                dev_opt_value = "cpu"
                subprocess.call([r, n_opt, n_opt_value, dev_opt, dev_opt_value], stdout=output, stderr=error, universal_newlines=True)
            elif program == OG:
                r = "./opencl"
                dev_opt_value = "gpu"
                subprocess.call([r, n_opt, n_opt_value, dev_opt, dev_opt_value], stdout=output, stderr=error, universal_newlines=True)

        with open("output.txt", 'r') as output:
            for line in output.read().split('\n'):
                if line != '':
                    if exetime_str in line:
                        exetime_l.append(float(line.replace(exetime_str, "")))
                    elif thr_str in line:
                        thr_l.append(float(line.replace(thr_str, "")))
                    else:
                        continue
    
    exetime = 0
    for et in exetime_l:
        exetime += et
    exetime /= ripetition

    thr = 0
    for t in thr_l:
        thr += t
    thr /= ripetition
     
    l.write(exetime_str + str(exetime))
    l.write('\n')
    l.write(thr_str + str(thr))
    l.write('\n')

    ti_sh.write(n, 0, n)
    th_sh.write(n, 0, n)

    if program == M:
        ti_sh.write(n, 1, exetime)
        th_sh.write(n, 1, thr)
    elif program == OC:
        ti_sh.write(n, 2, exetime)
        th_sh.write(n, 2, thr)
    elif program == OG:
        ti_sh.write(n, 3, exetime)
        th_sh.write(n, 3, thr)



def benchmark():
    
    with open("log.txt", "w") as log:
        book = Workbook()
        time_sh = book.add_sheet('TIME', True)
        thr_sh = book.add_sheet('THROUGHPUT', True)

        time_sh.write(0,0, "NUM_MOLECULES")
        time_sh.write(0,1, M)
        time_sh.write(0,2, OC)
        time_sh.write(0,3, OG)

        thr_sh.write(0,0, "NUM_MOLECULES")
        thr_sh.write(0,1, M)
        thr_sh.write(0,2, OC)
        thr_sh.write(0,3, OG)
        step = 1		
        n_mols = 1
        
        while n_mols < DB_DIMENSION:    
            run_program(M, n_mols, log, time_sh, thr_sh)
            log.write('\n')

            run_program(OC, n_mols,log, time_sh, thr_sh)
            log.write('\n')

            run_program(OG, n_mols, log, time_sh, thr_sh)

            log.write(delimitation)
            log.write('\n')

            if n_mols in range(1, 11):
                step = 1
            elif n_mols in range(1, 100):
                step = 10
            elif n_mols in range(100, 1000):
                step = 100
            elif n_mols > 1000:
                step = 250

            n_mols += step

        book.save('benchmark.xls')
        book.save(TemporaryFile())


if __name__ == "__main__":
    benchmark();
