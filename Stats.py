import pandas as pd
import numpy as np
import argparse, sys, os, re

def basic_metrics(data, val):
    max = np.max(data)
    min = np.min(data)
    mean = np.mean(data)
    std = np.std(data)
    median = np.median(data)
    print(f"Basic stats for {val}:\nMean = {mean}\nStd = {std}\nMedian = {median}\nMax = {max}\nMin = {min}\n")

def sum_harmo(exp, diff):
    res = 0
    for dif in diff:
        res += 1/(dif**exp)
    return res

def moy_harmo(exp, diff, size):
    sum_h = sum_harmo(exp, diff)
    res = (size/abs(sum_h))**1/exp
    return res

def writeCSV(data, out, type):
    #type = "size" OR "ram" OR "time"
    #print(data)
    with open(out, 'w') as out_file:
        out_file.write("value,type,threshold,"+type+",tool\n")
        for key in data:
            for name in data[key]:
                if('diff' in data[key][name]):
                    out_file.write(str(data[key][name]['diff']))
                    out_file.write(",")
                    out_file.write("error")
                    out_file.write(",")
                    out_file.write(key)
                    out_file.write(",")
                    out_file.write(str(data[key][name]['type']))
                    out_file.write(",")
                    out_file.write(name)
                    out_file.write("\n")
                else:
                    print(f"Missing differences for subsampling rate {key}, tool is {name}. Maybe comparisons did not go through ?")

def get_error(res_spsp, res_simka, out, type):
    print(f"Type of data is {type}")
    #READING SIMKA
    simka = pd.read_csv(res_simka, sep = ";", header = 0)
    simka = simka.drop(simka.columns[[0]], axis = 1)
    simka =  simka.applymap(lambda x: 1-x)
    simka = simka.to_numpy()
    #ltri_simka = simka[np.tril_indices_from(simka, k = -1)]

    files_spsp = []
    with open(res_spsp, 'r') as fof_spsp:
        line = fof_spsp.readline().strip()
        while line != "":
            files_spsp.append(line)
            line = fof_spsp.readline().strip()

    #COMPUTING DIFFERENCES FOR SPSP + SAVING IN DICT
    data = {}
    for i in range(len(files_spsp)):
        spsp = pd.read_csv(files_spsp[i], sep = ",", header = 0)
        spsp = spsp.to_numpy()
        #ltri_spsp = spsp[np.tril_indices_from(spsp, k = -1)]
        name = files_spsp[i].split("/")[-1]
        key = [s for s in re.findall(r'\d+', name)][0]
        tmp = "SuperSampler_decycling"
        if key in data:
            if tmp in data[key]:
                #MEAN DIFF
                #data[key][tmp]['diff'] = np.mean(abs(simka - spsp)/simka)
                #SUM DIFF
                data[key][tmp]['type'] = abs(np.mean(simka) - np.mean(spsp))
                #data[key][tmp]['diff'] = np.mean(abs(ltri_simka - ltri_spsp))
            else:
                data[key][tmp] = {}
                data[key][tmp]['type'] = abs(np.mean(simka) - np.mean(spsp))
        else:
            data[key] = {}
            data[key][tmp] = {}
            data[key][tmp]['type'] = abs(np.mean(simka) - np.mean(spsp))

    with open(out, 'w') as out_file:
        out_file.write("value,type,threshold,"+type+",tool\n")
        for key in data:
            for name in data[key]:
                out_file.write(str(data[key][name]['type']))
                out_file.write(",")
                out_file.write("error")
                out_file.write(",")
                out_file.write(key)
                out_file.write(",")
                out_file.write("0")
                out_file.write(",")
                out_file.write(name)
                out_file.write("\n")

def read_bench(fof, type):
    data = {}
    with open(fof, 'r') as file_of_file:
        name = file_of_file.readline().strip()
        while name != "":
            tmp = name.split("/")[-1]
            key = [s for s in re.findall(r'\d+', tmp)][0]
            parts = name.split("_")
            with open(name, 'r') as bench:
                skip = bench.readline();
                values = bench.readline().strip().split("\t")
                if len(parts) > 5:
                    tool = "SuperSampler_m"+parts[3]
                else:
                    tool = "sourmash"
                if key in data:
                    if tool in data[key]:
                        if type == "ram":
                            data[key][tool]['type'] = values[2]
                        else:
                            data[key][tool]['type'] = values[-1]
                    else:
                        data[key][tool] = {}
                        if type == "ram":
                            data[key][tool]['type'] = values[2]
                        else:
                            data[key][tool]['type'] = values[-1]
                else:
                    data[key] = {}
                    data[key][tool] = {}
                    if type == "ram":
                        data[key][tool]['type'] = values[2]
                    else:
                        data[key][tool]['type'] = values[-1]
            name = file_of_file.readline().strip()
    return data

            

def read_index_size(sub_sourmash, sub_spsp):
    data = {}
    with open(sub_sourmash, 'r') as fof_sub_sourmash:
        line = fof_sub_sourmash.readline().strip()
        while line != "":
            name = "sourmash_zipped"
            #SAVING NAMES FOR DISPLAY IN FIGURE
            tmp = line.split("/")[-1]
            key = [s for s in re.findall(r'\d+', tmp)][0]
            if key in data:
                data[key][name] = {}
                data[key][name]['type'] = os.stat(line).st_size/(1024*1024)
            else:
                data[key] = {}
                data[key][name] = {}
                data[key][name]['type'] = os.stat(line).st_size/(1024*1024)
            line = fof_sub_sourmash.readline().strip()
    with open(sub_spsp, 'r') as fof_sub_spsp:
        f_name = fof_sub_spsp.readline().strip()
        while f_name != "":
            size = 0
            #SAVING NAMES FOR DISPLAY IN FIGURE
            tmp = f_name.split("/")[-1]
            name = "SuperSampler_m" + [s for s in re.findall(r'\d+', tmp)][1]
            key = [s for s in re.findall(r'\d+', tmp)][0]
            if key in data:
                data[key][name] = {}
                data[key][name]['type'] = os.stat(f_name).st_size/(1024*1024)
            else:
                data[key] = {}
                data[key][name] = {}
                data[key][name]['type'] = os.stat(f_name).st_size/(1024*1024)
            f_name = fof_sub_spsp.readline().strip()
    return data

def compare_results(res_spsp, res_sourmash, res_simka, data, out, type):
    print(f"Type of data is {type}")
    #READING SIMKA
    simka = pd.read_csv(res_simka, sep = ";", header = 0)
    simka = simka.drop(simka.columns[[0]], axis = 1)
    simka =  simka.applymap(lambda x: 1-x)
    simka = simka.to_numpy()
    #ltri_simka = simka[np.tril_indices_from(simka, k = -1)]

    #GETTING EVERY CSV FILENAME
    files_spsp, files_sourmash = [], []
    with open(res_spsp, 'r') as fof_spsp:
        line = fof_spsp.readline().strip()
        while line != "":
            files_spsp.append(line)
            line = fof_spsp.readline().strip()

    with open(res_sourmash, 'r') as fof_sourmash:
        line = fof_sourmash.readline().strip()
        while line != "":
            files_sourmash.append(line)
            if type == "size":
                files_sourmash.append(line)
            line = fof_sourmash.readline().strip()

    #COMPUTING DIFFERENCES FOR SPSP + SAVING IN DICT
    for i in range(len(files_spsp)):
        spsp = pd.read_csv(files_spsp[i], sep = ",", header = 0)
        spsp = spsp.to_numpy()
        #ltri_spsp = spsp[np.tril_indices_from(spsp, k = -1)]
        name = files_spsp[i].split("/")[-1]
        key = [s for s in re.findall(r'\d+', name)][0]
        tmp = "SuperSampler_m" + [s for s in re.findall(r'\d+', name)][1]
        if key in data:
            if tmp in data[key]:
                #MEAN DIFF
                #data[key][tmp]['diff'] = np.mean(abs(simka - spsp)/simka)
                #SUM DIFF
                data[key][tmp]['diff'] = abs(np.mean(simka) - np.mean(spsp))
                #data[key][tmp]['diff'] = np.mean(abs(ltri_simka - ltri_spsp))
            else:
                print(f"should not happen, {tmp} not in dict[{key}]")
        else:
            print(f"should not happen, {key} not in dict.")

    #COMPUTING DIFFERENCES FOR SOURMASH + SAVING IN DICT
    for i in range(len(files_sourmash)):
        sourmash = pd.read_csv(files_sourmash[i], sep = ",", header = 0)
        sourmash = sourmash.to_numpy()
        #ltri_sm = sourmash[np.tril_indices_from(sourmash, k = -1)]
        name = files_sourmash[i].split("/")[-1]
        key = [s for s in re.findall(r'\d+', name)][0]
        tmp = "sourmash"
        if type == "size":
                tmp += "_zipped"
        if key in data:
            if tmp in data[key]:
                #MEAN DIFF
                #data[key][tmp]['diff'] = np.mean(abs(simka - sourmash)/simka)
                #SUM DIFF
                data[key][tmp]['diff'] = abs(np.mean(simka) - np.mean(sourmash))
                #data[key][tmp]['diff'] = np.mean(abs(ltri_simka - ltri_sm))
            else:
                print(f"should not happen, {tmp} not in dict[{key}]")
        else:
            print(f"should not happen, {key} not in dict.")

    #print(data)
    #norm_simka_sm = np.linalg.norm(diff_simka_sm,1)
    #norm_simka_spsp = np.linalg.norm(diff_simka_spsp,1)
    #norm_simka_sm2 = np.linalg.norm(diff_simka_sm,2)
    #norm_simka_spsp2 = np.linalg.norm(diff_simka_spsp,2)
    #print(f"Norm1 of differences between Simka and Sourmash: {norm_simka_sm}\nNorm of differences between Simka and Supersampler: {norm_simka_spsp}")
    #print(f"Norm2 of differences between Simka and Sourmash: {norm_simka_sm2}\nNorm of differences between Simka and Supersampler: {norm_simka_spsp2}\n")

    #basic_metrics(differences[0], "Differences Simka Sourmash")
    #basic_metrics(differences[1], "Differences Simka SuPerSamPler")
    writeCSV(data, out, type)
    #moy_harmo1 = moy_harmo(1, diff_simka_sm, len(diff_simka_sm))
    #moy_harmo2 = moy_harmo(2, diff_simka_sm, len(diff_simka_sm))
    #moy_harmo3 = moy_harmo(3, diff_simka_sm, len(diff_simka_sm))
    #print(f"Moyenne harmonique simka VS sourmash (1): {moy_harmo1}\nMoyenne harmonique simka VS sourmash (2): {moy_harmo2}\nMoyenne harmonique simka VS sourmash (3): {moy_harmo3}\n")

    #moy_harmo1 = moy_harmo(1, diff_simka_spsp, len(diff_simka_spsp))
    #moy_harmo2 = moy_harmo(2, diff_simka_spsp, len(diff_simka_spsp))
    #moy_harmo3 = moy_harmo(3, diff_simka_spsp, len(diff_simka_spsp))
    #print(f"Moyenne harmonique simka VS supersampler (1): {moy_harmo1}\nMoyenne harmonique simka VS supersampler (2): {moy_harmo2}\nMoyenne harmonique simka VS supersampler (3): {moy_harmo3}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Stats comparing sourmash and SPSP')
    parser.add_argument('spsp', help='file of file of results for SuPerSamPler')
    parser.add_argument('sourmash', help='file of file of results for sourmash')
    parser.add_argument('simka', help='results for SimKa')
    parser.add_argument('--subspsp', required=False, help='Subsampling file of file for SPSP (needed to output size)')
    parser.add_argument('--subsm', required=False, help='Subsampling file for sourmash (needed to output size)')
    parser.add_argument('-b', required=False, help="File of file for bench on comparisons (Needed for graph with ram and time).")
    parser.add_argument('-t', help="Type of graph wanted (size, ram, time).")
    parser.add_argument('-o', help='Out filename')
    args = parser.parse_args(sys.argv[1:])
    if args.t == "size": 
        data = read_index_size(args.subsm, args.subspsp)
    elif args.t == "ram" or args.t == 'time':
        data = read_bench(args.b, args.t)
    elif args.t == "error":
        get_error(args.spsp, args.simka, args.o, args.t)
        exit()
    else:
        sys.exit("INVALID VALUE FOR TYPE.")
        
    compare_results(args.spsp, args.sourmash, args.simka, data, args.o, args.t)

""""""
