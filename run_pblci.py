import subprocess
import sys
import re
import make_lncRNA_matrix
import pblci


# file to perform 5-fold CV correctly by creating a steve (heterogeneous matrix)
# from input of rna_to_disease_interactions file for each training sample


# method that calls the r script to create a disease-disease similarity matrix
def make_disease_sim(doid_set: set):
    # write list of do ids we want to create a disease sim matrix
    with open("target_do_list.txt", "w") as writer:
        writer.write(",".join(doid_set) + "\n")

    # make the disease sim matrix, writes it to d_sim_matrix.text
    subprocess.call(["/usr/local/bin/Rscript", "--vanilla",
                     "/Users/elizabethzhang/csci1810/csci2810_rproject/make_d_sim_matrix .R"])

    # read the disease sim matrix into steve
    hm = {}
    with open("d_sim_matrix.txt", "r") as reader:
        line = reader.readline().strip()
        do_id_list = re.findall('"([^"]*)"', line)
        num_dis = len(do_id_list)

        count = 0
        for raw_ln in reader:
            ln: str = raw_ln.strip()
            score_arr = ln.split(",")
            for i in range(0, num_dis):
                this_dis = do_id_list[count]
                if this_dis not in hm.keys():
                    hm[this_dis] = {}
                hm[this_dis][do_id_list[i]] = score_arr[i]
            count = count + 1

    return hm


# makes the lncRNA-lncRNA similarity matrix
def make_lncrna_sim(rna_to_doid: dict, partial_steve: dict) -> dict:
    make_lncRNA_matrix.run("d_sim_matrix.txt", rna_to_doid,
                           "lncrna_sim_matrix.txt")

    with open("lncrna_sim_matrix.txt", 'r') as reader1:
        line = reader1.readline().strip()
        list_of_lncRNAs = line.split(",")
        num_rnas = len(list_of_lncRNAs)

        count = 0
        for ln in reader1:
            score_arr = ln.split(",")
            for i in range(0, num_rnas):
                partial_steve[list_of_lncRNAs[count]][list_of_lncRNAs[i]] = score_arr[i]
            count = count + 1

    return partial_steve

# makes heterogenous matrix, represented by a nested hashmap
def make_network_dict(test_dict, interactions: dict):
    steve_interactions = read_interactions_dict(test_dict, interactions)
    doid_set: set = make_doid_set(interactions)
    steve_diseases = make_disease_sim(doid_set)
    steve_rnas = make_lncrna_sim(interactions, steve_interactions)
    final_steve = {**steve_diseases, **steve_rnas}
    return interactions, doid_set, final_steve


# makes heterogenous matrix, represented by a nested hashmap
def make_network_file(interactions: str):
    lrna_doid, steve_interactions = run_sample_from_file(interactions)
    doid_set: set = make_doid_set(lrna_doid)
    steve_diseases = make_disease_sim(doid_set)
    steve_rnas = make_lncrna_sim(lrna_doid, steve_interactions)
    final_steve = {**steve_diseases, **steve_rnas}
    return lrna_doid, doid_set, final_steve


# reads interactions file and makes a dictionary that maps
# lncs to the DOIDs they interact with
def read_interactions_file(interactions_file) -> (dict, dict):
    lncrna_to_doid: dict = {}
    steve_subset = {}

    with open(interactions_file, "r") as reader:
        for raw_line in reader:
            line: str = raw_line.strip()
            split: list = line.split(";")
            key: str = split[0]
            values: list = [x for x in split[1].split(",") if x]
            lncrna_to_doid[key] = values
            steve_subset[key] = {}
            for v in values:
                steve_subset[key][v] = 1.0

    return lncrna_to_doid, steve_subset


# makes doid set
def make_doid_set(rna_to_doid: dict) -> set:
    all_doid_list = []
    for disease_list in rna_to_doid.values():
        all_doid_list = all_doid_list + disease_list
    all_doid_list.sort()
    cancer_set = set(all_doid_list)
    return cancer_set


def read_interactions_dict(test_dict, rna_to_disease_dict: dict):
    steve_subset = {}

    for rna, dis_lst in rna_to_disease_dict.items():
        steve_subset[rna] = {}
        for d in dis_lst:
            if rna not in test_dict.keys():
                steve_subset[rna][d] = 1.0

    return steve_subset


def run_sample_from_dict(test_dict, rna_to_disease_dict: dict) -> list[tuple]:
    rna_to_dis_to_score: dict[tuple, float] = {}
    sample_dict, all_diseases, my_network = make_network_dict(test_dict, rna_to_disease_dict)

    for rna in sample_dict.keys():
        for disease in all_diseases:
            my_score = pblci.get_pblci_score(rna, disease, my_network)
            rna_to_dis_to_score[rna, disease] = my_score

    sorted_by_score = sorted(rna_to_dis_to_score, key=rna_to_dis_to_score.get,
                             reverse=True)
    return sorted_by_score


def run_sample_from_file(interactions_str: str) -> list[tuple]:
    rna_to_dis_to_score: dict[tuple, float] = {}
    sample_dict, all_diseases, my_network = make_network_file(interactions_str)

    for rna in sample_dict.keys():
        for disease in all_diseases:
            my_score = pblci.get_pblci_score(rna, disease, my_network)
            rna_to_dis_to_score[rna, disease] = my_score

    sorted_by_score = sorted(rna_to_dis_to_score, key=rna_to_dis_to_score.get)
    return sorted_by_score


# method to run file
def run_dict(test_dict, interactions: dict):
    print(run_sample_from_dict(test_dict, interactions))


# method to run file
def run_file(interactions: str):
    print(run_sample_from_file(interactions))


if __name__ == '__main__':
    run_file(sys.argv[1])
