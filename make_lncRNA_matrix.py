import sys
import re


# read in disease similarity matrix
def make_d_sim_matrix(fn: str):
    disease_matrix = []
    with open(fn, "r") as reader:
        line_1 = reader.readline().strip()
        do_id_list = re.findall('"([^"]*)"', line_1)
        for ln in reader:
            ln_as_list = [float(x) for x in ln.strip().split(",")]
            disease_matrix.append(ln_as_list)
    return do_id_list, disease_matrix


# read in lncRNA to disease group dictionary
def make_lncrna_to_d_dict(fn: str) -> dict[str, set[str]]:
    lncrna_to_d_dict: dict[str, set[str]] = {}
    with open(fn, "r") as reader:
        for ln in reader:
            split = ln.strip().split(";")
            key: str = split[0]
            values: set = set(split[1].split(",")[:-1])
            lncrna_to_d_dict[key] = values
    return lncrna_to_d_dict


# calculate similarity between two lncRNAs
def lncrna_sim_score(id_list: list, d_sim_matrix: list, lnc_to_d_group: dict,
                     lncrna_a: str, lncrna_b: str):
    def d_to_group_sim(d1: str, group_b: set[str]):
        max_score: float = float("-inf")
        for d2 in group_b:
            d1_index: int = id_list.index(d1)
            d2_index: int = id_list.index(d2)
            current_score = d_sim_matrix[d1_index][d2_index]
            if current_score > max_score:
                max_score = current_score
        return max_score

    # calculate similarity between two given disease groups
    def group_to_group_sim(group_a: set[str], group_b: set[str]):
        group_a_sim: float = 0
        group_b_sim: float = 0
        for d1 in group_a:
            group_a_sim += d_to_group_sim(d1, group_b)
        for d2 in group_b:
            group_b_sim += d_to_group_sim(d2, group_a)
        return group_a_sim, group_b_sim

    d_group_a: set[str] = lnc_to_d_group[lncrna_a]
    d_group_b: set[str] = lnc_to_d_group[lncrna_b]

    d_group_a_sim, d_group_b_sim = group_to_group_sim(d_group_a, d_group_b)
    lncrna_sim = (d_group_a_sim + d_group_b_sim) / \
                 (len(d_group_a) + len(d_group_b))

    return lncrna_sim


def make_lncrna_sim_matrix(id_list: list, d_sim_matrix: list,
                           lnc_to_d_group: dict):
    lncrna_list: list = list(lnc_to_d_group.keys())
    matrix_size: int = len(lncrna_list)
    lncrna_sim_matrix: list[list[float]] = [[float("-inf")] * matrix_size
                                            for _ in range(matrix_size)]

    for i in range(matrix_size):
        for j in range(matrix_size):
            lncrna_a = lncrna_list[i]
            lncrna_b = lncrna_list[j]
            sim_score = lncrna_sim_score(id_list, d_sim_matrix, lnc_to_d_group,
                                         lncrna_a, lncrna_b)
            lncrna_sim_matrix[i][j] = sim_score
    return lncrna_list, lncrna_sim_matrix


def write_matrix_to_file(fn: str, lncrna_list: list, lncrna_sim_matrix):
    with open(fn, "w") as writer:
        lncrna_header: str = ",".join(lncrna_list)
        writer.write(lncrna_header + "\n")
        for inner_list in lncrna_sim_matrix:
            inner_list_as_str = [str(x) for x in inner_list]
            row_as_str: str = ",".join(inner_list_as_str)
            writer.write(row_as_str + "\n")


def run(d_sim_matrix_file: str, lnc_to_disease: dict, rna_sim_file: str):
    id_list, d_sim_matrix = make_d_sim_matrix(d_sim_matrix_file)
    lncrna_list, lncrna_sim_matrix = make_lncrna_sim_matrix(id_list,
                                                            d_sim_matrix,
                                                            lnc_to_disease)
    write_matrix_to_file(rna_sim_file, lncrna_list, lncrna_sim_matrix)


# [path to d_sim_matrix] [path to lnc_to_do_group]
if __name__ == "__main__":
    print("")
