import sys


def read_data(samples_file: str, results_file: str):
    true_interactions = []  # list of list of tuples (lnc, doid)
    our_interactions = []

    list_to_insert1 = []
    with open(samples_file, "r") as sample_reader:
        for raw_ln in sample_reader:
            ln = raw_ln.strip()
            if ln.startswith("#"):
                if list_to_insert1:
                    true_interactions.append(list_to_insert1)
                list_to_insert1 = []
            else:
                split: list = ln.split(";")
                tuple_1: str = split[0]
                tuple_2_values: list = [x for x in split[1].split(",") if x]
                for tuple_2 in tuple_2_values:
                    list_to_insert1.append((tuple_1, tuple_2))
        if list_to_insert1:
            true_interactions.append(list_to_insert1)

    list_to_insert2 = []
    with open(results_file, "r") as result_reader:
        for raw_ln in result_reader:
            ln = raw_ln.strip()
            if ln.startswith("#"):
                if list_to_insert2:
                    our_interactions.append(list_to_insert2)
                list_to_insert2 = []
            else:
                split: list = ln.split(",")
                list_to_insert2.append((split[0], split[1]))
        if list_to_insert2:
            our_interactions.append(list_to_insert2)

    # need two columns: x values (FPR) =
    # percentage of samples whose ranks are lower than threshold
    # y values (TPR) =
    # percentage of samples whose ranks are higher than threshold

    return true_interactions, our_interactions


def calculate_TPR(inter_to_ranks_dict: dict, rank_threshold: int) -> float:
    dict_size: int = len(inter_to_ranks_dict)
    num_within_threshold: int = 0
    for interaction in inter_to_ranks_dict:
        if inter_to_ranks_dict[interaction] <= rank_threshold:
            num_within_threshold += 1
    return num_within_threshold / dict_size


def calculate_FPR(inter_to_ranks_dict: dict, rank_threshold: int) -> float:
    dict_size: int = len(inter_to_ranks_dict)
    num_within_threshold: int = 0
    for interaction in inter_to_ranks_dict:
        if inter_to_ranks_dict[interaction] > rank_threshold:
            num_within_threshold += 1
    return num_within_threshold / dict_size

#
# # FPR: the percentage of random interactions (not confirmed) above threshold
# def calculate_FPR(all_interactions: list, tru_interactions, rank_threshold: int) -> float:
#     dict_size: int = len(all_interactions)
#     num_within_threshold: int = 0
#
#     for interaction in all_interactions:
#         if inc in tru_interactions:
#             if inter_to_ranks_dict[interaction] <= rank_threshold:
#                 num_within_threshold += 1
#     return num_within_threshold / dict_size


def run(samples_file: str, results_file: str):
    tpr_per_thresh_per_sample: list = []
    fpr_per_thresh_per_sample: list = []
    positivity_rate_per_sample: list = []

    true_interactions, our_interactions = read_data(samples_file, results_file)

    for i in range(5):
        interactions_to_ranks: dict = {}
        sample_true_interactions = true_interactions[i]
        sample_predicted_interactions = our_interactions[i]
        positivity_rate_per_sample.append(len(sample_true_interactions)/len(sample_predicted_interactions))

        for rna_d_tuple in sample_true_interactions:
            rank: int = sample_predicted_interactions.index(rna_d_tuple)
            interactions_to_ranks[rna_d_tuple] = rank

        num_ranks: float = float(len(sample_predicted_interactions))
        thresholds = [round((x / 20) * num_ranks) for x in range(1, 21)]


        tpr_per_threshold = [calculate_TPR(interactions_to_ranks, x) for x in thresholds]
        tpr_per_thresh_per_sample.append(tpr_per_threshold)
        fpr_per_threshold = [calculate_FPR(interactions_to_ranks, x) for x in thresholds]
        fpr_per_thresh_per_sample.append(fpr_per_threshold)

        avg_positivity_rate: float = sum(positivity_rate_per_sample)/5

    print(tpr_per_thresh_per_sample)
    print(calculate_averages(tpr_per_thresh_per_sample))
    print(fpr_per_thresh_per_sample)
    print(calculate_averages(fpr_per_thresh_per_sample))
    print(avg_positivity_rate)


# given list of lists, condense to list of average of the first elements, second elems, etc. of each list
def calculate_averages(big_lst):
    build_averages_list = []

    if big_lst:
        for i in range(0, len(big_lst[0])):
            my_sum = 0.0
            for lst in big_lst:
                my_sum = my_sum + lst[i]

            avg = my_sum / len(big_lst)
            build_averages_list.append(avg)

    return build_averages_list


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])
    print("hello")
