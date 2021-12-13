import sys
import split_samples
import run_pblci


# file that performs 5-fold cross validation


# runs the results for one sample
def run_sample(test: dict, sample_dict: dict):
    rna_to_dis_to_score: dict[tuple, float] = {}
    # diseases_with_dupes: list = []
    # for disease_values in sample_dict.values():
    #     diseases_with_dupes = diseases_with_dupes + disease_values
    # all_diseases: set = set(diseases_with_dupes)
    #
    # my_network = make_big_matrix.read_matrices()
    sorted_by_score: list[tuple] = run_pblci.run_sample_from_dict(test, sample_dict)
    return sorted_by_score


# performs the cross-validation
def five_fold_cv(interactions_file: str, samples_file: str):
    samples: list[dict] = split_samples.run(interactions_file, samples_file)
    samples_to_scores = {}

    for i in range(5):
        test_sample: dict = samples[i]  # left out
        others: list[dict] = samples.copy()
        others.remove(test_sample)
        training_sample: dict = {}
        all_sample: dict = {}

        for sample in others:
            training_sample.update(sample)
        for s in samples:
            all_sample.update(s)

        # samples_dict =
        all_rankings: list[tuple] = run_sample(test_sample, all_sample)
        novel_rankings: list[tuple] = []
        for rna, doid in all_rankings:
            # if rna in training_sample and training_sample[rna] == doid:
            #     novel_rankings.append((rna, doid))
            if not (rna in training_sample and doid in training_sample[rna]):
                novel_rankings.append((rna, doid))
        samples_to_scores[i] = novel_rankings

    return samples_to_scores


def write_samples_to_scores(filename: str, sample_to_scores: dict):
    with open(filename, "w") as writer:
        for i in range(5):
            writer.write("#" + str(i) + "\n")
            rankings: list[tuple] = sample_to_scores[i]
            for rank in rankings:
                line_to_write: str = rank[0] + "," + rank[1] + "\n"
                writer.write(line_to_write)


# method to run file
def run(interactions_file: str, samples_file: str, scores_file: str):
    samples_to_scores = five_fold_cv(interactions_file, samples_file)
    write_samples_to_scores(scores_file, samples_to_scores)


# [interactions files] [samples file to write to] [path file to write to]
if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3])
