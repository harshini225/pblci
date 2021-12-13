import random
import sys


def read_big_interactions_file(filename: str):
    big_interactions_dict: dict = {}
    with open(filename, "r") as reader:
        for raw_line in reader:
            line: str = raw_line.strip()
            split: list = line.split(";")
            key: str = split[0]
            values: list = [x for x in split[1].split(",") if x]
            big_interactions_dict[key] = values
    return big_interactions_dict


def take_rand_sample(interactions_dict: dict, sample_size: int) -> dict:
    sample_dict: dict = {}
    rna_key_list: list = list(interactions_dict.keys())
    rna_sample: list = random.sample(rna_key_list, sample_size)
    for rna in rna_sample:
        disease_values: list = interactions_dict[rna]
        sample_dict[rna] = disease_values
    return sample_dict


def take_samples(interactions_dict: dict):
    interactions_left: dict = interactions_dict.copy()
    sample_size: int = 22
    samples: list[dict] = []
    for i in range(5):
        rand_sample: dict = take_rand_sample(interactions_left, sample_size)
        for rna in rand_sample.keys():
            del interactions_left[rna]
        samples.append(rand_sample)
        # last sample should have size = 21
        if i == 3:
            sample_size = 21
    return samples


def write_samples(filename: str, samples: list[dict]):
    with open(filename, "w") as writer:
        for i in range(5):
            sample = samples[i]
            # signal new sample with #
            writer.write("#" + str(i) + "\n")
            for rna in sample:
                values: list = sample[rna]
                values_as_str = ",".join(values)
                writer.write(rna + ";" + values_as_str + "\n")

    # for i in range(5):
    #     sample_to_check = all_samples[i]
    #     others = all_samples.copy()
    #     others.remove(sample_to_check)
    #     for something in sample_to_check:
    #         for sample in others:
    #             if something in sample:
    #                 print(str(something) + " is in " + str(sample))

    return samples


def run(interactions_file: str, samples_file: str):
    big_interactions_dict: dict = read_big_interactions_file(interactions_file)
    samples = take_samples(big_interactions_dict)
    write_samples(samples_file, samples)
    return samples


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])
