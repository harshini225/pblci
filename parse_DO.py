import re
import sys


def extract_disease_name(line: str) -> str:
    if line.startswith("# class"):
        disease_name_with_paren: str = re.findall("\\(.*\\)", line)[0]
        return disease_name_with_paren[1:-1]
    elif line.startswith("annotationassertion"):
        disease_name_with_paren: str = re.findall("\".*\"", line)[0]
        return disease_name_with_paren[1:-1]


def extract_id(line: str) -> str:
    obo_id: str = re.findall("obo:\\w*", line)[0]
    formatted_id = "DOID:" + obo_id[9:]
    return formatted_id


def find_do_ids(do_file: str, disease_list: set) -> dict:
    disease_dict: dict = {}
    with open(do_file, "r") as f:
        for raw_line in f:
            line: str = raw_line.strip().lower()
            if line.startswith("# class") or ("hasexactsynonym" in line) or \
                    ("hasrelatedsynonym" in line):
                disease_name: str = extract_disease_name(line)
                if disease_name in disease_list:
                    disease_id: str = extract_id(line)
                    if disease_id in disease_dict:
                        disease_dict[disease_id].append(disease_name)
                    else:
                        disease_dict[disease_id] = [disease_name]
    return disease_dict


def write_do_ids(filename: str, do_dict: dict):
    with open(filename, "w") as writer:
        writer.write(",".join(do_dict.keys()))


def read_disease_list(filename: str):
    disease_list: list = []
    with open(filename, "r") as f:
        for raw_line in f:
            line: str = raw_line.strip()
            disease_list.append(line)
    return disease_list


def write_do_to_disease(do_to_disease: dict, filename: str):
    with open(filename, "w") as f:
        for k in do_to_disease.keys():
            line_to_write: str = k
            for disease in do_to_disease[k]:
                line_to_write += "," + disease
            f.write(line_to_write + "\n")


def run(do_file: str, disease_):
    disease_list = read_disease_list(sys.argv[2])
    do_dict = find_do_ids(sys.argv[1], disease_list)
    write_do_ids(sys.argv[3], do_dict)
    write_do_to_disease(do_dict, sys.argv[4])


if __name__ == "__main__":
    # [DO file] [file with diseases] [file to write DO ids to] [file to write
    # dict to]
    run()
