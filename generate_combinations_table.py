import cyt_combinations as cytcomb

letters_dict = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7,
           'H': 8, 'I': 9, 'J': 10, 'K': 11, 'L': 12, 'M': 13, 'N': 14,
           'O': 15, 'P': 16, 'Q': 17, 'R': 18, 'S': 19, 'T': 20, 'U': 21,
           'V': 22, 'W': 23, 'X': 24, 'Y': 25, 'Z': 26}

inverse_letters_dict = {number: letter for letter, number in letters_dict.items()}


def generate_table(path_input):

    raw_data, variables_names, product_variables_names, variables_dict, labels, products_labels, smiles_dict, ld50_dict \
    = cytcomb.parsing_and_preparation_data(path_input)

    inverse_variables_dict = {variable: number for number, variable in variables_dict.items()}

    combinations, legend = cytcomb.generate_combinations(raw_data, variables_names, inverse_letters_dict)

    path_out = cytcomb.generate_outtable(path_input, raw_data, combinations, labels, products_labels, variables_dict, product_variables_names, letters_dict, inverse_variables_dict,)
    with open(path_out.split('.')[0] + 'out' + '.csv', "w", encoding="utf-8") as out_file:
        text = ''
        for variable in variables_names:
            text += variable + '-'
        text = text[:-1]
        first_string = print(f'Reaction name:, {text}', file=out_file)
        for name, values in legend.items():
            print(name, file=out_file)
            for name_val, val in values.items():
                print(f'{name_val},{list(val.keys())[0]},{list(val.values())[0]},{smiles_dict[name_val]},{ld50_dict[name_vale]}', file=out_file)

    number_of_combinations = len(combinations)

    return path_out, number_of_combinations
