from itertools import product

roles_of_molecules = ['starting materials', 'catalysts', 'reagents', 'solvents', 'products', 'byproducts']


def parsing_and_preparation_data(path_file):
    with open(path_file, "r", encoding="utf-8") as f:
        extension = path_file.split('.', 1)[1].lower()

        if extension == 'csv':
            lines = [line.replace('\n', '').split(';') for line in f.readlines()]
        elif extension == 'txt':
            lines = [line.replace('\n', '').split('\t') for line in f.readlines()]
        else:
            error_message = 'Error! The file extension "' + extension + '" is not allowed!'
            return error_message

        molecules = []
        generic_names = []
        cytotoxicity = []
        normal_cytotoxicity = []

        labels = []
        products_labels = []
        product_flag = 0

        cell = lines[0][1]

        # parsing variables
        if len(lines[1]) == 1:
            variables_names = []
        elif len(lines[1][1]) == 0:
            variables_names = []
        else:
            variables_names = lines[1][1].replace(' ', '').split(',')

        if len(lines[2]) == 1:
            product_variables_names = []
        elif len(lines[2][1]) == 0:
            product_variables_names = []
        else:
            product_variables_names = lines[2][1].replace(' ', '').split(',')

        # parsing data
        for count_line, line in enumerate(lines):
            name_string = line[0].lower()

            if count_line == 0 or count_line == 1 or count_line == 2 or count_line == 3 or count_line == 4:
                continue
            elif len(line) <= 1:
                continue
            elif name_string in roles_of_molecules:
                if name_string == 'products':
                    product_flag = 1
                continue
            else:
                # creating a list of labels. Needed for sorting
                label = get_original_label(line[1])

                labels = add_label(label, labels)
                if product_flag:
                    products_labels = add_label(label, products_labels)

                molecules.append(line[0])
                generic_names.append(line[1])
                cytotoxicity.append(line[4].replace(',', '.'))

                # calculation of normalized cytotoxicity
                if line[2].lower() == 'na' or line[3].lower() == 'na' or line[4].lower() == 'na':
                    NC = 'NA'
                else:
                    NC = 1000*float(line[3].replace(',', '.'))/(
                            float(line[2].replace(',', '.'))*float(line[4].replace(',', '.')))
                normal_cytotoxicity.append(NC)

        # creating a dictionary of variables
        variables_dict = {count_variable: variable_name for count_variable, variable_name in enumerate(variables_names)}

        # sorting data
        generic_names_sort, sort_indices = sort_names(generic_names, labels)
        molecules_sort, normal_cytotoxicity_sort, cytotoxicity_sort = sort_data_by_indices(molecules,
                                                                                           normal_cytotoxicity,
                                                                                           cytotoxicity, sort_indices)

        raw_data = [cell, molecules_sort, generic_names_sort, normal_cytotoxicity_sort, cytotoxicity_sort]

    return raw_data, variables_names, product_variables_names, variables_dict, labels, products_labels


def add_label(label, labels):
    if len(labels) == 0:
        labels.append(label)
    else:
        if label in labels:
            pass
        else:
            labels.append(label)

    return labels


def get_original_label(label_long):
    label = label_long.split('-')[0]
    original_label = ''

    for el in label:
        if el.isdigit():
            continue
        else:
            original_label += el

    return original_label


def get_index_by_label(label, names):
    for i in range(len(names)):
        if label == names[i]:
            index = i
        else:
            continue

    return index


def get_indices_products(names, product_labels):
    product_indices = []

    for label in product_labels:
        for i in range(len(names)):
            original_name = get_original_label(names[i][0])
            if len(names[i]) > 1 and label == original_name:
                product_indices.append(i)
            else:
                continue

    return product_indices


def add_molecule(label, molecules, generic_names, normal_cytotoxicity, cytotoxicity, raw_data):
    index = get_index_by_label(label, raw_data[2])

    labels = raw_data[2][index].split('-')

    molecules.append(raw_data[1][index])
    generic_names.append(labels[0])
    normal_cytotoxicity.append(raw_data[3][index])
    cytotoxicity.append(raw_data[4][index])

    return molecules, generic_names, normal_cytotoxicity, cytotoxicity


def get_product_combinations(variables_list, inverse_variables_dict, product_variables_names, letters_dict):
    numbers_variables = [inverse_variables_dict[variable] for variable in product_variables_names]

    variables_to_num = []

    for variable in variables_list:
        variables_to_num.append(letters_dict[variable])

    products_var_num = []
    for number in numbers_variables:
        products_var_num.append(['0', str(variables_to_num[number])])

    for i in range(len(products_var_num) - 1):
        if i == 0:
            products_combinations = ['-'.join(j) for j in list(product(products_var_num[i], products_var_num[i + 1]))]
        else:
            products_combinations = ['-'.join(j) for j in list(product(products_combinations, products_var_num[i + 1]))]

    products_combinations = ['-' + product for product in products_combinations]

    return products_combinations


def sort_names(names, labels):
    indices_labels = [i for i in range(len(names))]
    new_names, new_indices_labels = zip(*[(name, index) for name, index in sorted(zip(names, indices_labels))])

    sort_names = []
    sort_indices = []

    for label in labels:
        for i in range(len(new_names)):
            if label == get_original_label(new_names[i]):
                sort_names.append(new_names[i])
                sort_indices.append(new_indices_labels[i])

    return sort_names, sort_indices


def sort_data_by_indices(molecules, normal_cytotoxicity, cytotoxicity, sort_indices):
    molecules_sort = [molecules[index] for index in sort_indices]
    normal_cytotoxicity_sort = [normal_cytotoxicity[index] for index in sort_indices]
    cytotoxicity_sort = [cytotoxicity[index] for index in sort_indices]

    return molecules_sort, normal_cytotoxicity_sort, cytotoxicity_sort


def generate_outtable(path_input, raw_data, combinations, labels, products_labels, variables_dict, product_variables_names, letters_dict, inverse_variables_dict):
    combinations_list = [line.split('-') for line in combinations]

    roles_molecules = [line.split('-') for line in raw_data[2]]

    number_combinations = len(combinations)

    product_indices = get_indices_products(roles_molecules, products_labels)

    full_products = restore_full_products(product_indices, raw_data)

    path_outfile = path_input.split('.')[0] + '_comb.txt'

    with open(path_outfile, "w", encoding="utf-8") as out_file:
        print(raw_data[0], end='\n', file=out_file)
        print('***', end='\n', file=out_file)

        for i in range(number_combinations):
            molecules = []
            generic_names = []
            normal_cytotoxicity = []
            cytotoxicity = []

            # add constant molecules
            for role in roles_molecules:
                if len(role) == 1:
                    label = role[0]
                    add_molecule(label, molecules, generic_names, normal_cytotoxicity, cytotoxicity, raw_data)
                else:
                    continue

            # add variable reagents
            if len(variables_dict) != 0:
                for count_letter, letter in enumerate(combinations_list[i]):
                    label = variables_dict[count_letter] + '-' + str(letters_dict[letter])
                    full_label = find_full_label(label, raw_data)
                    add_molecule(full_label, molecules, generic_names, normal_cytotoxicity, cytotoxicity, raw_data)

            # search for possible product combinations
            if len(product_variables_names) != 0:
                product_combinations = get_product_combinations(combinations_list[i], inverse_variables_dict, product_variables_names, letters_dict)

                # add variable products
                for combination in product_combinations:
                    for product in full_products:
                        if combination in product[0]:
                            label = raw_data[2][product[1]]
                            add_molecule(label, molecules, generic_names, normal_cytotoxicity, cytotoxicity, raw_data)

            # sorting molecules
            generic_names_sort, sort_indices = sort_names(generic_names, labels)
            molecules_sort, normal_cytotoxicity_sort, cytotoxicity_sort = sort_data_by_indices(molecules,
                                                                                               normal_cytotoxicity,
                                                                                               cytotoxicity,
                                                                                               sort_indices)

            # writing data to a file
            print(combinations[i], end='\n', file=out_file)
            for j in range(len(molecules_sort)):
                print(molecules_sort[j], generic_names_sort[j], normal_cytotoxicity_sort[j],
                      cytotoxicity_sort[j], sep='\t', end='\n', file=out_file)
            print('***', end='\n', file=out_file)

    return path_outfile


def generate_combinations(raw_data, variables_names, inverse_letters_dict):
    if len(variables_names) == 0:
        combinations = inverse_letters_dict[1]
    else:
        combinations = []
        combinations_copy = []
        flag_start = 0
        flag_start_label = 0
        for count_line, line in enumerate(raw_data[2]):
            label = line.split('-')
            label_origin = label[0]

            if label_origin in variables_names:
                if flag_start == 0:
                    label_copy = label_origin
                    flag_start = 1

                if label_origin != label_copy:
                    label_copy = label_origin

                    if flag_start_label != 0:
                        combinations = combinations_copy
                        combinations_copy = []
                    flag_start_label = 1

                if len(label) <= 2:
                    if flag_start_label == 0:
                        combinations.append(label[1])
                    else:
                        for combination in combinations:
                            new_comb = combination + '-' + label[1]
                            combinations_copy.append(new_comb)
                else:
                    numbers_label = label[2].split(',')

                    for number in numbers_label:
                        indices = find_indices_of_element(combinations, number)

                        for index in indices:
                            new_comb = combinations[index] + '-' + label[1]
                            combinations_copy.append(new_comb)

        if len(combinations_copy) == 0:
            combinations_copy = combinations

        combinations = []

        # creating a list with variable values
        for combination in combinations_copy:
            combination_numbers = combination.split('-')

            combination_letters = [inverse_letters_dict[int(number)] for number in combination_numbers]
            combination_result = '-'.join(combination_letters)
            combinations.append(combination_result)

    role_dict = {}
    flag = ''
    for role in raw_data[2]:
        index = raw_data[2].index(role)
        true_role = role.split('-')[0]
        if true_role in variables_names:
            if flag == '':
                flag = true_role
                role_for_dict = {}
            if flag != true_role:
                role_dict[flag] = role_for_dict
                role_for_dict = {}
                flag = true_role
            true_letter = inverse_letters_dict[int(role.split('-')[1])]
            true_role_dict = {}
            true_role_dict[role] = true_letter
            role_for_dict[raw_data[1][index]] = true_role_dict
        else:
            pass
    return sorted(combinations), role_dict


def find_indices_of_element(combinations, number):
    indices = []

    for count_line, el in enumerate(combinations):
        el_number = el.split('-')[-1]
        if el_number == number:
            indices.append(count_line)

    return indices


def find_full_label(label, raw_data):
    for el in raw_data[2]:
        if label in el:
            full_label = el
            break

    return full_label


def restore_full_products(product_indices, raw_data):
    full_products = []

    for index in product_indices:
        variable_numbers = []
        label = raw_data[2][index]
        label_list = label.split('-')

        for count, el in enumerate(label_list):
            if count == 0:
                continue
            else:
                numbers = el.split('/')
                variable_numbers.append(numbers)

        for i in range(len(variable_numbers) - 1):
            if i == 0:
                products_combinations = ['-'.join(j) for j in
                                         list(product(variable_numbers[i], variable_numbers[i + 1]))]
            else:
                products_combinations = ['-'.join(j) for j in
                                         list(product(products_combinations, variable_numbers[i + 1]))]

        for combination in products_combinations:
            new_label = label_list[0] + '-' + combination
            full_products.append([new_label, index])

    return full_products

