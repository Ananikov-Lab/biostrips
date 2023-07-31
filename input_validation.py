import cyt_combinations as cytcomb

roles_of_molecules = ['starting materials', 'catalysts', 'reagents', 'solvents', 'products', 'byproducts']


def data_validation(path_file):
    try:
        with open(path_file, "r", encoding="utf-8") as f:
            extension = path_file.split('.', 1)[1].lower()

            if extension == 'csv':
                lines = [line.replace('\n', '').split(';') for line in f.readlines()]
            elif extension == 'txt':
                lines = [line.replace('\n', '').split('\t') for line in f.readlines()]
            else:
                error_message = 'Error! The file extension "' + extension + '" is not allowed!'
                return error_message
    except:
        error_message = 'Error! Make sure the file is encoded "utf-8"'
        return error_message

    error_message = 'Error! '

    labels = []
    general_labels = []
    product_flag = 0
    generic_names = {}

    for count_line, line in enumerate(lines):
        if count_line == 0:
            line_0 = line[0].lower().replace('\ufeff', '')
            if line_0 != 'cell':
                error_message += 'Add the string "cell" to the line 1!'
                return error_message
        elif count_line == 1:
            if len(line) > 1:
                variables_names = lines[1][1].replace(' ', '').split(',')
            else:
                variables_names = ''

            if line[0].lower() != 'variables':
                error_message += 'Add the string "variables" to the line 2!'
                return error_message
        elif count_line == 2:
            if len(line) > 1:
                product_variables_names = lines[2][1].replace(' ', '').split(',')
                number_product_variables = len(product_variables_names)
            else:
                product_variables_names = ''
                number_product_variables = 0

            if line[0].lower() != 'product variables':
                error_message += 'Add the string "product variables" to the line 3!'
                return error_message
        elif count_line == 3:  # table description
            continue
        elif line[0].lower() in roles_of_molecules:
            if line[0].lower() == 'products':
                product_flag = 1
                line_product_start = count_line + 1
            continue
        elif len(line) <= 1:
            continue
        else:
            repetition = get_repeating_element_index(line[1], generic_names)
            if repetition == 'new':
                generic_names[line[1]] = count_line + 1
            else:
                message = 'The substance notations in lines ' + str(repetition) + ' and ' + str(
                    count_line) + ' are the same!'
                error_message += message
                return error_message

            general_label = cytcomb.get_original_label(line[1])

            general_labels = cytcomb.add_label(general_label, general_labels)

            label_list = line[1].split('-')
            label = label_list[0]

            # comparison of the numbers of specified variables and variables in the file
            if product_flag and len(label_list) > 1:
                if len(label_list) != number_product_variables + 1:
                    message = 'Wrong product notation in string ' + str(count_line + 1) + '!'
                    error_message += message
                    return error_message
            else:
                if len(label_list) > 1:
                    if label in variables_names:
                        labels = cytcomb.add_label(label, labels)
                    else:
                        message = 'Wrong product notation in string ' + str(count_line + 1) + '!'
                        error_message += message
                        return error_message

            # checking if a string can be converted to float
            numerical_data = {'Mr': line[2], 'Mass': line[3], 'CC50': line[4]}

            for data in numerical_data.items():
                number_line = data[1].replace(',', '.')

                if data[1].lower() == "na":
                    continue
                else:
                    try:
                        number = float(number_line)
                        if number <= 0:
                            message = 'Element ' + str(data[0]) + ' in line ' + str(count_line + 1) + ' must be a positive number!'
                            error_message += message
                            return error_message
                    except:
                        message = 'Element ' + str(data[0]) + ' in line ' + str(count_line + 1) + ' is not a number!'
                        error_message += message
                        return error_message

    if len(variables_names) == 1 and len(variables_names[0]) == 0 or len(variables_names) == 0:
        pass
    else:
        if variables_names != labels:
            error_message += 'The variables specified in the string "variables" do not match the product variables in the file!'
            return error_message

    if product_flag == 0:
        error_message += 'Add the string "products" before the reaction products!'
        return error_message

    # validation of variable notation
    generic_names_sort = sort_names_dict(generic_names, general_labels)
    flag_start_test = 0
    previous_gen_label = ''
    number_variables_label = []
    number_variables_label_copy = []

    for count, el in enumerate(generic_names_sort.items()):
        el_split = el[0].split('-')
        if count == 0:
            el_name_copy = el_split[0]

            if len(el_split) > 1:
                number_variables_label.append(el_split[1])
            continue

        if el[1] >= line_product_start:
            break

        if el_split[0] != el_name_copy:
            flag_start_test = 1
            previous_gen_label = el_name_copy
            number_variables_label_copy = number_variables_label
            number_variables_label = []

        if len(el_split) > 1:
            if el_split[1] not in number_variables_label:
                number_variables_label.append(el_split[1])

        if len(el_split) > 2:
            numbers_labels = el_split[-1].split(',')
            for number in numbers_labels:
                if flag_start_test == 1 and number not in number_variables_label_copy:
                    message = 'Wrong variable notation in line ' + str(el[1])
                    message += '! Reagent ' + previous_gen_label + '-' + str(number) + ' is not found!'
                    error_message += message
                    return error_message

        el_name_copy = el_split[0]

    message_success = 'SUCCESS'
    return message_success


def get_repeating_element_index(label, names):
    repetition_found = 0

    for el in names.items():
        if label == el[0]:
            index = el[1] + 1
            repetition_found = 1

    if repetition_found == 0:
        return 'new'
    else:
        return index


def sort_names_dict(names, labels):
    new_names = dict(sorted(names.items()))
    sort_names = {}

    for label in labels:
        for name in new_names.items():
            if label == cytcomb.get_original_label(name[0]):
                sort_names[name[0]] = name[1]

    return sort_names
