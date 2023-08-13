import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
import numpy as np
import os

def read_data(path_file):
    with open(path_file, "r", encoding="utf-8") as f:
        extention = path_file.split('.', 1)[1].lower()

        if extention == 'csv':
            lines = [line.replace('\n', '').split(';') for line in f.readlines()]
        elif extention == 'txt':
            lines = [line.replace('\n', '').split('\t') for line in f.readlines()]
        else:
            error_message = 'Error! The file extension "' + extention + '" is not allowed!'
            return error_message

        len_lines = len(lines)

        reactions = []
        count_reaction = 0
        molecules = []
        generic_names = []
        normal_cytotoxicity = []
        cytotoxicity = []
        full_data_flag = 0

        for count_line, line in enumerate(lines):
            if count_line == 0:
                cell_name = line

            elif line[0] == '***':
                if count_line == 1:
                    reactions.append(cell_name)
                else:
                    count_reaction = -1

                    if full_data_flag == 0:
                        full_data_line = "enough data"
                    else:
                        full_data_line = "not enough data"

                    reactions.append([full_data_line, reaction_name, molecules, generic_names,
                                      np.array(normal_cytotoxicity), cytotoxicity])
                    molecules = []
                    generic_names = []
                    normal_cytotoxicity = []
                    cytotoxicity = []
                count_reaction = -1
            else:
                if count_reaction == 0:
                    reaction_name = line
                else:
                    if line[2].lower() == 'na' or line[3].lower() == 'na':
                        full_data_flag = 1
                    else:
                        for count_el, el in enumerate(line):
                            if count_el == 0:
                                molecules.append(el)
                            elif count_el == 1:
                                generic_names.append(el)
                            elif count_el == 2:
                                normal_cytotoxicity.append(float(el.replace(',', '.')))
                            else:
                                cytotoxicity.append(float(el.replace(',', '.')))

                if count_line == len_lines - 1:
                    if full_data_flag == 0:
                        full_data_line = "enough data"
                    else:
                        full_data_line = "not enough data"
                    reactions.append([full_data_line, reaction_name, molecules, generic_names,
                                      np.array(normal_cytotoxicity), cytotoxicity])

            count_reaction += 1

    return reactions


def fill_colors(data, category_colors, cytotoxity_scale):
    colors_list = []
    for el in data:
        el_color = 0
        for index, cytotoxity in enumerate(cytotoxity_scale):
            if el >= cytotoxity:
                el_color = category_colors[index]
            else:
                break

        colors_list.append(el_color)

    colors = np.array(colors_list)

    return colors


def cyt_chart(path_graph, full_data, reaction_name, molecules, generic_names, normal_cytotoxicity, biofactor, colors, formats, scale_coef):
    data_cum = normal_cytotoxicity.cumsum(axis=0)
    len_graph = data_cum[len(molecules) - 1]

    width = 25 * len_graph/scale_coef
    height = 2
    cyt_median = 1.2 * scale_coef * np.median(normal_cytotoxicity)

    fig, ax = plt.subplots(figsize=(width, height))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    count_decrease_cyt = 0
    cycle_len = len(normal_cytotoxicity)

    for i in range(cycle_len):
        widths = normal_cytotoxicity[i]
        starts = data_cum[i] - widths

        xcenters = starts + widths / 2

        ax.barh('cell', widths, left=starts, height=0.8, color=colors[i],
                edgecolor='black', linewidth=1)

        # color selection for reagent type
        if 'SM' == generic_names[i][:2]:
            text_color = 'darkblue'
        elif 'CT' == generic_names[i][:2]:
            text_color = 'black'
        elif 'R' == generic_names[i][:1]:
            text_color = 'maroon'
        elif 'S' == generic_names[i][:1]:
            text_color = 'maroon'
        elif 'P' == generic_names[i][:1]:
            text_color = 'red'
        elif 'BP' == generic_names[i][:2]:
            text_color = 'purple'
        else:
            text_color = 'black'

        # display of reagent names
        if normal_cytotoxicity[i] < cyt_median:
            flag_single_decrease_cyt = 0
            if count_decrease_cyt == 0:

                for k in range(i, cycle_len):
                    if normal_cytotoxicity[k] < cyt_median:
                        count_decrease_cyt += 1
                    else:
                        break

                if count_decrease_cyt == 1:
                    flag_single_decrease_cyt = 1

            height_text = -0.4 - 0.3 * count_decrease_cyt
            ax.annotate("", xy=(xcenters, -0.4), xytext=(xcenters, height_text), arrowprops=dict(arrowstyle="-"))

            if flag_single_decrease_cyt == 1:
                ax.text(xcenters, height_text, generic_names[i], ha='center', va='bottom', color=text_color,
                        fontsize=35)
                flag_single_decrease_cyt = 0
            else:
                ax.text(xcenters, height_text, generic_names[i], ha='left', va='bottom', color=text_color, fontsize=35)

            count_decrease_cyt -= 1

        else:
            height_text = -0.4
            ax.text(xcenters, height_text, generic_names[i], ha='center', va='bottom', color=text_color, fontsize=35)

    for spine in ax.spines:
        ax.spines[spine].set_visible(False)

    ax.text(0, 0, reaction_name + ' ', color='black', fontsize=50, horizontalalignment='right',
            verticalalignment='center')

    ax.text(0, 0.5, 'BF = ' + str(biofactor), color='black', fontsize=50, horizontalalignment='left',
            verticalalignment='top')

    if full_data == 'not enough data':
        ax.text(0, 1.0, 'not enough data!', color='red', fontsize=50, horizontalalignment='left',
                verticalalignment='top')

    for i in range(len(formats)):
        plt.savefig(path_graph[i] + '.' + formats[i], bbox_inches='tight')
    plt.clf()

    return fig, ax


def calc_cyt_metrics(generic_names, normal_cytotoxicity):
    cyt_sum_reagents = 0
    cyt_sum_products = 0
    cyt_sum_target_products = 0
    cyt_sum_const = 0

    for names, cyt in zip(generic_names, normal_cytotoxicity):

        if names[:2] == 'SM' or names[:1] == 'R':
            cyt_sum_reagents += cyt
        elif names[:1] == 'P' or names[:2] == 'BP':
            cyt_sum_products += cyt

            if names[:1] == 'P':
                cyt_sum_target_products += cyt
        else:
            cyt_sum_const += cyt

    CPi = cyt_sum_reagents + cyt_sum_const
    CPf = cyt_sum_products + cyt_sum_const
    CPf_rel = CPf - cyt_sum_target_products
    biofactor = (cyt_sum_products + cyt_sum_const) / (cyt_sum_reagents + cyt_sum_const)

    cyt_metrics = [round(biofactor, 2), round(CPi, 2), round(CPf, 2), round(CPf_rel, 2)]

    return cyt_metrics


def get_all_cytotoxity(reactions):
    cytotoxity_array = []

    for count, reaction in enumerate(reactions):
        if count == 0:
            continue
        else:
            for el in reaction[5]:
                if el not in cytotoxity_array:
                    cytotoxity_array.append(el)

    cytotoxity_array.sort()

    cytotox_percent = np.percentile(cytotoxity_array, 50)

    number_of_colors_1 = 0
    number_of_colors_2 = 0

    for el in cytotoxity_array:
        if el < cytotox_percent:
            number_of_colors_1 += 1
        else:
            number_of_colors_2 += 1

    return cytotoxity_array, number_of_colors_1, number_of_colors_2


def cyt_colormap(cytotoxity_scale, colors, cell, colormap, path_dirs, formats):
    if colormap == 'linear':
        cycle_len = len(cytotoxity_scale)

        fig, ax = plt.subplots(figsize=(np.sqrt(cycle_len), 1))

        ax.invert_yaxis()
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        for i in range(cycle_len):
            ax.barh('cell', 1, left=i, height=1, color=colors[i], linewidth=1)
            if i == 0 or i == cycle_len - 1 or i % 10 == 0:
                ax.text(i, 0.6, '|', ha='center', va='center_baseline', color='black')
                ax.text(i, 0.8, round(cytotoxity_scale[i], 1), ha='center', va='center_baseline', color='black')

        for spine in ax.spines:
            ax.spines[spine].set_visible(False)

        ax.text(len(cytotoxity_scale) / 2, -0.8, 'cell line ' + cell, color='black', fontsize=15,
                horizontalalignment='center',
                verticalalignment='top')
    else:
        cycle_len = len(cytotoxity_scale)

        fig, ax = plt.subplots(figsize=(cycle_len, 1))

        ax.invert_yaxis()
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        for i in range(cycle_len):
            ax.barh('cell', 1, left=i, height=1, color=colors[i],
                    edgecolor='black', linewidth=1)
            ax.text(i + 0.5, 0, round(cytotoxity_scale[i], 1), ha='center', va='center_baseline', color='black')

        for spine in ax.spines:
            ax.spines[spine].set_visible(False)

        ax.text(len(cytotoxity_scale) / 2, 0.6, 'cell line ' + cell, color='black', fontsize=15,
                horizontalalignment='center',
                verticalalignment='top')


    for i in range(len(formats)):
        path_graph = os.path.join(path_dirs[i], 'colormap.' + formats[i])
        plt.savefig(path_graph, bbox_inches='tight')
    plt.clf()

    return fig, ax


def choice_colormap(colormap, reactions):
    # parsing unique cytotoxicity values to build a color map
    cytotoxity_array, number_of_colors_1, number_of_colors_2 = get_all_cytotoxity(reactions)

    if colormap == 'percentile':
        # creating an array of colors to sample
        color_start = 0.0
        color_end_1 = 0.13
        color_end_2 = 0.3
        color_1 = plt.get_cmap('hsv')(np.linspace(color_start, color_end_1, number_of_colors_1))
        color_2 = plt.get_cmap('hsv')(np.linspace(color_end_1, color_end_2, number_of_colors_2))
        colors = np.vstack([color_1, color_2])
        cytotoxity_scale = cytotoxity_array

    if colormap == 'linear':
        number_of_colors = 100

        # creating an array of colors
        cytotoxity_scale_start = min(cytotoxity_array)
        cytotoxity_scale_end = max(cytotoxity_array)
        colors = plt.get_cmap('hsv')(np.linspace(0.0, 0.35, number_of_colors))
        cytotoxity_scale = np.linspace(cytotoxity_scale_start, cytotoxity_scale_end, number_of_colors)

    return colors, cytotoxity_scale


def find_top_combinations(cyt_metrics, metric, top_size):
    if metric == 'BF':
        index_metric = 0
    elif metric == 'CPi':
        index_metric = 1
    elif metric == 'CPf':
        index_metric = 2
    elif metric == 'CPf_rel':
        index_metric = 3

    names, metric_vector = zip(*[(metrics[0], metrics[1][index_metric]) for metrics in cyt_metrics.items()])

    indices_metrics = [i for i in range(len(metric_vector))]

    sort_indices_metrics = [index for metric_num, index in sorted(zip(metric_vector, indices_metrics))]

    top_combinations = []
    for i in range(top_size):
        index_sort = sort_indices_metrics[i]
        top_combinations.append([names[index_sort], metric_vector[index_sort]])

    return top_combinations


def calc_scaling_coef(reactions):
    NC_len_max = 0
    NC_len_min = 10000000

    for count, el in enumerate(reactions):
        if count == 0:
            continue
        else:
            NC_len_el = np.sum(el[4])

            if NC_len_el > NC_len_max:
                NC_len_max = NC_len_el

            if NC_len_el < NC_len_min:
                NC_len_min = NC_len_el

    scale_coef = np.sqrt(NC_len_min * NC_len_max)

    return scale_coef
