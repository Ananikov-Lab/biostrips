#!/home/nkolomoets/miniconda3/envs/biostrips/bin/python
from model import get_ld50
from flask import Flask, render_template, request, flash, redirect, url_for, send_from_directory
from werkzeug.utils import secure_filename
from zipfile import ZipFile
import os
import generate_chart as gench
import generate_combinations_table as gentab
import input_validation as inpval
from flask_wtf import FlaskForm
from wtforms import StringField, FormField, FieldList, SelectField, Form, DecimalField
from wtforms.validators import DataRequired, Optional, NumberRange
from flask_bootstrap import Bootstrap
import logging
from logging.handlers import RotatingFileHandler
import pickle as pkl
import uuid

path_data = 'data'
path_meta = 'metadata'
path_results = 'results'
path_examples = 'examples'
path_new_chart_json = os.path.join(path_meta, 'new_chart.json')
path_figures = os.path.join('static', 'figures')

ALLOWED_EXTENSIONS = ['txt', 'csv']

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = path_data

bootstrap = Bootstrap(app)

with open(os.path.join(path_meta, 'secret_key.txt'), 'r') as f:
    app.config['SECRET_KEY'] = f.read()

menu = [{"name": "Main", "url": "/"},
        {"name": "Intro to bio-Strips", "url": '/intro_to_biostrips'},
        {"name": "Build Charts", "url": "/create_chart"},
        {"name": "Previous runs", "url": "/send_check"},
        {"name": "Manual", "url": "/manual"},
        {"name": "About", "url": "/about"}]

colormap_list = [{'name': 'percentile'},
                 {'name': 'linear'}]

cyt_potentials_list = [{'name': 'BF'},
                       {'name': 'CPi'},
                       {'name': 'CPf'},
                       {'name': 'CPf_rel'}]

mesure = ['mol/kg or mol/L', 'g/kg or g/L']


class ReagentLineForm(Form):
    reagent_name = StringField("Reagent name: ", validators=[DataRequired()], description="Reagent name")
    reagent_role = StringField("Reagent role: ", validators=[DataRequired()], description="Reagent role")
    molar_mass = DecimalField("Molar mass: ", places=3,
                              validators=[DataRequired(), NumberRange(0.001)],
                              description="Molar mass")
    mass = DecimalField("Mass: ", places=10,
                        validators=[DataRequired(), NumberRange(0.0000000001)],
                        description="Mass")
    cc50 = StringField("CC50: ",
                        validators=[DataRequired()],
                        description="CC50")


class OneChartForm(FlaskForm):
    filename = StringField("Enter filename: ", validators=[DataRequired()], description="Filename")
    cell_name = StringField("Enter cell name: ", validators=[DataRequired()], description="Cell name")
    reagents_info = FieldList(FormField(ReagentLineForm), min_entries=1)
    products_info = FieldList(FormField(ReagentLineForm), min_entries=1)
    colormap = SelectField("Choose colormap: ")
    cyt_potential = SelectField("Choose cytotoxic potential: ")
    variables = StringField("Enter variables separated by commas: ", description="Variables")
    products_variables = StringField("Enter product variables separated by commas: ", description="Products variables")
    type_mesure = SelectField("Select variable dimension: ")

class CheckFile(FlaskForm):
    filename = StringField("Enter access code: ", description="Access code")


@app.route('/')
@app.route('/home')
def main():
    return render_template('main.html', menu=menu)


@app.route('/manual')
def manual():
    return render_template('manual.html', menu=menu)


@app.route('/intro_to_biostrips')
def intro():
    return render_template('intro_to_bio.html', menu=menu)


@app.route('/about')
def about():
    return render_template('About.html', menu=menu)


@app.route('/send_check', methods=['POST', 'GET'])
def check_file_in_system():
    form = CheckFile()
    if 'send_check' in request.form:
        filename = form.filename.data
        filename = filename + '.txt'
        if filename in os.listdir(path_data):
            return redirect(url_for('output_file', filename=filename.rsplit(".")[0]), 302)
        else:
            return render_template('send_check.html', menu=menu, form=form, text='This file was not found')
    elif 'new_exp' in request.form:
        return redirect(url_for('create_chart'), 302)
    return render_template('send_check.html', menu=menu, form=form, text='')


@app.route("/output_file/<filename>", methods=['POST', 'GET'])
def output_file(filename):
    with open(f'static/figures/{filename.rsplit(".")[0]}/top_combinations.txt', 'rt') as f:
        lines = f.readlines()
        filename = lines[0].strip()
        colormap = lines[1].strip()
        cyt_potential = lines[2].strip()
        number_of_combs = lines[3].strip()
    with open(f'static/figures/{filename.rsplit(".")[0]}/file_info.pkl', 'rb') as f:
        top_combinations = pkl.load(f)
    return render_template('output_file.html', menu=menu, number=number_of_combs, colormap=colormap, filename=filename,
                           top_combinations=top_combinations, cyt_potential=cyt_potential)


@app.route('/create_chart', methods=['POST', 'GET'])
def create_chart():
    filename = str(uuid.uuid4())
    form = OneChartForm()
    form.colormap.choices = [(el["name"], el["name"]) for el in colormap_list]
    form.cyt_potential.choices = [(el["name"], el["name"]) for el in cyt_potentials_list]
    form.type_mesure.choices = [(el, el) for el in mesure]

    if 'send_create' in request.form:
        meta_filename = filename + '.txt'
        if meta_filename in os.listdir(path_data):
            flash('Sorry, this name is taken. Try another filename')
            return redirect(url_for('create_chart'), 302)
        cell_name = form.cell_name.data
        reagents_info = form.reagents_info.data
        products_info = form.products_info.data
        variables = form.variables.data
        products_variables = form.products_variables.data
        mesure_variable = request.form.get('type_mesure').replace('"', '')
        save_chart_data(filename, cell_name, reagents_info, products_info, variables, products_variables)
        colormap_name = request.form.get('colormap').replace('"', '')
        cyt_potential_name = request.form.get('cyt_potential').replace('"', '')
        file_info = {'title': meta_filename, 'colormap': colormap_name,
                     'cyt_potential': cyt_potential_name,
                     'mesure_type': mesure_variable}
        
        path_table, number_of_combinations = calc_combinations(file_info)
        top_combinations = make_chart(file_info, path_table)

        with open(f'static/figures/{filename}/top_combinations.txt', 'wt') as f:
            text = f'{filename}\n{colormap_name}\n{cyt_potential_name}\n{number_of_combinations}\n'
            f.write(text)
        with open(f'static/figures/{filename}/file_info.pkl', 'wb') as f:
            pkl.dump(top_combinations, f)
        return redirect(url_for('output_file', filename=filename.rsplit(".")[0]), 302)
    elif 'send_upload' in request.form:
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            upload_filename = secure_filename(file.filename)
            extension_upload = upload_filename.rsplit('.')[-1]
            filename = filename + '.' + extension_upload
            if filename in os.listdir(path_data):
                flash('Sorry, this name is taken. Try another filename')
                return redirect(url_for('create_chart'))
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            colormap_name = request.form.get('colormap').replace('"', '')
            cyt_potential_name = request.form.get('cyt_potential').replace('"', '')
            mesure_variable = request.form.get('type_mesure').replace('"', '')
            file_info = {'title': filename, 'colormap': colormap_name,
                         'cyt_potential': cyt_potential_name,
                         'mesure_type': mesure_variable}
            path_table, number_of_combinations = calc_combinations(file_info)
            top_combinations = make_chart(file_info, path_table)
            with open(f'static/figures/{filename.rsplit(".")[0]}/top_combinations.txt', 'wt') as f:
                text = f'{filename.rsplit(".")[0]}\n{colormap_name}\n{cyt_potential_name}\n{number_of_combinations}'
                f.write(text)

            with open(f'static/figures/{filename.rsplit(".")[0]}/file_info.pkl', 'wb') as f:
                pkl.dump(top_combinations, f)
            return redirect(url_for('output_file', filename=filename.rsplit(".")[0]), 302)
        else:
            extensions = 'Invalid file extension! The file must have one of the following extensions: '
            len_extensions = len(ALLOWED_EXTENSIONS)
            for count, extension in enumerate(ALLOWED_EXTENSIONS):
                if count == len_extensions - 1:
                    extensions += extension
                else:
                    extensions += extension + ','

            flash(extensions)
    return render_template('create_chart.html', menu=menu, colormap_list=colormap_list,
                           cyt_potentials_list=cyt_potentials_list, form=form, type_mesure=mesure)


@app.route('/save_chart_data')
def save_chart_data(filename, cell_name, reagents_info, products_info, variables, products_variables):
    path_file = os.path.join(path_data, filename + '.txt')

    with open(path_file, "w", encoding='utf-8') as out_file:
        print("Cell", cell_name, sep='\t', end='\n', file=out_file)
        print("Variables", variables, sep='\t', end='\n', file=out_file)
        print("Product variables", products_variables, sep='\t', end='\n', file=out_file)
        print("Samples", "Abbreviation", "Mr, g*mol-1", "Mass, g", "CC50, mM", sep='\t', end='\n', file=out_file)
        print("Starting materials", end='\t', file=out_file)
        print("", end='\t', file=out_file)
        print("", end='\t', file=out_file)
        print("", end='\t', file=out_file)
        print("", end='\n', file=out_file)

        flag = "S"
        dict_output_reag = {
            "C": "Catalysts",
            "R": "Reagents",
            "S": "Solvents"
        }
        for el in reagents_info:
            if el["reagent_role"][0] != flag:
                flag = el["reagent_role"][0]
                print(dict_output_reag[flag], end='\t', file=out_file)
                print("", end='\t', file=out_file)
                print("", end='\t', file=out_file)
                print("", end='\t', file=out_file)
                print("", end='\n', file=out_file)
            print(el["reagent_name"], end='\t', file=out_file)
            print(el["reagent_role"], end='\t', file=out_file)
            print(el["molar_mass"], end='\t', file=out_file)
            print(el["mass"], end='\t', file=out_file)
            print(el["cc50"], end='\n', file=out_file)

        print("Products", end='\t', file=out_file)
        print("", end='\t', file=out_file)
        print("", end='\t', file=out_file)
        print("", end='\t', file=out_file)
        print("", end='\n', file=out_file)

        flag = "P"
        for el in products_info:
            if el["reagent_role"][0] != flag:
                flag = el["reagent_role"][0]
                print("Byproducts", end='\t', file=out_file)
                print("", end='\t', file=out_file)
                print("", end='\t', file=out_file)
                print("", end='\t', file=out_file)
                print("", end='\n', file=out_file)
            print(el["reagent_name"], end='\t', file=out_file)
            print(el["reagent_role"], end='\t', file=out_file)
            print(el["molar_mass"], end='\t', file=out_file)
            print(el["mass"], end='\t', file=out_file)
            print(el["cc50"], end='\n', file=out_file)

    return 0


@app.route('/data_validation')
def data_validation(metadata):
    filename = metadata['title']
    path_data = os.path.join('data', filename)

    message = inpval.data_validation(path_data)

    return message


@app.route('/calc_combinations')
def calc_combinations(metadata):
    filename = metadata['title']
    path_data = os.path.join('data', filename)
    mesure_variable = metadata['mesure_type']
    path_table, number_of_combinations = gentab.generate_table(path_data, mesure_variable)

    return path_table, number_of_combinations


@app.route('/make_chart')
def make_chart(metadata, path_table):
    filename = metadata['title']
    colormap = metadata['colormap']
    cyt_potential = metadata['cyt_potential']
    path_data = os.path.join('data', filename)

    top_combinations = gench.generate_charts(path_data, path_table, colormap, cyt_potential)

    return top_combinations


@app.route('/display_chart/<filename>')
def display_chart(filename):
    dir_name = filename.rsplit('.', 1)[0]
    path_graphs = os.path.join(path_figures, dir_name)

    colormap_path = 'figures/' + dir_name + '/colormap.png'
    graphs = [colormap_path]

    for el in os.listdir(path_graphs):
        if el == 'colormap.png':
            continue
        elif el.split('.')[-1] not in ('pkl', 'txt'):
            graphs.append('figures/' + dir_name + '/' + el)
    with open(f'static/figures/{filename.rsplit(".")[0]}/top_combinations.txt', 'rt') as f:
        lines = f.readlines()
        filename = lines[0].strip()
        colormap = lines[1].strip()
        cyt_potential = lines[2].strip()
        number_of_combs = lines[3].strip()
    with open(f'static/figures/{filename.rsplit(".")[0]}/file_info.pkl', 'rb') as f:
        top_combinations = pkl.load(f)
    return render_template('display_chart.html', menu=menu, number=number_of_combs, colormap=colormap,
                           filename=filename, top_combinations=top_combinations, cyt_potential=cyt_potential,
                           graphs=graphs)


@app.route("/download/<path:filename>")
def download(filename):
    dir_name = filename.rsplit('.', 1)[0]
    path_folder = os.path.join(path_results, dir_name)
    zip_path = path_folder + '.zip'
    zip_name = dir_name + '.zip'

    with ZipFile(zip_path, "w") as zip_arch:
        for dirpath, _, filenames in os.walk(path_folder):
            for filename_path in filenames:
                print(filename_path)
                if filename_path.split('_')[-1].split('.')[0] == 'comb':
                    filename_path_txt = filename_path.split('_')[0] + '.txt'
                    filename_path_legend = filename_path.split('_')[0] + '_combout.csv'
                    filename_path_csv = filename_path.split('_')[0] + '.csv'
                    if filename_path_csv in os.listdir('data'):
                        zip_arch.write(f'data/{filename_path_csv}')
                    else:
                        zip_arch.write(f'data/{filename_path_txt}')
                        zip_arch.write(f'data/{filename_path_legend}')
                else:
                    zip_arch.write(os.path.join(dirpath, filename_path))
    directory = os.path.join(app.root_path, path_results)
    return send_from_directory(directory, zip_name)


@app.route("/download_file/<path:filename>", methods=['GET', 'POST'])
def download_file(filename):
    directory = os.path.join(app.root_path, path_examples)

    return send_from_directory(directory, filename, as_attachment=True)


def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.errorhandler(404)
def pageNotFound(error):
    app.logger.error(error)
    return render_template('page404.html', title="Page 404", menu=menu), 404


@app.errorhandler(500)
def internal_error(exception):
    app.logger.error(exception)
    return render_template('page404.html', title="Page 404", menu=menu), 500


if __name__ == '__main__':
    log_file = 'flask.log'

    file_handler = RotatingFileHandler(log_file, maxBytes=1024 * 1024 * 100, backupCount=20)
    file_handler.setLevel(logging.ERROR)
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    app.logger.addHandler(file_handler)

    app.run(host="0.0.0.0", port=5001)

