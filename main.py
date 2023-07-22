# -*- coding: utf-8 -*-
from flask import Flask, render_template, request, flash, redirect, url_for, send_from_directory, session
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
        {"name": "About Us", "url": "/about"},
        {"name": "Titorial", "url": "/tutorial"},
        {"name": "Build Chart", "url": "/send_check"}]

colormap_list = [{'name': 'percentile'},
                 {'name': 'linear'}]

cyt_potentials_list = [{'name': 'BF'},
                       {'name': 'CPi'},
                       {'name': 'CPf'},
                       {'name': 'CPf_rel'}]


class ReagentLineForm(Form):
    reagent_name = StringField("Reagent name: ", validators=[DataRequired()], description="Reagent name")
    reagent_role = StringField("Reagent role: ", validators=[DataRequired()], description="Reagent role")
    molar_mass = DecimalField("Molar mass: ",
                              validators=[DataRequired(), NumberRange(0)],
                              description="Molar mass")
    mass = DecimalField("Mass: ",
                        validators=[DataRequired(), NumberRange(0)],
                        description="Mass")
    cc50 = DecimalField("CC50: ",
                        validators=[DataRequired(), NumberRange(0.00001)],
                        description="CC50")


class OneChartForm(FlaskForm):
    filename = StringField("Enter filename: ", validators=[DataRequired()], description="Filename")
    cell_name = StringField("Enter cell name: ", validators=[DataRequired()], description="Cell name")
    reagents_info = FieldList(FormField(ReagentLineForm), min_entries=1)
    products_info = FieldList(FormField(ReagentLineForm), min_entries=1)
    colormap = SelectField("Choose colormap: ")
    cyt_potential = SelectField("Choose cytotoxic potential: ")


class CheckFile(FlaskForm):
    filename = StringField("Enter filename: ", validators=[DataRequired()], description="Filename")


@app.route('/')
@app.route('/home')
def main():
    return render_template('main.html', menu=menu)


@app.route('/tutorial')
def manual():
    return render_template('manual.html', menu=menu)


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
            return redirect(url_for('create_chart', flag='False'), 302)
    return render_template('send_check.html', menu=menu, form=form)


@app.route("/output_file/<filename>", methods=['POST', 'GET'])
def output_file(filename):
    return render_template('output_file.html', menu=menu, filename=filename)


@app.route('/create_chart_<flag>', methods=['POST', 'GET'])
def create_chart(flag):
    if flag == 'True':
        file_info = {'title': session["filename"], 'colormap': session["colormap"],
                     'cyt_potential': session["cyt_potential"]}
    else:
        file_info = {'title': 'None', 'colormap': 'None',
                     'cyt_potential': 0}
        reset_session()

    form = OneChartForm()
    form.colormap.choices = [(el["name"], el["name"]) for el in colormap_list]
    form.cyt_potential.choices = [(el["name"], el["name"]) for el in cyt_potentials_list]

    if 'send_create' in request.form:
        filename = form.filename.data
        meta_filename = filename + '.txt'
        if meta_filename in os.listdir(path_data):
            flash('Sorry, this name is taken. Try another filename')
            return redirect(url_for('create_chart', flag='False'))
        cell_name = form.cell_name.data
        colormap = form.colormap.data
        cyt_potential = form.cyt_potential.data
        reagents_info = form.reagents_info.data
        products_info = form.products_info.data
        try:
            save_chart_data(filename, cell_name, reagents_info, products_info)

            session["filename"] = filename + '.txt'
            session["colormap"] = colormap
            session["cyt_potential"] = cyt_potential
            session["data_upload"] = 1
        except:
            flash('Sorry, this name is taken. You have entered incorrect data')
            return redirect(url_for('create_chart', flag='False'))
    elif "data_upload" in session:
        if 'send_upload' in request.form:
            session["data_upload"] = 1
    else:
        session["data_upload"] = 0
        session["graphs"] = []

    if session["data_upload"] == 1:
            file_info = {'title': session["filename"], 'colormap': session["colormap"],
                         'cyt_potential': session["cyt_potential"]}
            if 'gencombs' in request.form:
                message = data_validation(file_info)
                if message == 'SUCCESS':
                    path_table, number_of_combinations = calc_combinations(file_info)
                    session["path_table"] = path_table
                    session["number_of_combinations"] = number_of_combinations
                    top_combinations = make_chart(file_info, path_table)
                    session["top_combinations"] = top_combinations
                    session["combs_found"] = 1
                else:
                    flash(message)
            return render_template('create_chart.html', menu=menu, colormap_list=colormap_list,
                               cyt_potentials_list=cyt_potentials_list, json=file_info, session=session, form=form)
    return render_template('create_chart.html', menu=menu, colormap_list=colormap_list,
                           cyt_potentials_list=cyt_potentials_list, json=file_info, session=session, form=form)


@app.route('/save_chart_data')
def save_chart_data(filename, cell_name, reagents_info, products_info):
    path_file = os.path.join(path_data, filename + '.txt')

    with open(path_file, "w", encoding='utf-8') as out_file:
        print("Cell", cell_name, sep='\t', end='\n', file=out_file)
        print("Variables", end='\n', file=out_file)
        print("Product variables", end='\n', file=out_file)
        print("Samples", "Abbreviation", "Mr, g*mol-1", "Mass, g", "CC50, mM", sep='\t', end='\n', file=out_file)
        print("Starting materials", end='\n', file=out_file)

        for el in reagents_info:
            print(el["reagent_name"], end='\t', file=out_file)
            print(el["reagent_role"], end='\t', file=out_file)
            print(el["molar_mass"], end='\t', file=out_file)
            print(el["mass"], end='\t', file=out_file)
            print(el["cc50"], end='\n', file=out_file)

        print("Products", end='\n', file=out_file)

        for el in products_info:
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

    path_table, number_of_combinations = gentab.generate_table(path_data)

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
        else:
            graphs.append('figures/' + dir_name + '/' + el)

    return render_template('display_chart.html', menu=menu, graphs=graphs, filename=filename)


@app.route("/download/<path:filename>")
def download(filename):
    dir_name = filename.rsplit('.', 1)[0]
    path_folder = os.path.join(path_results, dir_name)
    zip_path = path_folder + '.zip'
    zip_name = dir_name + '.zip'

    with ZipFile(zip_path, "w") as zip_arch:
        for dirpath, _, filenames in os.walk(path_folder):
            for filename in filenames:
                zip_arch.write(os.path.join(dirpath, filename))

    directory = os.path.join(app.root_path, path_results)

    reset_session()

    return send_from_directory(directory, zip_name)


@app.route("/reset_session")
def reset_session():
    session["data_upload"] = 0
    session["graphs"] = []
    session["filename"] = ''
    session["colormap"] = ''
    session["cyt_potential"] = ''
    session["combs_found"] = 0
    session["path_table"] = ''
    session["charts_built"] = 0
    session["top_combinations"] = ''
    session["number_of_combinations"] = 0


@app.route("/download_file/<path:filename>", methods=['GET', 'POST'])
def download_file(filename):
    directory = os.path.join(app.root_path, path_examples)

    return send_from_directory(directory, filename, as_attachment=True)


def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/uploader', methods=['GET', 'POST'])
def uploader():
    file_info = {'title': 'NONE', 'colormap': 'NONE', 'cyt_potential': 'NONE'}
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
        if file.filename in os.listdir(path_data):
            flash('Sorry, this name is taken. Try another filename')
            return redirect(request.url)
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

            colormap_name = request.form.get('colormap').replace('"', '')
            cyt_potential_name = request.form.get('cyt_potential').replace('"', '')

            session["filename"] = filename
            session["colormap"] = colormap_name
            session["cyt_potential"] = cyt_potential_name
            session["data_upload"] = 1
        else:
            extensions = 'Invalid file extension! The file must have one of the following extensions: '
            len_extensions = len(ALLOWED_EXTENSIONS)
            for count, extension in enumerate(ALLOWED_EXTENSIONS):
                if count == len_extensions - 1:
                    extensions += extension
                else:
                    extensions += extension + ','

            flash(extensions)

            return redirect(url_for('create_chart', flag='False'), 302)
    return redirect(url_for('create_chart', flag='True'), 302)


@app.errorhandler(404)
def pageNotFound(error):
    return render_template('page404.html', title="Page 404", menu=menu), 404


@app.errorhandler(500)
def internal_error(exception):
    app.logger.error(exception)
    return render_template('500.html'), 500


if __name__ == '__main__':
    log_file = 'flask.log'

    file_handler = RotatingFileHandler(log_file, maxBytes=1024 * 1024 * 100, backupCount=20)
    file_handler.setLevel(logging.ERROR)
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    app.logger.addHandler(file_handler)

    app.run(host="0.0.0.0", port=8080)
