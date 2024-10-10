from flask import Flask, render_template, request, jsonify
import subprocess
import os

app = Flask(__name__)

DEFAULT_DIRECTORY = "/media/anegin97/DATA/DATA"

def execute_script(input_directory, script_name):
    full_path = os.path.join(DEFAULT_DIRECTORY, input_directory)

    if not os.path.isdir(full_path):
        return "Invalid directory path"

    script_dir = os.path.dirname(os.path.realpath(__file__))
    bash_script_path = os.path.join(script_dir, script_name)

    try:
        process = subprocess.Popen(['bash', bash_script_path, full_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        output = stdout.decode('utf-8') + '\n' + stderr.decode('utf-8')
        return output
    except Exception as e:
        return f"Error executing script: {str(e)}"

def get_folders(directory):
    folders = []
    for root, dirs, files in os.walk(directory):
        for dir_name in dirs:
            folders.append(os.path.relpath(os.path.join(root, dir_name), directory))
    return folders

@app.route('/', methods=['GET', 'POST'])
def index():
    output = None
    status = None
    folders = get_folders(DEFAULT_DIRECTORY)
    if request.method == 'POST':
        directory_name = request.form.get('directory_path', '')
        script_name = request.form.get('script_name', '')
        output = execute_script(directory_name, script_name)
        if output:
            status = 'Finished' if "Process finished" in output else 'Error'
    return render_template('analysis.html', folders=folders, output=output, status=status)

@app.route('/execute', methods=['POST'])
def execute_command():
    command = request.json['command']
    return execute_script(command)

@app.route('/samplesheet')
def samplesheet():
    return render_template('samplesheet.html')

@app.route('/visualization')
def visualization():
    return render_template('visualization.html')

@app.route('/database')
def database():
    return render_template('database.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=9999)



